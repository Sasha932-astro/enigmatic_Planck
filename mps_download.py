#!/usr/bin/env python3
"""
Download MPS-ATLAS stellar spectra from the Max Planck Digital Library (Edmond).

Dataset: "An extended MPS-ATLAS library of stellar model atmospheres and spectra"
         doi:10.17617/3.NJ56TR  (Witzke et al. 2021; Kostogryz et al. 2022)

The library ships as two ~9.5 GB zip files. Downloading them whole is usually
unnecessary: this script reads the remote zip's central directory with HTTP
range requests and pulls out ONLY the model files you ask for (a few hundred kB
each), so fetching a handful of temperatures takes seconds instead of hours.

USAGE (command line)
--------------------
    # one temperature
    python mps_download.py --teff 5800

    # several temperatures (the 100 K experiment pair)
    python mps_download.py --teff 5700 5800

    # every temperature on the grid, solar metallicity, logg 4.4
    python mps_download.py --all

    # disk-integrated fluxes instead of centre-to-limb intensities
    python mps_download.py --teff 5800 --kind flux

Files land in   mps_atlas/<set>/MH<mh>/teff<T>/logg<logg>/ ,
i.e. exactly the layout the official reader (model_spectra_data.py) expects.

USAGE (as a module)
-------------------
    from mps_download import fetch, read_clv
    fetch([4500, 5450, 5800])                 # download if missing
    wav, I, mu = read_clv(5800)               # wavelengths, I_nu(lam, mu), mu grid

GRID
----
    Teff  3500 ... 9000 K, step 100 K
    logg  3.0, 3.5, 4.0, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 5.0
    M/H   -5.0 ... 1.5
    set1  Grevesse & Sauval (1998) abundances, mixing length 1.25
    set2  Asplund et al. (2009) abundances, Viani et al. (2018) mixing length

UNITS
-----
    clv  : I_nu  [erg s^-1 cm^-2 Hz^-1 ster^-1]  on (n_wav, n_mu)
    flux : F_nu  [erg s^-1 cm^-2 Hz^-1] at 1 AU
    wavelengths in nm (1221 points, 9.1 - 160000 nm)
"""

import argparse
import io
import os
import sys
import urllib.request
import zipfile

# Edmond datafile ids for the two computation sets
DATAFILE_ID = {"set1": 199597, "set2": 199601}
API = "https://edmond.mpg.de/api/access/datafile/{}"
FILENAME = {"clv": "mpsa_intensity_spectra.dat",
            "flux": "mpsa_flux_spectra.dat",
            "model": "mpsa_model_atmosphere.dat"}


class _RemoteFile(io.RawIOBase):
    """Minimal seekable file-like object backed by HTTP range requests."""

    def __init__(self, url, timeout=60):
        self.pos = 0
        req = urllib.request.Request(url, method="HEAD")
        with urllib.request.urlopen(req, timeout=timeout) as r:
            self.size = int(r.headers["Content-Length"])
            self.url = r.url          # resolved URL: avoids re-following redirects

    def seek(self, off, whence=0):
        self.pos = {0: off, 1: self.pos + off, 2: self.size + off}[whence]
        return self.pos

    def tell(self):
        return self.pos

    def readable(self):
        return True

    def seekable(self):
        return True

    def read(self, n=-1):
        if n == -1:
            n = self.size - self.pos
        if n <= 0 or self.pos >= self.size:
            return b""
        end = min(self.pos + n, self.size) - 1
        req = urllib.request.Request(self.url, headers={"Range": f"bytes={self.pos}-{end}"})
        with urllib.request.urlopen(req, timeout=300) as r:
            data = r.read()
        self.pos += len(data)
        return data

    def readinto(self, b):
        data = self.read(len(b))
        b[:len(data)] = data
        return len(data)


def open_remote_zip(set_type="set1"):
    """Open the remote Edmond zip without downloading it."""
    url = API.format(DATAFILE_ID[set_type])
    return zipfile.ZipFile(io.BufferedReader(_RemoteFile(url), buffer_size=4 * 1024 * 1024))


def fetch(teffs, mh=0.0, logg=4.4, set_type="set1", kind="clv",
          outdir="mps_atlas", skip_existing=True, verbose=True):
    """Download the requested spectra. Returns the list of local paths."""
    teffs = [int(t) for t in teffs]
    want = {}
    for t in teffs:
        rel = os.path.join(set_type, f"MH{mh}", f"teff{t}", f"logg{logg}", FILENAME[kind])
        want[rel.replace(os.sep, "/")] = os.path.join(outdir, rel)

    todo = {k: v for k, v in want.items()
            if not (skip_existing and os.path.exists(v))}
    if not todo:
        if verbose:
            print(f"all {len(want)} file(s) already present")
        return list(want.values())

    if verbose:
        print(f"opening remote {set_type}.zip ...")
    z = open_remote_zip(set_type)
    names = set(z.namelist())

    paths = []
    for i, (rel, dest) in enumerate(sorted(todo.items()), 1):
        if rel not in names:
            raise FileNotFoundError(
                f"{rel} is not in the dataset - check teff/logg/mh against the grid")
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        with z.open(rel) as src, open(dest, "wb") as out:
            out.write(src.read())
        paths.append(dest)
        if verbose:
            print(f"  [{i}/{len(todo)}] {dest}")
    return list(want.values())


# ---------------------------------------------------------------- readers
def read_clv(teff, mh=0.0, logg=4.4, set_type="set1", outdir="mps_atlas", mu_subset=None):
    """Return (wav_nm, I_nu[n_wav, n_mu], mu). Downloads the file if missing.

    mu_subset: optional array of mu values; the nearest library columns are
    selected (the library's 24 mu values contain most common grids exactly).
    """
    import numpy as np
    path = fetch([teff], mh, logg, set_type, "clv", outdir, verbose=False)[0]
    data = np.loadtxt(path, skiprows=2)
    wav, I = data[:, 0], data[:, 1:]
    with open(path) as f:
        f.readline()
        mu = np.array(f.readline().split()[2:], dtype=float)
    if mu_subset is not None:
        idx = [int(np.argmin(abs(mu - m))) for m in mu_subset]
        I, mu = I[:, idx], mu[idx]
    return wav, I, mu


def read_flux(teff, mh=0.0, logg=4.4, set_type="set1", outdir="mps_atlas"):
    """Return (wav_nm, F_nu) disk-integrated flux at 1 AU. Downloads if missing."""
    import numpy as np
    path = fetch([teff], mh, logg, set_type, "flux", outdir, verbose=False)[0]
    return np.loadtxt(path, skiprows=1, unpack=True)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--teff", type=int, nargs="+", help="effective temperature(s) in K")
    p.add_argument("--all", action="store_true", help="every Teff from 3500 to 9000 K")
    p.add_argument("--mh", type=float, default=0.0, help="metallicity [M/H] (default 0.0)")
    p.add_argument("--logg", type=float, default=4.4, help="log g (default 4.4)")
    p.add_argument("--set", dest="set_type", default="set1", choices=["set1", "set2"])
    p.add_argument("--kind", default="clv", choices=["clv", "flux", "model"],
                   help="clv = centre-to-limb intensities (default), flux = disk-integrated")
    p.add_argument("--outdir", default="mps_atlas")
    a = p.parse_args(argv)

    teffs = list(range(3500, 9100, 100)) if a.all else a.teff
    if not teffs:
        p.error("give --teff T [T ...] or --all")
    fetch(teffs, a.mh, a.logg, a.set_type, a.kind, a.outdir)
    print("done")


if __name__ == "__main__":
    sys.exit(main())
