#!/usr/bin/env python3
"""Final real-spectra comparison for enigmatic_Planck.

Part 1
------
Rerun the canonical MPS-ATLAS 5800->5700 K experiment and the MURaM
G2->K0/K4/M0/M2 spot-surrogate experiments, then add the 100 K MPS-ATLAS
curve/reference to all class and mixture plots.

Part 2
------
Use the real G2V quiet/spot/penumbra/umbra disk-integrated spectra from
Edmond doi:10.17617/3.HS2EE6, datafile 252373.  The file is disk-integrated,
not mu-resolved, so the appropriate geometry-free filling-factor experiment is

    F(lambda,t) = F_q(lambda) + alpha(t) [F_feature(lambda)-F_q(lambda)]

which gives exactly

    a(lambda) = [(F_feature-F_q)/F_q]
                / [ integral(F_feature-F_q) dlambda / integral(F_q) dlambda ].

We score this response against Planck(5780 K) over 400--1600 nm and optimize
physical umbra/penumbra mixtures with non-negative fractions summing to one.

All spectra are downloaded/read directly from the published source datasets.
No digitization, spectral emulation, or interpolation between spectral classes
is used.
"""

from __future__ import annotations

import itertools
import json
import pathlib
import urllib.request

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize, minimize_scalar

import experiment_100K as base
import experiment_full_comparison as fc
import experiment_muram_spot_mixtures as ms

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "scenarios"
OUT.mkdir(exist_ok=True)
SUNSPOT_CACHE = HERE / "sunspot_data"
SUNSPOT_CACHE.mkdir(exist_ok=True)

LO, HI = 400.0, 1600.0
T_PLANCK = 5780.0
HS2EE6_G2_FILE_ID = 252373


def planck(wav):
    return base.G["planck_temperature_slope"](wav, T_PLANCK)


# ----------------------------------------------------------------------
# Part 1: rerun MPS + MURaM and remake previous plots with 100 K added.
# ----------------------------------------------------------------------
def rerun_previous_experiments():
    # Canonical MPS 100 K pair.
    wav_mps, i_mps_qs, i_mps_feat = fc._load_mps_100k_pair()
    r100 = fc._run_pair("MPS-ATLAS 5800->5700 K", wav_mps, i_mps_qs, i_mps_feat)

    # G2 hydro quiet star and all MURaM hydro spot surrogates.
    wav_q, i_q, teff_q, mh_q, logg_q = fc._load_muram_state("G2_MH_00", "hydro")
    single = {}
    spectra = {}
    metadata = {
        "quiet": {"star": "G2_MH_00 hydro", "teff_K": teff_q, "MH": mh_q, "logg": logg_q},
        "components": {},
    }
    for cls in ms.MIX_CLASSES:
        wav, I, teff, mh, logg = ms.load_hydro(cls)
        if not np.allclose(wav, wav_q):
            raise ValueError(f"MURaM wavelength grid mismatch for {cls}")
        spectra[cls] = I
        metadata["components"][cls] = {"teff_K": teff, "MH": mh, "logg": logg}
        if cls in ms.PRIMARY_CLASSES:
            single[cls] = fc._run_pair(f"MURaM G2->{cls} hydro", wav_q, i_q, I)

    component_spectra = [spectra[c] for c in ms.MIX_CLASSES]
    E_qs, delta_basis = ms.prepare_caseB_component_basis(wav_q, i_q, component_spectra)

    pure_fits = {}
    for j, cls in enumerate(ms.MIX_CLASSES):
        w = np.zeros(len(ms.MIX_CLASSES)); w[j] = 1.0
        pure_fits[cls] = ms.fit_caseB_weights(wav_q, E_qs, delta_basis, w)

    grid, curves, pairs = ms.pair_scan(wav_q, E_qs, delta_basis, ms.MIX_CLASSES)
    best_pair = pairs[0]

    rng = np.random.default_rng(20260811)
    multi = []
    for n in (3, 4, 5):
        for subset in itertools.combinations(ms.MIX_CLASSES, n):
            multi.append(ms.optimize_subset(wav_q, E_qs, delta_basis, ms.MIX_CLASSES, subset, rng))
    multi.sort(key=lambda x: x["rms"])
    best_multi = multi[0]

    return {
        "r100": r100,
        "wav": wav_q,
        "single": single,
        "pure_fits": pure_fits,
        "grid": grid,
        "curves": curves,
        "pairs": pairs,
        "best_pair": best_pair,
        "best_multi": best_multi,
        "metadata": metadata,
    }


def plot_classes_with_100(runs, case: str, ratio: bool, path: pathlib.Path):
    fig, ax = plt.subplots(figsize=(10, 6))
    r100 = runs["r100"]
    w = r100["wav"]
    m = r100[case]["valid"] & (w >= LO) & (w <= HI)
    y = r100[case]["slope"][m]
    if ratio:
        y = y / planck(w[m])
    ax.plot(w[m], y, lw=1.65,
            label=f"MPS-ATLAS 5800->5700 K (RMS {r100[f'{case}_rms']:.3f})")

    for cls in ms.PRIMARY_CLASSES:
        r = runs["single"][cls]
        w = r["wav"]
        m = r[case]["valid"] & (w >= LO) & (w <= HI)
        y = r[case]["slope"][m]
        if ratio:
            y = y / planck(w[m])
        ax.plot(w[m], y, lw=1.45,
                label=f"MURaM G2->{cls} hydro (RMS {r[f'{case}_rms']:.3f})")

    if ratio:
        ax.axhline(1.0, ls="--", lw=1.8, label="Planck 5780 K")
        ax.set_ylabel("a(lambda) / Planck(5780 K)")
        ax.set_title(f"MPS 100 K + MURaM spot surrogates, Case {case}: ratio to Planck")
    else:
        wp = np.linspace(LO, HI, 1201)
        ax.plot(wp, planck(wp), "--", lw=1.8, label="Planck 5780 K")
        ax.set_ylabel("a(lambda)")
        ax.set_title(f"MPS 100 K + MURaM spot surrogates, Case {case}")
    ax.set_xlabel("Wavelength (nm)")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8.5)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_pair_scan_with_100(runs, path):
    fig, ax = plt.subplots(figsize=(11, 6.5))
    for name, rms in runs["curves"].items():
        ax.plot(runs["grid"], rms, lw=1.2, label=name)
    ax.axhline(runs["r100"]["B_rms"], ls="--", lw=1.8,
               label=f"MPS-ATLAS 100 K Case B (RMS {runs['r100']['B_rms']:.3f})")
    ax.set_xlabel("Area fraction of first MURaM class in pair")
    ax.set_ylabel("RMS distance from Planck(5780 K)")
    ax.set_title("Disk-centre MURaM two-component mixtures + MPS 100 K reference")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_best_with_100(runs, path, ratio=False):
    wav = runs["wav"]
    m = (wav >= LO) & (wav <= HI)
    fig, ax = plt.subplots(figsize=(10, 6))

    r100 = runs["r100"]
    wm = r100["wav"]
    mm = r100["B"]["valid"] & (wm >= LO) & (wm <= HI)
    y = r100["B"]["slope"][mm]
    if ratio:
        y = y / planck(wm[mm])
    ax.plot(wm[mm], y, lw=1.6,
            label=f"MPS 100 K Case B (RMS {r100['B_rms']:.3f})")

    best_pure_name, best_pure = min(runs["pure_fits"].items(), key=lambda kv: kv[1]["rms"])
    entries = [
        (f"best pure: {best_pure_name} (RMS {best_pure['rms']:.3f})", best_pure),
        (f"best pair: {'+'.join(runs['best_pair']['classes'])} (RMS {runs['best_pair']['rms']:.3f})",
         runs["best_pair"]["fit"]),
        (f"best multi: {'+'.join(runs['best_multi']['classes'])} (RMS {runs['best_multi']['rms']:.3f})",
         runs["best_multi"]["fit"]),
    ]
    for label, fit in entries:
        y = fit["slope"][m]
        if ratio:
            y = y / planck(wav[m])
        ax.plot(wav[m], y, lw=1.55, label=label)

    if ratio:
        ax.axhline(1.0, ls="--", lw=1.8, label="Planck 5780 K")
        ax.set_ylabel("a(lambda) / Planck(5780 K)")
        ax.set_title("MPS 100 K + best MURaM mixtures, disk centre: ratio to Planck")
    else:
        ax.plot(wav[m], planck(wav[m]), "--", lw=1.8, label="Planck 5780 K")
        ax.set_ylabel("a(lambda)")
        ax.set_title("MPS 100 K + best MURaM mixtures, disk centre")
    ax.set_xlabel("Wavelength (nm)")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8.5)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


# ----------------------------------------------------------------------
# Part 2: G2V quiet / penumbra / umbra from doi:10.17617/3.HS2EE6.
# ----------------------------------------------------------------------
def download_hs2ee6_g2():
    dest = SUNSPOT_CACHE / "G_star_disk_integrated_flux.txt"
    if not dest.exists():
        url = f"https://edmond.mpg.de/api/access/datafile/{HS2EE6_G2_FILE_ID}"
        print(f"Downloading G_star_disk_integrated_flux.txt from Edmond datafile {HS2EE6_G2_FILE_ID} ...")
        urllib.request.urlretrieve(url, dest)
    return dest


def read_hs2ee6_g2(path):
    rows = []
    for line in pathlib.Path(path).read_text(errors="replace").splitlines():
        parts = line.replace(",", " ").split()
        if len(parts) < 5:
            continue
        try:
            vals = [float(parts[i]) for i in range(5)]
        except ValueError:
            continue
        rows.append(vals)
    a = np.asarray(rows, dtype=float)
    if a.ndim != 2 or a.shape[1] != 5 or a.shape[0] < 20:
        raise ValueError(f"Unexpected HS2EE6 G-star table shape {a.shape}")
    # wavelength, quiet, published spot, penumbra, umbra
    order = np.argsort(a[:, 0])
    a = a[order]
    return a[:, 0], a[:, 1], a[:, 2], a[:, 3], a[:, 4]


def flux_response(wav, quiet, feature):
    valid = np.isfinite(wav) & np.isfinite(quiet) & np.isfinite(feature) & (quiet > 0)
    qbol = np.trapezoid(quiet[valid], wav[valid])
    d = feature - quiet
    db = np.trapezoid(d[valid], wav[valid])
    if not np.isfinite(db) or abs(db) < 1e-30 * abs(qbol):
        raise ValueError("Bolometric contrast is zero/invalid")
    slope = np.full_like(wav, np.nan, dtype=float)
    slope[valid] = (d[valid] / quiet[valid]) / (db / qbol)
    m = valid & (wav >= LO) & (wav <= HI)
    rms = float(np.sqrt(np.mean((slope[m] - planck(wav[m])) ** 2)))
    bestT = base.best_fit_planck_T(wav, slope, valid)
    return {"slope": slope, "valid": valid, "rms": rms, "bestT": bestT,
            "qbol": float(qbol), "dbol": float(db)}


def hs2ee6_experiment():
    path = download_hs2ee6_g2()
    wav, quiet, published_spot, penumbra, umbra = read_hs2ee6_g2(path)

    fits = {
        "penumbra": flux_response(wav, quiet, penumbra),
        "umbra": flux_response(wav, quiet, umbra),
        "published spot": flux_response(wav, quiet, published_spot),
    }

    def blend(fu):
        return fu * umbra + (1.0 - fu) * penumbra

    def obj(fu):
        return flux_response(wav, quiet, blend(float(fu)))["rms"]

    scan_f = np.linspace(0.0, 1.0, 1001)
    scan_rms = np.array([obj(f) for f in scan_f])
    k = int(np.argmin(scan_rms))
    lo = scan_f[max(0, k - 5)]
    hi = scan_f[min(scan_f.size - 1, k + 5)]
    opt = minimize_scalar(obj, bounds=(lo, hi), method="bounded", options={"xatol": 1e-10})
    fu = float(opt.x)
    best_feature = blend(fu)
    best = flux_response(wav, quiet, best_feature)
    fits["optimized U+P"] = best

    # Infer whether the dataset's published 'spot' column is itself a fixed U/P mixture.
    # Use the spectrum directly, not the Planck score, for this diagnostic.
    scale = np.maximum(np.abs(published_spot), 1e-300)
    def spot_match(f):
        model = blend(float(f))
        return float(np.sqrt(np.mean(((model - published_spot) / scale) ** 2)))
    sm = minimize_scalar(spot_match, bounds=(0.0, 1.0), method="bounded")
    published_fu = float(sm.x)
    published_match_rms = float(sm.fun)

    # Plots.
    m = (wav >= LO) & (wav <= HI)
    fig, ax = plt.subplots(figsize=(10, 6))
    for name in ("penumbra", "umbra", "published spot", "optimized U+P"):
        f = fits[name]
        label = name
        if name == "optimized U+P":
            label += f" ({fu:.3f} U + {1-fu:.3f} P, RMS {f['rms']:.3f})"
        else:
            label += f" (RMS {f['rms']:.3f})"
        ax.plot(wav[m], f["slope"][m], lw=1.55, label=label)
    ax.plot(wav[m], planck(wav[m]), "--", lw=1.8, label="Planck 5780 K")
    ax.set(xlabel="Wavelength (nm)", ylabel="a(lambda)",
           title="HS2EE6 G2V disk-integrated spot spectra: geometry-free response")
    ax.grid(alpha=0.25); ax.legend(fontsize=8.5)
    fig.tight_layout(); fig.savefig(OUT / "_hs2ee6_gspot_components_a.png", dpi=160); plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 6))
    for name in ("penumbra", "umbra", "published spot", "optimized U+P"):
        f = fits[name]
        label = name
        if name == "optimized U+P":
            label += f" ({fu:.3f} U + {1-fu:.3f} P, RMS {f['rms']:.3f})"
        else:
            label += f" (RMS {f['rms']:.3f})"
        ax.plot(wav[m], f["slope"][m] / planck(wav[m]), lw=1.55, label=label)
    ax.axhline(1.0, ls="--", lw=1.8, label="Planck 5780 K")
    ax.set(xlabel="Wavelength (nm)", ylabel="a(lambda) / Planck(5780 K)",
           title="HS2EE6 G2V spot spectra: ratio to Planck")
    ax.grid(alpha=0.25); ax.legend(fontsize=8.5)
    fig.tight_layout(); fig.savefig(OUT / "_hs2ee6_gspot_components_ratio.png", dpi=160); plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.plot(scan_f, scan_rms, lw=1.6)
    ax.axvline(fu, ls="--", lw=1.4,
               label=f"best f_umbra={fu:.4f}, RMS={best['rms']:.4f}")
    ax.set(xlabel="Umbra area fraction in umbra+penumbra blend",
           ylabel="RMS distance from Planck(5780 K)",
           title="HS2EE6 G2V: optimize the umbra/penumbra area ratio")
    ax.grid(alpha=0.25); ax.legend()
    fig.tight_layout(); fig.savefig(OUT / "_hs2ee6_gspot_fraction_scan.png", dpi=160); plt.close(fig)

    result = {
        "source": "doi:10.17617/3.HS2EE6",
        "datafile_id": HS2EE6_G2_FILE_ID,
        "file": "G_star_disk_integrated_flux.txt",
        "wavelength_min_nm": float(wav.min()),
        "wavelength_max_nm": float(wav.max()),
        "n_wavelength": int(wav.size),
        "normalization": "disk-integrated filling-factor response; bolometric integral over full supplied wavelength grid",
        "fits": {name: {"rms": f["rms"], "bestT_K": f["bestT"],
                         "qbol": f["qbol"], "dbol": f["dbol"]}
                 for name, f in fits.items()},
        "optimized_umbra_fraction": fu,
        "optimized_penumbra_fraction": 1.0 - fu,
        "published_spot_implied_umbra_fraction": published_fu,
        "published_spot_relative_spectral_match_rms": published_match_rms,
    }

    np.savez_compressed(
        OUT / "_hs2ee6_gspot_arrays.npz",
        wav=wav, quiet=quiet, published_spot=published_spot,
        penumbra=penumbra, umbra=umbra, optimized_feature=best_feature,
        penumbra_a=fits["penumbra"]["slope"],
        umbra_a=fits["umbra"]["slope"],
        published_spot_a=fits["published spot"]["slope"],
        optimized_a=best["slope"],
        scan_f_umbra=scan_f, scan_rms=scan_rms,
    )
    return result


def main():
    print("Final comparison using only real source spectra\n")
    runs = rerun_previous_experiments()

    plot_classes_with_100(runs, "A", False, OUT / "_with100K_classes_caseA_a.png")
    plot_classes_with_100(runs, "A", True, OUT / "_with100K_classes_caseA_ratio.png")
    plot_classes_with_100(runs, "B", False, OUT / "_with100K_classes_caseB_a.png")
    plot_classes_with_100(runs, "B", True, OUT / "_with100K_classes_caseB_ratio.png")
    plot_pair_scan_with_100(runs, OUT / "_with100K_mixture_pair_scan.png")
    plot_best_with_100(runs, OUT / "_with100K_mixture_best_a.png", ratio=False)
    plot_best_with_100(runs, OUT / "_with100K_mixture_best_ratio.png", ratio=True)

    hs = hs2ee6_experiment()

    summary = {
        "mps_100K": {
            "caseA_rms": runs["r100"]["A_rms"], "caseB_rms": runs["r100"]["B_rms"],
            "caseA_bestT_K": runs["r100"]["A_bestT"], "caseB_bestT_K": runs["r100"]["B_bestT"],
        },
        "muram_best_pair": {
            "classes": runs["best_pair"]["classes"],
            "weights": runs["best_pair"]["weights"],
            "rms": runs["best_pair"]["rms"],
            "bestT_K": runs["best_pair"]["bestT"],
        },
        "hs2ee6": hs,
    }
    with open(OUT / "_final_sunspot_comparison_results.json", "w") as f:
        json.dump(summary, f, indent=2)

    print("MPS 100 K: A RMS={:.5f}, B RMS={:.5f}".format(runs["r100"]["A_rms"], runs["r100"]["B_rms"]))
    print("MURaM best class-mixture: {} weights={} RMS={:.5f}".format(
        "+".join(runs["best_pair"]["classes"]), runs["best_pair"]["weights"], runs["best_pair"]["rms"]))
    print("\nHS2EE6 G2V disk-integrated spectra:")
    for name, f in hs["fits"].items():
        print(f"  {name:<16s} RMS={f['rms']:.5f}  bestT={f['bestT_K']:.0f} K")
    print("  optimized U/P: f_umbra={:.6f}, f_penumbra={:.6f}".format(
        hs["optimized_umbra_fraction"], hs["optimized_penumbra_fraction"]))
    print("  published spot ~= f_umbra={:.6f}, relative spectral mismatch={:.3e}".format(
        hs["published_spot_implied_umbra_fraction"], hs["published_spot_relative_spectral_match_rms"]))
    print(f"\nSaved final plots and arrays to {OUT}")


if __name__ == "__main__":
    main()
