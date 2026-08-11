#!/usr/bin/env python3
"""Reproducible comparison of the canonical MPS-ATLAS 100 K experiment
with MURaM/MPS-ATLAS G2->K0 hydro and SSD spectral pairs.

Everything is recomputed from source spectra. No curve digitization is used.

Outputs (in scenarios/):
  _compare_caseA_a.png
  _compare_caseA_ratio.png
  _compare_caseB_a.png
  _compare_caseB_ratio.png
  _compare_full_arrays.npz
  _compare_full_results.json

MURaM source dataset:
  doi:10.17617/3.FBTIYY
  G2_MH_00.h5, Edmond datafile 344883
  K0_MH_00.h5, Edmond datafile 344889
"""

from __future__ import annotations

import json
import pathlib
import urllib.request

import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import experiment_100K as base

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "scenarios"
OUT.mkdir(exist_ok=True)
MURAM_CACHE = HERE / "muram_data"
MURAM_CACHE.mkdir(exist_ok=True)

LO, HI = 400.0, 1600.0
T_PLANCK = 5780.0
RADIUS = 10.0
N_STEPS = 120

MURAM_FILES = {
    "G2_MH_00": (344883, "G2_MH_00.h5"),
    "K0_MH_00": (344889, "K0_MH_00.h5"),
}


def _download_edmond(datafile_id: int, filename: str) -> pathlib.Path:
    dest = MURAM_CACHE / filename
    if not dest.exists():
        url = f"https://edmond.mpg.de/api/access/datafile/{datafile_id}"
        print(f"Downloading {filename} from Edmond datafile {datafile_id} ...")
        urllib.request.urlretrieve(url, dest)
    return dest


def _load_muram_state(star: str, magnetization: str):
    """Return wavelength, intensities on project MU_GRID, Teff.

    The MURaM library supplies mu = 0.1..1.0.  The project grid additionally
    contains mu=0.05.  For that single near-limb point we use constant
    extrapolation from mu=0.1, equivalent to clipping mu<0.1 to the nearest
    available library angle.  All other project mu nodes are exact library
    nodes.
    """
    datafile_id, filename = MURAM_FILES[star]
    path = _download_edmond(datafile_id, filename)
    with h5py.File(path, "r") as h:
        g = h[star]
        wav = g["wavelengths"][:].astype(float)
        mu_src = g["mu"][:].astype(float)                # 0.1 ... 1.0
        mags = [x.decode() if isinstance(x, bytes) else str(x)
                for x in g["magnetizations"][:]]
        if magnetization not in mags:
            raise ValueError(f"{magnetization} not in {star}: {mags}")
        k = mags.index(magnetization)
        intensity_src = g["spectra"][k].astype(float)    # (lambda, mu increasing)
        teff = float(g["teff"][k])
        mh = float(g.attrs["MH"])
        logg = float(g.attrs["logg"])

    # model.ipynb expects descending mu in exactly base.MU order.
    assert np.allclose(mu_src[::-1], base.MU[:-1]), (
        f"Unexpected MURaM mu grid {mu_src}; project grid is {base.MU}")
    intensity_desc = intensity_src[:, ::-1]               # 1.0 ... 0.1
    intensity = np.column_stack([intensity_desc, intensity_desc[:, -1]])
    assert intensity.shape[1] == base.MU.size
    return wav, intensity, teff, mh, logg


def _load_mps_100k_pair():
    wav_qs, i_qs = base.clv_si(5800)
    wav_feat, i_feat = base.clv_si(5700)
    if not np.allclose(wav_qs, wav_feat):
        raise ValueError("MPS-ATLAS 5800/5700 wavelength grids differ")
    return wav_qs, i_qs, i_feat


def _run_pair(name: str, wav: np.ndarray, i_qs: np.ndarray, i_feat: np.ndarray):
    """Run exactly the same Case A and Case B geometry as experiment_100K.py."""
    phi, lon, dphi, dlon = base.G["build_sphere_grid"](n_lat=180, n_lon=360)
    mu0, vis0 = base.G["compute_mu"](phi, lon, B0=0.0, lambda0=0.0)
    dOm0 = base.G["pixel_solid_angle"](phi, mu0, vis0, dphi, dlon)
    E_qs = (base.G["interp_I_vs_mu"](i_qs, mu0, mu_grid=base.MU)
            * dOm0.ravel()[None, :]).sum(axis=1)

    # Case A: normal rotating transit.
    dist = base.G["build_distribution"](
        phi, lon,
        [{"type": "umbra", "shape": "circle", "lat": 0.0,
          "lon": -100.0, "radius": RADIUS}],
    )
    _, E_A, _ = base.G["simulate_rotation_multi"](
        wav, i_qs, {"umbra": i_feat}, dist,
        phi, lon, dphi, dlon, E_qs,
        n_steps=N_STEPS, rotation_period_days=27.0, B0=0.0,
    )
    fitA = base.G["fit_ssi_vs_tsi"](wav, E_A, E_qs)

    # Case B: feature grows/shrinks at disk centre, no rotation.
    radii = RADIUS * np.sin(np.pi * np.arange(N_STEPS) / (N_STEPS - 1))
    E_B = np.empty((N_STEPS, wav.size))
    for i, r in enumerate(radii):
        if r < 0.5:
            E_B[i] = E_qs
            continue
        mask = base.G["make_circular_patch_mask"](
            phi, lon, lat0_deg=0.0, lon0_deg=0.0, radius_deg=r)
        E_B[i] = E_qs + base.G["feature_delta_irradiance"](
            i_qs, i_feat, mu0, vis0, dOm0, mask)
    fitB = base.G["fit_ssi_vs_tsi"](wav, E_B, E_qs)

    out = {"name": name, "wav": wav, "E_qs": E_qs, "A": fitA, "B": fitB}
    for case in ("A", "B"):
        fit = out[case]
        m = fit["valid"] & (wav >= LO) & (wav <= HI)
        out[f"{case}_rms"] = base.distance_from_planck(
            wav, fit["slope"], fit["valid"], T=T_PLANCK)
        out[f"{case}_bestT"] = base.best_fit_planck_T(
            wav, fit["slope"], fit["valid"])
        out[f"{case}_corr"] = float(np.nanmean(fit["corr"][m]))
    return out


def _planck(wav):
    return base.G["planck_temperature_slope"](wav, T_PLANCK)


def _plot_absolute(results, case: str, path: pathlib.Path):
    fig, ax = plt.subplots(figsize=(10, 6))
    for r in results:
        wav = r["wav"]
        m = r[case]["valid"] & (wav >= LO) & (wav <= HI)
        ax.plot(wav[m], r[case]["slope"][m], lw=1.6,
                label=f"{r['name']}  (RMS {r[f'{case}_rms']:.3f})")
    wavp = np.linspace(LO, HI, 1201)
    ax.plot(wavp, _planck(wavp), "--", lw=1.8, label="Planck 5780 K")
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("a(lambda)")
    ax.set_title(f"Case {case}: direct comparison from rerun spectra")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _plot_ratio(results, case: str, path: pathlib.Path):
    fig, ax = plt.subplots(figsize=(10, 6))
    for r in results:
        wav = r["wav"]
        m = r[case]["valid"] & (wav >= LO) & (wav <= HI)
        ratio = r[case]["slope"][m] / _planck(wav[m])
        ax.plot(wav[m], ratio, lw=1.6,
                label=f"{r['name']}  (RMS {r[f'{case}_rms']:.3f})")
    ax.axhline(1.0, ls="--", lw=1.8, label="Planck 5780 K")
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("a(lambda) / Planck(5780 K)")
    ax.set_title(f"Case {case}: ratio to Planck from rerun spectra")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main():
    print("Fresh full comparison: no digitization\n")

    # 1. Canonical MPS-ATLAS 100 K pair.
    wav_mps, mps_qs, mps_feat = _load_mps_100k_pair()
    r100 = _run_pair("MPS-ATLAS 5800->5700 K", wav_mps, mps_qs, mps_feat)

    # 2. MURaM/MPS-ATLAS hydro pair, G2 quiet -> K0 cooler feature.
    wav_g_h, g_h, teff_g_h, mh_g, logg_g = _load_muram_state("G2_MH_00", "hydro")
    wav_k_h, k_h, teff_k_h, mh_k, logg_k = _load_muram_state("K0_MH_00", "hydro")
    if not np.allclose(wav_g_h, wav_k_h):
        raise ValueError("MURaM G2/K0 hydro wavelength grids differ")
    rh = _run_pair("MURaM G2->K0 hydro", wav_g_h, g_h, k_h)

    # 3. MURaM/MPS-ATLAS SSD pair, G2 SSD quiet -> K0 SSD feature.
    wav_g_s, g_s, teff_g_s, _, _ = _load_muram_state("G2_MH_00", "ssd")
    wav_k_s, k_s, teff_k_s, _, _ = _load_muram_state("K0_MH_00", "ssd")
    if not np.allclose(wav_g_s, wav_k_s):
        raise ValueError("MURaM G2/K0 SSD wavelength grids differ")
    rs = _run_pair("MURaM G2->K0 SSD", wav_g_s, g_s, k_s)

    results = [r100, rh, rs]

    # Numerical report.
    metadata = {
        "MPS-ATLAS 100K": {"quiet_teff_K": 5800.0, "feature_teff_K": 5700.0},
        "MURaM hydro": {
            "quiet": "G2_MH_00 hydro", "quiet_teff_K": teff_g_h,
            "feature": "K0_MH_00 hydro", "feature_teff_K": teff_k_h,
            "MH": mh_g, "quiet_logg": logg_g, "feature_logg": logg_k,
        },
        "MURaM SSD": {
            "quiet": "G2_MH_00 ssd", "quiet_teff_K": teff_g_s,
            "feature": "K0_MH_00 ssd", "feature_teff_K": teff_k_s,
            "MH": mh_g, "quiet_logg": logg_g, "feature_logg": logg_k,
        },
        "MURaM_mu_note": "library mu=0.1..1.0; project mu=0.05 is constant-extrapolated from mu=0.1",
        "metrics": {},
    }

    print("model                                Case A RMS   Case B RMS   A bestT   B bestT   A corr    B corr")
    print("-" * 105)
    for r in results:
        metadata["metrics"][r["name"]] = {
            "caseA_rms": r["A_rms"], "caseB_rms": r["B_rms"],
            "caseA_bestT_K": r["A_bestT"], "caseB_bestT_K": r["B_bestT"],
            "caseA_mean_corr": r["A_corr"], "caseB_mean_corr": r["B_corr"],
        }
        print(f"{r['name']:<36s} {r['A_rms']:10.4f} {r['B_rms']:12.4f} "
              f"{r['A_bestT']:9.0f} {r['B_bestT']:9.0f} {r['A_corr']:9.4f} {r['B_corr']:9.4f}")

    with open(OUT / "_compare_full_results.json", "w") as f:
        json.dump(metadata, f, indent=2)

    # Raw arrays, so every plotted point is preserved.
    np.savez_compressed(
        OUT / "_compare_full_arrays.npz",
        wav_mps=r100["wav"],
        mps_A_slope=r100["A"]["slope"], mps_A_corr=r100["A"]["corr"], mps_A_valid=r100["A"]["valid"],
        mps_B_slope=r100["B"]["slope"], mps_B_corr=r100["B"]["corr"], mps_B_valid=r100["B"]["valid"],
        wav_muram=rh["wav"],
        hydro_A_slope=rh["A"]["slope"], hydro_A_corr=rh["A"]["corr"], hydro_A_valid=rh["A"]["valid"],
        hydro_B_slope=rh["B"]["slope"], hydro_B_corr=rh["B"]["corr"], hydro_B_valid=rh["B"]["valid"],
        ssd_A_slope=rs["A"]["slope"], ssd_A_corr=rs["A"]["corr"], ssd_A_valid=rs["A"]["valid"],
        ssd_B_slope=rs["B"]["slope"], ssd_B_corr=rs["B"]["corr"], ssd_B_valid=rs["B"]["valid"],
    )

    _plot_absolute(results, "A", OUT / "_compare_caseA_a.png")
    _plot_ratio(results, "A", OUT / "_compare_caseA_ratio.png")
    _plot_absolute(results, "B", OUT / "_compare_caseB_a.png")
    _plot_ratio(results, "B", OUT / "_compare_caseB_ratio.png")

    print(f"\nSaved raw arrays and comparison plots to {OUT}")


if __name__ == "__main__":
    main()
