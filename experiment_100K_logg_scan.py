#!/usr/bin/env python3
"""100 K MPS-ATLAS experiment with lower gravity in the cooler component.

Quiet model is fixed at Teff=5800 K, logg=4.4.  The cooler model is fixed at
Teff=5700 K and scanned over logg = 4.4, 4.3, 4.2, 4.0, 3.5, 3.0.

The calculation repeats the two cases of experiment_100K.py:
  A: a circular feature rotates across the disk.
  B: the feature grows/shrinks at disk centre (geometry-free limit).

All spectra are from one MPS-ATLAS family: set1, [M/H]=0.0.
Outputs are written to scenarios/_100K_logg_scan_*.
"""
import json
import pathlib

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar

from mps_download import fetch

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "scenarios"
OUT.mkdir(exist_ok=True)

T_QS = 5800
LOGG_QS = 4.4
T_FEATURE = 5700
FEATURE_LOGGS = [4.4, 4.3, 4.2, 4.0, 3.5, 3.0]
T_PLANCK = 5780.0
LO, HI = 400.0, 1600.0
RADIUS = 10.0

FUNCS = ("read_mu_intensity_txt", "inu_cgs_to_ilambda_si_per_nm",
         "build_sphere_grid", "compute_mu", "pixel_solid_angle",
         "interp_I_vs_mu", "make_circular_patch_mask", "build_distribution",
         "feature_delta_irradiance", "simulate_rotation_multi",
         "fit_ssi_vs_tsi", "planck_temperature_slope")
CONSTS = ("MU_GRID =", "C_CGS =")
G = {"np": np}
for cell in json.loads((HERE / "model.ipynb").read_text())["cells"]:
    src = "".join(cell["source"])
    if cell["cell_type"] != "code":
        continue
    if not (any(f"def {n}" in src for n in FUNCS) or any(c in src for c in CONSTS)):
        continue
    src = src.split("# ---------- example")[0]
    exec(compile(src, "<model.ipynb>", "exec"), G)
MU = G["MU_GRID"]


def clv_si(teff, logg):
    path = fetch([teff], mh=0.0, logg=logg, set_type="set1", kind="clv", verbose=False)[0]
    data = np.loadtxt(path, skiprows=2)
    wav, I = data[:, 0], data[:, 1:]
    with open(path) as f:
        f.readline()
        mu24 = np.array(f.readline().split()[2:], dtype=float)
    idx = [int(np.argmin(abs(mu24 - m))) for m in MU]
    assert np.allclose(mu24[idx], MU)
    return wav, G["inu_cgs_to_ilambda_si_per_nm"](I[:, idx], wav)


def band_mask(wav, fit):
    return fit["valid"] & (wav >= LO) & (wav <= HI)


def distance_from_planck(wav, fit):
    m = band_mask(wav, fit)
    ap = G["planck_temperature_slope"](wav[m], T_PLANCK)
    return float(np.sqrt(np.mean((fit["slope"][m] - ap) ** 2)))


def best_fit_planck_T(wav, fit):
    m = band_mask(wav, fit)
    r = minimize_scalar(
        lambda T: np.sqrt(np.mean((fit["slope"][m] - G["planck_temperature_slope"](wav[m], T)) ** 2)),
        bounds=(3500, 9000), method="bounded")
    return float(r.x)


def main():
    print(f"QS: {T_QS} K logg={LOGG_QS}; feature: {T_FEATURE} K; feature logg scan {FEATURE_LOGGS}")
    wav, i_qs = clv_si(T_QS, LOGG_QS)

    phi, lon, dphi, dlon = G["build_sphere_grid"](n_lat=180, n_lon=360)
    mu0, vis0 = G["compute_mu"](phi, lon, B0=0.0, lambda0=0.0)
    dOm0 = G["pixel_solid_angle"](phi, mu0, vis0, dphi, dlon)
    E_qs = (G["interp_I_vs_mu"](i_qs, mu0, mu_grid=MU) * dOm0.ravel()[None, :]).sum(axis=1)

    dist = G["build_distribution"](
        phi, lon,
        [{"type": "umbra", "shape": "circle", "lat": 0.0, "lon": -100.0, "radius": RADIUS}])

    n = 120
    radii = RADIUS * np.sin(np.pi * np.arange(n) / (n - 1))
    masks_B = []
    for r in radii:
        if r < 0.5:
            masks_B.append(None)
        else:
            masks_B.append(G["make_circular_patch_mask"](
                phi, lon, lat0_deg=0.0, lon0_deg=0.0, radius_deg=r))

    rows = []
    fitsA, fitsB = {}, {}
    for logg in FEATURE_LOGGS:
        print(f"\nfeature logg={logg}")
        w2, i_feat = clv_si(T_FEATURE, logg)
        assert np.allclose(w2, wav)

        _, E_A, _ = G["simulate_rotation_multi"](
            wav, i_qs, {"umbra": i_feat}, dist, phi, lon, dphi, dlon, E_qs,
            n_steps=120, rotation_period_days=27.0, B0=0.0)
        fitA = G["fit_ssi_vs_tsi"](wav, E_A, E_qs)
        fitsA[logg] = fitA

        E_B = np.empty((n, wav.size))
        for i, mask in enumerate(masks_B):
            if mask is None:
                E_B[i] = E_qs
            else:
                E_B[i] = E_qs + G["feature_delta_irradiance"](
                    i_qs, i_feat, mu0, vis0, dOm0, mask)
        fitB = G["fit_ssi_vs_tsi"](wav, E_B, E_qs)
        fitsB[logg] = fitB

        row = {"feature_logg": logg}
        for tag, fit in (("A", fitA), ("B", fitB)):
            m = band_mask(wav, fit)
            row[f"case{tag}_rms_planck5780"] = distance_from_planck(wav, fit)
            row[f"case{tag}_best_planck_T"] = best_fit_planck_T(wav, fit)
            row[f"case{tag}_mean_corr"] = float(np.nanmean(fit["corr"][m]))
        rows.append(row)
        print(json.dumps(row, indent=2))

    (OUT / "_100K_logg_scan_results.json").write_text(json.dumps(rows, indent=2))
    with open(OUT / "_100K_logg_scan_results.csv", "w") as f:
        keys = list(rows[0].keys())
        f.write(",".join(keys) + "\n")
        for r in rows:
            f.write(",".join(str(r[k]) for k in keys) + "\n")

    aP = G["planck_temperature_slope"](wav, T_PLANCK)
    for tag, fits in (("A", fitsA), ("B", fitsB)):
        fig, ax = plt.subplots(figsize=(10.5, 6.4))
        for logg in FEATURE_LOGGS:
            fit = fits[logg]
            m = band_mask(wav, fit)
            rms = next(r[f"case{tag}_rms_planck5780"] for r in rows if r["feature_logg"] == logg)
            ax.plot(wav[m], fit["slope"][m] / aP[m], lw=1.5,
                    label=fr"5700 K, $\log g$={logg:g}  (RMS={rms:.3f})")
        ax.axhline(1.0, color="black", ls="--", lw=1.5)
        ax.set(xlim=(LO, HI), xlabel="Wavelength (nm)",
               ylabel=fr"$a(\lambda)/a_{{\rm Planck}}(5780\,{{\rm K}})$",
               title=f"100 K MPS-ATLAS experiment, Case {tag}: lower gravity in cooler model")
        ax.grid(alpha=.25)
        ax.legend(fontsize=8, ncol=2)
        fig.tight_layout()
        fig.savefig(OUT / f"_100K_logg_scan_case{tag}_ratio.png", dpi=170)
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    gvals = np.array(FEATURE_LOGGS)
    for tag, marker in (("A", "o"), ("B", "s")):
        y = [r[f"case{tag}_rms_planck5780"] for r in rows]
        ax.plot(gvals, y, marker=marker, lw=1.7, label=f"Case {tag}")
    ax.invert_xaxis()
    ax.set(xlabel=r"Cooler component $\log g$", ylabel="RMS distance from Planck(5780)",
           title="Does lower gravity make the 100 K response more Planckian?")
    ax.grid(alpha=.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT / "_100K_logg_scan_rms.png", dpi=170)
    plt.close(fig)

    # compact selected-wavelength ratios for interpretation
    sample_nm = [500, 700, 1000, 1300, 1500, 1600]
    sample = {}
    for tag, fits in (("A", fitsA), ("B", fitsB)):
        sample[tag] = {}
        for logg, fit in fits.items():
            vals = {}
            for lam in sample_nm:
                j = int(np.argmin(abs(wav - lam)))
                vals[str(lam)] = float(fit["slope"][j] / aP[j]) if fit["valid"][j] else None
            sample[tag][str(logg)] = vals
    (OUT / "_100K_logg_scan_samples.json").write_text(json.dumps(sample, indent=2))

    print("\nFinal table")
    for r in rows:
        print(r)


if __name__ == "__main__":
    main()
