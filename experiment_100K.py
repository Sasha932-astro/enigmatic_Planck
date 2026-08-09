#!/usr/bin/env python3
"""
The 100 K experiment: is the deviation from Planck(5780 K) caused by geometry?

A feature only 100 K cooler than the quiet Sun is placed on the disk, using a
SELF-CONSISTENT MPS-ATLAS pair (quiet Sun = 5800 K, feature = 5700 K).  Two
scenarios are compared:

  Case A  the feature ROTATES across the disk (normal transit; the contrast is
          sampled along a whole mu-path, with foreshortening and limb darkening)
  Case B  the feature GROWS AND SHRINKS AT DISK CENTRE and does not rotate
          (the contrast is sampled at mu ~ 1 only; geometry is switched off)

If the deviation of a(lambda) from the Planck 5780 K curve were caused by
foreshortening / limb darkening, Case B would be much closer to Planck.

RESULT (this is what the script reproduces):
    Case A   distance from Planck(5780) = 0.112
    Case B   distance from Planck(5780) = 0.123   (correlation +1.0000)
Both miss the Planck curve by the same amount, so geometry is NOT the cause -
the mismatch is in the model spectra themselves (line blanketing: molecular
bands amplify the violet, the deep-forming infrared is damped).

Run:      python experiment_100K.py
Outputs:  scenarios/_100K_caseA.png, _100K_caseB.png, _100K_overlay.png
          (MPS-ATLAS spectra are downloaded automatically if missing)
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

T_QS, T_FEATURE = 5800, 5700          # MPS-ATLAS grid temperatures, 100 K apart
T_PLANCK = 5780.0                     # the curve the Sun is observed to follow
LO, HI = 400.0, 1600.0                # analysis band (nm)
RADIUS = 10.0                         # feature radius (deg)

# ---- load the physics from the notebook (single source of truth) -------
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
    # only cells that DEFINE what we need - never cells that merely use it,
    # so loading the notebook stays free of side effects
    if not (any(f"def {n}" in src for n in FUNCS) or any(c in src for c in CONSTS)):
        continue
    src = src.split("# ---------- example")[0]      # drop demo blocks that do file I/O
    exec(compile(src, "<model.ipynb>", "exec"), G)
for fn in ("build_sphere_grid", "fit_ssi_vs_tsi", "planck_temperature_slope"):
    assert fn in G, f"{fn} not found in model.ipynb"
MU = G["MU_GRID"]


def clv_si(teff):
    """MPS-ATLAS intensities on the project mu grid, converted to SI per nm."""
    path = fetch([teff], verbose=False)[0]
    data = np.loadtxt(path, skiprows=2)
    wav, I = data[:, 0], data[:, 1:]
    with open(path) as f:
        f.readline()
        mu24 = np.array(f.readline().split()[2:], dtype=float)
    idx = [int(np.argmin(abs(mu24 - m))) for m in MU]
    assert np.allclose(mu24[idx], MU), "project mu grid is not a subset of the library grid"
    return wav, G["inu_cgs_to_ilambda_si_per_nm"](I[:, idx], wav)


def distance_from_planck(wav, a, valid, T=T_PLANCK):
    """RMS distance of a(lambda) from the Planck curve at T, over the band."""
    m = valid & (wav >= LO) & (wav <= HI)
    return float(np.sqrt(np.mean((a[m] - G["planck_temperature_slope"](wav[m], T)) ** 2)))


def best_fit_planck_T(wav, a, valid):
    m = valid & (wav >= LO) & (wav <= HI)
    r = minimize_scalar(
        lambda T: np.sqrt(np.mean((a[m] - G["planck_temperature_slope"](wav[m], T)) ** 2)),
        bounds=(3500, 9000), method="bounded")
    return float(r.x)


def main():
    print(f"quiet Sun = MPS-ATLAS {T_QS} K, feature = {T_FEATURE} K (100 K cooler)\n")
    wav, i_qs = clv_si(T_QS)
    _, i_feat = clv_si(T_FEATURE)

    phi, lon, dphi, dlon = G["build_sphere_grid"](n_lat=180, n_lon=360)
    mu0, vis0 = G["compute_mu"](phi, lon, B0=0.0, lambda0=0.0)
    dOm0 = G["pixel_solid_angle"](phi, mu0, vis0, dphi, dlon)
    E_qs = (G["interp_I_vs_mu"](i_qs, mu0, mu_grid=MU) * dOm0.ravel()[None, :]).sum(axis=1)

    # ---- Case A: the feature rotates across the disk -------------------
    dist = G["build_distribution"](phi, lon, [{"type": "umbra", "shape": "circle",
                                               "lat": 0.0, "lon": -100.0, "radius": RADIUS}])
    _, E_A, _ = G["simulate_rotation_multi"](wav, i_qs, {"umbra": i_feat}, dist,
                                             phi, lon, dphi, dlon, E_qs,
                                             n_steps=120, rotation_period_days=27.0, B0=0.0)
    fitA = G["fit_ssi_vs_tsi"](wav, E_A, E_qs)

    # ---- Case B: the feature grows/shrinks at disk centre, no rotation --
    n = 120
    radii = RADIUS * np.sin(np.pi * np.arange(n) / (n - 1))
    E_B = np.empty((n, wav.size))
    for i, r in enumerate(radii):
        if r < 0.5:
            E_B[i] = E_qs
            continue
        mask = G["make_circular_patch_mask"](phi, lon, lat0_deg=0.0, lon0_deg=0.0, radius_deg=r)
        E_B[i] = E_qs + G["feature_delta_irradiance"](i_qs, i_feat, mu0, vis0, dOm0, mask)
    fitB = G["fit_ssi_vs_tsi"](wav, E_B, E_qs)

    # ---- report --------------------------------------------------------
    results = []
    for tag, label, fit in (("A", "rotating across the disk", fitA),
                            ("B", "emerging at disk centre (no rotation)", fitB)):
        d = distance_from_planck(wav, fit["slope"], fit["valid"])
        m = fit["valid"] & (wav >= LO) & (wav <= HI)
        results.append((tag, label, fit, d))
        print(f"Case {tag}  {label}")
        print(f"          distance from Planck({T_PLANCK:.0f} K) = {d:.3f}")
        print(f"          best-fit Planck T = {best_fit_planck_T(wav, fit['slope'], fit['valid']):.0f} K"
              f"   mean SSI-TSI correlation = {fit['corr'][m].mean():+.4f}\n")
    print("Both cases miss the Planck(5780) curve by the same amount:")
    print("geometry is not the cause - the mismatch lives in the spectra.")

    # ---- figures -------------------------------------------------------
    aP = G["planck_temperature_slope"](wav, T_PLANCK)
    for tag, label, fit, d in results:
        m = fit["valid"] & (wav >= LO) & (wav <= HI)
        fig, ax = plt.subplots(1, 2, figsize=(13, 5))
        ax[0].plot(wav[m], fit["slope"][m], lw=1.6, label=f"100 K cooler feature, {label}")
        ax[0].plot(wav[m], aP[m], "r--", lw=1.8, label=f"Planck {T_PLANCK:.0f} K")
        ax[0].set(xlabel="Wavelength (nm)", ylabel="a(λ) (dimensionless)", title=f"Case {tag}: a(λ)")
        ax[0].legend(fontsize=9); ax[0].grid(alpha=.3)
        ax[1].plot(wav[m], fit["slope"][m] / aP[m], lw=1.6)
        ax[1].axhline(1.0, color="red", ls="--", lw=1.8)
        ax[1].set(xlabel="Wavelength (nm)", ylabel="a(λ) / Planck(5780 K)", ylim=(0.6, 1.4),
                  title=f"ratio to Planck(5780)   distance = {d:.3f}")
        ax[1].grid(alpha=.3)
        fig.tight_layout(); fig.savefig(OUT / f"_100K_case{tag}.png", dpi=140); plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 6))
    for (tag, label, fit, d), c in zip(results, ("C0", "C3")):
        m = fit["valid"] & (wav >= LO) & (wav <= HI)
        ax.plot(wav[m], fit["slope"][m] / aP[m], color=c, lw=1.8,
                label=f"Case {tag}: {label}  (distance {d:.3f})")
    ax.axhline(1.0, color="red", ls="--", lw=2, label=f"Planck {T_PLANCK:.0f} K")
    ax.set(xlabel="Wavelength (nm)", ylabel="a(λ) / Planck(5780 K)", ylim=(0.6, 1.4),
           title="Does geometry cause the deviation from Planck(5780 K)?  No.")
    ax.legend(fontsize=9); ax.grid(alpha=.3)
    fig.tight_layout(); fig.savefig(OUT / "_100K_overlay.png", dpi=140); plt.close(fig)
    print(f"\nfigures written to {OUT}/_100K_*.png")


if __name__ == "__main__":
    main()
