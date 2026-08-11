#!/usr/bin/env python3
"""MURaM spot-temperature surrogates and disk-centre mixture scan.

Uses only real spectra from the MURaM/MPS-ATLAS Magnetic LD Library
(Edmond doi:10.17617/3.FBTIYY).

Scientific setup
----------------
Quiet star: G2_MH_00, hydro.
Single-component spot surrogates: K0, K4, M0, M2 (plus M4 as an exploratory
cool endpoint), all hydro spectra from the same MURaM library.

For each single-component spot we run the same Case A (rotating transit) and
Case B (grows/shrinks at disk centre) geometry used by experiment_100K.py and
experiment_full_comparison.py.

Mixtures are investigated only for Case B, as requested. A mixture is a fixed
sub-area composition inside the active region:

    I_mix(lambda, mu) = sum_j w_j I_j(lambda, mu),
    w_j >= 0, sum_j w_j = 1.

All subareas grow/shrink together, so the weights are directly analogous to a
fixed umbra/penumbra area ratio. Because the irradiance perturbation is linear
in intensity, the mixture time series is computed exactly from precomputed
component perturbations; no spectral interpolation between spectral classes,
curve digitization, or surrogate fitting is used.

Outputs in scenarios/:
  _muram_classes_caseA_a.png
  _muram_classes_caseA_ratio.png
  _muram_classes_caseB_a.png
  _muram_classes_caseB_ratio.png
  _muram_mixture_pair_scan.png
  _muram_mixture_best_a.png
  _muram_mixture_best_ratio.png
  _muram_mixture_results.json
  _muram_mixture_results.txt
  _muram_mixture_arrays.npz
"""

from __future__ import annotations

import itertools
import json
import pathlib

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import minimize, minimize_scalar

import experiment_100K as base
import experiment_full_comparison as fc

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "scenarios"
OUT.mkdir(exist_ok=True)

LO, HI = 400.0, 1600.0
T_PLANCK = 5780.0
RADIUS = 10.0
N_STEPS = 120

# Edmond datafile IDs from doi:10.17617/3.FBTIYY, Version 1.0.
fc.MURAM_FILES.update({
    "K4_MH_00": (344881, "K4_MH_00.h5"),
    "M0_MH_00": (344891, "M0_MH_00.h5"),
    "M2_MH_00": (344886, "M2_MH_00.h5"),
    "M4_MH_00": (344890, "M4_MH_00.h5"),
})

CLASS_TO_STAR = {
    "K0": "K0_MH_00",
    "K4": "K4_MH_00",
    "M0": "M0_MH_00",
    "M2": "M2_MH_00",
    "M4": "M4_MH_00",
}
PRIMARY_CLASSES = ["K0", "K4", "M0", "M2"]
MIX_CLASSES = ["K0", "K4", "M0", "M2", "M4"]


def planck(wav):
    return base.G["planck_temperature_slope"](wav, T_PLANCK)


def load_hydro(cls: str):
    return fc._load_muram_state(CLASS_TO_STAR[cls], "hydro")


def plot_class_results(results, case: str, ratio: bool, path: pathlib.Path):
    fig, ax = plt.subplots(figsize=(10, 6))
    for cls in PRIMARY_CLASSES:
        r = results[cls]
        wav = r["wav"]
        m = r[case]["valid"] & (wav >= LO) & (wav <= HI)
        y = r[case]["slope"][m]
        if ratio:
            y = y / planck(wav[m])
        ax.plot(wav[m], y, lw=1.5,
                label=f"G2 -> {cls} hydro (RMS {r[f'{case}_rms']:.3f})")
    if ratio:
        ax.axhline(1.0, ls="--", lw=1.8, label="Planck 5780 K")
        ax.set_ylabel("a(lambda) / Planck(5780 K)")
        ax.set_title(f"MURaM spot-temperature surrogates, Case {case}: ratio to Planck")
    else:
        w = np.linspace(LO, HI, 1201)
        ax.plot(w, planck(w), "--", lw=1.8, label="Planck 5780 K")
        ax.set_ylabel("a(lambda)")
        ax.set_title(f"MURaM spot-temperature surrogates, Case {case}")
    ax.set_xlabel("Wavelength (nm)")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def prepare_caseB_component_basis(wav, i_qs, component_spectra):
    """Precompute exact disk-centre perturbation time series for each component.

    Returns
    -------
    E_qs : (L,)
    delta_basis : (C,T,L), where each entry is the exact irradiance
                  perturbation for a pure component at that time/radius.
    """
    phi, lon, dphi, dlon = base.G["build_sphere_grid"](n_lat=180, n_lon=360)
    mu0, vis0 = base.G["compute_mu"](phi, lon, B0=0.0, lambda0=0.0)
    dOm0 = base.G["pixel_solid_angle"](phi, mu0, vis0, dphi, dlon)
    E_qs = (base.G["interp_I_vs_mu"](i_qs, mu0, mu_grid=base.MU)
            * dOm0.ravel()[None, :]).sum(axis=1)

    # Only pixels inside the largest disk-centre patch can ever contribute.
    max_mask = base.G["make_circular_patch_mask"](
        phi, lon, lat0_deg=0.0, lon0_deg=0.0, radius_deg=RADIUS)
    active = max_mask & vis0
    mu_active = mu0[active]
    dOm_active = dOm0[active]

    Iq = base.G["interp_I_vs_mu"](i_qs, mu_active, mu_grid=base.MU)

    radii = RADIUS * np.sin(np.pi * np.arange(N_STEPS) / (N_STEPS - 1))
    submasks = []
    for r in radii:
        if r < 0.5:
            submasks.append(np.zeros(mu_active.size, dtype=bool))
        else:
            full = base.G["make_circular_patch_mask"](
                phi, lon, lat0_deg=0.0, lon0_deg=0.0, radius_deg=float(r))
            submasks.append(full[active])

    C = len(component_spectra)
    delta_basis = np.zeros((C, N_STEPS, wav.size), dtype=float)
    for j, Ifull in enumerate(component_spectra):
        If = base.G["interp_I_vs_mu"](Ifull, mu_active, mu_grid=base.MU)
        weighted = (If - Iq) * dOm_active[None, :]
        for t, sm in enumerate(submasks):
            if np.any(sm):
                delta_basis[j, t] = weighted[:, sm].sum(axis=1)
    return E_qs, delta_basis


def fit_caseB_weights(wav, E_qs, delta_basis, weights):
    """Exact Case-B regression for a convex spectral-class mixture."""
    w = np.asarray(weights, dtype=float)
    D = np.tensordot(w, delta_basis, axes=(0, 0))  # (T,L)

    TSI_qs = np.trapezoid(E_qs, wav)
    x = np.trapezoid(D, wav, axis=1) / TSI_qs
    valid = E_qs > (1e-30 * np.nanmax(E_qs))

    Y = np.full_like(D, np.nan)
    Y[:, valid] = D[:, valid] / E_qs[None, valid]

    xc = x - x.mean()
    x_var = np.mean(xc ** 2)
    slope = np.full(wav.size, np.nan)
    corr = np.full(wav.size, np.nan)
    if x_var > 0:
        Yv = Y[:, valid]
        yc = Yv - Yv.mean(axis=0)[None, :]
        cov = np.mean(xc[:, None] * yc, axis=0)
        slope[valid] = cov / x_var
        with np.errstate(divide="ignore", invalid="ignore"):
            corr_v = cov / (np.sqrt(x_var) * Yv.std(axis=0))
        corr[valid] = np.where(np.isfinite(corr_v), corr_v, np.nan)

    m = valid & (wav >= LO) & (wav <= HI)
    rms = float(np.sqrt(np.mean((slope[m] - planck(wav[m])) ** 2)))
    bestT = base.best_fit_planck_T(wav, slope, valid)
    mean_corr = float(np.nanmean(corr[m]))
    return {
        "weights": w,
        "slope": slope,
        "corr": corr,
        "valid": valid,
        "rms": rms,
        "bestT": bestT,
        "mean_corr": mean_corr,
    }


def pair_scan(wav, E_qs, delta_basis, classes):
    """Dense exact scan plus local refinement for every class pair."""
    out = []
    grid = np.linspace(0.0, 1.0, 501)
    curves = {}

    for ia, ib in itertools.combinations(range(len(classes)), 2):
        ca, cb = classes[ia], classes[ib]
        rms_grid = []
        for fa in grid:
            w = np.zeros(len(classes))
            w[ia] = fa
            w[ib] = 1.0 - fa
            rms_grid.append(fit_caseB_weights(wav, E_qs, delta_basis, w)["rms"])
        rms_grid = np.asarray(rms_grid)
        curves[f"{ca}+{cb}"] = rms_grid
        k = int(np.argmin(rms_grid))
        lo = grid[max(0, k - 4)]
        hi = grid[min(grid.size - 1, k + 4)]

        def obj(fa):
            w = np.zeros(len(classes))
            w[ia] = fa
            w[ib] = 1.0 - fa
            return fit_caseB_weights(wav, E_qs, delta_basis, w)["rms"]

        opt = minimize_scalar(obj, bounds=(lo, hi), method="bounded",
                              options={"xatol": 1e-7})
        fa = float(opt.x)
        w = np.zeros(len(classes))
        w[ia] = fa
        w[ib] = 1.0 - fa
        fit = fit_caseB_weights(wav, E_qs, delta_basis, w)
        out.append({
            "classes": [ca, cb],
            "weights": [float(fa), float(1.0 - fa)],
            "rms": fit["rms"],
            "bestT": fit["bestT"],
            "mean_corr": fit["mean_corr"],
            "fit": fit,
        })
    out.sort(key=lambda x: x["rms"])
    return grid, curves, out


def optimize_subset(wav, E_qs, delta_basis, classes, subset, rng):
    """Robust constrained optimization for one 3+ component subset."""
    inds = np.array([classes.index(c) for c in subset], dtype=int)
    n = len(inds)

    def full_weights(local):
        w = np.zeros(len(classes))
        w[inds] = local
        return w

    def obj(local):
        return fit_caseB_weights(wav, E_qs, delta_basis, full_weights(local))["rms"]

    starts = [np.full(n, 1.0 / n)]
    # Real convex mixtures, sampled broadly.  These are only starting points;
    # every reported result is subsequently optimized exactly.
    for _ in range(80):
        starts.append(rng.dirichlet(np.ones(n)))
    starts.sort(key=obj)

    best = None
    cons = {"type": "eq", "fun": lambda z: np.sum(z) - 1.0}
    bounds = [(0.0, 1.0)] * n
    for x0 in starts[:8]:
        res = minimize(obj, x0, method="SLSQP", bounds=bounds,
                       constraints=cons,
                       options={"ftol": 1e-12, "maxiter": 300})
        z = np.clip(res.x, 0.0, 1.0)
        z /= z.sum()
        fit = fit_caseB_weights(wav, E_qs, delta_basis, full_weights(z))
        if best is None or fit["rms"] < best["rms"]:
            best = {
                "classes": list(subset),
                "weights": [float(x) for x in z],
                "rms": fit["rms"],
                "bestT": fit["bestT"],
                "mean_corr": fit["mean_corr"],
                "fit": fit,
                "success": bool(res.success),
            }
    return best


def plot_pair_scan(grid, curves, path):
    fig, ax = plt.subplots(figsize=(11, 6.5))
    for name, rms in curves.items():
        a, b = name.split("+")
        ax.plot(grid, rms, lw=1.2, label=f"{a}+{b}")
    ax.set_xlabel("Area fraction of first class in pair")
    ax.set_ylabel("RMS distance from Planck(5780 K)")
    ax.set_title("Disk-centre MURaM two-component spot mixtures")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_best_mixtures(wav, pure_fits, best_pair, best_multi, path, ratio=False):
    fig, ax = plt.subplots(figsize=(10, 6))
    m = (wav >= LO) & (wav <= HI)

    # Best pure class among mixture candidates.
    best_pure_name, best_pure = min(pure_fits.items(), key=lambda kv: kv[1]["rms"])
    entries = [
        (f"best pure: {best_pure_name} (RMS {best_pure['rms']:.3f})", best_pure),
        (f"best pair: {'+'.join(best_pair['classes'])} (RMS {best_pair['rms']:.3f})", best_pair["fit"]),
        (f"best multi: {'+'.join(best_multi['classes'])} (RMS {best_multi['rms']:.3f})", best_multi["fit"]),
    ]
    for label, fit in entries:
        y = fit["slope"][m]
        if ratio:
            y = y / planck(wav[m])
        ax.plot(wav[m], y, lw=1.6, label=label)

    if ratio:
        ax.axhline(1.0, ls="--", lw=1.8, label="Planck 5780 K")
        ax.set_ylabel("a(lambda) / Planck(5780 K)")
        ax.set_title("Best physical MURaM spot mixtures, disk centre: ratio to Planck")
    else:
        ax.plot(wav[m], planck(wav[m]), "--", lw=1.8, label="Planck 5780 K")
        ax.set_ylabel("a(lambda)")
        ax.set_title("Best physical MURaM spot mixtures, disk centre")
    ax.set_xlabel("Wavelength (nm)")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main():
    print("MURaM class + mixture investigation using real Edmond HDF5 spectra\n")

    # Quiet G2 hydro and all cooler hydro states.
    wav_g, i_g, teff_g, mh_g, logg_g = fc._load_muram_state("G2_MH_00", "hydro")
    states = {}
    meta = {}
    for cls in MIX_CLASSES:
        wav, inten, teff, mh, logg = load_hydro(cls)
        if not np.allclose(wav, wav_g):
            raise ValueError(f"wavelength grid differs for {cls}")
        states[cls] = inten
        meta[cls] = {"teff_K": teff, "MH": mh, "logg": logg}

    print(f"quiet: G2 hydro, Teff={teff_g:.1f} K, [M/H]={mh_g:+.1f}, logg={logg_g:.3f}")
    for cls in MIX_CLASSES:
        print(f"{cls:>2s} hydro: Teff={meta[cls]['teff_K']:.1f} K, logg={meta[cls]['logg']:.3f}")
    print()

    # Full Case A/B runs for requested classes, plus K0 reference.
    class_results = {}
    for cls in PRIMARY_CLASSES:
        print(f"Running full G2 -> {cls} hydro experiment ...")
        class_results[cls] = fc._run_pair(
            f"MURaM G2->{cls} hydro", wav_g, i_g, states[cls])

    print("\nSingle-component spot surrogates:")
    print("class   Teff[K]   Case A RMS   Case B RMS   A bestT   B bestT")
    print("-" * 67)
    for cls in PRIMARY_CLASSES:
        r = class_results[cls]
        print(f"{cls:>4s} {meta[cls]['teff_K']:9.1f} {r['A_rms']:12.4f} {r['B_rms']:12.4f} "
              f"{r['A_bestT']:9.0f} {r['B_bestT']:9.0f}")

    plot_class_results(class_results, "A", False, OUT / "_muram_classes_caseA_a.png")
    plot_class_results(class_results, "A", True, OUT / "_muram_classes_caseA_ratio.png")
    plot_class_results(class_results, "B", False, OUT / "_muram_classes_caseB_a.png")
    plot_class_results(class_results, "B", True, OUT / "_muram_classes_caseB_ratio.png")

    # Exact disk-centre basis for mixtures.
    component_spectra = [states[c] for c in MIX_CLASSES]
    E_qs, delta_basis = prepare_caseB_component_basis(wav_g, i_g, component_spectra)

    # Pure-class disk-centre fits on the same mixture basis, including M4.
    pure_fits = {}
    for j, cls in enumerate(MIX_CLASSES):
        w = np.zeros(len(MIX_CLASSES)); w[j] = 1.0
        pure_fits[cls] = fit_caseB_weights(wav_g, E_qs, delta_basis, w)

    # All physical pair mixtures.
    grid, pair_curves, pairs = pair_scan(wav_g, E_qs, delta_basis, MIX_CLASSES)
    best_pair = pairs[0]

    # All 3-, 4-, and 5-component subsets.  Multiple randomized starts protect
    # against a local constrained minimum.
    rng = np.random.default_rng(20260811)
    multis = []
    for n in (3, 4, 5):
        for subset in itertools.combinations(MIX_CLASSES, n):
            multis.append(optimize_subset(
                wav_g, E_qs, delta_basis, MIX_CLASSES, subset, rng))
    multis.sort(key=lambda x: x["rms"])
    best_multi = multis[0]

    plot_pair_scan(grid, pair_curves, OUT / "_muram_mixture_pair_scan.png")
    plot_best_mixtures(wav_g, pure_fits, best_pair, best_multi,
                       OUT / "_muram_mixture_best_a.png", ratio=False)
    plot_best_mixtures(wav_g, pure_fits, best_pair, best_multi,
                       OUT / "_muram_mixture_best_ratio.png", ratio=True)

    print("\nBest pairwise mixtures (disk centre):")
    print("pair          weights(first,second)      RMS      bestT")
    print("-" * 62)
    for p in pairs:
        print(f"{p['classes'][0]}+{p['classes'][1]:<7s}  "
              f"({p['weights'][0]:.4f}, {p['weights'][1]:.4f})   "
              f"{p['rms']:.5f}   {p['bestT']:.0f}")

    print("\nBest 3+ component solutions:")
    for x in multis[:10]:
        weights = ", ".join(f"{c}:{w:.4f}" for c, w in zip(x['classes'], x['weights']))
        print(f"RMS={x['rms']:.5f}  T={x['bestT']:.0f} K  {weights}")

    best_pure_name, best_pure = min(pure_fits.items(), key=lambda kv: kv[1]["rms"])
    print("\nGlobal summary:")
    print(f"best pure  : {best_pure_name}, RMS={best_pure['rms']:.5f}, bestT={best_pure['bestT']:.0f} K")
    print(f"best pair  : {'+'.join(best_pair['classes'])}, RMS={best_pair['rms']:.5f}, "
          f"weights={best_pair['weights']}, bestT={best_pair['bestT']:.0f} K")
    print(f"best multi : {'+'.join(best_multi['classes'])}, RMS={best_multi['rms']:.5f}, "
          f"weights={best_multi['weights']}, bestT={best_multi['bestT']:.0f} K")

    # Machine-readable results, omitting bulky slope arrays from JSON.
    results_json = {
        "dataset": "doi:10.17617/3.FBTIYY",
        "quiet": {"class": "G2", "state": "hydro", "teff_K": teff_g,
                  "MH": mh_g, "logg": logg_g},
        "components": meta,
        "analysis_band_nm": [LO, HI],
        "planck_reference_K": T_PLANCK,
        "mixture_definition": "convex area-weighted sum of MURaM hydro intensities; fixed weights; Case B disk centre only",
        "single_component": {
            cls: {
                "caseA_rms": class_results[cls]["A_rms"],
                "caseB_rms": class_results[cls]["B_rms"],
                "caseA_bestT_K": class_results[cls]["A_bestT"],
                "caseB_bestT_K": class_results[cls]["B_bestT"],
                "caseA_mean_corr": class_results[cls]["A_corr"],
                "caseB_mean_corr": class_results[cls]["B_corr"],
            } for cls in PRIMARY_CLASSES
        },
        "pure_disk_center": {
            cls: {"rms": f["rms"], "bestT_K": f["bestT"], "mean_corr": f["mean_corr"]}
            for cls, f in pure_fits.items()
        },
        "pairs": [
            {k: p[k] for k in ("classes", "weights", "rms", "bestT", "mean_corr")}
            for p in pairs
        ],
        "multicomponent": [
            {k: x[k] for k in ("classes", "weights", "rms", "bestT", "mean_corr", "success")}
            for x in multis
        ],
        "best_pure": {"class": best_pure_name, "rms": best_pure["rms"],
                      "bestT_K": best_pure["bestT"]},
        "best_pair": {k: best_pair[k] for k in ("classes", "weights", "rms", "bestT", "mean_corr")},
        "best_multicomponent": {k: best_multi[k] for k in ("classes", "weights", "rms", "bestT", "mean_corr", "success")},
    }
    with open(OUT / "_muram_mixture_results.json", "w") as f:
        json.dump(results_json, f, indent=2)

    # Preserve all curves and the dense pair scans.
    arrays = {
        "wav_nm": wav_g,
        "E_qs": E_qs,
        "planck_5780": planck(wav_g),
        "pair_fraction_grid": grid,
        "best_pair_slope": best_pair["fit"]["slope"],
        "best_multi_slope": best_multi["fit"]["slope"],
    }
    for cls, f in pure_fits.items():
        arrays[f"pure_{cls}_slope"] = f["slope"]
    for name, curve in pair_curves.items():
        arrays[f"pairscan_{name.replace('+', '_')}"] = curve
    for cls in PRIMARY_CLASSES:
        arrays[f"{cls}_caseA_slope"] = class_results[cls]["A"]["slope"]
        arrays[f"{cls}_caseB_slope"] = class_results[cls]["B"]["slope"]
    np.savez_compressed(OUT / "_muram_mixture_arrays.npz", **arrays)

    print(f"\nSaved figures, arrays, and JSON to {OUT}")


if __name__ == "__main__":
    main()
