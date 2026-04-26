# enigmatic_Planck — project context

Modeling **SSI/TSI variability** caused by solar rotation transporting distributions of magnetic features (umbrae, penumbrae, faculae) across the visible disk. The main scientific deliverable is the regression slope **a(λ) = (ΔSSI/E_qs(λ)) / (ΔTSI/TSI_qs)** for each feature-distribution scenario.

## Authoritative source
`model.ipynb` is the canonical source of all physics/plot code. `run_all.py` and `build_pptx.py` re-use functions from the notebook (no duplicate code) — any new function added to the notebook needs to be whitelisted in `run_all.py`'s `_is_def_cell()`.

## Files
- `qs.txt`, `umbra.txt`, `penumbra.txt`, `faculae.txt` — input I_ν(λ, μ) tables (CGS per Hz, 1221 wavelengths × 11 μ values, 9–160 000 nm). Below ~10 nm the values are the sentinel `1e-99` — handled by `eps_rel` in `fit_ssi_vs_tsi`.
- `scenarios/` — output folders `01_..09_*` from the 9-scenario batch + analysis figures (`_planck_*.png`, `_overlay_*.png`, `_comparison_*.png`).
- `scenarios/_per_feature_cache.npz` — cached per-feature ΔE(t, λ) arrays for the **(lat=0°, lon=−100°, R=10°)** reference geometry. Read by cells 30–32 to evaluate any (w_u, w_p, w_f) mixture without re-running rotation.
- `requirements.txt` — `numpy scipy matplotlib jupyterlab Pillow python-pptx`.

## Conventions (do not change silently)
- **Rotation direction**: features move left → right on the disk (ω = −2π/P). Implemented in `simulate_rotation_multi` and `animate_scenario`.
- **Feature priority** when masks overlap: umbra > penumbra > faculae > QS. Enforced in `build_distribution`; downstream code assumes the three masks are disjoint.
- **Spot at east limb at t=0**: place at lon = −90° − R (so lon = −100° for R=10°). Used in scenarios 01–03.
- **Junk-λ masking**: `fit_ssi_vs_tsi` flags λ where E_qs(λ) < eps_rel · max(E_qs) (default `eps_rel=1e-30`). Always filter by `result["valid"]` before plotting slope/corr.
- **Units**: intensities converted to SI per nm once (`inu_cgs_to_ilambda_si_per_nm`). Relative variations plotted in **ppm** (axis label "Relative variation (ppm)"). a(λ) is dimensionless, λ in nm.
- **Smoothing**: when overplotting per-scenario a(λ), Gaussian-smooth in wavelength space (FWHM = 25 nm) — helper `_smooth_wav` is defined in cell 28.
- **Grid defaults**: `n_lat=180, n_lon=360` (1°×1°), `n_steps=180`, `rotation_period_days=27.0`, `B0=0.0`.

## Where we are (current session: 2026-04-26)
Cells 26–32 layer a Planck-δT baseline analysis on top of the existing pipeline:

- **Cell 26** — `planck_temperature_slope(wav_nm, T_eff)`: analytic a(λ) for a Planck-spectrum Sun under uniform δT.
  Derivation: ΔTSI/TSI = 4ΔT/T (Stefan–Boltzmann) and ΔSSI/SSI = (∂ln B/∂T)·ΔT (Planck) ⇒ **a(λ) = (x/4)/(1−exp(−x))**, x = hc/(λ k T). Pearson r ≡ 1 by construction.
- **Cells 27–28** — comparison + overlay plots vs scenarios 01–07.
- **Cell 29** — cached per-feature ΔE arrays at (0°, −100°, R=10°).
- **Cell 30** — mixture optimization at fixed T = 5780 K (scipy L-BFGS-B, w_k ≥ 0). Best (u:p:f) = **0.10 : 0.56 : 0.34**, RMS = **0.078** in 300–1800 nm.
- **Cell 31** — joint fit (w, T_eff). Best T = **5374 K** (~400 K cooler than nominal), weights u:p:f = **0.07 : 0.62 : 0.31**, RMS = **0.043** (45 % improvement).
- **Cell 32** — `planck_T_from_slope` pointwise Planck inversion. T(λ) of the best mixture rises from ~5500 K (visible) to ~10–12 kK (IR); concentric AR goes the *opposite* direction (down toward ~2 kK in the IR). Faculae alone are intrinsically non-Planckian (a < 1/4 in the IR ⇒ inversion fails).

## Open questions for next session
1. **Reference geometry alternatives** — pick one and re-run the cache+optimization pipeline (cells 29–32):
   - **(b-i) Concentric AR co-located** at (0°, −100°): umbra 2°, penumbra 5°, faculae 10°. Each component runs as its own single-feature scenario.
   - **(b-ii) Spot inside brightening halo**: pure umbra (no penumbra) inside a faculae halo — umbra 2° + faculae 10° at (0°, −100°).
2. **IR mismatch diagnostic** (proposed two-panel figure):
   - per-feature ΔSSI(λ, t)/ΔTSI(t) at peak rotation phase
   - μ-dependent contrast (I_feat − I_qs)/I_qs at μ=1 vs μ=0.3
   to test whether faculae's μ-dependent contrast crossover drives the IR residual.

## How to resume on a fresh machine
```bash
git clone https://github.com/Sasha932-astro/enigmatic_Planck.git
cd enigmatic_Planck
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
jupyter lab model.ipynb
```
Then run cells **0, 4, 6, 7, 8, 10, 12, 14, 26** to load definitions; jump to cell 30 (or 31) for mixture experiments. Skip cell 29 — `_per_feature_cache.npz` is committed.

## Working style
- **Step-by-step**: after each routine, run a short test, show numerical/physical sanity check, then pause for confirmation ("All agreed", "proceed") before the next step. Don't batch multiple design changes into one reply.
- **Approve before expensive compute**: enumerate scenarios + parameters and wait for explicit approval ("yes", "go ahead") before multi-minute runs or disk writes.
- **Concise responses**, direct answers; assume the user has read the code. Keep verdicts short (1–2 sentences) with file:line references.
- **Always label ppm** on relative-variation axes; nm on wavelength; "dimensionless" for a(λ).
- **Distributions in natural language**: user describes ("a 5° spot at 15°N with a faculae halo"); translate to the `[{"type": ..., "lat": ..., ...}, ...]` spec list and show before running.
