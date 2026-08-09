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
- **Cells 33–35** (session 2026-06-30) — **umbra spot-distribution sensitivity of a(λ)** over 400–1600 nm vs Planck-δT(5780 K). `umbra_a_of_lambda`, `planck_band_metrics`, `scatter_umbra_spots`, `build_umbra_scenarios`, `run_umbra_distribution_study` (cell 35 guarded by `RUN_UMBRA_DISTRIBUTION_STUDY=False`). **Finding: a(λ) is nearly independent of how umbrae are distributed** — best-fit apparent T_eff = **5246–5285 K** across all realistic |lat|≤35° distributions; only latitude is a weak lever (T_eff 5270→5247 K, RMS 0.110→0.134 over 0→35°); radius/number/size/longitude are negligible (a(λ) is intensive). **Caveat: this apparent ~5260 K is NOT the solar 5780 K**, so spot-distribution robustness does not by itself explain the observed 5780 K Planck match — that requires the facula balance. Figures: `scenarios/_umbra_dist_{master,latitude,invariances,realistic}.png`.

- **Cells 36–38** (session 2026-06-30) — **composite-spot (umbra-core + penumbra-annulus) distribution study** vs Planck-δT(5780 K), 400–1600 nm. `composite_spot_specs`, `a_of_lambda_composite`, `scatter_positions`, `build_composite_scenarios`, `run_composite_spot_study` (cell 38 guarded by `RUN_COMPOSITE_SPOT_STUDY=False`). A spot = `dict(lat,lon,r_p,f_u)`, r_u=r_p·√f_u. **Findings:** apparent T_eff = **5270 K (pure umbra) → ~5515 K (pure penumbra)** — a CEILING ~265 K short of 5780; deviation from 5780 **minimized at the real sunspot ratio F≈0.2** (RMS 0.077); Phase 2 shows a(λ) depends only on the **global** umbral area fraction F (heterogeneous/random/size-correlated/bimodal collapse onto the uniform curve within ~10–20 K), not on per-spot heterogeneity or spatial distribution. **Dark features alone cannot reach 5780 K.** Figures: `scenarios/_spot_dist_{curves,lever,robust,globalF}.png`.
- **Cells 39–41** (session 2026-06-30) — **three-component model: composite spots + faculae** (AR-linked halos & diffuse independent network) vs Planck-δT(5780 K), 400–1600 nm. `full_specs`, `a_of_lambda_full`, `network_faculae`, `run_faculae_study` (cell 41 guarded by `RUN_FACULAE_STUDY=False`). A spot may carry a facular halo via `r_f>r_p`. **Findings:** (1) independent diffuse-network faculae raise apparent T_eff and cross 5780 K only at **extreme coverage G=A_fac/A_spot≈25–29**; as T→5780 the RMS-to-Planck doubles (0.08→0.14) and **no config reaches T∈[5750,5810] with RMS<0.10**; (2) the **IR a(λ) collapses to ratio ≈0.3** beyond ~1.3 µm even when best-fit T≈5800 (faculae contrast crossover); (3) **co-located** AR facular halos CANCEL spots (T↓ to ~4400 K, fit degrades) — only diffuse separated network helps; (4) faculae **destroy distribution-robustness**: apparent-T scatter grows from ~40 K (spots) to ~1000 K (spots+faculae, G~8) across random placements; (5) facula *latitude* at matched G is weak (~50 K). **Verdict: no spot+facula configuration yields a robust, Planck-clean a(λ) at 5780 K** — the observed match likely needs the real facular IR/μ-contrast or cycle-scale averaging. Figures: `scenarios/_fac_{coverage,colocated,robust,IR}.png`.

- **Cells 42–44** (session 2026-07-03) — **MPS-ATLAS RE spectra** (`mps_atlas/`, Edmond doi:10.17617/3.NJ56TR, set1, MH0.0, logg4.4, all 56 Teff; same 1221-λ grid as qs.txt, MU_GRID an exact subset of its 24 μ). `mps_clv`, `mps_penumbra_5450`, `run_single_feature_arrays`, `re_delta_t_curve`, guarded `RUN_RE_PAIR_SCAN`. **Findings:** (1) **Model I vs Model II** (project files vs self-consistent grid at SATIRE temps QS5800/P5450/U4500): apparent T agrees to ~40 K (5270→5312 umbra, 5496→5532 penumbra), a(λ) ratios within ~2% — **library choice is not a lever**; (2) **full (T_u,T_p,F) pair scan** (3500–5700 K, 1518 combos, composition validated ≤10 K vs exact runs): **max apparent T = 5612 K — no RE pair lands on Planck(5780)**; best shape match rms 0.060; single-T curve U-shaped (min 5290 K at ~4200 K); (3) **RE-δT limit** (true δT of the RE atmosphere, disk-integrated 5700/5800/5900): best-fit Planck T = **5560 K**, rms 0.081, Stefan–Boltzmann check passes to 0.06% — **even a genuine RE temperature perturbation is not Planck(5780)-like; the real Sun's a(λ) is MORE Planckian than a self-consistent RE response**. Figures: `scenarios/_modelcmp_{umbra,penumbra}.png`, `scenarios/_pairscan_{maps,RElimit}.png`.

## Open questions for next session
0. **Why is the Sun more Planckian than RE models?** The RE-δT limit fits Planck at 5560 K, not 5780 K — line blanketing redistributes the δT response. Candidates: NLTE/chromospheric contributions, the facular (FAL-P-like) component compensating the line-blanketing deficit, or wavelength-band selection effects in the observations. This is now the sharpest form of the project's central question.
1. **Root-cause the IR breakdown** — examine the facular input spectra's center-to-limb (μ) contrast and the IR crossover directly (proposed μ=1 vs μ=0.3 contrast diagnostic); test cycle-scale/ensemble averaging (many rotations) vs the single-rotation snapshot; check whether the observed 5780 K match is dominated by a network regime this geometry under-samples.
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
