# enigmatic_Planck — project primer (read this first)

**Purpose of this file:** a self-contained handoff so an assistant (or a new collaborator)
can start working on this code immediately, with no prior conversation. Read it top to
bottom once; it tells you what the project does, how the code is laid out, the conventions
you must not break, what has already been found, and what is open.

---

## 1. What this project does

We model **spectral and total solar irradiance (SSI / TSI) variability** caused by **solar
rotation** carrying distributions of magnetic features (umbrae, penumbrae, faculae) across
the visible solar disk. As features rotate on/off the disk, the disk-integrated irradiance
varies with time. We measure how that variation is shared between total irradiance (TSI)
and per-wavelength irradiance (SSI).

**Main scientific deliverable** — for each feature-distribution scenario, the
per-wavelength regression slope

> **a(λ) = (ΔSSI(λ)/E_qs(λ)) / (ΔTSI/TSI_qs)**   — dimensionless

i.e. how strongly the relative SSI at wavelength λ responds to a relative change in TSI,
obtained by linear regression of the SSI time series on the TSI time series, wavelength by
wavelength. We also report the Pearson correlation r(λ) of that regression.

A second analysis layer asks: **what would a(λ) look like if the Sun were a Planck blackbody
with a small uniform temperature change δT**, and which *mixture* of umbra/penumbra/faculae
best reproduces that Planck-like a(λ)?

---

## 2. Input data

Four ASCII tables of specific intensity I(λ, μ), one per surface component:

| file | component |
|------|-----------|
| `qs.txt`       | quiet Sun (reference) |
| `umbra.txt`    | sunspot umbra (dark) |
| `penumbra.txt` | sunspot penumbra |
| `faculae.txt`  | facula / bright network |

- Each table is **1221 wavelengths × 11 μ values**, wavelengths 9–160 000 nm.
- Intensities are **CGS per Hz** (I_ν). Below ~10 nm the entries are the sentinel `1e-99`
  (no real data) — handled downstream by junk-λ masking, never delete them.
- μ = cos(heliocentric angle); μ=1 is disk center, μ→0 is the limb.

---

## 3. File inventory

- **`model.ipynb`** — the **canonical source of all physics and plotting code** (33 code
  cells). Everything else re-uses functions from it; do not duplicate logic elsewhere.
- `run_all.py` — batch driver: loads the notebook's function definitions by `exec`-ing the
  relevant cells, then runs all 9 scenarios and writes their figures. If you add a new
  function to the notebook that the batch needs, add its name to `_is_def_cell()`'s
  whitelist. *(Known snag: `run_all.py` currently fails to import because a loaded cell
  references `planck_temperature_slope`, which is not in the whitelist — the notebook itself
  runs fine; only the standalone batch script trips. Fix the whitelist before relying on it.)*
- `build_pptx.py` — assembles `enigmatic_Planck_scenarios.pptx` from the scenario figures.
- `scenarios/` — output folders `01_…09_*`, each with `regression_big.png` (main result),
  `diagnostic.png`, `animation.gif`, `regression_coefficients.txt`. Plus comparison figures
  `_comparison_all.png`, `_comparison_350_1300.png`.
- `scenarios/_per_feature_cache.npz` — cached per-feature ΔE(t, λ) arrays for the reference
  geometry **(lat=0°, lon=−100°, R=10°)**. Lets cells 30–32 evaluate any (w_u, w_p, w_f)
  mixture without re-running rotation. Committed, so cell 29 can be skipped.
- `requirements.txt` — `numpy scipy matplotlib jupyterlab Pillow python-pptx`.
- `CLAUDE.md` — the working-notes / conventions file (this primer is an expanded version).

---

## 4. Notebook structure (cell map)

**Core pipeline (load these to have all functions defined): cells 0, 4, 6, 7, 8, 10, 12, 14, 26.**

| cell | role |
|------|------|
| 0  | `read_mu_intensity_txt(filename)` → (wav, I[λ,μ]); parses the intensity tables |
| 4  | grid + geometry + unit conversion: `build_sphere_grid`, `compute_mu`, `pixel_solid_angle`, `inu_cgs_to_ilambda_si_per_nm` (I_ν CGS → I_λ SI per nm), `interp_I_vs_mu`, `disk_irradiance_at_1au_from_file`. Also physical constants (`C_CGS`, `R_SUN_M`, `AU_M`, `MU_GRID`). |
| 6  | `make_patch_mask`, `feature_delta_irradiance`, `simulate_rotation` (single-feature) |
| 7  | `make_circular_patch_mask(phi, lon, lat0_deg, lon0_deg, radius_deg)` |
| 8  | `build_distribution(phi, lon, specs)` → dict of disjoint boolean masks per feature type |
| 10 | `simulate_rotation_multi(...)` → time series E_t(λ) as features rotate across the disk |
| 12 | `fit_ssi_vs_tsi(wav_nm, E_t, E_qs, eps_rel=1e-30)` → slope a(λ), intercept, corr r(λ), valid mask |
| 14 | `run_scenario(specs, ...)` (end-to-end one-liner) + `plot_diagnostics(result, ...)` |
| 16 | `animate_scenario(result, wavelengths_nm=(400,800,1200), ...)` → rotating-disk GIF |
| 26 | `planck_temperature_slope(wav_nm, T_eff=5772.0)` → analytic Planck-δT baseline a(λ) |
| 27–28 | Planck baseline comparison + overlay plots; `_smooth_wav` (Gaussian FWHM=25 nm) |
| 29 | builds `_per_feature_cache.npz` (slow; skip — file is committed) |
| 30 | mixture optimization at fixed T=5780 K (scipy L-BFGS-B, w_k ≥ 0) |
| 31 | joint fit of (w_u, w_p, w_f, T_eff) |
| 32 | `planck_T_from_slope` — pointwise Planck inversion T(λ) from a(λ) |

### Key signatures
```python
run_scenario(specs, filenames=None, n_lat=180, n_lon=360, n_steps=180,
             rotation_period_days=27.0, B0=0.0, t_days=None, eps_rel=1e-30)
#   -> result dict with: wav_nm, t_days, E_qs, E_t, TSI_rel, SSI_rel,
#      slope (a(λ)), intercept, corr (r(λ)), valid, distribution

fit_ssi_vs_tsi(wav_nm, E_t, E_qs, eps_rel=1e-30)
simulate_rotation_multi(wav_nm, i_qs_lam, intensities, distribution,
                        phi, lon, dphi, dlon, E_qs,
                        t_days=None, n_steps=180, rotation_period_days=27.0, B0=0.0)
planck_temperature_slope(wav_nm, T_eff=5772.0)   # a(λ) = (x/4)/(1−exp(−x)), x = hc/(λkT)
```

### Distribution spec format
A scenario is a list of feature dicts, e.g.
```python
specs = [
    {"type": "faculae",  "shape": "circle", "lat": 0, "lon": -100, "radius": 10},
    {"type": "penumbra", "shape": "circle", "lat": 0, "lon": -100, "radius": 5},
    {"type": "umbra",    "shape": "circle", "lat": 0, "lon": -100, "radius": 2},
]
result = run_scenario(specs)
```
`type` ∈ {umbra, penumbra, faculae}; `lat`/`lon` in degrees; `radius` in degrees.

---

## 5. Conventions — do NOT change silently

- **Rotation direction:** features move **left → right** across the disk (ω = −2π/P).
- **Feature priority when masks overlap:** umbra > penumbra > faculae > quiet Sun. Enforced
  in `build_distribution`; downstream code assumes the three feature masks are disjoint.
- **Spot at east limb at t=0:** place at lon = −90° − R (so lon = −100° for R=10°).
- **Junk-λ masking:** `fit_ssi_vs_tsi` flags λ where E_qs(λ) < eps_rel·max(E_qs). Always
  filter by `result["valid"]` before plotting slope/corr.
- **Units:** intensities converted to SI per nm once via `inu_cgs_to_ilambda_si_per_nm`.
  Relative variations are plotted in **ppm** (axis label "Relative variation (ppm)").
  a(λ) is **dimensionless**; λ in **nm**.
- **Smoothing:** when overplotting per-scenario a(λ), Gaussian-smooth in wavelength
  (FWHM = 25 nm) via `_smooth_wav` (cell 28).
- **Grid defaults:** `n_lat=180, n_lon=360` (1°×1°), `n_steps=180`,
  `rotation_period_days=27.0`, `B0=0.0`.

---

## 6. The 9 standard scenarios (in `run_all.py`)

1. `01_single_umbra` — one 10° umbra at the equator, behind east limb at t=0
2. `02_single_penumbra` — same, penumbra
3. `03_single_faculae` — same, faculae
4. `04_multi_umbrae` — four scattered umbrae
5. `05_multi_penumbrae` — same positions, penumbrae
6. `06_multi_faculae` — same positions, faculae
7. `07_concentric_AR` — umbra(2°) ⊂ penumbra(5°) ⊂ faculae(10°) at 15°N
8. `08_multi_concentric` — three concentric active regions
9. `09_activity_belts` — six concentric ARs in ±15° belts

**Qualitative result:** pure-spot scenarios flip a(λ) sign in the mid-UV; faculae-dominated
and mixed scenarios keep a(λ) > 0 across the visible.

---

## 7. Results so far (Planck-mixture layer)

- Planck-δT baseline: **a(λ) = (x/4)/(1−exp(−x))**, x = hc/(λkT); Pearson r ≡ 1 by construction.
- Mixture at fixed T=5780 K (cell 30): best u:p:f = **0.10 : 0.56 : 0.34**, RMS = **0.078**
  over 300–1800 nm.
- Joint fit (cell 31): best T = **5374 K** (~400 K cooler than nominal), u:p:f =
  **0.07 : 0.62 : 0.31**, RMS = **0.043** (≈45% better than the fixed-T fit).
- Pointwise inversion (cell 32): T(λ) of the best mixture rises ~5500 K (visible) → 10–12 kK
  (IR); a concentric active region goes the opposite way (down toward ~2 kK in the IR);
  faculae alone are intrinsically non-Planckian in the IR (a < 1/4 ⇒ inversion fails).

---

## 8. Open questions / next steps

1. **Alternative reference geometries** — re-run the cache + optimization (cells 29–32) for:
   - concentric AR co-located at (0°, −100°): umbra 2°, penumbra 5°, faculae 10°; or
   - spot-inside-halo: pure umbra 2° + faculae 10° (no penumbra) at (0°, −100°).
2. **IR-mismatch diagnostic** (two-panel figure): per-feature ΔSSI(λ,t)/ΔTSI(t) at peak
   rotation phase, and μ-dependent contrast (I_feat − I_qs)/I_qs at μ=1 vs μ=0.3 — to test
   whether faculae's μ-dependent contrast crossover drives the IR residual.

---

## 9. How to run

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
jupyter lab model.ipynb
```
Run cells **0, 4, 6, 7, 8, 10, 12, 14, 26** to load all definitions, then jump to cell 30
or 31 for mixture experiments. Skip cell 29 (`_per_feature_cache.npz` is included). To
regenerate all scenario figures + the PPTX from the command line, fix the `run_all.py`
whitelist snag noted in §3, then `python run_all.py && python build_pptx.py`.

---

## 10. Working style (how the original author likes to collaborate)

- Step-by-step: after each routine, run a short test, show a numerical/physical sanity
  check, then pause for confirmation before the next step.
- Get explicit approval before expensive multi-minute runs or disk writes.
- Concise, direct answers with `file:line` references; assume the reader has read the code.
- Always label axes: **ppm** on relative-variation axes, **nm** on wavelength,
  "dimensionless" for a(λ).
- Describe distributions in natural language, then translate to the `specs` list and show it
  before running.
