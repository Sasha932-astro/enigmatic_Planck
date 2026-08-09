# enigmatic_Planck — handoff

**Read this file first. It is a complete brief: the question, the answer so far, and how to run everything.**

## The question

The Sun's spectral irradiance responds to magnetic activity with a wavelength-dependent
sensitivity

    a(λ) = (ΔSSI/E_qs) / (ΔTSI/TSI_qs)         [dimensionless]

obtained by regressing SSI(λ,t) on TSI(t) as magnetic features cross the disk.
**Observations show the Sun follows the Planck-δT curve at its true effective temperature,
T = 5780 K**, i.e. it behaves as if the whole star simply changed temperature:

    a_Planck(λ; T) = (x/4) / (1 − e^−x),    x = hc / (λ k T)

That is physically "wrong" — the variability comes from a few localized dark and bright
features, not a uniform temperature change — yet it matches. **Why?**

## The answer so far: no model reproduces it

Every configuration is scored by ONE number: the RMS distance of its a(λ) from the
Planck **5780 K** curve over 400–1600 nm (0 = perfect match).

| model | distance from Planck(5780) |
|---|---|
| perfect match | 0 |
| **best of ~2000 configurations tried** (optimized MPS-ATLAS umbra+penumbra) | **0.060** |
| penumbra, best of 27 spatial layouts | 0.077 |
| umbra+penumbra spot at the real sunspot ratio, best of 64 | 0.077 |
| spots + faculae, best of 30 | 0.079 |
| umbra, best of 34 layouts | 0.095 |
| spots evolving, lifetimes 1–27 d over 5 rotations | 0.119 |
| spot+plage populations with different lifetimes, best | 0.128 |
| faculae only, best of 29 layouts | 0.132 |

**Nothing gets below ~0.06.** All models miss the curve the *same way*: a few percent too
high around 0.7–1.1 µm, far too low beyond ~1.3 µm.

What has been ruled out as the cause:

- **Spatial distribution** — a(λ) is intensive: latitude, longitude, size, number, and
  clustering of dark features barely change it (apparent-T spread ≈ 40–110 K).
  (Faculae *are* placement-sensitive, but they never improve the Planck match.)
- **Spectral library** — swapping SATIRE-family spectra for a self-consistent MPS-ATLAS
  radiative-equilibrium set at the same temperatures changes a(λ) by ~2 %.
- **Feature temperatures** — a scan of all 1518 (T_umbra, T_penumbra, area-fraction)
  combinations from 3500–5700 K never reaches distance < 0.05.
- **Geometry** — see `experiment_100K.py` below: removing rotation, foreshortening and
  limb darkening entirely does *not* improve the match.
- **Time evolution** — spot lifetimes from 1 day to 1 month, random emergence over
  5 rotations, and spot/plage populations with different lifetimes all fail too.
  (Differential lifetimes *can* sweep the best-fit temperature through 5780 K, but there
  the Planck *shape* is worse, the SSI–TSI correlation collapses to ≈ 0.86, and the result
  scatters by ±400 K between realizations — the Sun shows none of that.)

**Conclusion.** The mismatch is in the temperature response of the model atmospheres
themselves (line blanketing: molecular bands amplify the violet, the deep-forming infrared
is damped), not in how features are arranged in space or time. Candidates for the missing
physics: NLTE / chromospheric contributions, facular contrast physics beyond FAL-P-derived
tables, or a band/timescale selection effect in the observations.

## The decisive experiment (`experiment_100K.py`)

A feature only **100 K cooler** than the quiet Sun, using a self-consistent MPS-ATLAS pair
(quiet Sun 5800 K, feature 5700 K), in two scenarios:

- **Case A** — the feature rotates across the disk (normal transit).
- **Case B** — the feature grows and shrinks **at disk centre and does not rotate**, so the
  contrast is sampled at μ ≈ 1 only: no foreshortening, no limb darkening, no transit.

    Case A   distance from Planck(5780) = 0.112
    Case B   distance from Planck(5780) = 0.123    (SSI–TSI correlation +1.0000)

Case B is a perfectly clean, geometry-free measurement — and it misses the Planck curve
just as badly. **Geometry is exonerated; the spectra are responsible.**

Run it with `python experiment_100K.py` (≈ 1 min) → `scenarios/_100K_*.png`.

## Downloading MPS-ATLAS spectra (`mps_download.py`)

The MPS-ATLAS library (Edmond, doi:10.17617/3.NJ56TR; Witzke et al. 2021, Kostogryz et al.
2022) ships as two ~9.5 GB zips. This script reads the remote zip directory with HTTP range
requests and extracts **only the models you ask for** (a few hundred kB each, seconds):

```bash
python mps_download.py --teff 5700 5800        # the 100 K pair
python mps_download.py --all                   # all 56 temperatures, logg 4.4, [M/H]=0
python mps_download.py --teff 4500 --logg 4.5 --set set2 --kind flux
```

```python
from mps_download import fetch, read_clv, read_flux
wav, I, mu = read_clv(5800)        # I_nu(λ, μ) in CGS per Hz; downloads if missing
```

Grid: Teff 3500–9000 K (100 K steps), logg 3.0–5.0, [M/H] −5.0…+1.5, `set1` (Grevesse &
Sauval 1998, mixing length 1.25) or `set2` (Asplund 2009, Viani 2018).
Files land in `mps_atlas/<set>/MH<mh>/teff<T>/logg<logg>/`, the layout the official reader
(`mps_atlas/model_spectra_data.py`) expects. Solar-metallicity logg 4.4 spectra for all 56
temperatures are already in the repo.

## Repository map

| file | what it is |
|---|---|
| `model.ipynb` | **canonical source of all physics** — grid, μ-interpolation, rotation, regression, Planck curve, plus the study runners (guarded by `RUN_*` flags) |
| `experiment_100K.py` | the geometry test above (standalone) |
| `mps_download.py` | MPS-ATLAS downloader (standalone) |
| `qs.txt`, `umbra.txt`, `penumbra.txt`, `faculae.txt` | SATIRE-family input spectra I_ν(λ,μ), 1221 λ × 11 μ, CGS per Hz |
| `mps_atlas/` | MPS-ATLAS library + official reader |
| `scenarios/` | output figures |
| `run_all.py` | batch scenario runner + PPTX builder |
| `CLAUDE.md` | detailed session-by-session log of every study and its numbers |

Conventions that matter: features move left→right (ω = −2π/P); mask priority
umbra > penumbra > faculae; grid 180×360, 180 steps, P = 27 d, B₀ = 0; always filter by
`result["valid"]` before plotting; relative variations in ppm, λ in nm, a(λ) dimensionless.
**Contrast must always be computed within one spectral family** — never a spot from one
library against a quiet Sun from another.

## Setup

```bash
pip install -r requirements.txt
python experiment_100K.py
```

## Open questions

1. **Why is the Sun more Planckian than the models?** Even a genuine temperature
   perturbation of a self-consistent radiative-equilibrium atmosphere misses the 5780 K
   curve. What physics is missing — NLTE, chromosphere, real facular contrast?
2. **Facular μ-contrast**: dissect the facular spectra's centre-to-limb behaviour
   (μ = 1 vs μ = 0.3) and the infrared crossover that drives the collapse beyond 1.3 µm.
3. **Testable prediction**: models with mixed spot/plage lifetimes predict a *noisy*
   SSI–TSI regression (r ≈ 0.86) and ±400 K epoch-to-epoch scatter. If observations show a
   clean, stable relation across epochs with different activity mixes, this whole class of
   models is excluded and the spectra are the culprit.
