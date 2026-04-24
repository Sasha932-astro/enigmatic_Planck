"""
Run all scenarios: per-scenario folder with regression.txt, diagnostic.png,
regression_big.png, animation.gif. Then assemble a single PPTX.

Core physics functions are loaded directly from model.ipynb so the notebook
remains the single source of truth.
"""
from __future__ import annotations
import json
import pathlib
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")  # headless
import matplotlib.pyplot as plt

HERE = pathlib.Path(__file__).resolve().parent
NB_PATH = HERE / "model.ipynb"
OUT = HERE / "scenarios"
OUT.mkdir(exist_ok=True)

# --- Load function definitions from the notebook -----------------------
nb = json.loads(NB_PATH.read_text())
_g = globals()

def _is_def_cell(cell):
    src = "".join(cell["source"])
    # Keep only cells that define our core functions (and constants in cell 4).
    names = ("read_mu_intensity_txt", "build_sphere_grid", "make_patch_mask",
             "make_circular_patch_mask", "build_distribution",
             "simulate_rotation_multi", "fit_ssi_vs_tsi",
             "run_scenario", "animate_scenario", "C_CGS =")
    return cell["cell_type"] == "code" and any(n in src for n in names)

for cell in nb["cells"]:
    if not _is_def_cell(cell):
        continue
    src = "".join(cell["source"])
    # Cell 4 ends with an "# example" block that executes file I/O. Strip it.
    marker = "# ---------- example"
    if marker in src:
        src = src.split(marker)[0]
    exec(compile(src, "<notebook>", "exec"), _g)

# Sanity check
for fn in ("read_mu_intensity_txt", "run_scenario", "animate_scenario"):
    assert fn in _g, f"{fn} not loaded from notebook"

os.chdir(HERE)  # so "qs.txt" etc. resolve for run_scenario


# --- Helpers: save-variants of the plots -------------------------------
def save_regression_txt(path: pathlib.Path, result) -> None:
    wav = result["wav_nm"]; a = result["slope"]; r = result["corr"]; valid = result["valid"]
    with path.open("w") as f:
        f.write("# wavelength_nm   slope_a   correlation_r\n")
        for i in np.where(valid)[0]:
            f.write(f"{wav[i]:12.4f}  {a[i]:+.8e}  {r[i]:+.8e}\n")


def save_regression_big(path: pathlib.Path, result, wav_range=(200, 2000)) -> None:
    """A clean, large 2-panel figure of a(λ) and r(λ) — the 'main result'."""
    wav = result["wav_nm"]; valid = result["valid"]
    mask = (wav >= wav_range[0]) & (wav <= wav_range[1]) & valid
    fig, ax = plt.subplots(2, 1, figsize=(11, 8), sharex=True)
    ax[0].plot(wav[mask], result["slope"][mask], color="C3", lw=1.5)
    ax[0].axhline(0, color="grey", lw=0.6)
    ax[0].set_ylabel("slope a(λ)  [ΔSSI/E$_{qs}$ per ΔTSI/TSI$_{qs}$]", fontsize=12)
    ax[0].set_title("Per-wavelength regression of SSI on TSI", fontsize=14)
    ax[0].grid(alpha=0.3)
    ax[1].plot(wav[mask], result["corr"][mask], color="C2", lw=1.5)
    ax[1].axhline(0, color="grey", lw=0.6)
    ax[1].set_xlabel("Wavelength (nm)", fontsize=12)
    ax[1].set_ylabel("Pearson correlation r(λ)", fontsize=12)
    ax[1].set_ylim(-1.05, 1.05)
    ax[1].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def save_diagnostic(path: pathlib.Path, result,
                    wavelengths=(400, 600, 800, 1000, 1200, 1600),
                    wav_range=(200, 2000)) -> None:
    """Same as plot_diagnostics in the notebook but saved to file."""
    wav_nm = result["wav_nm"]
    t_days = result["t_days"]
    x      = result["TSI_rel"]
    Y      = result["SSI_rel"]
    valid  = result["valid"]

    ks  = [int(np.argmin(np.abs(wav_nm - w))) for w in wavelengths]
    mask = (wav_nm >= wav_range[0]) & (wav_nm <= wav_range[1]) & valid

    n_scat = len(ks); n_cols = 3
    n_rows_scat = int(np.ceil(n_scat / n_cols))
    fig, ax = plt.subplots(1 + n_rows_scat, n_cols,
                           figsize=(4.2 * n_cols, 3.2 * (1 + n_rows_scat)))
    if ax.ndim == 1:
        ax = ax[None, :]

    ax[0, 0].plot(t_days, x * 1e6, color="black")
    ax[0, 0].axhline(0, color="grey", lw=0.5)
    ax[0, 0].set_xlabel("Time (days)")
    ax[0, 0].set_ylabel("ΔTSI / TSI$_{qs}$ (ppm)")
    ax[0, 0].set_title("TSI vs time")

    ax[0, 1].plot(wav_nm[mask], result["slope"][mask], color="C3")
    ax[0, 1].axhline(0, color="grey", lw=0.5)
    for w in wavelengths: ax[0, 1].axvline(w, color="grey", lw=0.4, ls=":")
    ax[0, 1].set_xlabel("Wavelength (nm)")
    ax[0, 1].set_ylabel("slope a(λ)")
    ax[0, 1].set_title("Regression slope")

    ax[0, 2].plot(wav_nm[mask], result["corr"][mask], color="C2")
    ax[0, 2].axhline(0, color="grey", lw=0.5)
    for w in wavelengths: ax[0, 2].axvline(w, color="grey", lw=0.4, ls=":")
    ax[0, 2].set_xlabel("Wavelength (nm)")
    ax[0, 2].set_ylabel("correlation r(λ)")
    ax[0, 2].set_ylim(-1.05, 1.05)
    ax[0, 2].set_title("Pearson correlation")

    for idx, k in enumerate(ks):
        row = 1 + idx // n_cols
        col = idx % n_cols
        lam = wav_nm[k]
        a_k = result["slope"][k]; b_k = result["intercept"][k]; r_k = result["corr"][k]
        ax[row, col].scatter(x * 1e6, Y[:, k] * 1e6, c=t_days, cmap="viridis", s=10)
        xline = np.linspace(x.min(), x.max(), 200)
        ax[row, col].plot(xline * 1e6, (a_k * xline + b_k) * 1e6,
                          color="red", lw=1.3,
                          label=f"a={a_k:.3g}\nr={r_k:+.4f}")
        ax[row, col].axhline(0, color="grey", lw=0.5)
        ax[row, col].axvline(0, color="grey", lw=0.5)
        ax[row, col].set_xlabel("ΔTSI/TSI$_{qs}$ (ppm)")
        ax[row, col].set_ylabel("ΔSSI/E$_{qs}$ (ppm)")
        ax[row, col].set_title(f"{lam:.1f} nm")
        ax[row, col].legend(loc="best", fontsize=8)

    for idx in range(n_scat, n_rows_scat * n_cols):
        row = 1 + idx // n_cols
        col = idx % n_cols
        ax[row, col].axis("off")

    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


def save_animation(path: pathlib.Path, result) -> None:
    anim = animate_scenario(
        result,
        wavelengths_nm=(400, 800, 1200),
        show_tsi=True,
        n_img=180, fps=15,
    )
    anim.save(str(path), writer="pillow", fps=15)


# --- Scenarios ----------------------------------------------------------
def concentric(lat, lon, r_umbra, r_penumbra, r_faculae):
    return [
        {"type": "faculae",  "shape": "circle", "lat": lat, "lon": lon, "radius": r_faculae},
        {"type": "penumbra", "shape": "circle", "lat": lat, "lon": lon, "radius": r_penumbra},
        {"type": "umbra",    "shape": "circle", "lat": lat, "lon": lon, "radius": r_umbra},
    ]


SCENARIOS = [
    dict(name="01_single_umbra",
         description="One 10° circular umbra at the equator; placed just behind the east limb at t=0.",
         specs=[{"type": "umbra", "shape": "circle", "lat": 0, "lon": -100, "radius": 10}]),

    dict(name="02_single_penumbra",
         description="One 10° circular penumbra at the equator; placed just behind the east limb at t=0.",
         specs=[{"type": "penumbra", "shape": "circle", "lat": 0, "lon": -100, "radius": 10}]),

    dict(name="03_single_faculae",
         description="One 10° circular faculae patch at the equator; just behind the east limb at t=0.",
         specs=[{"type": "faculae", "shape": "circle", "lat": 0, "lon": -100, "radius": 10}]),

    dict(name="04_multi_umbrae",
         description="Four scattered circular umbrae at varying latitudes and longitudes.",
         specs=[
             {"type": "umbra", "shape": "circle", "lat":  +20, "lon":  -60, "radius": 5},
             {"type": "umbra", "shape": "circle", "lat":  -15, "lon":   30, "radius": 4},
             {"type": "umbra", "shape": "circle", "lat":   +5, "lon": -150, "radius": 6},
             {"type": "umbra", "shape": "circle", "lat":  -25, "lon":  120, "radius": 3},
         ]),

    dict(name="05_multi_penumbrae",
         description="Four scattered circular penumbrae (same positions as multi_umbrae).",
         specs=[
             {"type": "penumbra", "shape": "circle", "lat":  +20, "lon":  -60, "radius": 5},
             {"type": "penumbra", "shape": "circle", "lat":  -15, "lon":   30, "radius": 4},
             {"type": "penumbra", "shape": "circle", "lat":   +5, "lon": -150, "radius": 6},
             {"type": "penumbra", "shape": "circle", "lat":  -25, "lon":  120, "radius": 3},
         ]),

    dict(name="06_multi_faculae",
         description="Four scattered circular faculae patches (same positions as multi_umbrae).",
         specs=[
             {"type": "faculae", "shape": "circle", "lat":  +20, "lon":  -60, "radius": 5},
             {"type": "faculae", "shape": "circle", "lat":  -15, "lon":   30, "radius": 4},
             {"type": "faculae", "shape": "circle", "lat":   +5, "lon": -150, "radius": 6},
             {"type": "faculae", "shape": "circle", "lat":  -25, "lon":  120, "radius": 3},
         ]),

    dict(name="07_concentric_AR",
         description="One concentric active region at (15°N, -90°): umbra(2°)⊂penumbra(5°)⊂faculae(10°).",
         specs=concentric(lat=15, lon=-90, r_umbra=2, r_penumbra=5, r_faculae=10)),

    dict(name="08_multi_concentric",
         description="Three concentric active regions at different positions; realistic spot+faculae complexes.",
         specs=(
             concentric(lat=+15, lon=-60, r_umbra=2, r_penumbra=5, r_faculae=10) +
             concentric(lat=-10, lon= 45, r_umbra=3, r_penumbra=6, r_faculae=12) +
             concentric(lat=+20, lon=120, r_umbra=2, r_penumbra=4, r_faculae= 8)
         )),

    dict(name="09_activity_belts",
         description="Butterfly-like two-belt activity at ±15°: 6 concentric active regions distributed in longitude.",
         specs=(
             concentric(lat=+15, lon=-150, r_umbra=2, r_penumbra=4, r_faculae= 8) +
             concentric(lat=-15, lon=-100, r_umbra=2, r_penumbra=5, r_faculae= 9) +
             concentric(lat=+12, lon= -40, r_umbra=3, r_penumbra=6, r_faculae=10) +
             concentric(lat=-18, lon=  20, r_umbra=2, r_penumbra=4, r_faculae= 8) +
             concentric(lat=+17, lon=  80, r_umbra=3, r_penumbra=5, r_faculae=10) +
             concentric(lat=-13, lon= 150, r_umbra=2, r_penumbra=4, r_faculae= 9)
         )),
]


# --- Run loop -----------------------------------------------------------
def run_one(sc):
    folder = OUT / sc["name"]
    folder.mkdir(exist_ok=True)
    print(f"\n=== {sc['name']} ===")
    print(f"    {sc['description']}")
    result = run_scenario(sc["specs"], n_steps=120, rotation_period_days=27.0, B0=0.0)

    counts = {k: int(v.sum()) for k, v in result["distribution"].items()}
    TSI_t  = np.trapezoid(result["E_t"], result["wav_nm"], axis=1)
    TSI_qs = np.trapezoid(result["E_qs"], result["wav_nm"])
    excursion = (TSI_t.max() - TSI_qs) * 1e6 / TSI_qs, (TSI_t.min() - TSI_qs) * 1e6 / TSI_qs
    print(f"    pixels: {counts}")
    print(f"    TSI excursion: +{excursion[0]:.1f} / {excursion[1]:+.1f} ppm")

    save_regression_txt(folder / "regression_coefficients.txt", result)
    save_diagnostic(folder / "diagnostic.png", result)
    save_regression_big(folder / "regression_big.png", result)
    save_animation(folder / "animation.gif", result)
    return result


if __name__ == "__main__":
    all_results = {}
    for sc in SCENARIOS:
        all_results[sc["name"]] = run_one(sc)

    # Save comparison overlays of a(λ) across all scenarios
    def _save_comparison(wav_lo, wav_hi, out_path, title_suffix=""):
        fig, ax = plt.subplots(2, 1, figsize=(12, 9), sharex=True)
        for sc in SCENARIOS:
            r = all_results[sc["name"]]
            wav = r["wav_nm"]; valid = r["valid"]
            mask = (wav >= wav_lo) & (wav <= wav_hi) & valid
            ax[0].plot(wav[mask], r["slope"][mask], lw=1.0, label=sc["name"].replace("_", " "))
            ax[1].plot(wav[mask], r["corr"][mask],  lw=1.0)
        ax[0].axhline(0, color="grey", lw=0.5)
        ax[0].set_ylabel("slope a(λ)", fontsize=12)
        ax[0].set_title(f"Comparison across scenarios{title_suffix}", fontsize=14)
        ax[0].legend(loc="best", fontsize=8, ncol=2)
        ax[0].grid(alpha=0.3)
        ax[1].axhline(0, color="grey", lw=0.5)
        ax[1].set_xlabel("Wavelength (nm)", fontsize=12)
        ax[1].set_ylabel("correlation r(λ)", fontsize=12)
        ax[1].set_ylim(-1.05, 1.05)
        ax[1].grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(out_path, dpi=150)
        plt.close(fig)
        print(f"Saved -> {out_path}")

    _save_comparison(200, 2000, OUT / "_comparison_all.png",       "")
    _save_comparison(350, 1300, OUT / "_comparison_350_1300.png",  " — zoomed (350–1300 nm)")
    print(f"\nAll scenario folders at: {OUT}")
