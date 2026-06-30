# Getting started (macOS) — run this project out of the box

This kit contains the **source code + data + a `requirements.txt`**, but **not** a Python
environment. You build the environment once (two commands, ~1 minute). After that it runs
in Jupyter, VS Code, or alongside Claude Code.

Prerequisite: **Python 3.10+**. Check with `python3 --version`. If missing, install from
<https://www.python.org/downloads/> or `brew install python`.

---

## Step 0 — unzip

Double-click the ZIP (or `unzip enigmatic_Planck_for_colleague.zip`). You get a folder with
`model.ipynb`, the four `*.txt` data files, `run_all.py`, `requirements.txt`,
`PROJECT_PRIMER.md`, `CLAUDE.md`, and a `scenarios/` folder. `cd` into that folder.

## Step 1 — create the environment (once)

```bash
cd /path/to/enigmatic_Planck      # the unzipped folder
python3 -m venv .venv             # create a private environment
source .venv/bin/activate         # activate it (prompt shows (.venv))
pip install -r requirements.txt   # install numpy/scipy/matplotlib/jupyterlab/Pillow/python-pptx
```

That's the whole setup. The `.venv` folder now holds everything.

---

## Run it in **Jupyter Lab**

```bash
source .venv/bin/activate         # if not already active
jupyter lab model.ipynb
```
In the notebook, **Run** cells **0, 4, 6, 7, 8, 10, 12, 14, 26** first (these define every
function), then jump to cell **30** or **31** for the mixture experiments. **Skip cell 29** —
its output (`scenarios/_per_feature_cache.npz`) is already included.

Quick end-to-end test in a fresh cell:
```python
result = run_scenario([{"type": "umbra", "shape": "circle",
                        "lat": 0, "lon": -100, "radius": 10}])
plot_diagnostics(result)
```

---

## Run it in **VS Code**

1. Install the **Python** and **Jupyter** extensions (Microsoft).
2. `File → Open Folder…` → choose the unzipped folder.
3. Open `model.ipynb`. Top-right **Select Kernel** → **Python Environments** → pick
   `.venv` (the one you created in Step 1).
4. Run cells the same order as above (0, 4, 6, 7, 8, 10, 12, 14, 26, then 30/31).

---

## Use it with **Claude Code** (macOS desktop app)

1. Open the Claude Code desktop app.
2. Open the **unzipped project folder** as the working directory (so Claude can see the files).
3. Claude Code automatically reads `CLAUDE.md` on startup. To get fully oriented, type as your
   **first message**:

   > Read `PROJECT_PRIMER.md` and `CLAUDE.md`, then give me a 5-line summary of what this
   > project does, the current results, and the open questions. Don't change anything yet.

4. To let Claude run the code, tell it once:

   > The environment is in `./.venv`. Activate it before running Python
   > (`source .venv/bin/activate`). `model.ipynb` is the source of truth; load cells
   > 0, 4, 6, 7, 8, 10, 12, 14, 26 to define all functions.

That's all the context Claude needs — `PROJECT_PRIMER.md` is fully self-contained, so no
prior-conversation history is required.

---

## Working from GitHub + staying in sync (collaboration loop)

The repo lives on GitHub. To get it and keep your local copy current:

```bash
git clone https://github.com/Sasha932-astro/enigmatic_Planck.git
cd enigmatic_Planck
git checkout claude/setup-project-Gj4Ui    # the active working branch
```

**One-time, in every clone** — install the notebook output-stripping filter so that *running*
cells locally never creates spurious git changes or merge conflicts:

```bash
source .venv/bin/activate     # or your env
pip install nbstripout
nbstripout --install          # activates the filter defined in .gitattributes for THIS clone
```

After that the loop is painless:

- **To get the latest changes:** `git pull` → in Jupyter, **File → Reload Notebook from Disk**
  (or restart the kernel). Because of nbstripout, your locally-run outputs won't block the pull.
- **The committed notebook is code-only by design.** Your plots still appear when you run cells;
  they just aren't stored in git (they regenerate, and key figures are shared inline / in
  `scenarios/`).

## Notes / gotchas

- The notebook in this kit has its **outputs stripped** (to keep the download small). They
  regenerate when you run the cells — nothing is missing.
- `run_all.py` (the batch script that regenerates all 9 scenario figures) currently has a
  small import snag — a leftover reference to `planck_temperature_slope`. The **notebook
  itself runs fine**; only that standalone script needs a one-line whitelist fix (see
  `PROJECT_PRIMER.md` §3) before `python run_all.py` will work.
- Always work inside the activated `.venv` so packages resolve correctly.
