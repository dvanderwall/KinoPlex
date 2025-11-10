# ATTENTION USERS AND MANUSCRIPT REVIEWERS
The KinoPlex Website is hosted on PhosphositePlus by CST at https://kinoplex.phosphosite.org and code can be found at https://github.com/dvanderwall/KinoPlex_CST_WebApp

# Structural Library

Tools for building a structural-context library from PDB/mmCIF structures and running the end-to-end pipeline.

## Folder layout

```
Structural_Library/
├─ main/
│  ├─ build_structure_context_library.py
│  ├─ build_structure_context_library_update_v2.py
│  └─ run_structural_context_pipeline.py
├─ Structure_Depot/          # put your .cif/.pdb (optionally .gz) here
├─ temp_results/             # scratch/intermediate outputs
├─ runner_temp/              # temporary working dir
├─ results.csv               # example final results file
├─ requirements.txt
└─ setup.py
```

---

## 1) Installation

> Requires Python 3.10–3.12.

### Create & activate a virtual environment



### Install dependencies
```bash
pip install --upgrade pip
pip install -r requirements.txt
```

The minimal `requirements.txt` expected by the scripts is:
```text
numpy>=1.24,<3
pandas>=2.1,<3
scipy>=1.11,<2
biopython>=1.83,<2
tqdm>=4.66,<5
psutil>=5.9,<6
```

---

## 2) Prepare inputs

Place your input structures in `Structure_Depot/` (several structures included as toy examples).  
Accepted formats: `.cif`, `.pdb`, and their gzipped equivalents (`.cif.gz`, `.pdb.gz`).  
You can organize them in subfolders if you like.

---

## 3) Running from the command line

Each script supports `-h/--help` to print the exact arguments. Run this first to confirm names:
```bash
python main/run_structural_context_pipeline.py -h
python main/build_structure_context_library.py -h
python main/build_structure_context_library_update_v2.py -h
```

Below are **typical** invocations. If an option name differs on your machine, use `-h` and substitute (e.g., some scripts use `--input` instead of `--structures`, or `--workers` instead of `--n-procs`).

### A) End-to-end pipeline (recommended)
Runs the full process on everything in `Structure_Depot/` and writes intermediates to `temp_results/` and a summary CSV.

```bash
python main/run_structural_context_pipeline.py   --structures Structure_Depot   --out temp_results   --results results.csv   --n-procs 8
```

Common optional flags you may see:
- `--resume` : skip structures already processed in `temp_results/`
- `--verbose` or `-v` : more logging
- `--seed 42` : make any randomized steps reproducible

### B) Build the library from scratch
```bash
python main/build_structure_context_library.py   --structures Structure_Depot   --out temp_results   --results results.csv   --n-procs 8
```

### C) Incremental update (add new/changed structures only)
```bash
python main/build_structure_context_library_update_v2.py   --structures Structure_Depot   --out temp_results   --results results.csv   --n-procs 8   --resume
```

> Tip: If your script prefers `--workers`, use that in place of `--n-procs`.

---

## 4) Outputs

- **Intermediate artifacts**: written under `temp_results/` (and/or `runner_temp/`).
- **Final table**: `results.csv` (path can be changed with `--results`).

Progress bars (from `tqdm`) will show per-file processing; CPU usage can be controlled with `--n-procs/--workers`.

---

## 5) Troubleshooting

- `ModuleNotFoundError: Bio` → `pip install biopython` (already in `requirements.txt`).
- `ImportError: scipy.spatial.cKDTree` → ensure SciPy installed and Python version is supported.
- Process gets killed / runs out of memory → lower `--n-procs` (try 2–4), free disk in `temp_results/`.
- “unrecognized arguments” → check `-h` and adjust flag names (`--input` vs `--structures`, `--output-dir` vs `--out`, etc.).
- Gzipped files not read → ensure file extension ends with `.gz` and that the script’s input glob includes `*.gz`.

---

## 6) Reproducing an example run

```bash
# from the project root
source .venv/bin/activate            # or Windows Activate.ps1
python main/run_structural_context_pipeline.py   --structures Structure_Depot   --out temp_results   --results results.csv   --n-procs 8 --resume -v
```

When it finishes, open `results.csv` to inspect the aggregated results.

---


