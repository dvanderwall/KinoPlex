#!/usr/bin/env python3
"""
run_structural_context_pipeline.py

One-stop CLI to produce:
  1) Structural-context atlas (build_structure_context_library.py)
  2) Disorder/boundary panel (script2.py)
  3) A merged CSV keyed by (UniProtID, ChainID, ResidueNumber, ResidueType, Site)

Requirements: both source scripts must be present in the same directory:
 - build_structure_context_library.py
 - script2.py

Example:
  python run_structural_context_pipeline.py Complete_AF_Proteome \
      -o combined_features.csv \
      -t cif \
      --fallback-pdb \
      -p 8 -b 16 -w 7 -m 0 \
      -d ./temp_results \
      --keep-intermediate

Use --only-atlas or --only-panel if you want to run a single stage.
"""

import argparse
import os
import sys
import subprocess
import pandas as pd

HERE = os.path.abspath(os.path.dirname(__file__))

def exe():
    return sys.executable or "python"

def run(cmd):
    print("[RUN]", " ".join(cmd))
    res = subprocess.run(cmd, stdout=sys.stdout, stderr=sys.stderr, check=True)
    return res.returncode

def main():
    ap = argparse.ArgumentParser(description="Run structural atlas + disorder panel and merge outputs")
    ap.add_argument("structure_dir", help="Directory containing structure files (.cif/.pdb/.gz)")
    ap.add_argument("-o", "--output", default="combined_structural_features.csv", help="Final merged CSV filename")
    ap.add_argument("-d", "--temp-dir", default="./temp_results", help="Directory for intermediate outputs")
    ap.add_argument("-t", "--file-types", default="cif", help="Comma-separated file types to include (e.g., cif,pdb)")
    ap.add_argument("-f", "--fallback-pdb", action="store_true", help="If set, include PDB when CIF not present")
    ap.add_argument("-p", "--num-processes", type=int, default=0, help="Worker processes for script2 (0=auto)")
    ap.add_argument("-b", "--batch-size", type=int, default=16, help="Batch size for script2 file submissions")
    ap.add_argument("-w", "--window", type=int, default=7, help="Half-window size for script2 (±window)")
    ap.add_argument("-m", "--max-files", type=int, default=0, help="Max files to process (0=all)")
    ap.add_argument("-q", "--quiet", action="store_true", help="Suppress per-file output from script2")
    ap.add_argument("--only-atlas", action="store_true", help="Only run the atlas (build_structure_context_library.py)")
    ap.add_argument("--only-panel", action="store_true", help="Only run the panel (script2.py)")
    ap.add_argument("--keep-intermediate", action="store_true", help="Keep intermediate CSVs")
    args = ap.parse_args()

    os.makedirs(args.temp_dir, exist_ok=True)

    atlas_csv = os.path.join(args.temp_dir, "atlas_features.csv")
    panel_csv = os.path.join(args.temp_dir, "panel_features.csv")

    build_py = os.path.join(HERE, "build_structure_context_library.py")
    panel_py = os.path.join(HERE, "build_structure_context_library_update.py")

    if not os.path.exists(build_py):
        print(f"ERROR: Missing {build_py}")
        sys.exit(1)
    if not os.path.exists(panel_py):
        print(f"ERROR: Missing {panel_py}")
        sys.exit(1)

    if not args.only_panel:
        # Run atlas builder
        cmd_build = [
            exe(), build_py,
            args.structure_dir,
            "--output", atlas_csv,
            "--temp-dir", args.temp_dir,
            "--file-types", args.file_types,
            "--max-files", str(args.max_files),
        ]
        if args.fallback_pdb:
            cmd_build.append("--fallback-pdb")
        # --summary is optional if present in the build script
        try:
            run(cmd_build + ["--summary"])
        except subprocess.CalledProcessError:
            # retry without --summary (in case the script lacks it)
            run(cmd_build)

    if not args.only_atlas:
        # Run panel generator
        cmd_panel = [
            exe(), panel_py,
            args.structure_dir,
            "-o", panel_csv,
            "-d", args.temp_dir,
            "-t", args.file_types,
            "-m", str(args.max_files),
            "-p", str(args.num_processes),
            "-b", str(args.batch_size),
            "-w", str(args.window),
        ]
        if args.fallback_pdb:
            cmd_panel.append("-f")
        if args.quiet:
            cmd_panel.append("-q")
        run(cmd_panel)

    # Merge outputs if we have both
    have_atlas = os.path.exists(atlas_csv)
    have_panel = os.path.exists(panel_csv)

    if args.only_panel and not have_panel:
        print("Panel did not produce output (expected at {})".format(panel_csv))
        sys.exit(2)

    if args.only_atlas and not have_atlas:
        print("Atlas did not produce output (expected at {})".format(atlas_csv))
        sys.exit(2)

    if have_atlas and have_panel:
        left = pd.read_csv(atlas_csv)
        right = pd.read_csv(panel_csv)

        # Standardize key set and merge
        key_cols = ["UniProtID", "ChainID", "ResidueNumber", "ResidueType", "Site"]
        for col in key_cols:
            if col not in left.columns or col not in right.columns:
                print(f"ERROR: Missing key column '{col}' in one of the outputs.")
                print("Atlas columns:", list(left.columns)[:20], "...")
                print("Panel columns:", list(right.columns)[:20], "...")
                sys.exit(3)

        combined = pd.merge(left, right, on=key_cols, how="outer", suffixes=("_atlas", "_panel"))
        combined.to_csv(args.output, index=False)
        print(f"[OK] Wrote merged CSV: {args.output}")
        if not args.keep_intermediate:
            try:
                os.remove(atlas_csv)
                os.remove(panel_csv)
            except OSError:
                pass
    elif have_atlas:
        # Atlas only → rename to requested output
        os.replace(atlas_csv, args.output)
        print(f"[OK] Wrote atlas-only CSV: {args.output}")
    elif have_panel:
        os.replace(panel_csv, args.output)
        print(f"[OK] Wrote panel-only CSV: {args.output}")
    else:
        print("ERROR: No outputs produced.")
        sys.exit(4)

if __name__ == "__main__":
    main()
