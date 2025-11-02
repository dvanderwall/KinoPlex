#!/usr/bin/env python3
"""
Run Script: executes two STY feature builders with the *same* args, then merges outputs.

- Script 1 (context): build_structure_context_library.py
- Script 2 (disorder): build_structure_context_library_update_2.py
  (If your file is named differently, use --script2 to point to it.)

Join keys: (UniProtID, ChainID, ResidueNumber)
Result: a single CSV with all columns (outer join), with duplicate names suffixed.

Example:
  python run_sty_dual.py ./Complete_AF_Proteome \
    -t cif -p 8 -b 10 -m 0 --fallback-pdb \
    -w 7 \
    -o out/structure_context_disorder_merged.csv
"""
import os
import sys
import argparse
import subprocess
import pandas as pd
from datetime import datetime

def run(cmd, cwd=None):
    print(">>", " ".join(cmd))
    res = subprocess.run(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    print(res.stdout, end="")
    if res.returncode != 0:
        raise SystemExit(f"Command failed (exit {res.returncode}): {' '.join(cmd)}")

def main():
    ap = argparse.ArgumentParser(description="Run both STY builders with shared args, then merge CSVs.")
    # Common, forwarded to both scripts
    ap.add_argument("structure_dir", help="Directory with structures (e.g., Complete_AF_Proteome)")
    ap.add_argument("-t", "--file-types", default="cif", help="Comma-separated types (e.g., cif,pdb)")
    ap.add_argument("-m", "--max-files", type=int, default=0, help="Max files (0=all)")
    ap.add_argument("-p", "--num-processes", type=int, default=0, help="Workers (0=auto)")
    ap.add_argument("-b", "--batch-size", type=int, default=10, help="Batch size per submit")
    ap.add_argument("-f", "--fallback-pdb", action="store_true", help="Also allow PDB if CIF missing")
    # Only script2 needs this; we accept it as common and forward only to script2
    ap.add_argument("-w", "--window", type=int, default=7, help="±N window for disorder features (script2 only)")
    # Runner outputs/paths
    ap.add_argument("-o", "--out", default="sty_context_disorder_merged.csv", help="Merged output CSV")
    ap.add_argument("--work", default="./runner_temp", help="Working directory for intermediates")
    # Script paths (override if your filenames differ)
    ap.add_argument("--script1", default="build_structure_context_library.py",
                    help="Path to script 1 (context/geometry)")
    ap.add_argument("--script2", default="build_structure_context_library_update_2.py",
                    help="Path to script 2 (disorder update)")
    # Keep intermediates?
    ap.add_argument("--keep-intermediate", action="store_true", help="Keep per-script CSVs/temp dirs")
    args = ap.parse_args()

    os.makedirs(args.work, exist_ok=True)
    s1_temp = os.path.join(args.work, "s1_temp")
    s2_temp = os.path.join(args.work, "s2_temp")
    os.makedirs(s1_temp, exist_ok=True)
    os.makedirs(s2_temp, exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    s1_out = os.path.join(args.work, f"s1_context_{timestamp}.csv")
    s2_out = os.path.join(args.work, f"s2_disorder_{timestamp}.csv")

    # --- Run Script 1 (context) ---
    s1_cmd = [
        sys.executable, args.script1,
        args.structure_dir,
        "-o", s1_out,
        "-t", args.file_types,
        "-m", str(args.max_files),
        "-p", str(args.num_processes),
        "-b", str(args.batch_size),
        "-d", s1_temp,
    ]
    if args.fallback_pdb:
        s1_cmd.append("--fallback-pdb")
    run(s1_cmd)

    # --- Run Script 2 (disorder) ---
    s2_cmd = [
        sys.executable, args.script2,
        args.structure_dir,
        "-o", s2_out,
        "-t", args.file_types,
        "-m", str(args.max_files),
        "-p", str(args.num_processes),
        "-b", str(max(args.batch_size, 1)),  # script2 default 16; any positive is fine
        "-d", s2_temp,
        "-w", str(args.window),
    ]
    if args.fallback_pdb:
        s2_cmd.append("--fallback-pdb")
    run(s2_cmd)

    # --- Merge ---
    print(f"\nMerging:\n  S1: {s1_out}\n  S2: {s2_out}\n=> {args.out}")
    df1 = pd.read_csv(s1_out)
    df2 = pd.read_csv(s2_out)

    join_keys = ["UniProtID", "ChainID", "ResidueNumber"]
    for k in join_keys:
        if k not in df1.columns or k not in df2.columns:
            raise SystemExit(f"Missing join key '{k}' in one of the inputs.")

    merged = df1.merge(df2, on=join_keys, how="outer", suffixes=("_ctx", "_dis"))
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    merged.to_csv(args.out, index=False)

    print(f"\nDone. Rows: s1={len(df1):,}, s2={len(df2):,}, merged={len(merged):,}")
    print(f"Merged CSV written to: {args.out}")

    if not args.keep_intermediate:
        # best-effort cleanup of per-run CSVs (keeps temp dirs for each script’s own batching)
        try:
            os.remove(s1_out)
            os.remove(s2_out)
        except Exception:
            pass

if __name__ == "__main__":
    main()
