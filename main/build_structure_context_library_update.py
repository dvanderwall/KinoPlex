#!/usr/bin/env python3
"""
build_structure_context_library_update.py
(Updated 'script2' — disorder/boundary panel for STY sites)

Computes a compact panel of sequence-window features around Ser/Thr/Tyr residues:
- pLDDT run-lengths below thresholds (50/70), signed distance to OD boundaries
- AUC-under (threshold - pLDDT)+, total variation, mean pLDDT
- Uversky hydropathy–charge classifier (H̄ vs R)
- Composition/patterning: NCPR, κ (simplified), AA entropy, Frac_PG, Frac_BulkyHydro
- Structural IDR surrogates: CA radius of gyration, compaction ratios, coil fraction
- Low-confidence contact fraction (CA–CA <= 8 Å with either pLDDT < 70)
- MoRF-likeness in disordered context (hydrophobic 5-mers)

Key choices (documented, not sacred):
- Window half-size W=7 (±W) by default, configurable with -w/--window
- Contacts use CA–CA distances (Å). Cutoff = 8.0 Å
- CA coordinates from the parsed structure; missing coords are skipped

CLI:
  python build_structure_context_library_update.py STRUCTURE_DIR \
    -o sty_disorder_update.csv -d ./temp_results -t cif -f -p 8 -b 16 -m 0 -w 7 -q

Outputs:
  One CSV row per STY site with identifiers:
  UniProtID, ChainID, ResidueNumber, ResidueType (S/T/Y), Site
  and all metrics with suffixes _pm{W} as applicable.

This script is chain-aware and robust to PDB/mmCIF. AlphaFold pLDDT is parsed
from per-atom B-factors and averaged to per-residue means keyed by (chain, resseq, icode).
"""

import os
import sys
import math
import gzip
import time
import signal
import argparse
from collections import defaultdict, Counter
from functools import lru_cache

import numpy as np
import pandas as pd

from Bio.PDB import PDBParser, MMCIFParser, MMCIF2Dict
from Bio.PDB.Polypeptide import is_aa, three_to_index
from Bio.PDB.vectors import Vector

# ----------------------------
# Utilities
# ----------------------------

AA3_TO1 = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G",
    "HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S",
    "THR":"T","TRP":"W","TYR":"Y","VAL":"V","SEC":"U","PYL":"O"
}
STY = {"SER":"S","THR":"T","TYR":"Y"}

FILVWY = set(list("FILVWY"))
CHARGE = { # per residue net at pH ~7 (simplified)
    "D":-1, "E":-1, "K":+1, "R":+1, "H":0
}
KYTE_DOOLITTLE = {
    "I":4.5,"V":4.2,"L":3.8,"F":2.8,"C":2.5,"M":1.9,"A":1.8,"G":-0.4,"T":-0.7,"S":-0.8,
    "W":-0.9,"Y":-1.3,"P":-1.6,"H":-3.2,"E":-3.5,"Q":-3.5,"D":-3.5,"N":-3.5,"K":-3.9,"R":-4.5
}

def safe_mean(x):
    x = [v for v in x if v is not None and not np.isnan(v)]
    return float(np.mean(x)) if x else np.nan

def fasta_like_seq(residues):
    out = []
    for r in residues:
        resname = r.get_resname()
        out.append(AA3_TO1.get(resname, "X"))
    return "".join(out)

def residue_id_tuple(res):
    # (chain, resseq, icode)
    ch = res.get_parent().get_id()
    rid = res.get_id()
    return (ch, int(rid[1]), (rid[2] or " ").strip() if isinstance(rid[2], str) else " ")

def uniprot_from_name(fname):
    # For AF naming like AF-P12345-F1-model_v4.cif
    base = os.path.basename(fname)
    if "AF-" in base and "-F" in base:
        try:
            return base.split("AF-")[1].split("-F")[0]
        except Exception:
            return os.path.splitext(base)[0]
    return os.path.splitext(base)[0]

# ----------------------------
# File discovery
# ----------------------------

def discover_files(structure_dir, exts=("cif",), fallback_pdb=False, max_files=0):
    exts = tuple([e.lower().strip().lstrip(".") for e in exts if e.strip()])
    out = []
    for root, _, files in os.walk(structure_dir):
        for f in files:
            fl = f.lower()
            ok = False
            for e in exts:
                if fl.endswith("." + e) or fl.endswith("." + e + ".gz"):
                    ok = True
                    break
            if not ok and fallback_pdb:
                if fl.endswith(".pdb") or fl.endswith(".pdb.gz"):
                    ok = True
            if not ok:
                continue
            out.append(os.path.join(root, f))
            if max_files and len(out) >= max_files:
                return out
    return out

# ----------------------------
# pLDDT parsing (chain-aware)
# ----------------------------

def parse_plddt_map(structure_path):
    """
    Return dict keyed by (chain_id, resseq, icode) -> mean B-factor (pLDDT surrogate).
    Supports .cif/.pdb and their .gz variants; uses text parsing for speed/robustness.
    """
    opener = gzip.open if structure_path.lower().endswith(".gz") else open
    plddt = {}
    sums = defaultdict(float)
    counts = defaultdict(int)

    with opener(structure_path, "rt", encoding="utf-8", errors="ignore") as fh:
        is_pdb = any(structure_path.lower().endswith(x) for x in (".pdb",".pdb.gz"))
        is_cif = any(structure_path.lower().endswith(x) for x in (".cif",".cif.gz"))
        for line in fh:
            if is_pdb:
                if not line.startswith("ATOM"):
                    continue
                chain_id = (line[21].strip() or " ")
                try:
                    resseq = int(line[22:26])
                except ValueError:
                    continue
                icode = (line[26].strip() or " ")
                try:
                    b = float(line[60:66])
                except ValueError:
                    continue
                key = (chain_id, resseq, icode)
                sums[key] += b
                counts[key] += 1
            else:
                # mmCIF is trickier to parse by hand; try MMCIF2Dict as a fallback after this loop
                pass

    if sums:
        for k in sums:
            plddt[k] = sums[k] / max(1, counts[k])

    # If not populated (e.g., CIF path), try MMCIF2Dict
    if not plddt:
        try:
            d = MMCIF2Dict(structure_path)
            # Use label_asym_id (chain), label_seq_id (resseq), label_alt_id? mmCIF typically no icode
            chains = d.get("_atom_site.label_asym_id", [])
            seqids = d.get("_atom_site.label_seq_id", [])
            bvals  = d.get("_atom_site.B_iso_or_equiv", [])
            if isinstance(bvals, str):  # single atom case
                chains, seqids, bvals = [chains], [seqids], [bvals]
            for c, s, b in zip(chains, seqids, bvals):
                try:
                    resseq = int(s)
                    b = float(b)
                except Exception:
                    continue
                key = (c.strip() or "A", resseq, " ")
                sums[key] += b
                counts[key] += 1
            for k in sums:
                plddt[k] = sums[k] / max(1, counts[k])
        except Exception:
            pass

    return plddt

# ----------------------------
# Secondary structure via mmCIF _struct_conf (optional)
# ----------------------------

def parse_secondary_structure_map(structure_path):
    """
    Return dict keyed by (chain_id, resseq, icode) -> {'class': one of HELIX/STRAND/TURN/COIL}
    For PDB or CIF without struct_conf, returns {} (treated as unknown → COIL in some metrics).
    """
    out = {}
    if not any(structure_path.lower().endswith(x) for x in (".cif",".cif.gz")):
        return out
    try:
        d = MMCIF2Dict(structure_path)
        beg_chain = d.get("_struct_conf.beg_auth_asym_id", d.get("_struct_conf.beg_label_asym_id", []))
        end_chain = d.get("_struct_conf.end_auth_asym_id", d.get("_struct_conf.end_label_asym_id", []))
        beg_seq   = d.get("_struct_conf.beg_auth_seq_id", d.get("_struct_conf.beg_label_seq_id", []))
        end_seq   = d.get("_struct_conf.end_auth_seq_id", d.get("_struct_conf.end_label_seq_id", []))
        conf_type = d.get("_struct_conf.conf_type_id", [])
        if isinstance(conf_type, str):
            beg_chain, end_chain, beg_seq, end_seq, conf_type = [beg_chain],[end_chain],[beg_seq],[end_seq],[conf_type]
        for bc, ec, bs, es, ct in zip(beg_chain, end_chain, beg_seq, end_seq, conf_type):
            try:
                bsi = int(bs); esi = int(es)
            except Exception:
                continue
            ch = (bc or ec or "A").strip()
            lab = str(ct or "").upper()
            cls = "COIL"
            if "HELX" in lab or "HELIX" in lab:
                cls = "HELIX"
            elif "STRN" in lab or "SHEET" in lab or "BETA" in lab or "STRAND" in lab:
                cls = "STRAND"
            elif "TURN" in lab or "BEND" in lab:
                cls = "TURN"
            for i in range(bsi, esi+1):
                out[(ch, i, " ")] = {"class": cls}
    except Exception:
        pass
    return out

# ----------------------------
# Structure parsing (coordinates)
# ----------------------------

def load_structure(structure_path):
    parser = None
    if structure_path.lower().endswith((".cif",".cif.gz")):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    with (gzip.open(structure_path, "rt") if structure_path.lower().endswith(".gz") else open(structure_path, "rt")) as fh:
        # Bio.PDB parsers accept file paths; to avoid tempfile, pass path directly
        structure = parser.get_structure("model", structure_path)
    return structure

# ----------------------------
# Feature helpers (panel)
# ----------------------------

def window_indices(k, n, W):
    a = max(0, k - W); b = min(n, k + W + 1)
    return list(range(a,b))

def run_len_at_threshold(vec, thr, center_idx):
    """
    Return contiguous run length (< thr) that includes center index.
    vec: list of floats (pLDDT); NaN treated as >=thr (i.e., not disordered).
    """
    if center_idx<0 or center_idx>=len(vec): return 0
    def lt(i): 
        v = vec[i]
        return (not np.isnan(v)) and (v < thr)
    if not lt(center_idx):
        return 0
    L=1
    i=center_idx-1
    while i>=0 and lt(i):
        L+=1; i-=1
    i=center_idx+1
    while i<len(vec) and lt(i):
        L+=1; i+=1
    return L

def dist_to_boundary(vec, thr, center_idx):
    """
    Signed distance (in residues) from center to nearest order↔disorder boundary,
    negative inside 'disordered' (<thr), positive inside 'ordered' (>=thr).
    NaN treated as ordered.
    """
    n = len(vec)
    def state(i):
        v = vec[i]
        return 1 if (not np.isnan(v) and v < thr) else 0  # 1=disordered, 0=ordered
    s0 = state(center_idx)
    # search outward for index where state changes between i and i+1
    best = None
    for d in range(0, n):
        left = center_idx - d
        right = center_idx + d
        if left>0 and state(left) != state(left-1):
            best = d; break
        if right+1<n and state(right) != state(right+1):
            best = d; break
    if best is None:
        best = n  # no boundary in window
    return -best if s0==1 else best

def auc_under(vec, thr):
    v = np.array([0.0 if np.isnan(x) else max(0.0, thr - x) for x in vec], dtype=float)
    return float(np.sum(v))

def total_variation(vec):
    v = np.array([np.nan if x is None else x for x in vec], dtype=float)
    # interpolate NaNs with nearest valid to avoid propagating NaN through TV
    if np.all(np.isnan(v)): 
        return np.nan
    idx = np.where(~np.isnan(v))[0]
    if idx.size==0: return np.nan
    for i in range(len(v)):
        if np.isnan(v[i]):
            j = idx[np.argmin(np.abs(idx - i))]
            v[i] = v[j]
    return float(np.sum(np.abs(np.diff(v))))

def mean_plddt(vec):
    v = np.array([x for x in vec if not np.isnan(x)], dtype=float)
    return float(np.mean(v)) if v.size else np.nan

def charges_of(seq_chars):
    return [CHARGE.get(a,0) for a in seq_chars]

def ncp_r(seq_chars):
    c = charges_of(seq_chars)
    return float(sum(c))/max(1,len(c))

def kappa_simplified(seq_chars):
    """
    A simple, bounded [0,1] charge blockiness metric:
    κ = (number of adjacent pairs with same nonzero charge sign) / (number of adjacent pairs with any nonzero charge)
    """
    c = charges_of(seq_chars)
    same = 0; denom=0
    for i in range(len(c)-1):
        if c[i]==0 or c[i+1]==0: 
            continue
        denom += 1
        if (c[i]>0 and c[i+1]>0) or (c[i]<0 and c[i+1]<0):
            same += 1
    return float(same)/denom if denom>0 else 0.0

def aa_entropy(seq_chars):
    counts = Counter(seq_chars)
    n = float(len(seq_chars))
    if n==0: return np.nan
    H = 0.0
    for a,k in counts.items():
        p = k/n
        H -= p*math.log(p,2)
    return H

def radius_of_gyration(coords):
    # coords: list of (x,y,z), may include None entries; ignore missing
    pts = np.array([c for c in coords if c is not None], dtype=float)
    if pts.size==0:
        return np.nan
    center = np.mean(pts, axis=0)
    diff = pts - center
    return float(np.sqrt(np.mean(np.sum(diff*diff, axis=1))))

def morf_like_features(seq_chars, mean_plddt_window, plddt_thr=70, k=5):
    """
    Identify hydrophobic micro-clusters (>=3 FILVWY in any k-mer) in windows with mean pLDDT < plddt_thr.
    Return (cluster_count, max_cluster_run_len).
    """
    if np.isnan(mean_plddt_window) or mean_plddt_window >= plddt_thr:
        return 0, 0
    n = len(seq_chars)
    flagged = [False]*n
    for i in range(0, n-k+1):
        block = seq_chars[i:i+k]
        if sum(1 for a in block if a in FILVWY) >= 3:
            for j in range(i, i+k):
                flagged[j] = True
    # clusterization over flagged
    cnt = 0; maxlen = 0
    i=0
    while i<n:
        if not flagged[i]:
            i+=1; continue
        j=i
        while j<n and flagged[j]:
            j+=1
        run_len = j-i
        cnt += 1
        maxlen = max(maxlen, run_len)
        i=j
    return cnt, maxlen

def coil_fraction(ss_classes):
    if not ss_classes: return np.nan
    n = len(ss_classes); return float(sum(1 for s in ss_classes if s=="COIL"))/n

def ca_distance(a, b):
    if a is None or b is None: return np.nan
    d = a - b
    return float(np.sqrt(np.dot(d,d)))

# ----------------------------
# Core processing
# ----------------------------

def process_structure(structure_path, window=7, file_types=("cif",), fallback_pdb=False, quiet=False):
    """
    Returns a list of dicts (rows). Each row is one STY site with all metrics.
    """
    t0 = time.time()
    try:
        structure = load_structure(structure_path)
    except Exception as e:
        if not quiet:
            print(f"[WARN] Failed to parse structure: {structure_path}: {e}")
        return []
    plddt_map = parse_plddt_map(structure_path)
    ss_map = parse_secondary_structure_map(structure_path)
    uniprot = uniprot_from_name(structure_path)

    rows = []
    for model in structure:
        for chain in model:
            chain_id = chain.get_id()
            # sequence of residues for indexing
            residues = [r for r in chain if is_aa(r, standard=True)]
            if not residues:
                continue
            seq = [AA3_TO1.get(r.get_resname(),"X") for r in residues]
            # CA coords per residue
            ca_coords = []
            keys = []
            for r in residues:
                key = residue_id_tuple(r)
                keys.append(key)
                ca = r["CA"].get_vector().get_array() if "CA" in r else None
                ca_coords.append(ca)

            n = len(residues)
            for k, r in enumerate(residues):
                if r.get_resname() not in STY:
                    continue
                aa1 = STY[r.get_resname()]
                rid = keys[k][1]

                # build window
                idxs = window_indices(k, n, window)
                window_seq = [seq[i] for i in idxs]
                window_ca = [ca_coords[i] for i in idxs]
                window_keys = [keys[i] for i in idxs]
                window_plddt = [plddt_map.get(keys[i], np.nan) for i in idxs]
                window_ss = [ss_map.get(keys[i], {"class":"COIL"}).get("class","COIL") for i in idxs]

                # pLDDT-centric metrics
                center = idxs.index(k)
                rl50 = run_len_at_threshold(window_plddt, 50.0, center)
                rl70 = run_len_at_threshold(window_plddt, 70.0, center)
                dOD50 = dist_to_boundary(window_plddt, 50.0, center)
                dOD70 = dist_to_boundary(window_plddt, 70.0, center)
                auc50 = auc_under(window_plddt, 50.0)
                auc70 = auc_under(window_plddt, 70.0)
                tv = total_variation(window_plddt)
                mp = mean_plddt(window_plddt)

                # Uversky (H vs R)
                H_vals = [KYTE_DOOLITTLE.get(a, 0.0) for a in window_seq]
                Hbar = float(np.mean(H_vals)) if H_vals else np.nan
                R = ncp_r(window_seq)
                uversky_flag = 1 if (not np.isnan(Hbar) and Hbar < 1.151*R + 0.660) else 0

                # Composition/patterning
                entropy = aa_entropy(window_seq)
                lcr_flag = 1 if (not np.isnan(entropy) and entropy < 2.2) else 0
                frac_pg = float(sum(1 for a in window_seq if a in ("P","G")))/len(window_seq) if window_seq else np.nan
                frac_bulky = float(sum(1 for a in window_seq if a in FILVWY))/len(window_seq) if window_seq else np.nan
                kappa = kappa_simplified(window_seq)
                ncpr = ncp_r(window_seq)

                # CA radius of gyration & compaction surrogates
                rg = radius_of_gyration(window_ca)
                nW = len([c for c in window_ca if c is not None])
                comp_norm = (rg / (nW**(1/3))) if (not np.isnan(rg) and nW>0) else np.nan
                comp_fold = (rg / (2.2*(nW**(1/3)))) if (not np.isnan(rg) and nW>0) else np.nan

                # Coil fraction
                coil_frac = coil_fraction(window_ss)

                # Low-confidence contact fraction (CA–CA <=8 Å, either pLDDT < 70)
                pairs = []
                lowconf_contacts = 0
                total_contacts = 0
                for i in range(len(window_ca)):
                    for j in range(i+1, len(window_ca)):
                        ci = window_ca[i]; cj = window_ca[j]
                        if ci is None or cj is None: 
                            continue
                        d = ca_distance(ci, cj)
                        if d <= 8.0:
                            total_contacts += 1
                            pli = window_plddt[i]; plj = window_plddt[j]
                            if (not np.isnan(pli) and pli < 70.0) or (not np.isnan(plj) and plj < 70.0):
                                lowconf_contacts += 1
                lowconf_frac = (lowconf_contacts/total_contacts) if total_contacts>0 else np.nan

                # MoRF-likeness
                cl_cnt, cl_max = morf_like_features(window_seq, mp, plddt_thr=70, k=5)

                row = {
                    "UniProtID": uniprot,
                    "ChainID": chain_id,
                    "ResidueNumber": rid,
                    "ResidueType": aa1,
                    "Site": f"{uniprot}_{rid}",
                    "LowPLDDT_RunLen50": rl50,
                    "LowPLDDT_RunLen70": rl70,
                    "DistTo_OD_Boundary50": dOD50,
                    "DistTo_OD_Boundary70": dOD70,
                    "pLDDT_AUC50_pm{W}".format(W=window): auc50,
                    "pLDDT_AUC70_pm{W}".format(W=window): auc70,
                    "pLDDT_TV_pm{W}".format(W=window): tv,
                    "Mean_pLDDT_pm{W}".format(W=window): mp,
                    "Uversky_H_pm{W}".format(W=window): Hbar,
                    "Uversky_R_pm{W}".format(W=window): R,
                    "Uversky_IDP_pm{W}".format(W=window): uversky_flag,
                    "NCPR_pm{W}".format(W=window): ncpr,
                    "Kappa_pm{W}".format(W=window): kappa,
                    "AA_Entropy_pm{W}".format(W=window): entropy,
                    "LCR_flag_pm{W}".format(W=window): lcr_flag,
                    "Frac_PG_pm{W}".format(W=window): frac_pg,
                    "Frac_BulkyHydro_pm{W}".format(W=window): frac_bulky,
                    "Rg_CA_pm{W}".format(W=window): rg,
                    "CompactionNorm_pm{W}".format(W=window): comp_norm,
                    "CompactionRatio_Folded_pm{W}".format(W=window): comp_fold,
                    "CoilFrac_pm{W}".format(W=window): coil_frac,
                    "LowConf_ContactFrac_pm{W}".format(W=window): lowconf_frac,
                    "MoRFlike_ClusterCount_pm{W}".format(W=window): cl_cnt,
                    "MoRFlike_MaxClusterLen_pm{W}".format(W=window): cl_max,
                }
                rows.append(row)

    t1 = time.time()
    return rows

# ----------------------------
# CLI and batch execution
# ----------------------------

def main():
    ap = argparse.ArgumentParser(description="Disorder/boundary panel for STY sites (updated script2)")
    ap.add_argument("structure_dir", nargs="?", default="Complete_AF_Proteome", help="Directory of structures (.cif/.pdb/.gz)")
    ap.add_argument("-o", "--output", default="sty_disorder_update.csv", help="Output CSV filename")
    ap.add_argument("-d", "--temp-dir", default="./temp_results", help="Directory for temp CSVs")
    ap.add_argument("-t", "--file-types", default="cif", help="Comma-separated list: e.g., cif,pdb")
    ap.add_argument("-f", "--fallback-pdb", action="store_true", help="If set, also include PDB files when CIF not present")
    ap.add_argument("-p", "--num-processes", type=int, default=0, help="Worker processes (0=auto, 1=single-thread)")
    ap.add_argument("-b", "--batch-size", type=int, default=16, help="Files per batch submitted")
    ap.add_argument("-m", "--max-files", type=int, default=0, help="Max number of files to process (0=all)")
    ap.add_argument("-w", "--window", type=int, default=7, help="Half-window size (±W)")
    ap.add_argument("-q", "--quiet", action="store_true", help="Suppress per-file messages")
    args = ap.parse_args()

    os.makedirs(args.temp_dir, exist_ok=True)

    exts = [x.strip() for x in args.file_types.split(",") if x.strip()]
    files = discover_files(args.structure_dir, exts=exts, fallback_pdb=args.fallback_pdb, max_files=args.max_files)
    if not files:
        print("No structure files found. Check --file-types and directory.")
        sys.exit(1)

    # Batch processing with optional multiprocessing
    batches = [files[i:i+args.batch_size] for i in range(0, len(files), args.batch_size)]
    tmp_paths = []

    def process_batch(batch, batch_idx):
        rows_all = []
        for path in batch:
            if not args.quiet:
                print(f"[{batch_idx}] Processing {os.path.basename(path)}")
            try:
                rows = process_structure(path, window=args.window, file_types=exts, fallback_pdb=args.fallback_pdb, quiet=args.quiet)
                rows_all.extend(rows)
            except Exception as e:
                if not args.quiet:
                    print(f"[WARN] {path}: {e}")
        out = pd.DataFrame(rows_all)
        tmp = os.path.join(args.temp_dir, f"panel_batch_{batch_idx:04d}.csv")
        out.to_csv(tmp, index=False)
        return tmp

    if args.num_processes and args.num_processes > 1:
        import multiprocessing as mp
        try:
            mp.set_start_method("spawn")
        except RuntimeError:
            pass
        with mp.Pool(processes=args.num_processes) as pool:
            for i, tmp in enumerate(pool.starmap(process_batch, [(b, bi) for bi, b in enumerate(batches)])):
                tmp_paths.append(tmp)
    else:
        for bi, b in enumerate(batches):
            tmp_paths.append(process_batch(b, bi))

    # Concatenate and write output
    df_all = pd.concat((pd.read_csv(p) for p in tmp_paths), ignore_index=True) if tmp_paths else pd.DataFrame()
    df_all.to_csv(args.output, index=False)
    print(f"[OK] Wrote {args.output} with {len(df_all)} rows.")
    # Cleanup temp CSVs
    for p in tmp_paths:
        try:
            os.remove(p)
        except OSError:
            pass

if __name__ == "__main__":
    # Handle Ctrl+C gracefully
    def signal_handler(sig, frame):
        print("Received termination signal. Exiting.")
        sys.exit(0)
    signal.signal(signal.SIGINT, signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)
    main()
