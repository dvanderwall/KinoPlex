
# Write the structure_library_update.py script implementing disorder features (1)-(5)
# -*- coding: utf-8 -*-
"""
structure_library_update.py

Compute a compact panel of *disorder- and boundary*-centric features around
Ser/Thr/Tyr (STY) residues from AlphaFold structures — as a standalone updater.
This script DOES NOT merge into your existing feature CSV; it writes its own file.

Features implemented (per STY site):
(1) pLDDT run-lengths, boundary proximity, AUC-under-thresholds, total variation
(2) Uversky charge–hydropathy classifier in ±N window (mean hydropathy, NCPR, flag)
(3) Charge patterning & composition: NCPR, κ (simplified), Pro/Gly fraction,
    bulky-hydrophobe fraction, AA Shannon entropy, low-complexity flag
(4) Structural “IDR-ness” from model: local CA radius of gyration (Rg),
    compaction ratios vs N^(1/3) and a folded baseline, coil fraction (if CIF SS),
    low-confidence contact fraction (CA–CA ≤ 8Å where either pLDDT < 70)
(5) Simple MoRF-likeness: count and maximum length of hydrophobic micro-clusters
    (≥3 of FILVWY in any 5-res sliding window) within ±N, in a disordered context

Keys in output:
- UniProtID, ChainID, ResidueNumber, ResidueType (S/T/Y), Site (UniProt_#)
- All metrics above with suffixes “_pm{N}” for the ±N window where applicable.
- No merging occurs here; you can join by (UniProtID, ChainID, ResidueNumber).

Usage example:
  python structure_library_update.py ./Complete_AF_Proteome \
      -o sty_disorder_update.csv -p 8 -t cif -w 7 --max-files 0

Notes:
- Optimized for AlphaFold mmCIF; PDB works too (uses B-factor for pLDDT).
- Secondary structure (coil fraction) uses mmCIF _struct_conf if present.
- Parallel processing writes per-file temp CSVs then concatenates them.

Author: generated for your workflow
"""

import os
import re
import sys
import gzip
import math
import time
import argparse
import traceback
from collections import defaultdict, Counter

import numpy as np
import pandas as pd

from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.Polypeptide import is_aa, three_to_index

# ---------- Utility: three_to_one (no Biopython dependency on three_to_one) ----------
def three_to_one(res_name: str) -> str:
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    try:
        idx = three_to_index(res_name)
        return amino_acids[idx]
    except Exception:
        return "X"


# ---------- pLDDT parsing (from B-factor) ----------
def parse_plddt_from_structure_file(structure_file: str):
    """
    Extract residue-wise pLDDT ~ average B-factor per label_seq_id / resseq.
    Works for .cif/.cif.gz (preferred) and .pdb/.pdb.gz.
    Returns: dict {res_seq_number(int) -> mean_Bfactor(float)}
    """
    plddt = {}
    is_gz = structure_file.endswith(".gz")
    is_cif = structure_file.endswith(".cif") or structure_file.endswith(".cif.gz")

    if is_cif:
        try:
            if is_gz:
                with gzip.open(structure_file, "rt") as f:
                    d = MMCIF2Dict(f)
            else:
                d = MMCIF2Dict(structure_file)
            if "_atom_site.B_iso_or_equiv" in d and "_atom_site.label_seq_id" in d:
                b = [float(x) for x in d["_atom_site.B_iso_or_equiv"]]
                s = [int(x) if x != "?" else None for x in d["_atom_site.label_seq_id"]]
                sums = defaultdict(float)
                counts = defaultdict(int)
                for bi, si in zip(b, s):
                    if si is None: 
                        continue
                    sums[si] += bi
                    counts[si] += 1
                for k in sums:
                    plddt[k] = sums[k] / counts[k]
                return plddt
        except Exception:
            # fall back to text scan below
            pass

    # Fallback: scan ATOM records (PDB or CIF-as-text)
    opener = gzip.open if is_gz else open
    with opener(structure_file, "rt") as f:
        for line in f:
            if line.startswith("ATOM"):
                try:
                    res_id = int(line[22:26].strip())
                    b_factor = float(line[60:66].strip())
                    plddt[res_id] = b_factor  # last-atom wins (ok for AF)
                except Exception:
                    continue
    return plddt


# ---------- Secondary structure from mmCIF (optional) ----------
def extract_secondary_structure(structure_file: str):
    """
    Returns: dict {(auth_asym_id, auth_seq_id_int) -> conf_type_id_str}
             e.g., "HELX_RH_PP_P", "STRN", "TURN_TY1_P", etc.
    If absent or non-CIF, returns {}.
    """
    ss = {}
    if not (structure_file.endswith(".cif") or structure_file.endswith(".cif.gz")):
        return ss
    try:
        if structure_file.endswith(".gz"):
            with gzip.open(structure_file, "rt") as f:
                d = MMCIF2Dict(f)
        else:
            d = MMCIF2Dict(structure_file)

        if "_struct_conf.conf_type_id" not in d:
            return ss

        beg_chain = d.get("_struct_conf.beg_auth_asym_id", [])
        beg_res   = d.get("_struct_conf.beg_auth_seq_id", [])
        end_chain = d.get("_struct_conf.end_auth_asym_id", [])
        end_res   = d.get("_struct_conf.end_auth_seq_id", [])
        types     = d.get("_struct_conf.conf_type_id", [])

        n = min(len(beg_chain), len(beg_res), len(end_chain), len(end_res), len(types))
        # Convert residue ids to int when possible
        for i in range(n):
            try:
                rs = int(beg_res[i]); reid = int(end_res[i])
            except Exception:
                continue
            c = beg_chain[i]
            t = types[i]
            for r in range(rs, reid + 1):
                ss[(c, r)] = t
    except Exception:
        pass
    return ss


# ---------- Chain data helpers ----------
def build_chain_data(chain):
    seq = []
    seq_ids = []
    residues = []
    ca_coords = []
    res_index_map = {}  # auth_seq_id -> index in arrays

    for res in chain:
        if not is_aa(res, standard=True):
            continue
        aa = three_to_one(res.resname)
        seq.append(aa)
        seq_ids.append(res.id[1])  # auth seq id
        residues.append(res)
        if "CA" in res:
            ca_coords.append(np.asarray(res["CA"].coord, dtype=float))
        else:
            ca_coords.append(np.array([np.nan, np.nan, np.nan], dtype=float))

    seq_str = "".join(seq)
    for i, rid in enumerate(seq_ids):
        res_index_map[rid] = i

    return {
        "sequence": seq_str,
        "seq_ids": seq_ids,
        "residues": residues,
        "ca_coords": np.vstack(ca_coords) if ca_coords else np.empty((0, 3)),
        "res_index": res_index_map,
    }


# ---------- Hydropathy & charge dictionaries ----------
KYTE_DOOLITTLE = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
    'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
    'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
}

# Charge at ~pH 7 (simplified; histidine treated 0 or +0.1 as you prefer)
CHARGE = {aa: 0 for aa in "ACDEFGHIKLMNPQRSTVWY"}
for aa in "KR": CHARGE[aa] = +1
for aa in "DE": CHARGE[aa] = -1
CHARGE['H'] = +0  # keep neutral in this simplified scheme

BULKY_HYDRO = set("FILVWY")

# ---------- Windows and helpers ----------
def window_indices(center_idx, n, length):
    start = max(0, center_idx - n)
    end   = min(length, center_idx + n + 1)
    return start, end

def shannon_entropy_window(seq_window: str):
    if not seq_window:
        return None
    counts = Counter(seq_window)
    total = len(seq_window)
    probs = [c / total for c in counts.values()]
    # avoid log(0)
    entropy = -sum(p * math.log(p, 2) for p in probs if p > 0)
    return entropy

def ncp_r(seq_window: str):
    if not seq_window:
        return None
    q = sum(CHARGE.get(a, 0) for a in seq_window)
    return q / len(seq_window)

def mean_hydropathy(seq_window: str):
    if not seq_window:
        return None
    vals = [KYTE_DOOLITTLE.get(a, 0) for a in seq_window]
    return float(np.mean(vals)) if len(vals) else None

def kappa_simplified(seq_window: str):
    """
    A simple, bounded [0,1] blockiness score for charge patterning.
    Use σ_i ∈ {+1,-1,0} for K/R (+1), D/E (-1), others 0.
    For all pairs (i<j) with |σ_i|+|σ_j|>0:
      s = σ_i*σ_j ∈ {-1,0,+1}; weight = (j-i).
    κ' = (Σ s*weight) / (Σ |s|*weight) ∈ [-1,1]; map to [0,1] => (κ'+1)/2.
    1.0 = perfectly segregated like-charges; 0.0 = alternating opposite charges.
    """
    if not seq_window:
        return None
    sig = []
    for a in seq_window:
        if a in ("K","R"):
            sig.append(1)
        elif a in ("D","E"):
            sig.append(-1)
        else:
            sig.append(0)
    n = len(sig)
    num = 0.0
    den = 0.0
    for i in range(n):
        si = sig[i]
        if si == 0: 
            continue
        for j in range(i+1, n):
            sj = sig[j]
            if sj == 0:
                continue
            w = (j - i)
            s = si * sj
            num += s * w
            den += abs(s) * w
    if den == 0:
        return None
    return (num / den + 1.0) / 2.0

def radius_of_gyration(coords: np.ndarray):
    """coords: (m,3) array, NaN rows ignored"""
    if coords.size == 0:
        return None
    mask = ~np.isnan(coords).any(axis=1)
    pts = coords[mask]
    if pts.shape[0] < 2:
        return None
    com = pts.mean(axis=0)
    dif = pts - com
    rg2 = np.mean(np.sum(dif * dif, axis=1))
    return float(np.sqrt(rg2))

def compaction_ratios(rg: float, n: int):
    """
    Return two compaction ratios:
      - Rg / (n**(1/3)) [dimensionless normalization]
      - Rg / (2.2 * n**(1/3)) heuristic folded baseline (Å)
    """
    if rg is None or n <= 0:
        return None, None
    base = (n ** (1/3))
    dimless = rg / base if base > 0 else None
    folded = rg / (2.2 * base) if base > 0 else None
    return dimless, folded

def low_conf_contact_fraction(plddt_vec, ca_coords, idx_start, idx_end, thresh=70.0, cutoff=8.0):
    """
    Fraction of CA–CA contacts in [idx_start, idx_end) for which either residue has pLDDT<thresh.
    """
    sub = np.arange(idx_start, idx_end, dtype=int)
    # mask out residues without coords or pLDDT
    ok = []
    for i in sub:
        if i >= len(ca_coords): 
            continue
        if np.isnan(ca_coords[i]).any():
            continue
        ok.append(i)
    ok = np.array(ok, dtype=int)
    if ok.size < 2:
        return None
    pts = ca_coords[ok]
    # Pairwise distances (O(m^2) but m<=~15 window => fine)
    dif = pts[:,None,:] - pts[None,:,:]
    dist = np.sqrt(np.sum(dif*dif, axis=2))
    # consider upper triangle
    m = len(ok)
    contact_pairs = []
    lowconf_pairs = 0
    total_pairs = 0
    for a in range(m):
        for b in range(a+1, m):
            if dist[a,b] <= cutoff:
                total_pairs += 1
                ia = ok[a]; ib = ok[b]
                pa = plddt_vec[ia]; pb = plddt_vec[ib]
                if (pa is not None and pa < thresh) or (pb is not None and pb < thresh):
                    lowconf_pairs += 1
    if total_pairs == 0:
        return None
    return lowconf_pairs / total_pairs

def morf_like_features(seq_window: str, mean_plddt_window: float, plddt_threshold=70.0):
    """
    Simple MoRF-likeness: within a generally disordered window (mean pLDDT < threshold),
    count hydrophobic micro-clusters: any 5-res window with ≥3 of FILVWY.
    Return count and max contiguous length of flagged positions.
    """
    if not seq_window:
        return 0, 0
    n = len(seq_window)
    flags = [False]*n
    for i in range(0, n-4):
        chunk = seq_window[i:i+5]
        if sum(aa in BULKY_HYDRO for aa in chunk) >= 3:
            for k in range(i, i+5):
                flags[k] = True
    # disordered context check
    if mean_plddt_window is not None and mean_plddt_window >= plddt_threshold:
        # if ordered on average, down-weight: return 0 clusters
        return 0, 0
    # count clusters and max run
    count = 0
    maxrun = 0
    run = 0
    in_block = False
    for f in flags:
        if f:
            run += 1
            if not in_block:
                count += 1
                in_block = True
        else:
            maxrun = max(maxrun, run)
            run = 0
            in_block = False
    maxrun = max(maxrun, run)
    return count, maxrun


# ---------- Per-residue feature computation (for a single chain) ----------
def compute_features_for_chain(chain, chain_id, plddt_map, ss_map, uniprot_id, window=7):
    cd = build_chain_data(chain)
    seq = cd["sequence"]
    seq_ids = cd["seq_ids"]
    residues = cd["residues"]
    ca = cd["ca_coords"]
    nres = len(seq)

    # pLDDT vector aligned to seq positions
    plddt_vec = [plddt_map.get(rid, None) for rid in seq_ids]

    # SS to coil/non-coil vector (if available)
    coil = [None]*nres
    if ss_map:
        for i, rid in enumerate(seq_ids):
            t = ss_map.get((chain_id, rid), None)
            if t is None:
                coil[i] = None
            else:
                # treat anything not STRN/HELX as coil-ish
                if "STRN" in t or "SHEET" in t or "HELX" in t:
                    coil[i] = 0
                else:
                    coil[i] = 1

    results = []
    for i, res in enumerate(residues):
        aa = seq[i]
        if aa not in ("S","T","Y"):
            continue  # only STY
        rid = seq_ids[i]
        p_center = plddt_vec[i]

        # window bounds
        s, e = window_indices(i, window, nres)
        seq_win = seq[s:e]
        plddt_win = [p for p in plddt_vec[s:e] if p is not None]
        mean_plddt_win = float(np.mean(plddt_win)) if plddt_win else None

        # (1) pLDDT run-lengths & boundary metrics
        def run_len_at_threshold(thr):
            if p_center is None:
                return None
            # bool vector of low-threshold across full chain
            low = [(p is not None and p < thr) for p in plddt_vec]
            if p_center >= thr:
                # If center is not low, run-length of low at site is 0
                # But we define "run containing site" only if site is low.
                return 0
            # extend from i outwards while low
            L = 1
            j = i-1
            while j >= 0 and low[j]:
                L += 1; j -= 1
            j = i+1
            while j < nres and low[j]:
                L += 1; j += 1
            return L

        def dist_to_boundary(thr):
            if p_center is None:
                return None
            center_is_low = (p_center < thr)
            dist = None
            # search outward until classification flips
            # on either side, in sequence units
            # limit to chain bounds
            for d in range(1, nres):
                left = i-d
                right = i+d
                hit = False
                if left >= 0 and plddt_vec[left] is not None:
                    if (plddt_vec[left] < thr) != center_is_low:
                        dist = d; hit = True
                if right < nres and plddt_vec[right] is not None and not hit:
                    if (plddt_vec[right] < thr) != center_is_low:
                        dist = d; hit = True
                if hit:
                    break
            return dist

        def auc_under(thr):
            # sum max(0, thr - p) over window
            if not plddt_win:
                return None
            diffs = [max(0.0, thr - p) for p in plddt_vec[s:e] if p is not None]
            return float(np.sum(diffs)) if len(diffs) else 0.0

        def total_variation():
            vals = [p for p in plddt_vec[s:e] if p is not None]
            if len(vals) < 2:
                return None
            return float(np.sum(np.abs(np.diff(vals))))

        rl50 = run_len_at_threshold(50.0)
        rl70 = run_len_at_threshold(70.0)
        dOD50 = dist_to_boundary(50.0)
        dOD70 = dist_to_boundary(70.0)
        auc50 = auc_under(50.0)
        auc70 = auc_under(70.0)
        tv7  = total_variation()

        # (2) Uversky charge–hydropathy (window)
        Hbar = mean_hydropathy(seq_win)
        R = ncp_r(seq_win)
        U_flag = None
        if Hbar is not None and R is not None:
            # classic boundary: H < 1.151 R + 0.660 => IDP-like
            U_flag = 1 if (Hbar < 1.151*R + 0.660) else 0

        # (3) Patterning & composition
        NCPR = R
        KAP  = kappa_simplified(seq_win)
        frac_PG = (seq_win.count("P") + seq_win.count("G")) / len(seq_win) if seq_win else None
        frac_bulky = sum(1 for a in seq_win if a in BULKY_HYDRO) / len(seq_win) if seq_win else None
        aa_entropy = shannon_entropy_window(seq_win)
        LCR_flag = None
        if aa_entropy is not None:
            # heuristic: entropy < 2.2 bits => low complexity-like
            LCR_flag = 1 if aa_entropy < 2.2 else 0

        # (4) IDR-ness structural
        rg = radius_of_gyration(ca[s:e]) if (e-s) >= 3 else None
        comp_norm, comp_fold = compaction_ratios(rg, (e - s))
        coil_frac = None
        if coil and any(c is not None for c in coil[s:e]):
            sub = [c for c in coil[s:e] if c is not None]
            coil_frac = float(np.mean(sub)) if len(sub) else None
        lc_contact_frac = low_conf_contact_fraction(plddt_vec, ca, s, e, thresh=70.0, cutoff=8.0)

        # (5) MoRF-likeness
        morf_count, morf_maxlen = morf_like_features(seq_win, mean_plddt_win, plddt_threshold=70.0)

        # Assemble record
        rec = {
            "UniProtID": uniprot_id,
            "ChainID": chain_id,
            "ResidueNumber": rid,
            "ResidueType": aa,
            "Site": f"{uniprot_id}_{rid}",

            # (1)
            "LowPLDDT_RunLen50": rl50,
            "LowPLDDT_RunLen70": rl70,
            "DistTo_OD_Boundary50": dOD50,
            "DistTo_OD_Boundary70": dOD70,
            "pLDDT_AUC50_pm{}".format(window): auc50,
            "pLDDT_AUC70_pm{}".format(window): auc70,
            "pLDDT_TV_pm{}".format(window): tv7,
            "Mean_pLDDT_pm{}".format(window): mean_plddt_win,

            # (2)
            "Uversky_H_pm{}".format(window): Hbar,
            "Uversky_R_pm{}".format(window): R,
            "Uversky_IDPflag_pm{}".format(window): U_flag,

            # (3)
            "NCPR_pm{}".format(window): NCPR,
            "Kappa_pm{}".format(window): KAP,
            "Frac_PG_pm{}".format(window): frac_PG,
            "Frac_BulkyHydro_pm{}".format(window): frac_bulky,
            "AA_Entropy_pm{}".format(window): aa_entropy,
            "LCR_flag_pm{}".format(window): LCR_flag,

            # (4)
            "Rg_CA_pm{}".format(window): rg,
            "CompactionNorm_pm{}".format(window): comp_norm,
            "CompactionRatio_folded_pm{}".format(window): comp_fold,
            "CoilFrac_pm{}".format(window): coil_frac,
            "LowConf_ContactFrac_pm{}".format(window): lc_contact_frac,

            # (5)
            "MoRFlike_ClusterCount_pm{}".format(window): morf_count,
            "MoRFlike_MaxClusterLen_pm{}".format(window): morf_maxlen,
        }
        results.append(rec)
    return results


# ---------- Per-file processing ----------
def process_structure_file(structure_file, file_idx, total_files, temp_dir, window=7):
    filename = os.path.basename(structure_file)
    # UniProt ID from AlphaFold naming "AF-<ID>-F1-model_v*.cif"
    if "AF-" in filename and "-F" in filename:
        uniprot_id = filename.split("AF-")[1].split("-F")[0]
    else:
        uniprot_id = os.path.splitext(filename)[0]

    print(f"[{file_idx}/{total_files}] {filename}")

    # Choose parser
    parser = MMCIFParser(QUIET=True) if (filename.endswith(".cif") or filename.endswith(".cif.gz")) else PDBParser(QUIET=True)

    # Extract SS map (optional)
    ss_map = extract_secondary_structure(structure_file)

    # Parse structure (handle .gz temporarily)
    try:
        if structure_file.endswith(".gz"):
            with gzip.open(structure_file, "rt") as f:
                temp_path = os.path.join(temp_dir, f"__tmp_{file_idx}.cif" if filename.endswith(".cif.gz") else f"__tmp_{file_idx}.pdb")
                with open(temp_path, "w") as out:
                    out.write(f.read())
            structure = parser.get_structure(uniprot_id, temp_path)
            try:
                os.remove(temp_path)
            except Exception:
                pass
        else:
            structure = parser.get_structure(uniprot_id, structure_file)
    except Exception as e:
        print(f"  ! Parse error: {e}")
        return 0

    # pLDDT
    plddt_map = parse_plddt_from_structure_file(structure_file)

    # Process chains
    records = []
    try:
        model = structure[0]
    except Exception:
        model = None

    if model is None:
        return 0

    for chain in model:
        chain_id = chain.id
        try:
            recs = compute_features_for_chain(chain, chain_id, plddt_map, ss_map, uniprot_id, window=window)
            records.extend(recs)
        except Exception as e:
            print(f"  ! Chain {chain_id} error: {e}")
            traceback.print_exc()

    # Write temp CSV
    if records:
        df = pd.DataFrame.from_records(records)
        temp_file = os.path.join(temp_dir, f"temp_{file_idx}.csv")
        df.to_csv(temp_file, index=False)
        return len(records)
    return 0


# ---------- Parallel batch driver ----------
def discover_files(structure_dir, exts=("cif",), fallback_pdb=False, max_files=0):
    # collect by UniProt ID and prefer CIF
    candidates = []
    for fn in os.listdir(structure_dir):
        path = os.path.join(structure_dir, fn)
        if not os.path.isfile(path):
            continue
        # extension filter
        ok = any(fn.endswith(ext) or fn.endswith(ext + ".gz") for ext in exts)
        if fallback_pdb and not ok:
            ok = any(fn.endswith(x) or fn.endswith(x + ".gz") for x in ("pdb",))
        if not ok:
            continue
        # Must look like AlphaFold or at least contain a plausible residue table
        candidates.append(path)

    # If >0, cap
    if max_files and max_files > 0:
        candidates = sorted(candidates)[:max_files]
    else:
        candidates = sorted(candidates)

    return candidates

def main():
    ap = argparse.ArgumentParser(description="Compute STY disorder-centric update features (no merge).")
    ap.add_argument("structure_dir", nargs="?", default="Complete_AF_Proteome", help="Directory of structures (.cif/.pdb/.gz)")
    ap.add_argument("-o", "--output", default="sty_disorder_update.csv", help="Output CSV filename")
    ap.add_argument("-d", "--temp-dir", default="./temp_results", help="Directory for temp CSVs")
    ap.add_argument("-t", "--file-types", default="cif", help="Comma-separated list of file types to include (e.g., cif,pdb)")
    ap.add_argument("-f", "--fallback-pdb", action="store_true", help="If set, also include PDB files when CIF not present")
    ap.add_argument("-p", "--num-processes", type=int, default=0, help="Number of worker processes (0=auto)")
    ap.add_argument("-b", "--batch-size", type=int, default=16, help="Files per batch submitted to the pool")
    ap.add_argument("-m", "--max-files", type=int, default=0, help="Max number of files to process (0=all)")
    ap.add_argument("-w", "--window", type=int, default=7, help="Half-window size for ±N window metrics (default 7)")

    args = ap.parse_args()

    os.makedirs(args.temp_dir, exist_ok=True)

    exts = tuple(x.strip() for x in args.file_types.split(",") if x.strip())
    files = discover_files(args.structure_dir, exts=exts, fallback_pdb=args.fallback_pdb, max_files=args.max_files)
    if not files:
        print("No structure files found. Check --file-types and directory.")
        sys.exit(1)

    total = len(files)
    print(f"Found {total} file(s). Output => {args.output}")

    # Parallel processing
    from concurrent.futures import ProcessPoolExecutor, as_completed
    import multiprocessing as mp

    nproc = args.num_processes if args.num_processes and args.num_processes > 0 else max(1, mp.cpu_count() - 1)
    submitted = 0
    completed = 0
    total_rows = 0

    def submit_batch(executor, batch_files, start_idx):
        futs = []
        for k, path in enumerate(batch_files):
            fut = executor.submit(process_structure_file, path, start_idx + k + 1, total, args.temp_dir, args.window)
            futs.append(fut)
        return futs

    with ProcessPoolExecutor(max_workers=nproc) as ex:
        idx = 0
        while idx < total:
            batch = files[idx : idx + args.batch_size]
            futs = submit_batch(ex, batch, idx)
            for fut in as_completed(futs):
                try:
                    res = fut.result()
                    completed += 1
                    if res:
                        total_rows += res
                except Exception as e:
                    print(f"  ! Worker error: {e}")
            idx += args.batch_size

    # Concatenate temp CSVs
    temp_csvs = [os.path.join(args.temp_dir, f) for f in os.listdir(args.temp_dir) if f.startswith("temp_") and f.endswith(".csv")]
    if not temp_csvs:
        print("No results to write.")
        sys.exit(0)

    dfs = []
    for f in sorted(temp_csvs):
        try:
            dfs.append(pd.read_csv(f))
        except Exception as e:
            print(f"  ! Skipping {f}: {e}")
    if not dfs:
        print("No readable temp CSVs.")
        sys.exit(0)

    out_df = pd.concat(dfs, ignore_index=True)
    out_df.to_csv(args.output, index=False)
    print(f"Wrote {len(out_df):,} rows to {args.output}")

if __name__ == "__main__":
    main()

