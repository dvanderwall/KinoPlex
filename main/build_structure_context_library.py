import pandas as pd
import numpy as np
from Bio.PDB import PDBParser, Selection, Vector, MMCIFParser, MMCIF2Dict
import argparse
import warnings
from math import degrees
import re
import time
import gzip
from Bio.PDB.Polypeptide import three_to_index, is_aa
from scipy.spatial import cKDTree  # Fast neighbor search
import os
from Bio.PDB.HSExposure import HSExposureCA, HSExposureCB
from Bio.PDB.vectors import calc_dihedral, calc_angle
from collections import defaultdict
import traceback
import random
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
import psutil
import signal
import sys
from tqdm import tqdm


# Add this function to replace three_to_one
def three_to_one(res_name):
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    try:
        index = three_to_index(res_name)
        return amino_acids[index]
    except:
        return "X"


# Suppress warnings
warnings.filterwarnings("ignore")

# Define amino acid properties
aa_properties = {
    'A': {'hydrophobic': True, 'polar': False, 'charged': False, 'charge': 0, 'size': 'small', 'vdw_volume': 67,
          'neg_charge': 0, 'pos_charge': 0, 'aromatic': False, 'aliphatic': True, 'sulfur': False, 'hydroxyl': False},
    'R': {'hydrophobic': False, 'polar': True, 'charged': True, 'charge': 1, 'size': 'large', 'vdw_volume': 148,
          'neg_charge': 0, 'pos_charge': 1, 'aromatic': False, 'aliphatic': False, 'sulfur': False, 'hydroxyl': False},
    'N': {'hydrophobic': False, 'polar': True, 'charged': False, 'charge': 0, 'size': 'medium', 'vdw_volume': 96,
          'neg_charge': 0, 'pos_charge': 0, 'aromatic': False, 'aliphatic': False, 'sulfur': False, 'hydroxyl': False},
    'D': {'hydrophobic': False, 'polar': True, 'charged': True, 'charge': -1, 'size': 'medium', 'vdw_volume': 91,
          'neg_charge': 1, 'pos_charge': 0, 'aromatic': False, 'aliphatic': False, 'sulfur': False, 'hydroxyl': False},
    'C': {'hydrophobic': True, 'polar': False, 'charged': False, 'charge': 0, 'size': 'medium', 'vdw_volume': 86,
          'neg_charge': 0, 'pos_charge': 0, 'aromatic': False, 'aliphatic': False, 'sulfur': True, 'hydroxyl': False},
    'Q': {'hydrophobic': False, 'polar': True, 'charged': False, 'charge': 0, 'size': 'medium', 'vdw_volume': 114,
          'neg_charge': 0, 'pos_charge': 0, 'aromatic': False, 'aliphatic': False, 'sulfur': False, 'hydroxyl': False},
    'E': {'hydrophobic': False, 'polar': True, 'charged': True, 'charge': -1, 'size': 'medium', 'vdw_volume': 109,
          'neg_charge': 1, 'pos_charge': 0, 'aromatic': False, 'aliphatic': False, 'sulfur': False, 'hydroxyl': False},
    'G': {'hydrophobic': True, 'polar': False, 'charged': False, 'charge': 0, 'size': 'small', 'vdw_volume': 48,
          'neg_charge': 0, 'pos_charge': 0, 'aromatic': False, 'aliphatic': False, 'sulfur': False, 'hydroxyl': False},
    'H': {'hydrophobic': False, 'polar': True, 'charged': True, 'charge': 0.5, 'size': 'medium', 'vdw_volume': 118,
          'neg_charge': 0, 'pos_charge': 0.5, 'aromatic': True, 'aliphatic': False, 'sulfur': False, 'hydroxyl': False},
    'I': {'hydrophobic': True, 'polar': False, 'charged': False, 'charge': 0, 'size': 'large', 'vdw_volume': 124,
          'neg_charge': 0, 'pos_charge': 0, 'aromatic': False, 'aliphatic': True, 'sulfur': False, 'hydroxyl': False},
    'L': {'hydrophobic': True, 'polar': False, 'charged': False, 'charge': 0, 'size': 'large', 'vdw_volume': 124,
          'neg_charge': 0, 'pos_charge': 0, 'aromatic': False, 'aliphatic': True, 'sulfur': False, 'hydroxyl': False},
    'K': {'hydrophobic': False, 'polar': True, 'charged': True, 'charge': 1, 'size': 'large', 'vdw_volume': 135,
          'neg_charge': 0, 'pos_charge': 1, 'aromatic': False, 'aliphatic': False, 'sulfur': False, 'hydroxyl': False},
    'M': {'hydrophobic': True, 'polar': False, 'charged': False, 'charge': 0, 'size': 'large', 'vdw_volume': 124,
          'neg_charge': 0, 'pos_charge': 0, 'aromatic': False, 'aliphatic': False, 'sulfur': True, 'hydroxyl': False},
    'F': {'hydrophobic': True, 'polar': False, 'charged': False, 'charge': 0, 'size': 'large', 'vdw_volume': 135,
          'neg_charge': 0, 'pos_charge': 0, 'aromatic': True, 'aliphatic': False, 'sulfur': False, 'hydroxyl': False},
    'P': {'hydrophobic': True, 'polar': False, 'charged': False, 'charge': 0, 'size': 'medium', 'vdw_volume': 90,
          'neg_charge': 0, 'pos_charge': 0, 'aromatic': False, 'aliphatic': False, 'sulfur': False, 'hydroxyl': False},
    'S': {'hydrophobic': False, 'polar': True, 'charged': False, 'charge': 0, 'size': 'small', 'vdw_volume': 73,
          'neg_charge': 0, 'pos_charge': 0, 'aromatic': False, 'aliphatic': False, 'sulfur': False, 'hydroxyl': True},
    'T': {'hydrophobic': False, 'polar': True, 'charged': False, 'charge': 0, 'size': 'medium', 'vdw_volume': 93,
          'neg_charge': 0, 'pos_charge': 0, 'aromatic': False, 'aliphatic': False, 'sulfur': False, 'hydroxyl': True},
    'W': {'hydrophobic': True, 'polar': True, 'charged': False, 'charge': 0, 'size': 'large', 'vdw_volume': 163,
          'neg_charge': 0, 'pos_charge': 0, 'aromatic': True, 'aliphatic': False, 'sulfur': False, 'hydroxyl': False},
    'Y': {'hydrophobic': False, 'polar': True, 'charged': False, 'charge': 0, 'size': 'large', 'vdw_volume': 141,
          'neg_charge': 0, 'pos_charge': 0, 'aromatic': True, 'aliphatic': False, 'sulfur': False, 'hydroxyl': True},
    'V': {'hydrophobic': True, 'polar': False, 'charged': False, 'charge': 0, 'size': 'medium', 'vdw_volume': 105,
          'neg_charge': 0, 'pos_charge': 0, 'aromatic': False, 'aliphatic': True, 'sulfur': False, 'hydroxyl': False},
}

# Define standard VDW radii for atoms
VDW_RADII = {
    'C': 1.7, 'N': 1.55, 'O': 1.52, 'S': 1.8,
    'P': 1.8, 'H': 1.2, 'F': 1.47, 'CL': 1.75,
    'BR': 1.85, 'I': 1.98, 'CA': 2.0, 'ZN': 1.39,
    'FE': 1.32, 'MG': 1.73
}
DEFAULT_RADIUS = 1.7  # Default for unknown atoms
WATER_RADIUS = 1.4  # Radius of water probe

# Precompute some common values
STY_RESIDUES = set(['SER', 'THR', 'TYR'])
STY_SINGLE_LETTER = set(['S', 'T', 'Y'])
MAX_DISTANCE_SQUARED = 16  # 4Å squared


def parse_plddt_from_structure_file(structure_file):
    """Extract pLDDT scores from AlphaFold structure files (B-factor column)
    Works with both .pdb, .cif and their gzipped versions"""
    plddt_dict = {}

    # Check if file is gzipped
    is_gzipped = structure_file.endswith('.gz')

    # Check if file is CIF
    is_cif = structure_file.endswith('.cif') or structure_file.endswith('.cif.gz')

    # For CIF files, try to extract B-factors using MMCIF2Dict
    if is_cif:
        try:
            if is_gzipped:
                with gzip.open(structure_file, 'rt') as f:
                    mmcif_dict = MMCIF2Dict.MMCIF2Dict(f)
            else:
                mmcif_dict = MMCIF2Dict.MMCIF2Dict(structure_file)

            # Check if the necessary data is in the CIF dictionary
            if '_atom_site.B_iso_or_equiv' in mmcif_dict and '_atom_site.label_seq_id' in mmcif_dict:
                b_factors = mmcif_dict['_atom_site.B_iso_or_equiv']
                seq_ids = mmcif_dict['_atom_site.label_seq_id']

                # Convert values to appropriate types
                b_factors = [float(b) for b in b_factors]
                seq_ids = [int(id) if id != '?' else None for id in seq_ids]

                # Create a dictionary mapping residue IDs to average B-factors
                b_factor_sums = {}
                b_factor_counts = {}

                for i in range(len(seq_ids)):
                    if seq_ids[i] is not None:
                        if seq_ids[i] not in b_factor_sums:
                            b_factor_sums[seq_ids[i]] = 0
                            b_factor_counts[seq_ids[i]] = 0
                        b_factor_sums[seq_ids[i]] += b_factors[i]
                        b_factor_counts[seq_ids[i]] += 1

                # Calculate average B-factor for each residue
                for res_id in b_factor_sums:
                    plddt_dict[res_id] = b_factor_sums[res_id] / b_factor_counts[res_id]

                return plddt_dict
        except Exception as e:
            print(f"  Warning: Could not extract B-factors from CIF dictionary: {e}")
            print("  Falling back to text parsing...")

    # Fall back to text parsing for both PDB and CIF
    # Open file accordingly
    if is_gzipped:
        opener = gzip.open(structure_file, 'rt')
    else:
        opener = open(structure_file, 'r')

    with opener as f:
        for line in f:
            if line.startswith('ATOM'):
                res_id = int(line[22:26].strip())
                b_factor = float(line[60:66].strip())
                plddt_dict[res_id] = b_factor

    return plddt_dict


def find_hydrogen_bonds_vectorized(hydroxyl_atom, neighbor_residues, max_distance=3.5):
    """Vectorized calculation of hydrogen bonds"""
    if not hydroxyl_atom:
        return 0

    max_distance_squared = max_distance * max_distance

    # Extract all oxygen and nitrogen atoms from neighbor residues
    potential_partners = []

    for res in neighbor_residues:
        # Skip the residue itself
        if res is hydroxyl_atom.get_parent():
            continue

        for atom in res:
            # Only consider oxygen and nitrogen atoms as potential partners
            if atom.element in ['O', 'N']:
                potential_partners.append(atom.coord)

    if not potential_partners:
        return 0

    # Convert to numpy arrays
    potential_partners = np.array(potential_partners)
    hydroxyl_coord = np.array(hydroxyl_atom.coord)

    # Calculate squared distances
    vectors = potential_partners - hydroxyl_coord
    distances_sq = np.sum(vectors * vectors, axis=1)

    # Count potential hydrogen bonds
    h_bonds = np.sum(distances_sq <= max_distance_squared)

    return h_bonds


def calculate_hydroxyl_exposure_vectorized(hydroxyl_atom, neighbor_residues):
    """Vectorized calculation of hydroxyl exposure"""
    if not hydroxyl_atom or not neighbor_residues:
        return None, 0

    # Extract all atom coordinates from neighbor residues
    neighbor_atoms = []
    backbone_atoms = []

    for res in neighbor_residues:
        for atom in res:
            # Skip if this is from the same residue as hydroxyl_atom
            if atom.get_parent() is hydroxyl_atom.get_parent():
                continue

            neighbor_atoms.append(atom.coord)

            # Check if it's a backbone atom
            if atom.name in ['N', 'CA', 'C', 'O']:
                backbone_atoms.append(atom.coord)

    if not neighbor_atoms:
        return 1.0, 0  # Fully exposed if no neighbors

    # Convert to numpy arrays
    neighbor_atoms = np.array(neighbor_atoms)
    hydroxyl_coord = np.array(hydroxyl_atom.coord)

    # Calculate squared distances to hydroxyl
    # Use broadcasting to calculate all distances at once
    vectors = neighbor_atoms - hydroxyl_coord
    distances_sq = np.sum(vectors * vectors, axis=1)

    # Count atoms within cutoff distance
    nearby_atoms = np.sum(distances_sq <= MAX_DISTANCE_SQUARED)

    # Calculate backbone contacts if we have backbone atoms
    if backbone_atoms:
        backbone_atoms = np.array(backbone_atoms)
        backbone_vectors = backbone_atoms - hydroxyl_coord
        backbone_distances_sq = np.sum(backbone_vectors * backbone_vectors, axis=1)
        nearby_backbone = np.sum(backbone_distances_sq <= MAX_DISTANCE_SQUARED)
    else:
        nearby_backbone = 0

    # Calculate exposure
    exposure = 1.0 / (1.0 + nearby_atoms / 10.0)

    return exposure, nearby_backbone


def find_neighbors_fast(chain_data, residue, radius=10.0):
    """Use KD-tree for fast neighbor finding"""
    # Get CA atom coordinates for the central residue
    try:
        ca_atom = residue['CA']
    except KeyError:
        return []

    ca_coord = ca_atom.coord

    # Query KD-tree for neighbors
    if chain_data['kd_tree'] is None:
        return []

    # Find all atoms within radius
    neighbor_indices = chain_data['kd_tree'].query_ball_point(ca_coord, radius)

    # Get unique residues
    neighbor_residues = set()
    for idx in neighbor_indices:
        if idx < len(chain_data['atom_to_residue']):
            neighbor_residues.add(chain_data['atom_to_residue'][idx])

    # Remove the center residue from neighbors
    if residue in neighbor_residues:
        neighbor_residues.remove(residue)

    return list(neighbor_residues)


def build_chain_data(chain):
    """Pre-compute and cache chain-level data for reuse"""
    chain_data = {}

    # Build complete chain sequence once
    chain_seq = ''
    seq_ids = []
    residue_map = {}  # Map residue ID to position in chain

    for i, res in enumerate(chain):
        if is_aa(res):
            chain_seq += three_to_one(res.resname)
            seq_ids.append(res.id[1])
            residue_map[res.id[1]] = i

    chain_data['sequence'] = chain_seq
    chain_data['seq_ids'] = seq_ids
    chain_data['residue_map'] = residue_map

    # Collect atom coordinates for KD-tree
    atom_coords = []
    atom_to_residue = []  # Map atom index to residue

    # Also collect STY residues
    sty_residues = []

    for res in chain:
        if is_aa(res):
            if res.resname in STY_RESIDUES:
                sty_residues.append(res)

            for atom in res:
                atom_coords.append(atom.coord)
                atom_to_residue.append(res)

    chain_data['atom_coords'] = np.array(atom_coords)
    chain_data['atom_to_residue'] = atom_to_residue
    chain_data['sty_residues'] = sty_residues

    # Build KD-tree for fast neighbor lookup if we have atoms
    if atom_coords:
        chain_data['kd_tree'] = cKDTree(chain_data['atom_coords'])
    else:
        chain_data['kd_tree'] = None

    return chain_data


def calc_dihedral(v1, v2, v3, v4):
    """Calculate dihedral angle between 4 vectors"""
    try:
        ab = v1 - v2
        cb = v3 - v2
        db = v4 - v3

        u = cb.cross(ab)
        v = db.cross(cb)

        u.normalize()
        v.normalize()
        w = u.cross(v)

        angle = np.arctan2(w * cb, u * v)
        return angle
    except:
        return None


def get_sequence_motif(chain, residue, window=7):
    """Get the -7:+7 sequence motif around a residue"""
    res_id = residue.id[1]
    chain_seq = ''
    seq_ids = []

    # Build ordered sequence and track residue IDs
    for res in chain:
        if is_aa(res):
            chain_seq += three_to_one(res.resname)
            seq_ids.append(res.id[1])

    # Find the position in the sequence
    if res_id not in seq_ids:
        return "X" * (2 * window + 1)

    seq_pos = seq_ids.index(res_id)

    # Extract motif
    start = max(0, seq_pos - window)
    end = min(len(chain_seq), seq_pos + window + 1)

    # Pad with X if needed
    prefix = "X" * max(0, window - seq_pos)
    suffix = "X" * max(0, window - (len(chain_seq) - seq_pos - 1))

    motif = prefix + chain_seq[start:end] + suffix

    # Ensure motif is the correct length
    if len(motif) > 2 * window + 1:
        motif = motif[:2 * window + 1]
    elif len(motif) < 2 * window + 1:
        motif = motif + "X" * (2 * window + 1 - len(motif))

    return motif


def analyze_pocket_properties(neighbors):
    """Analyze biochemical properties of residues in the pocket"""
    if not neighbors:
        return {
            'hydrophobic_count': 0,
            'polar_count': 0,
            'positive_count': 0,
            'negative_count': 0,
            'net_charge': 0,
            'small_count': 0,
            'medium_count': 0,
            'large_count': 0,
            'hydrophobic_ratio': 0,
            'charged_ratio': 0,
            'aromatic_count': 0,
            'aliphatic_count': 0,
            'sulfur_count': 0
        }

    pocket_aa = [three_to_one(res.resname) for res in neighbors if is_aa(res)]
    total_count = len(pocket_aa)

    if total_count == 0:
        return {
            'hydrophobic_count': 0,
            'polar_count': 0,
            'positive_count': 0,
            'negative_count': 0,
            'net_charge': 0,
            'small_count': 0,
            'medium_count': 0,
            'large_count': 0,
            'hydrophobic_ratio': 0,
            'charged_ratio': 0,
            'aromatic_count': 0,
            'aliphatic_count': 0,
            'sulfur_count': 0
        }

    hydrophobic_count = sum(1 for aa in pocket_aa if aa in aa_properties and aa_properties[aa]['hydrophobic'])
    polar_count = sum(1 for aa in pocket_aa if aa in aa_properties and aa_properties[aa]['polar'])
    positive_count = sum(1 for aa in pocket_aa if
                         aa in aa_properties and aa_properties[aa]['charged'] and aa_properties[aa]['charge'] > 0)
    negative_count = sum(1 for aa in pocket_aa if
                         aa in aa_properties and aa_properties[aa]['charged'] and aa_properties[aa]['charge'] < 0)
    aromatic_count = sum(1 for aa in pocket_aa if aa in aa_properties and aa_properties[aa]['aromatic'])
    aliphatic_count = sum(1 for aa in pocket_aa if aa in aa_properties and aa_properties[aa]['aliphatic'])
    sulfur_count = sum(1 for aa in pocket_aa if aa in aa_properties and aa_properties[aa]['sulfur'])

    properties = {
        'hydrophobic_count': hydrophobic_count,
        'polar_count': polar_count,
        'positive_count': positive_count,
        'negative_count': negative_count,
        'net_charge': sum(aa_properties[aa]['charge'] for aa in pocket_aa if aa in aa_properties),
        'small_count': sum(1 for aa in pocket_aa if aa in aa_properties and aa_properties[aa]['size'] == 'small'),
        'medium_count': sum(1 for aa in pocket_aa if aa in aa_properties and aa_properties[aa]['size'] == 'medium'),
        'large_count': sum(1 for aa in pocket_aa if aa in aa_properties and aa_properties[aa]['size'] == 'large'),
        'hydrophobic_ratio': hydrophobic_count / total_count if total_count > 0 else 0,
        'charged_ratio': (positive_count + negative_count) / total_count if total_count > 0 else 0,
        'aromatic_count': aromatic_count,
        'aliphatic_count': aliphatic_count,
        'sulfur_count': sulfur_count
    }

    return properties


def identify_domain_and_motifs(sequence, uniprot_id, res_pos):
    """Identify known domains and motifs around the STY site"""
    # Pre-compile regex patterns for better performance
    motifs = {
        # Generic kinase motifs
        'PKA_consensus': re.compile(r'R[RK].S'),
        'PKC_consensus': re.compile(r'[ST].[RK]'),
        'CK2_consensus': re.compile(r'S..E'),
        'CK1_consensus': re.compile(r'S..[ST]'),
        'proline_directed': re.compile(r'[ST]P'),
        'GSK3_consensus': re.compile(r'S...S'),
        'acidic_directed': re.compile(r'S.[DE]'),
        'basic_directed': re.compile(r'S.[RK]'),
        'CDK_consensus': re.compile(r'[ST]P.[RK]'),
        'MAPK_consensus': re.compile(r'P.[ST]P'),
        'ATM_ATR_consensus': re.compile(r'[ST]Q'),
    }

    found_motifs = []
    window_size = 10
    start = max(0, res_pos - window_size)
    end = min(len(sequence), res_pos + window_size + 1)

    window_seq = sequence[start:end]
    centered_pos = min(res_pos, window_size)

    for motif_name, pattern in motifs.items():
        matches = pattern.finditer(window_seq)
        for match in matches:
            # Check if the match includes our STY site
            match_start = match.start()
            match_end = match.end()
            if match_start <= centered_pos < match_end:
                found_motifs.append(motif_name)

    return ', '.join(found_motifs) if found_motifs else 'None'


def extract_secondary_structure(structure_file):
    """Extract secondary structure information from a CIF file with optimized performance"""
    ss_dict = {}  # Map of residue_id -> secondary structure type

    # Early return for non-CIF files (unchanged, this is already efficient)
    if not (structure_file.endswith('.cif') or structure_file.endswith('.cif.gz')):
        return ss_dict

    try:
        # Handle gzipped files (no change needed here, this is standard practice)
        if structure_file.endswith('.gz'):
            with gzip.open(structure_file, 'rt') as f:
                mmcif_dict = MMCIF2Dict.MMCIF2Dict(f)
        else:
            mmcif_dict = MMCIF2Dict.MMCIF2Dict(structure_file)

        # Check if struct_conf is in the CIF file
        if '_struct_conf.conf_type_id' not in mmcif_dict:
            return ss_dict  # Early return if no secondary structure data

        # Get the secondary structure data
        # Use dict.get() with default empty list to avoid potential KeyErrors
        begin_chain_ids = mmcif_dict.get('_struct_conf.beg_auth_asym_id', [])
        begin_res_ids = mmcif_dict.get('_struct_conf.beg_auth_seq_id', [])
        end_chain_ids = mmcif_dict.get('_struct_conf.end_auth_asym_id', [])
        end_res_ids = mmcif_dict.get('_struct_conf.end_auth_seq_id', [])
        ss_types = mmcif_dict.get('_struct_conf.conf_type_id', [])

        # Skip if any required data is missing
        if not (begin_chain_ids and begin_res_ids and end_chain_ids and end_res_ids and ss_types):
            return ss_dict

        # Pre-calculate length for performance
        num_entries = len(ss_types)

        # Pre-convert all residue IDs to integers at once rather than repeatedly
        try:
            begin_res_ids = [int(x) for x in begin_res_ids]
            end_res_ids = [int(x) for x in end_res_ids]
        except ValueError:
            # Handle potential non-integer values
            begin_res_ids_clean = []
            end_res_ids_clean = []
            for i in range(num_entries):
                try:
                    begin_res_ids_clean.append(int(begin_res_ids[i]))
                    end_res_ids_clean.append(int(end_res_ids[i]))
                except ValueError:
                    # Skip entries with non-integer residue IDs
                    continue
            begin_res_ids = begin_res_ids_clean
            end_res_ids = end_res_ids_clean

            # Update num_entries in case we filtered out some values
            num_entries = min(len(begin_res_ids), len(end_res_ids), len(ss_types))

        # Map each residue to its secondary structure
        # Use zip to iterate over all lists simultaneously
        for i, (chain_id, start_res, end_res, ss_type) in enumerate(
                zip(begin_chain_ids[:num_entries], begin_res_ids[:num_entries],
                    end_res_ids[:num_entries], ss_types[:num_entries])):
            # Add each residue in the range to the dictionary
            # Use a dictionary comprehension for better performance
            ss_dict.update({(chain_id, res_id): ss_type for res_id in range(start_res, end_res + 1)})

    except Exception as e:
        print(f"Error extracting secondary structure: {e}")

    return ss_dict


def calculate_hse_fast(model, residues, chain_data, radius=13.0):
    """
    Fast HSE calculation that only processes specific residues using the existing
    find_neighbors_fast function for better performance.

    Args:
        model: BioPython Model object
        residues: List of Bio.PDB Residue objects to calculate HSE for
        chain_data: Pre-computed chain data with KD-tree
        radius: Sphere radius in Å (default: 13.0)

    Returns:
        dict: Mapping from residue IDs to HSE dictionaries
    """
    results = {}

    # Process each requested residue
    for residue in residues:
        res_id = residue.id[1]

        # Skip non-amino acids
        if not is_aa(residue):
            results[res_id] = {}
            continue

        # Get CA and CB atoms (if they exist)
        ca_atom = residue['CA'] if 'CA' in residue else None
        cb_atom = residue['CB'] if 'CB' in residue else None

        # Skip residues without required atoms
        if not ca_atom:
            results[res_id] = {}
            continue

        # Find neighbors using our fast function
        neighbors = find_neighbors_fast(chain_data, residue, radius)

        # Initialize HSE data
        hse = {}

        # Define up and down regions based on CA-CB vector (or pseudo-CB)
        if cb_atom:
            # Use actual CB atom if available
            ca_cb_vector = cb_atom.coord - ca_atom.coord
        else:
            # Calculate pseudo-CB vector based on backbone
            n_atom = residue['N'] if 'N' in residue else None
            c_atom = residue['C'] if 'C' in residue else None

            if n_atom and c_atom:
                # Pseudo-CB direction for Glycine (perpendicular to peptide plane)
                ca_n_vector = n_atom.coord - ca_atom.coord
                ca_c_vector = c_atom.coord - ca_atom.coord
                ca_cb_vector = np.cross(ca_n_vector, ca_c_vector)
                # Normalize vector
                norm = np.sqrt(np.sum(ca_cb_vector * ca_cb_vector))
                if norm > 0:
                    ca_cb_vector = ca_cb_vector / norm
                else:
                    ca_cb_vector = np.array([0, 0, 0])
            else:
                ca_cb_vector = np.array([0, 0, 0])

        # Count neighbors in up and down half-spheres for CA-based HSE
        up_count = 0
        down_count = 0

        for neighbor_res in neighbors:
            # Skip if it's the same residue
            if neighbor_res is residue:
                continue

            # Count contributions from each atom in the neighbor residue
            for atom in neighbor_res:
                # Calculate vector from CA to neighbor atom
                ca_neighbor_vector = atom.coord - ca_atom.coord

                # Check if this atom is within the radius
                dist_sq = np.sum(ca_neighbor_vector * ca_neighbor_vector)
                if dist_sq <= radius * radius:
                    # Dot product determines if neighbor is in up or down half-sphere
                    if np.dot(ca_cb_vector, ca_neighbor_vector) >= 0:
                        up_count += 1
                    else:
                        down_count += 1

        # Store CA-based HSE results
        hse["HSE_CA_U"] = up_count
        hse["HSE_CA_D"] = down_count
        hse["HSE_CA_RATIO"] = up_count / down_count if down_count > 0 else None

        # Calculate CB-based HSE if CB exists
        if cb_atom:
            cb_up_count = 0
            cb_down_count = 0

            for neighbor_res in neighbors:
                # Skip if it's the same residue
                if neighbor_res is residue:
                    continue

                # Count contributions from each atom in the neighbor residue
                for atom in neighbor_res:
                    # Calculate vector from CB to neighbor atom
                    cb_neighbor_vector = atom.coord - cb_atom.coord

                    # Check if this atom is within the radius
                    dist_sq = np.sum(cb_neighbor_vector * cb_neighbor_vector)
                    if dist_sq <= radius * radius:
                        # Dot product determines if neighbor is in up or down half-sphere
                        if np.dot(ca_cb_vector, cb_neighbor_vector) >= 0:
                            cb_up_count += 1
                        else:
                            cb_down_count += 1

            # Store CB-based HSE results
            hse["HSE_CB_U"] = cb_up_count
            hse["HSE_CB_D"] = cb_down_count
            hse["HSE_CB_RATIO"] = cb_up_count / cb_down_count if cb_down_count > 0 else None

        # Store results for this residue
        results[res_id] = hse

    return results


def cache_all_secondary_structures(model, ss_dict):
    """Cache all secondary structure residues in the model for faster lookup

    Args:
        model: BioPython Model object
        ss_dict: Dictionary mapping (chain_id, res_id) to secondary structure

    Returns:
        dict: Dictionary with keys for each secondary structure type, each containing:
            - 'residues': List of (chain_id, res_id) tuples
            - 'ca_coords': Dictionary mapping (chain_id, res_id) to CA coordinates
    """
    if not ss_dict:
        return {}

    # Define secondary structure categories
    ss_categories = {
        'helix': ['HELX_LH_PP_P', 'HELX_RH_PP_P', 'HELX_RH_3T_P', 'HELX_RH_PI_P', 'HELX_LH_P'],
        'strand': ['STRN', 'SHEET'],
        'turn': ['TURN_TY1_P', 'TURN_TY2_P', 'TURN_TY1_PM', 'TURN_TY2_PM', 'TURN_TY3_P', 'BEND'],
        'coil': ['COIL']
    }

    # Initialize result dictionary
    ss_data = {ss_type: {'residues': [], 'ca_coords': {}} for ss_type in ss_categories}

    # Create a mapping from specific SS codes to category
    ss_type_to_category = {}
    for category, ss_types in ss_categories.items():
        for ss_type in ss_types:
            ss_type_to_category[ss_type] = category

    # Process all residues in the SS dictionary
    for (c_id, r_id), ss_type in ss_dict.items():
        # Determine the secondary structure category
        category = None
        for cat, patterns in ss_categories.items():
            if any(pattern in ss_type for pattern in patterns):
                category = cat
                break

        # If no category matched, use the raw SS type
        if category is None:
            category = 'other'
            if 'other' not in ss_data:
                ss_data['other'] = {'residues': [], 'ca_coords': {}}

        # Add the residue to the appropriate category
        ss_data[category]['residues'].append((c_id, r_id))

        # Try to get CA coordinates
        if c_id in model:
            chain = model[c_id]
            try:
                # BioPython residue selector
                res = chain[(' ', r_id, ' ')]
                if 'CA' in res:
                    ss_data[category]['ca_coords'][(c_id, r_id)] = res['CA'].coord
            except (KeyError, Exception):
                continue

    return ss_data


def calculate_secondary_structure_distances(residue, seq_ids, ss_data):
    """Vectorized calculation of distances to nearest secondary structures of each type

    Args:
        residue: The target residue
        seq_ids: List of residue IDs in sequence order
        ss_data: Dictionary of secondary structure data from cache_all_secondary_structures

    Returns:
        dict: Dictionary mapping secondary structure types to tuples of
              (sequence distance, spatial distance in Ångstroms)
    """
    # Get residue info
    chain_id = residue.get_parent().id
    res_id = residue.id[1]

    # Get CA atom of the residue
    try:
        ca_atom = residue['CA']
        ca_coord = np.array(ca_atom.coord)
    except KeyError:
        return {ss_type: (None, None) for ss_type in ss_data}

    # Get target sequence position
    target_seq_pos = None
    if res_id in seq_ids:
        target_seq_pos = seq_ids.index(res_id)

    # Initialize results dictionary
    distances = {}

    # Calculate distances for each secondary structure type
    for ss_type, data in ss_data.items():
        # Get residues and coordinates for this SS type
        ss_residues = data['residues']
        ss_ca_coords = data['ca_coords']

        if not ss_residues:
            distances[ss_type] = (None, None)
            continue

        # Initialize distances
        min_seq_distance = None
        min_spatial_distance = None

        # Calculate sequence distance (if in same chain)
        if target_seq_pos is not None:
            # Find SS elements in the same chain
            same_chain_ss = [(c, r) for (c, r) in ss_residues if c == chain_id]

            if same_chain_ss:
                # Extract residue IDs and find those that are in seq_ids
                ss_res_ids = [r for (_, r) in same_chain_ss]
                valid_ss_res_ids = [r for r in ss_res_ids if r in seq_ids]

                if valid_ss_res_ids:
                    # Calculate sequence distances vectorized
                    ss_seq_pos = [seq_ids.index(r) for r in valid_ss_res_ids]
                    seq_distances = np.abs(np.array(ss_seq_pos) - target_seq_pos)
                    min_seq_distance = np.min(seq_distances)

        # Calculate spatial distance
        all_ss_coords = []

        for key in ss_ca_coords:
            if key != (chain_id, res_id):  # Skip the residue itself
                all_ss_coords.append(ss_ca_coords[key])

        if all_ss_coords:
            # Convert to numpy array for vectorized operations
            all_ss_coords = np.array(all_ss_coords)

            # Calculate distances vectorized
            vectors = all_ss_coords - ca_coord
            distances_sq = np.sum(vectors * vectors, axis=1)
            min_spatial_distance = np.sqrt(np.min(distances_sq))

        distances[ss_type] = (min_seq_distance, min_spatial_distance)

    return distances


# ============= NEW FEATURE IMPLEMENTATIONS =============
def calculate_sasa_vectorized(residue, neighbors, probe_radius=1.4, n_points=100):
    """Calculate approximate SASA using a vectorized ray-casting approach

    Args:
        residue: The residue to calculate SASA for
        neighbors: List of neighboring residues
        probe_radius: Radius of water probe (default 1.4 Å)
        n_points: Number of points to use on unit sphere (higher = more accurate, slower)

    Returns:
        total_sasa: Approximate SASA in Å²
        hydroxyl_sasa: SASA specific to hydroxyl oxygen (for STY residues)
    """

    # Get VDW radius for atoms
    def get_vdw_radius(atom):
        element = atom.element.upper()
        return VDW_RADII.get(element, DEFAULT_RADIUS)

    # Generate points on unit sphere (evenly distributed)
    # Using Fibonacci sphere algorithm for uniform distribution
    phi = (1 + np.sqrt(5)) / 2  # Golden ratio
    indices = np.arange(n_points)
    y = 1 - (indices / (n_points - 1)) * 2  # y goes from 1 to -1
    radius = np.sqrt(1 - y * y)  # Radius at each height

    theta = 2 * np.pi * indices / phi  # Golden angle increment

    x = np.cos(theta) * radius
    z = np.sin(theta) * radius

    # Combine into unit vectors
    unit_sphere = np.column_stack((x, y, z))

    # Get hydroxyl oxygen for STY residues
    hydroxyl_atom = None
    if residue.resname == 'SER' and 'OG' in residue:
        hydroxyl_atom = residue['OG']
    elif residue.resname == 'THR' and 'OG1' in residue:
        hydroxyl_atom = residue['OG1']
    elif residue.resname == 'TYR' and 'OH' in residue:
        hydroxyl_atom = residue['OH']

    # Extract neighbor atom positions and radii
    neighbor_positions = []
    neighbor_radii = []

    # Add atoms from neighboring residues
    for res in neighbors:
        for atom in res:
            # Skip if this is from the central residue
            if atom.get_parent() is residue:
                continue

            neighbor_positions.append(atom.coord)
            neighbor_radii.append(get_vdw_radius(atom) + probe_radius)

    # Convert to numpy arrays
    neighbor_positions = np.array(neighbor_positions) if neighbor_positions else np.empty((0, 3))
    neighbor_radii = np.array(neighbor_radii) if neighbor_radii else np.empty(0)

    # Process each atom in the residue
    total_area = 0
    hydroxyl_area = 0

    for atom in residue:
        # Get VDW radius
        atom_radius = get_vdw_radius(atom)
        atom_coord = atom.coord

        # Expand sphere points to actual position
        expanded_radius = atom_radius + probe_radius
        test_points = unit_sphere * expanded_radius + atom_coord

        # Area of each point on expanded surface
        point_area = 4 * np.pi * expanded_radius ** 2 / n_points

        # Test each ray for overlap with neighbors
        if len(neighbor_positions) > 0:
            # For each test point, compute distance to all neighbor centers
            # Shape: (n_points, n_neighbors)
            diff = test_points[:, np.newaxis, :] - neighbor_positions[np.newaxis, :, :]
            sq_dists = np.sum(diff ** 2, axis=2)

            # Compare squared distances to squared neighbor radii
            # If any squared distance is less than squared radius, point is blocked
            sq_radii = neighbor_radii ** 2
            blocked_mask = np.any(sq_dists <= sq_radii, axis=1)

            # Count accessible points
            accessible_points = n_points - np.sum(blocked_mask)
            accessible_area = accessible_points * point_area
        else:
            # No neighbors, all points are accessible
            accessible_area = 4 * np.pi * expanded_radius ** 2

        total_area += accessible_area

        # Track hydroxyl SASA separately
        if atom is hydroxyl_atom:
            hydroxyl_area = accessible_area

    return total_area, hydroxyl_area


def calculate_residue_depth_multi(residue, model, thresholds=[15], cutoffs=[10.0]):
    """Calculate minimum distance from residue atoms to protein surface with multiple parameters

    Args:
        residue: The target residue
        model: BioPython Model object
        thresholds: List of surface thresholds to try (default: [10, 15, 20])
        cutoffs: List of surface cutoffs to try (default: [5.0, 10.0, 15.0])

    Returns:
        dict: Dictionary of depth metrics with multiple parameter combinations
        min_depth: Minimum depth using default parameters (threshold=15, cutoff=5.0)
        avg_depth: Average depth using default parameters
        ca_depth: CA atom depth using default parameters
        hydroxyl_depth: Hydroxyl depth using default parameters
    """
    # Set up results dictionary
    results = {}

    # Track default parameter results for backward compatibility
    default_min_depth = None
    default_avg_depth = None
    default_ca_depth = None
    default_hydroxyl_depth = None

    # Collect all atom coordinates (once for all parameter combinations)
    all_atoms = []
    for chain in model:
        for r in chain:
            if is_aa(r):
                for atom in r:
                    all_atoms.append(atom.coord)

    if not all_atoms:
        return results, None, None, None, None

    # Convert to numpy array
    all_atoms = np.array(all_atoms)

    # Build KDTree for efficient neighbor search (once for all combinations)
    atom_tree = cKDTree(all_atoms)

    # Collect residue atoms (once for all parameter combinations)
    residue_atoms = []
    ca_index = None
    hydroxyl_index = None

    for i, atom in enumerate(residue):
        residue_atoms.append(atom.coord)

        # Track CA and hydroxyl indices
        if atom.name == 'CA':
            ca_index = i
        elif ((residue.resname == 'SER' and atom.name == 'OG') or
              (residue.resname == 'THR' and atom.name == 'OG1') or
              (residue.resname == 'TYR' and atom.name == 'OH')):
            hydroxyl_index = i

    if not residue_atoms:
        return results, None, None, None, None

    # Convert to numpy array
    residue_atoms = np.array(residue_atoms)

    # Calculate for each parameter combination
    for threshold in thresholds:
        for cutoff in cutoffs:
            # Key suffix for this parameter combination
            param_key = f"t{threshold}_c{int(cutoff)}"

            # Check for cached surface atoms with these parameters
            cache_key = f"_cached_surface_atoms_{param_key}"
            if hasattr(model, cache_key):
                surface_atoms = getattr(model, cache_key)
            else:
                # Find surface atoms (those with fewer neighbors)
                surface_atoms = []

                # Use a more efficient approach with query_ball_point
                for i, atom_coord in enumerate(all_atoms):
                    # Find all atoms within cutoff
                    indices = atom_tree.query_ball_point(atom_coord, cutoff)

                    # Count neighbors (excluding self)
                    neighbor_count = len(indices) - 1

                    # Add to surface atoms if below threshold
                    if neighbor_count < threshold:
                        surface_atoms.append(atom_coord)

                # If we have too few surface atoms, just use all atoms
                if len(surface_atoms) < max(5, len(all_atoms) * 0.01):  # At least 1% of atoms or 5 atoms
                    surface_atoms = all_atoms

                # Convert to numpy array
                surface_atoms = np.array(surface_atoms)

                # Cache the result for future calls
                setattr(model, cache_key, surface_atoms)

            # Check for cached surface tree with these parameters
            tree_cache_key = f"_cached_surface_tree_{param_key}"
            if hasattr(model, tree_cache_key):
                surface_tree = getattr(model, tree_cache_key)
            else:
                surface_tree = cKDTree(surface_atoms)
                # Cache the tree
                setattr(model, tree_cache_key, surface_tree)

            # Batch query to get all depths at once
            distances, _ = surface_tree.query(residue_atoms)

            # Calculate depths for this parameter combination
            min_depth = np.min(distances) if len(distances) > 0 else None
            avg_depth = np.mean(distances) if len(distances) > 0 else None
            ca_depth = distances[ca_index] if ca_index is not None else None
            hydroxyl_depth = distances[hydroxyl_index] if hydroxyl_index is not None else None

            # Store results for this parameter combination
            results[f"min_depth_{param_key}"] = min_depth
            results[f"avg_depth_{param_key}"] = avg_depth
            results[f"ca_depth_{param_key}"] = ca_depth
            results[f"hydroxyl_depth_{param_key}"] = hydroxyl_depth

            # Store default parameter results (threshold=15, cutoff=5.0)
            if threshold == 15 and cutoff == 5.0:
                default_min_depth = min_depth
                default_avg_depth = avg_depth
                default_ca_depth = ca_depth
                default_hydroxyl_depth = hydroxyl_depth

    # If default parameters weren't in the provided thresholds/cutoffs, use the first combination
    if default_min_depth is None and thresholds and cutoffs:
        first_key = f"t{thresholds[0]}_c{int(cutoffs[0])}"
        default_min_depth = results.get(f"min_depth_{first_key}")
        default_avg_depth = results.get(f"avg_depth_{first_key}")
        default_ca_depth = results.get(f"ca_depth_{first_key}")
        default_hydroxyl_depth = results.get(f"hydroxyl_depth_{first_key}")

    # Return both the complete results dictionary and the 4 specific values expected by the calling code
    return results, default_min_depth, default_avg_depth, default_ca_depth, default_hydroxyl_depth


def calculate_side_chain_angles(residue):
    """Calculate side chain dihedral and angles for STY residues

    Args:
        residue: The target residue (SER, THR, or TYR)

    Returns:
        dict: Dictionary of angles in degrees
    """
    # Initialize results
    angles = {
        'chi1': None,  # N-CA-CB-OG/OG1/CG
        'chi2': None,  # CA-CB-CG-CD1 (for TYR)
        'hydroxyl_angle': None,  # CA-CB-OG angle (for SER/THR)
        'ca_cb_og_dist': None,  # Distance from CA to hydroxyl
    }

    # Check if we have necessary atoms
    required_atoms = ['N', 'CA', 'CB']
    if not all(atom in residue for atom in required_atoms):
        return angles

    # Get atom coordinates
    n_coord = residue['N'].coord
    ca_coord = residue['CA'].coord
    cb_coord = residue['CB'].coord

    # Additional atoms depending on residue type
    if residue.resname == 'SER' and 'OG' in residue:
        og_coord = residue['OG'].coord

        # Calculate chi1 (N-CA-CB-OG)
        try:
            chi1 = calc_dihedral(Vector(*n_coord), Vector(*ca_coord),
                                 Vector(*cb_coord), Vector(*og_coord))
            angles['chi1'] = degrees(chi1)
        except Exception:
            pass

        # Calculate hydroxyl angle (CA-CB-OG)
        try:
            hydroxyl_angle = calc_angle(Vector(*ca_coord), Vector(*cb_coord), Vector(*og_coord))
            angles['hydroxyl_angle'] = degrees(hydroxyl_angle)
        except Exception:
            pass

        # Calculate distance from CA to hydroxyl
        try:
            ca_og_dist = np.sqrt(np.sum((np.array(ca_coord) - np.array(og_coord)) ** 2))
            angles['ca_cb_og_dist'] = ca_og_dist
        except Exception:
            pass

    elif residue.resname == 'THR' and 'OG1' in residue:
        og1_coord = residue['OG1'].coord

        # Calculate chi1 (N-CA-CB-OG1)
        try:
            chi1 = calc_dihedral(Vector(*n_coord), Vector(*ca_coord),
                                 Vector(*cb_coord), Vector(*og1_coord))
            angles['chi1'] = degrees(chi1)
        except Exception:
            pass

        # Calculate hydroxyl angle (CA-CB-OG1)
        try:
            hydroxyl_angle = calc_angle(Vector(*ca_coord), Vector(*cb_coord), Vector(*og1_coord))
            angles['hydroxyl_angle'] = degrees(hydroxyl_angle)
        except Exception:
            pass

        # Calculate distance from CA to hydroxyl
        try:
            ca_og_dist = np.sqrt(np.sum((np.array(ca_coord) - np.array(og1_coord)) ** 2))
            angles['ca_cb_og_dist'] = ca_og_dist
        except Exception:
            pass

    elif residue.resname == 'TYR':
        if 'CG' in residue:
            cg_coord = residue['CG'].coord

            # Calculate chi1 (N-CA-CB-CG)
            try:
                chi1 = calc_dihedral(Vector(*n_coord), Vector(*ca_coord),
                                     Vector(*cb_coord), Vector(*cg_coord))
                angles['chi1'] = degrees(chi1)
            except Exception:
                pass

            # Calculate chi2 if we have CD1
            if 'CD1' in residue:
                cd1_coord = residue['CD1'].coord
                try:
                    chi2 = calc_dihedral(Vector(*ca_coord), Vector(*cb_coord),
                                         Vector(*cg_coord), Vector(*cd1_coord))
                    angles['chi2'] = degrees(chi2)
                except Exception:
                    pass

            # Calculate hydroxyl angle and distance if we have OH
            if 'OH' in residue:
                oh_coord = residue['OH'].coord
                try:
                    ca_oh_dist = np.sqrt(np.sum((np.array(ca_coord) - np.array(oh_coord)) ** 2))
                    angles['ca_cb_og_dist'] = ca_oh_dist
                except Exception:
                    pass

    return angles


def calculate_packing_density(residue, neighbors, shell_radii=[5.0, 8.0, 10.0, 12.0, 15.0]):
    """Calculate atom packing density in concentric shells around residue

    Args:
        residue: The target residue
        neighbors: List of neighboring residues
        shell_radii: List of shell radii to calculate density for

    Returns:
        dict: Dictionary mapping shell radius to atom density
    """
    # Get CA atom as center
    try:
        ca_atom = residue['CA']
        ca_coord = np.array(ca_atom.coord)
    except KeyError:
        return {}

    # Also calculate packing around hydroxyl if it exists
    hydroxyl_atom = None
    if residue.resname == 'SER' and 'OG' in residue:
        hydroxyl_atom = residue['OG']
    elif residue.resname == 'THR' and 'OG1' in residue:
        hydroxyl_atom = residue['OG1']
    elif residue.resname == 'TYR' and 'OH' in residue:
        hydroxyl_atom = residue['OH']

    hydroxyl_coord = np.array(hydroxyl_atom.coord) if hydroxyl_atom else None

    # Collect neighbor atoms
    neighbor_atoms = []
    for neighbor in neighbors:
        for atom in neighbor:
            neighbor_atoms.append(atom.coord)

    # Add residue's own atoms
    for atom in residue:
        neighbor_atoms.append(atom.coord)

    if not neighbor_atoms:
        return {f"packing_density_{r}A": 0 for r in shell_radii}

    # Convert to numpy array
    neighbor_atoms = np.array(neighbor_atoms)

    # Calculate densities for each shell
    densities = {}
    prev_radius = 0

    # CA-centered packing
    for radius in shell_radii:
        # Calculate distances to CA atom
        distances = np.sqrt(np.sum((neighbor_atoms - ca_coord) ** 2, axis=1))

        # Count atoms in current shell
        atoms_in_shell = np.sum((distances > prev_radius) & (distances <= radius))

        # Calculate shell volume (spherical shell)
        shell_volume = (4 / 3) * np.pi * (radius ** 3 - prev_radius ** 3)

        # Calculate density
        density = atoms_in_shell / shell_volume if shell_volume > 0 else 0
        densities[f"packing_density_{radius}A"] = density

        prev_radius = radius

    # Calculate hydroxyl packing if available
    if hydroxyl_coord is not None:
        prev_radius = 0
        for radius in shell_radii:
            # Calculate distances to hydroxyl oxygen
            distances = np.sqrt(np.sum((neighbor_atoms - hydroxyl_coord) ** 2, axis=1))

            # Count atoms in current shell
            atoms_in_shell = np.sum((distances > prev_radius) & (distances <= radius))

            # Calculate shell volume
            shell_volume = (4 / 3) * np.pi * (radius ** 3 - prev_radius ** 3)

            # Calculate density
            density = atoms_in_shell / shell_volume if shell_volume > 0 else 0
            densities[f"hydroxyl_packing_{radius}A"] = density

            prev_radius = radius

    return densities


def calculate_contact_order(residue, chain_data, neighbors, max_dist=8.0):
    """Calculate contact order for the target residue

    Args:
        residue: The target residue
        chain_data: Pre-computed chain data
        neighbors: List of neighboring residues
        max_dist: Maximum distance for contacts

    Returns:
        dict: Dictionary of contact order metrics
    """
    # Get residue ID and sequence position
    res_id = residue.id[1]
    if res_id not in chain_data['seq_ids']:
        return {
            'contact_order': None,
            'long_range_contacts': 0,
            'medium_range_contacts': 0,
            'short_range_contacts': 0,
            'local_contacts': 0
        }

    seq_pos = chain_data['seq_ids'].index(res_id)

    # Get CA atom coordinates
    try:
        ca_atom = residue['CA']
    except KeyError:
        return {
            'contact_order': None,
            'long_range_contacts': 0,
            'medium_range_contacts': 0,
            'short_range_contacts': 0,
            'local_contacts': 0
        }

    ca_coord = np.array(ca_atom.coord)

    # Track contacts by sequence separation
    long_range = 0  # |i-j| >= 12
    medium_range = 0  # 6 <= |i-j| < 12
    short_range = 0  # 2 <= |i-j| < 6
    local_contacts = 0  # |i-j| = 1

    # Total sequence separation
    total_separation = 0
    total_contacts = 0

    # Threshold distance squared (for faster comparison)
    max_dist_squared = max_dist * max_dist

    # Check each neighbor
    for neighbor in neighbors:
        # Only consider residues in the same chain
        if neighbor.get_parent() != residue.get_parent():
            continue

        # Get neighbor's sequence position
        neighbor_id = neighbor.id[1]
        if neighbor_id not in chain_data['seq_ids']:
            continue

        neighbor_pos = chain_data['seq_ids'].index(neighbor_id)
        seq_separation = abs(seq_pos - neighbor_pos)

        # Skip if it's the same residue
        if seq_separation == 0:
            continue

        # Check spatial proximity (using CA atoms)
        try:
            neighbor_ca = neighbor['CA']
        except KeyError:
            continue

        # Calculate squared distance (faster than square root)
        dist_sq = np.sum((ca_coord - neighbor_ca.coord) ** 2)

        # Count as contact if within distance threshold
        if dist_sq <= max_dist_squared:
            total_contacts += 1
            total_separation += seq_separation

            # Categorize contact by sequence separation
            if seq_separation >= 12:
                long_range += 1
            elif seq_separation >= 6:
                medium_range += 1
            elif seq_separation >= 2:
                short_range += 1
            else:  # seq_separation == 1
                local_contacts += 1

    # Calculate contact order
    contact_order = total_separation / total_contacts if total_contacts > 0 else None

    return {
        'contact_order': contact_order,
        'long_range_contacts': long_range,
        'medium_range_contacts': medium_range,
        'short_range_contacts': short_range,
        'local_contacts': local_contacts
    }


def calculate_hydrogen_bond_energy(residue, neighbors, max_dist=3.5):
    """Calculate approximate hydrogen bond energies

    Args:
        residue: The target residue
        neighbors: List of neighboring residues
        max_dist: Maximum distance for hydrogen bond detection

    Returns:
        dict: Dictionary of hydrogen bond energy metrics
    """
    # Get the hydroxyl oxygen atom
    hydroxyl_atom = None
    if residue.resname == 'SER' and 'OG' in residue:
        hydroxyl_atom = residue['OG']
    elif residue.resname == 'THR' and 'OG1' in residue:
        hydroxyl_atom = residue['OG1']
    elif residue.resname == 'TYR' and 'OH' in residue:
        hydroxyl_atom = residue['OH']

    if not hydroxyl_atom:
        return {
            'hbond_energy': None,
            'hbond_count': 0,
            'strongest_hbond': None,
            'backbone_hbonds': 0,
            'sidechain_hbonds': 0
        }

    # Get hydroxyl coordinates
    hydroxyl_coord = np.array(hydroxyl_atom.coord)

    # Identify potential hydrogen bond partners (O and N atoms)
    h_bonds = []
    backbone_hbonds = 0
    sidechain_hbonds = 0

    for neighbor in neighbors:
        # Skip the residue itself
        if neighbor is residue:
            continue

        for atom in neighbor:
            # Only consider oxygen and nitrogen atoms as potential partners
            if atom.element not in ['O', 'N']:
                continue

            # Calculate distance
            dist = np.sqrt(np.sum((hydroxyl_coord - atom.coord) ** 2))

            # Consider as hydrogen bond if within distance cutoff
            if dist <= max_dist:
                # Calculate energy using a simple distance-dependent model
                # E = -5 * (1/d - 1/3.5) kcal/mol
                energy = -5.0 * (1 / dist - 1 / 3.5)
                h_bonds.append((dist, energy, atom))

                # Track backbone vs sidechain hydrogen bonds
                if atom.name in ['N', 'O', 'C', 'CA']:
                    backbone_hbonds += 1
                else:
                    sidechain_hbonds += 1

    # Calculate total energy and find strongest bond
    total_energy = sum(e for _, e, _ in h_bonds)
    strongest_hbond = min(h_bonds)[1] if h_bonds else None

    return {
        'hbond_energy': total_energy,
        'hbond_count': len(h_bonds),
        'strongest_hbond': strongest_hbond,
        'backbone_hbonds': backbone_hbonds,
        'sidechain_hbonds': sidechain_hbonds
    }


def calculate_b_factor_statistics(residue, neighbors, plddt_dict=None):
    """Calculate extended B-factor statistics for flexibility analysis

    Args:
        residue: The target residue
        neighbors: List of neighboring residues
        plddt_dict: Optional dictionary of pLDDT values (AlphaFold)

    Returns:
        dict: Dictionary of B-factor statistics
    """
    # Get the residue ID
    res_id = residue.id[1]

    # Get B-factors for the STY residue atoms
    sty_b_factors = [atom.bfactor for atom in residue]
    hydroxyl_b_factor = None

    # Get the hydroxyl atom's B-factor
    if residue.resname == 'SER' and 'OG' in residue:
        hydroxyl_b_factor = residue['OG'].bfactor
    elif residue.resname == 'THR' and 'OG1' in residue:
        hydroxyl_b_factor = residue['OG1'].bfactor
    elif residue.resname == 'TYR' and 'OH' in residue:
        hydroxyl_b_factor = residue['OH'].bfactor

    # Calculate statistics
    sty_avg_b = np.mean(sty_b_factors) if sty_b_factors else None
    sty_std_b = np.std(sty_b_factors) if len(sty_b_factors) > 1 else None
    sty_min_b = min(sty_b_factors) if sty_b_factors else None
    sty_max_b = max(sty_b_factors) if sty_b_factors else None

    # Get B-factors for the neighbor residues
    neighbor_b_factors = []

    for neighbor in neighbors:
        if neighbor is not residue:  # Skip the STY residue itself
            neighbor_b_factors.extend([atom.bfactor for atom in neighbor])

    # Calculate neighbor statistics
    neighbor_avg_b = np.mean(neighbor_b_factors) if neighbor_b_factors else None
    neighbor_std_b = np.std(neighbor_b_factors) if len(neighbor_b_factors) > 1 else None

    # Calculate relative flexibility (ratio of residue B-factor to neighbors)
    relative_flexibility = None
    if sty_avg_b is not None and neighbor_avg_b is not None and neighbor_avg_b > 0:
        relative_flexibility = sty_avg_b / neighbor_avg_b

    # Calculate normalized B-factor (Z-score relative to neighbors)
    normalized_b = None
    if sty_avg_b is not None and neighbor_avg_b is not None and neighbor_std_b is not None and neighbor_std_b > 0:
        normalized_b = (sty_avg_b - neighbor_avg_b) / neighbor_std_b

    # Compare with pLDDT if available
    plddt = plddt_dict.get(res_id) if plddt_dict else None

    return {
        'avg_bfactor': sty_avg_b,
        'std_bfactor': sty_std_b,
        'min_bfactor': sty_min_b,
        'max_bfactor': sty_max_b,
        'hydroxyl_bfactor': hydroxyl_b_factor,
        'neighbor_avg_bfactor': neighbor_avg_b,
        'neighbor_std_bfactor': neighbor_std_b,
        'relative_flexibility': relative_flexibility,
        'normalized_bfactor': normalized_b,
        'plddt': plddt
    }


def analyze_pocket_shape(residue, neighbors):
    """Analyze the shape of the pocket around STY hydroxyl

    Args:
        residue: The target residue
        neighbors: List of neighboring residues

    Returns:
        dict: Dictionary of pocket shape metrics
    """
    # Get the hydroxyl oxygen atom
    hydroxyl_atom = None
    if residue.resname == 'SER' and 'OG' in residue:
        hydroxyl_atom = residue['OG']
    elif residue.resname == 'THR' and 'OG1' in residue:
        hydroxyl_atom = residue['OG1']
    elif residue.resname == 'TYR' and 'OH' in residue:
        hydroxyl_atom = residue['OH']

    if not hydroxyl_atom:
        return {
            'pocket_width_x': None,
            'pocket_width_y': None,
            'pocket_width_z': None,
            'pocket_volume': None,
            'pocket_asymmetry': None,
            'pocket_atom_count': 0
        }

    # Get hydroxyl coordinates
    hydroxyl_coord = np.array(hydroxyl_atom.coord)

    # Collect coordinates of atoms within pocket (10Å from hydroxyl)
    pocket_atoms = []
    for neighbor in neighbors:
        for atom in neighbor:
            # Calculate distance to hydroxyl
            dist = np.sqrt(np.sum((hydroxyl_coord - atom.coord) ** 2))

            # Consider part of pocket if within 10Å
            if dist <= 10.0:
                pocket_atoms.append(atom.coord)

    # Add the residue's own atoms
    for atom in residue:
        if atom != hydroxyl_atom:  # Exclude hydroxyl itself
            dist = np.sqrt(np.sum((hydroxyl_coord - atom.coord) ** 2))
            if dist <= 10.0:
                pocket_atoms.append(atom.coord)

    # Check if we have enough atoms
    if len(pocket_atoms) < 4:
        return {
            'pocket_width_x': None,
            'pocket_width_y': None,
            'pocket_width_z': None,
            'pocket_volume': None,
            'pocket_asymmetry': None,
            'pocket_atom_count': len(pocket_atoms)
        }

    # Convert to numpy array
    pocket_atoms = np.array(pocket_atoms)

    # Center the pocket at the hydroxyl
    centered_atoms = pocket_atoms - hydroxyl_coord

    # Calculate PCA to find principal axes of the pocket
    try:
        # Calculate covariance matrix
        cov_matrix = np.cov(centered_atoms, rowvar=False)

        # Calculate eigenvalues and eigenvectors
        eigenvalues, eigenvectors = np.linalg.eigh(cov_matrix)

        # Sort eigenvalues in descending order
        idx = eigenvalues.argsort()[::-1]
        eigenvalues = eigenvalues[idx]
        eigenvectors = eigenvectors[:, idx]

        # Transform points to principal component space
        transformed = np.dot(centered_atoms, eigenvectors)

        # Calculate extents along each principal axis
        extents = np.max(transformed, axis=0) - np.min(transformed, axis=0)

        # Estimate volume as product of extents
        volume = np.prod(extents)

        # Calculate asymmetry as ratio of largest to smallest extent
        asymmetry = extents[0] / extents[2] if extents[2] > 0 else None

        return {
            'pocket_width_x': extents[0],
            'pocket_width_y': extents[1],
            'pocket_width_z': extents[2],
            'pocket_volume': volume,
            'pocket_asymmetry': asymmetry,
            'pocket_atom_count': len(pocket_atoms)
        }
    except (np.linalg.LinAlgError, ZeroDivisionError):
        # Handle cases where PCA fails
        return {
            'pocket_width_x': None,
            'pocket_width_y': None,
            'pocket_width_z': None,
            'pocket_volume': None,
            'pocket_asymmetry': None,
            'pocket_atom_count': len(pocket_atoms)
        }


def calculate_local_structural_features(residue, chain_data):
    """Calculate local structural features around the STY site

    Args:
        residue: The target residue
        chain_data: Pre-computed chain data

    Returns:
        dict: Dictionary of local structural features
    """
    # Get residue ID and sequence position
    res_id = residue.id[1]
    if res_id not in chain_data['seq_ids']:
        return {
            'local_rmsd': None,
            'phi': None,
            'psi': None,
            'omega': None,
            'ramachandran_region': None
        }

    seq_pos = chain_data['seq_ids'].index(res_id)

    # Get CA atom
    try:
        ca_atom = residue['CA']
    except KeyError:
        return {
            'local_rmsd': None,
            'phi': None,
            'psi': None,
            'omega': None,
            'ramachandran_region': None
        }

    # Calculate backbone torsion angles
    phi, psi, omega = None, None, None

    # Get chain
    chain = residue.get_parent()

    # Find previous and next residues
    prev_res, next_res = None, None

    # Check previous residue
    if seq_pos > 0:
        prev_id = chain_data['seq_ids'][seq_pos - 1]
        for res in chain:
            if res.id[1] == prev_id:
                prev_res = res
                break

    # Check next residue
    if seq_pos < len(chain_data['seq_ids']) - 1:
        next_id = chain_data['seq_ids'][seq_pos + 1]
        for res in chain:
            if res.id[1] == next_id:
                next_res = res
                break

    # Calculate phi (C-1, N, CA, C)
    if prev_res and 'C' in prev_res and 'N' in residue and 'CA' in residue and 'C' in residue:
        try:
            phi_rad = calc_dihedral(
                Vector(*prev_res['C'].coord),
                Vector(*residue['N'].coord),
                Vector(*residue['CA'].coord),
                Vector(*residue['C'].coord)
            )
            phi = degrees(phi_rad)
        except Exception:
            phi = None

    # Calculate psi (N, CA, C, N+1)
    if 'N' in residue and 'CA' in residue and 'C' in residue and next_res and 'N' in next_res:
        try:
            psi_rad = calc_dihedral(
                Vector(*residue['N'].coord),
                Vector(*residue['CA'].coord),
                Vector(*residue['C'].coord),
                Vector(*next_res['N'].coord)
            )
            psi = degrees(psi_rad)
        except Exception:
            psi = None

    # Calculate omega (CA, C, N+1, CA+1)
    if 'CA' in residue and 'C' in residue and next_res and 'N' in next_res and 'CA' in next_res:
        try:
            omega_rad = calc_dihedral(
                Vector(*residue['CA'].coord),
                Vector(*residue['C'].coord),
                Vector(*next_res['N'].coord),
                Vector(*next_res['CA'].coord)
            )
            omega = degrees(omega_rad)
        except Exception:
            omega = None

    # Determine Ramachandran region
    ramachandran_region = None
    if phi is not None and psi is not None:
        if (-180 <= phi <= -30) and (-100 <= psi <= 45):
            ramachandran_region = 'alpha'
        elif (phi <= -30) and (psi >= 45):
            ramachandran_region = 'beta'
        elif (-180 <= phi <= -30) and (-180 <= psi <= -100):
            ramachandran_region = 'beta'
        elif (-30 <= phi <= 180) and (-180 <= psi <= 180):
            ramachandran_region = 'left'
        else:
            ramachandran_region = 'other'

    # Calculate local RMSD of CA atoms (±3 residues)
    local_rmsd = None
    window_size = 3

    # Get local CA coordinates
    ca_coords = []

    # Check window around residue
    for i in range(max(0, seq_pos - window_size), min(len(chain_data['seq_ids']), seq_pos + window_size + 1)):
        window_id = chain_data['seq_ids'][i]

        # Find residue
        for res in chain:
            if res.id[1] == window_id and 'CA' in res:
                ca_coords.append(res['CA'].coord)
                break

    # Need at least 3 CA atoms to calculate RMSD
    if len(ca_coords) >= 3:
        # Calculate centroid
        ca_coords = np.array(ca_coords)
        centroid = np.mean(ca_coords, axis=0)

        # Calculate RMSD
        squared_diffs = np.sum((ca_coords - centroid) ** 2, axis=1)
        local_rmsd = np.sqrt(np.mean(squared_diffs))

    return {
        'local_rmsd': local_rmsd,
        'phi': phi,
        'psi': psi,
        'omega': omega,
        'ramachandran_region': ramachandran_region
    }


def calculate_hydrophobic_feature(residue, neighbors, radius=6.0):
    """Calculate hydrophobic environment features

    Args:
        residue: The target residue
        neighbors: List of neighboring residues
        radius: Radius for hydrophobic patch detection

    Returns:
        dict: Dictionary of hydrophobic environment metrics
    """
    # Get the hydroxyl oxygen atom
    hydroxyl_atom = None
    if residue.resname == 'SER' and 'OG' in residue:
        hydroxyl_atom = residue['OG']
    elif residue.resname == 'THR' and 'OG1' in residue:
        hydroxyl_atom = residue['OG1']
    elif residue.resname == 'TYR' and 'OH' in residue:
        hydroxyl_atom = residue['OH']

    if not hydroxyl_atom:
        return {
            'hydrophobic_moment': None,
            'hydrophobic_patch_size': 0,
            'hydrophobicity_variance': None
        }

    # Get hydroxyl coordinates
    hydroxyl_coord = np.array(hydroxyl_atom.coord)

    # Collect coordinates of residues within radius
    hydrophobic_residues = []
    hydrophobic_vectors = []
    residue_coords = []
    hydrophobicity_values = []

    # Hydrophobicity scale (Kyte-Doolittle)
    hydrophobicity = {
        'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
        'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
        'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
        'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2,
    }

    for neighbor in neighbors:
        # Skip the residue itself
        if neighbor is residue:
            continue

        # Get one-letter code
        one_letter = three_to_one(neighbor.resname)

        # Skip non-standard amino acids
        if one_letter == 'X':
            continue

        # Get CA atom
        if 'CA' not in neighbor:
            continue

        ca_coord = np.array(neighbor['CA'].coord)

        # Calculate distance to hydroxyl
        dist = np.sqrt(np.sum((hydroxyl_coord - ca_coord) ** 2))

        # Consider if within radius
        if dist <= radius:
            # Get hydrophobicity value
            h_value = hydrophobicity.get(one_letter, 0)
            hydrophobicity_values.append(h_value)

            # Calculate vector from hydroxyl to CA
            vector = ca_coord - hydroxyl_coord

            # Normalize vector
            norm = np.sqrt(np.sum(vector ** 2))
            if norm > 0:
                vector = vector / norm

            # Scale vector by hydrophobicity
            vector = vector * h_value

            # Store residue info
            residue_coords.append(ca_coord)
            hydrophobic_vectors.append(vector)

            # Track hydrophobic residues
            if h_value > 0:
                hydrophobic_residues.append(neighbor)

    # Count hydrophobic patch size
    hydrophobic_patch_size = len(hydrophobic_residues)

    # Calculate hydrophobicity variance
    hydrophobicity_variance = np.var(hydrophobicity_values) if len(hydrophobicity_values) > 0 else None

    # Calculate hydrophobic moment (sum of hydrophobicity vectors)
    if len(hydrophobic_vectors) > 0:
        # Sum vectors
        moment = np.sum(hydrophobic_vectors, axis=0)

        # Calculate magnitude
        hydrophobic_moment = np.sqrt(np.sum(moment ** 2))
    else:
        hydrophobic_moment = None

    return {
        'hydrophobic_moment': hydrophobic_moment,
        'hydrophobic_patch_size': hydrophobic_patch_size,
        'hydrophobicity_variance': hydrophobicity_variance
    }


def calculate_electrostatic_features(residue, neighbors, radius=10.0):
    """Calculate electrostatic environment features

    Args:
        residue: The target residue
        neighbors: List of neighboring residues
        radius: Radius for electrostatic calculations

    Returns:
        dict: Dictionary of electrostatic environment metrics
    """
    # Get the hydroxyl oxygen atom
    hydroxyl_atom = None
    if residue.resname == 'SER' and 'OG' in residue:
        hydroxyl_atom = residue['OG']
    elif residue.resname == 'THR' and 'OG1' in residue:
        hydroxyl_atom = residue['OG1']
    elif residue.resname == 'TYR' and 'OH' in residue:
        hydroxyl_atom = residue['OH']

    if not hydroxyl_atom:
        return {
            'net_charge': 0,
            'electrostatic_potential': 0,
            'charge_density': 0,
            'pos_charge_count': 0,
            'neg_charge_count': 0
        }

    # Get hydroxyl coordinates
    hydroxyl_coord = np.array(hydroxyl_atom.coord)

    # Count charged residues near hydroxyl
    pos_charges = []  # (position, charge)
    neg_charges = []  # (position, charge)

    for neighbor in neighbors:
        # Skip if not amino acid
        if not is_aa(neighbor):
            continue

        # Get one-letter code
        one_letter = three_to_one(neighbor.resname)

        # Skip non-standard amino acids
        if one_letter == 'X':
            continue

        # Get charge properties
        if one_letter in aa_properties:
            charge = aa_properties[one_letter]['charge']

            # Skip uncharged residues
            if charge == 0:
                continue

            # Get CA atom
            if 'CA' not in neighbor:
                continue

            ca_coord = np.array(neighbor['CA'].coord)

            # Calculate distance to hydroxyl
            dist = np.sqrt(np.sum((hydroxyl_coord - ca_coord) ** 2))

            # Consider if within radius
            if dist <= radius:
                if charge > 0:
                    pos_charges.append((ca_coord, charge))
                else:
                    neg_charges.append((ca_coord, charge))

    # Calculate net charge
    net_charge = sum(c for _, c in pos_charges) + sum(c for _, c in neg_charges)

    # Calculate electrostatic potential at hydroxyl
    potential = 0
    for pos, charge in pos_charges + neg_charges:
        # Calculate distance
        dist = np.sqrt(np.sum((hydroxyl_coord - pos) ** 2))

        # Skip if too close (avoid division by zero)
        if dist < 0.1:
            continue

        # Simple electrostatic potential: q/r
        potential += charge / dist

    # Calculate charge density (charges per unit volume)
    volume = (4 / 3) * np.pi * radius ** 3
    charge_density = (len(pos_charges) + len(neg_charges)) / volume

    return {
        'net_charge': net_charge,
        'electrostatic_potential': potential,
        'charge_density': charge_density,
        'pos_charge_count': len(pos_charges),
        'neg_charge_count': len(neg_charges)
    }


def calculate_aromatic_features(residue, neighbors, radius=8.0):
    """Calculate aromatic environment features

    Args:
        residue: The target residue
        neighbors: List of neighboring residues
        radius: Radius for aromatic interactions

    Returns:
        dict: Dictionary of aromatic environment metrics
    """
    # Get the hydroxyl oxygen atom
    hydroxyl_atom = None
    if residue.resname == 'SER' and 'OG' in residue:
        hydroxyl_atom = residue['OG']
    elif residue.resname == 'THR' and 'OG1' in residue:
        hydroxyl_atom = residue['OG1']
    elif residue.resname == 'TYR' and 'OH' in residue:
        hydroxyl_atom = residue['OH']

    if not hydroxyl_atom:
        return {
            'aromatic_count': 0,
            'nearest_aromatic_dist': None,
            'aromatic_stacking': 0,
            'pi_cation': 0
        }

    # Get hydroxyl coordinates
    hydroxyl_coord = np.array(hydroxyl_atom.coord)

    # Find aromatic residues
    aromatic_residues = []  # (residue, distance)
    aromatic_distances = []
    aromatic_stacking = 0  # Count potential aromatic stacking interactions
    pi_cation = 0  # Count potential pi-cation interactions

    # First, identify aromatic and charged residues
    for neighbor in neighbors:
        # Skip if not amino acid
        if not is_aa(neighbor):
            continue

        # Get one-letter code
        one_letter = three_to_one(neighbor.resname)

        # Skip non-standard amino acids
        if one_letter == 'X':
            continue

        # Get CA atom
        if 'CA' not in neighbor:
            continue

        ca_coord = np.array(neighbor['CA'].coord)

        # Calculate distance to hydroxyl
        dist = np.sqrt(np.sum((hydroxyl_coord - ca_coord) ** 2))

        # Consider if within radius
        if dist <= radius:
            # Check if aromatic
            if one_letter in aa_properties and aa_properties[one_letter]['aromatic']:
                aromatic_residues.append((neighbor, dist))
                aromatic_distances.append(dist)

            # Check for charged residues (potential pi-cation)
            if one_letter in aa_properties and aa_properties[one_letter]['charged'] and aa_properties[one_letter][
                'charge'] > 0:
                # Look for nearby aromatic residues
                for other_res, other_dist in aromatic_residues:
                    # Check if this charged residue is close to an aromatic residue
                    other_ca = np.array(other_res['CA'].coord)
                    charge_to_aromatic = np.sqrt(np.sum((ca_coord - other_ca) ** 2))

                    # Count if within 6Å
                    if charge_to_aromatic <= 6.0:
                        pi_cation += 1

    # Check for aromatic stacking
    for i in range(len(aromatic_residues)):
        for j in range(i + 1, len(aromatic_residues)):
            res1, _ = aromatic_residues[i]
            res2, _ = aromatic_residues[j]

            # Get CA atoms
            ca1 = np.array(res1['CA'].coord)
            ca2 = np.array(res2['CA'].coord)

            # Calculate distance
            dist = np.sqrt(np.sum((ca1 - ca2) ** 2))

            # Count if within stacking distance
            if dist <= 7.0:
                aromatic_stacking += 1

    # Calculate nearest aromatic distance
    nearest_aromatic_dist = min(aromatic_distances) if aromatic_distances else None

    return {
        'aromatic_count': len(aromatic_residues),
        'nearest_aromatic_dist': nearest_aromatic_dist,
        'aromatic_stacking': aromatic_stacking,
        'pi_cation': pi_cation
    }


def calculate_pocket_solvent_features(residue, neighbors):
    """Calculate solvent-related pocket features with optimized performance"""
    # Get the hydroxyl oxygen atom
    hydroxyl_atom = None
    if residue.resname == 'SER' and 'OG' in residue:
        hydroxyl_atom = residue['OG']
    elif residue.resname == 'THR' and 'OG1' in residue:
        hydroxyl_atom = residue['OG1']
    elif residue.resname == 'TYR' and 'OH' in residue:
        hydroxyl_atom = residue['OH']

    if not hydroxyl_atom:
        return {
            'polar_ratio': 0,
            'hydrophobic_ratio': 0,
            'water_channel_score': 0,
            'polar_atom_count': 0
        }

    # Get hydroxyl coordinates
    hydroxyl_coord = np.array(hydroxyl_atom.coord)

    # Collect all atom coordinates at once
    all_atom_coords = []
    is_polar = []
    hydrophobic_residues = 0

    for neighbor in neighbors:
        if not is_aa(neighbor):
            continue

        # Check if hydrophobic
        one_letter = three_to_one(neighbor.resname)
        if one_letter in aa_properties and aa_properties[one_letter]['hydrophobic']:
            if 'CA' in neighbor:
                ca_dist = np.sqrt(np.sum((hydroxyl_coord - np.array(neighbor['CA'].coord)) ** 2))
                if ca_dist <= 8.0:
                    hydrophobic_residues += 1

        # Collect atoms
        for atom in neighbor:
            all_atom_coords.append(atom.coord)
            is_polar.append(atom.element in ['O', 'N'])

    # Vectorized calculations
    if not all_atom_coords:
        return {
            'polar_ratio': 0,
            'hydrophobic_ratio': 0,
            'water_channel_score': 0,
            'polar_atom_count': 0
        }

    all_atom_coords = np.array(all_atom_coords)
    is_polar = np.array(is_polar)

    # Calculate distances
    vectors = all_atom_coords - hydroxyl_coord
    distances_sq = np.sum(vectors * vectors, axis=1)
    distances = np.sqrt(distances_sq)

    # Find pocket atoms
    pocket_mask = distances <= 6.0
    pocket_atoms = np.sum(pocket_mask)
    polar_atoms = np.sum(is_polar & pocket_mask)

    # Calculate ratios
    polar_ratio = polar_atoms / pocket_atoms if pocket_atoms > 0 else 0
    hydrophobic_ratio = hydrophobic_residues / len(neighbors) if neighbors else 0

    # Use KD-tree for efficient water channel detection
    atom_tree = cKDTree(all_atom_coords)

    # Use fewer vectors (4×4 grid instead of 6×6)
    n_theta, n_phi = 4, 4
    channel_vectors = []

    for i in range(n_theta):
        theta = i * np.pi / n_theta
        for j in range(n_phi):
            phi = j * 2 * np.pi / n_phi
            x = np.sin(theta) * np.cos(phi)
            y = np.sin(theta) * np.sin(phi)
            z = np.cos(theta)
            channel_vectors.append(np.array([x, y, z]))

    # Check channels using more efficient algorithm
    water_channels = 0
    water_radius = 1.8  # Typical atom radius + water radius
    check_distances = np.linspace(2.0, 8.0, 5)  # Check 5 points instead of 10

    for vector in channel_vectors:
        blocked = False
        # Calculate all positions at once
        positions = hydroxyl_coord + np.outer(check_distances, vector)

        # Check each position using KD-tree (much faster than nested loops)
        for pos in positions:
            dist, _ = atom_tree.query(pos)
            if dist < water_radius:
                blocked = True
                break

        if not blocked:
            water_channels += 1

    water_channel_score = water_channels / len(channel_vectors)

    return {
        'polar_ratio': polar_ratio,
        'hydrophobic_ratio': hydrophobic_ratio,
        'water_channel_score': water_channel_score,
        'polar_atom_count': polar_atoms
    }


def calculate_domain_boundary_features(residue, chain_data, window_size=7):
    """Calculate features related to domain boundaries and structure discontinuities

    Args:
        residue: The target residue
        chain_data: Pre-computed chain data
        window_size: Window size for local calculations

    Returns:
        dict: Dictionary of domain boundary features
    """
    # Get residue ID and position
    res_id = residue.id[1]

    if res_id not in chain_data['seq_ids']:
        return {
            'sequence_conservation': None,
            'structural_discontinuity': None,
            'hinge_score': None,
            'domain_boundary_score': None
        }

    seq_pos = chain_data['seq_ids'].index(res_id)

    # Get chain
    chain = residue.get_parent()

    # Calculate sequence conservation (based on amino acid similarity in window)
    window_start = max(0, seq_pos - window_size)
    window_end = min(len(chain_data['seq_ids']), seq_pos + window_size + 1)

    window_residues = []
    for i in range(window_start, window_end):
        curr_id = chain_data['seq_ids'][i]
        for res in chain:
            if res.id[1] == curr_id:
                window_residues.append(res)
                break

    # Calculate similarity scores within window
    similarities = []

    for i in range(len(window_residues)):
        for j in range(i + 1, len(window_residues)):
            res1 = window_residues[i]
            res2 = window_residues[j]

            # Get one-letter codes
            aa1 = three_to_one(res1.resname)
            aa2 = three_to_one(res2.resname)

            # Skip non-standard residues
            if aa1 == 'X' or aa2 == 'X':
                continue

            # Calculate similarity score
            similarity = 0

            # Check hydrophobicity
            if aa_properties[aa1]['hydrophobic'] == aa_properties[aa2]['hydrophobic']:
                similarity += 1

            # Check charge
            if aa_properties[aa1]['charged'] == aa_properties[aa2]['charged']:
                similarity += 1

            # Check size
            if aa_properties[aa1]['size'] == aa_properties[aa2]['size']:
                similarity += 1

            # Normalize to 0-1
            similarity /= 3.0
            similarities.append(similarity)

    # Calculate mean similarity
    sequence_conservation = np.mean(similarities) if similarities else None

    # Calculate structural discontinuity
    # Look for shifts in CA positions between residues
    ca_distances = []

    for i in range(len(window_residues) - 1):
        res1 = window_residues[i]
        res2 = window_residues[i + 1]

        # Check if both have CA atoms
        if 'CA' in res1 and 'CA' in res2:
            ca1 = np.array(res1['CA'].coord)
            ca2 = np.array(res2['CA'].coord)

            # Calculate CA-CA distance
            dist = np.sqrt(np.sum((ca1 - ca2) ** 2))
            ca_distances.append(dist)

    # Calculate structural discontinuity (standard deviation of distances)
    structural_discontinuity = np.std(ca_distances) if len(ca_distances) > 1 else None

    # Calculate hinge score (variation in angles between consecutive residues)
    angles = []

    for i in range(len(window_residues) - 2):
        res1 = window_residues[i]
        res2 = window_residues[i + 1]
        res3 = window_residues[i + 2]

        # Check if all have CA atoms
        if 'CA' in res1 and 'CA' in res2 and 'CA' in res3:
            ca1 = np.array(res1['CA'].coord)
            ca2 = np.array(res2['CA'].coord)
            ca3 = np.array(res3['CA'].coord)

            # Calculate vectors
            v1 = ca2 - ca1
            v2 = ca3 - ca2

            # Calculate angle
            dot_product = np.dot(v1, v2)
            norm1 = np.sqrt(np.sum(v1 ** 2))
            norm2 = np.sqrt(np.sum(v2 ** 2))

            if norm1 > 0 and norm2 > 0:
                cos_angle = dot_product / (norm1 * norm2)
                # Clamp to avoid numerical errors
                cos_angle = max(-1.0, min(1.0, cos_angle))
                angle = np.arccos(cos_angle)
                angles.append(angle)

    # Calculate hinge score (standard deviation of angles)
    hinge_score = np.std(angles) if len(angles) > 1 else None

    # Combine metrics for domain boundary score
    if structural_discontinuity is not None and hinge_score is not None:
        domain_boundary_score = (structural_discontinuity + hinge_score) / 2.0
    else:
        domain_boundary_score = structural_discontinuity if structural_discontinuity is not None else hinge_score

    return {
        'sequence_conservation': sequence_conservation,
        'structural_discontinuity': structural_discontinuity,
        'hinge_score': hinge_score,
        'domain_boundary_score': domain_boundary_score
    }


def calculate_bottleneck_metrics(residue, chain_data, neighbors, radius=10.0):
    """Calculate structural bottleneck metrics using residue interaction network analysis

    Args:
        residue: The target residue
        chain_data: Pre-computed chain data
        neighbors: List of neighboring residues
        radius: Contact radius for network building

    Returns:
        dict: Dictionary of bottleneck metrics
    """
    # Get residue ID and sequence position
    res_id = residue.id[1]
    if res_id not in chain_data['seq_ids']:
        return {
            'bottleneck_score': None,
            'edge_betweenness': None,
            'clustering_coefficient': None,
            'topological_centrality': None
        }

    # Build a simplified residue interaction network
    # Using CA-CA distances for contacts
    network_nodes = set()  # Residue IDs
    network_edges = set()  # Pairs of residue IDs

    # Add residues to network
    for res in chain_data.get('all_residues', []):
        if is_aa(res) and 'CA' in res:
            network_nodes.add(res.id[1])

    # Add edges (contacts) to network
    for i, res1 in enumerate(chain_data.get('all_residues', [])):
        if not (is_aa(res1) and 'CA' in res1):
            continue

        res1_id = res1.id[1]
        res1_ca = np.array(res1['CA'].coord)

        for j in range(i + 1, len(chain_data.get('all_residues', []))):
            res2 = chain_data.get('all_residues', [])[j]
            if not (is_aa(res2) and 'CA' in res2):
                continue

            res2_id = res2.id[1]
            res2_ca = np.array(res2['CA'].coord)

            # Calculate CA-CA distance
            dist = np.sqrt(np.sum((res1_ca - res2_ca) ** 2))

            # Add edge if within contact radius
            if dist <= radius:
                network_edges.add((min(res1_id, res2_id), max(res1_id, res2_id)))

    # Calculate bottleneck metrics
    # 1. Edge betweenness (simplified approximation)
    paths_through_residue = 0
    total_paths = 0

    # Sample a subset of paths for efficiency (can't compute all paths for large proteins)
    sample_nodes = random.sample(list(network_nodes), min(50, len(network_nodes)))

    for source in sample_nodes:
        for target in sample_nodes:
            if source == target or source == res_id or target == res_id:
                continue

            # Simple BFS for shortest path
            visited = {source}
            queue = [(source, [source])]
            found_paths = []

            while queue:
                node, path = queue.pop(0)

                if node == target:
                    found_paths.append(path)
                    continue

                # Find neighbors in network
                for e in network_edges:
                    if e[0] == node:
                        neighbor = e[1]
                        if neighbor not in visited:
                            visited.add(neighbor)
                            queue.append((neighbor, path + [neighbor]))
                    elif e[1] == node:
                        neighbor = e[0]
                        if neighbor not in visited:
                            visited.add(neighbor)
                            queue.append((neighbor, path + [neighbor]))

            # Check if paths go through our residue
            total_paths += len(found_paths)
            for path in found_paths:
                if res_id in path:
                    paths_through_residue += 1

    edge_betweenness = paths_through_residue / total_paths if total_paths > 0 else None

    # 2. Clustering coefficient
    # Find all neighbors of the target residue in the network
    target_neighbors = set()
    for e in network_edges:
        if e[0] == res_id:
            target_neighbors.add(e[1])
        elif e[1] == res_id:
            target_neighbors.add(e[0])

    # Count connections between neighbors
    neighbor_connections = 0
    possible_connections = len(target_neighbors) * (len(target_neighbors) - 1) / 2

    for i, n1 in enumerate(target_neighbors):
        for n2 in list(target_neighbors)[i + 1:]:
            if (min(n1, n2), max(n1, n2)) in network_edges:
                neighbor_connections += 1

    clustering_coefficient = neighbor_connections / possible_connections if possible_connections > 0 else None

    # 3. Topological centrality (simplified)
    # Distance to protein centroid
    try:
        ca_atom = residue['CA']
        ca_coord = np.array(ca_atom.coord)

        # Calculate protein centroid
        all_ca_coords = []
        for res in chain_data.get('all_residues', []):
            if is_aa(res) and 'CA' in res:
                all_ca_coords.append(res['CA'].coord)

        if all_ca_coords:
            centroid = np.mean(np.array(all_ca_coords), axis=0)

            # Distance to centroid
            dist_to_centroid = np.sqrt(np.sum((ca_coord - centroid) ** 2))

            # Normalize by maximum distance in protein
            max_dist = 0
            for coord in all_ca_coords:
                dist = np.sqrt(np.sum((np.array(coord) - centroid) ** 2))
                max_dist = max(max_dist, dist)

            topological_centrality = 1.0 - (dist_to_centroid / max_dist) if max_dist > 0 else None
        else:
            topological_centrality = None
    except KeyError:
        topological_centrality = None

    # Calculate bottleneck score (combined metric)
    bottleneck_score = None
    if edge_betweenness is not None and clustering_coefficient is not None:
        bottleneck_score = edge_betweenness * (1.0 - clustering_coefficient)

    return {
        'bottleneck_score': bottleneck_score,
        'edge_betweenness': edge_betweenness,
        'clustering_coefficient': clustering_coefficient,
        'topological_centrality': topological_centrality
    }


def calculate_cavity_features(hydroxyl_atom, neighbors, grid_spacing=1.0, radius=6.0):
    """Calculate cavity features around STY hydroxyl group using a grid-based approach

    Args:
        hydroxyl_atom: The hydroxyl oxygen atom of the STY residue
        neighbors: List of neighboring residues
        grid_spacing: Grid resolution in Ångstroms (higher = faster but less accurate)
        radius: Maximum distance from hydroxyl to consider for cavity in Ångstroms

    Returns:
        dict: Dictionary of cavity metrics
    """
    if not hydroxyl_atom:
        return {
            'cavity_volume': 0.0,
            'pocket_openness': 0.0,
            'cavity_depth': 0.0,
            'cavity_polarity': 0.0,
            'cavity_sphericity': 0.0
        }

    # Get hydroxyl center
    hydroxyl_coord = np.array(hydroxyl_atom.coord)

    # Extract neighbor atom coordinates and VDW radii
    atom_coords = []
    atom_radii = []
    atom_elements = []

    for res in neighbors:
        for atom in res:
            # Skip if this is the hydroxyl atom itself
            if atom == hydroxyl_atom:
                continue

            atom_coords.append(atom.coord)
            element = atom.element.upper()
            atom_radii.append(VDW_RADII.get(element, DEFAULT_RADIUS) + WATER_RADIUS)
            atom_elements.append(element)

    # Handle the case with no neighbors
    if not atom_coords:
        sphere_volume = (4 / 3) * np.pi * radius ** 3
        return {
            'cavity_volume': sphere_volume,
            'pocket_openness': 1.0,
            'cavity_depth': radius,
            'cavity_polarity': 0.0,
            'cavity_sphericity': 1.0
        }

    # Convert to numpy arrays
    atom_coords = np.array(atom_coords)
    atom_radii = np.array(atom_radii)
    atom_elements = np.array(atom_elements)

    # Build a KD-tree for neighbor atoms for faster distance calculations
    atom_tree = cKDTree(atom_coords)

    # Generate fibonacci spiral points - more efficient spherical sampling
    n_points = min(int((4 / 3) * np.pi * (radius ** 3) / (grid_spacing ** 3)), 8000)

    # Generate fibonacci spiral points
    indices = np.arange(n_points)
    phi = (1 + np.sqrt(5)) / 2  # Golden ratio

    # Create points on unit sphere
    z = 1 - (2 * indices) / (n_points - 1)  # Linear from 1 to -1
    radius_xy = np.sqrt(1 - z * z)  # Radius in xy-plane

    theta = 2 * np.pi * indices / phi  # Golden angle increment

    x = radius_xy * np.cos(theta)
    y = radius_xy * np.sin(theta)

    # Create points with random radii to fill sphere
    r = radius * np.cbrt(np.random.random(n_points))  # Cube root for uniform volume distribution

    # Scale to desired radius and position
    grid_points = np.column_stack([x, y, z]) * r.reshape(-1, 1)
    grid_coords = grid_points + hydroxyl_coord

    # Use KD-tree for efficient distance calculations
    is_occupied = np.zeros(len(grid_coords), dtype=bool)

    # For each grid point, find the nearest atom
    dists, idx = atom_tree.query(grid_coords)

    # For each point, check if it's inside any atom
    is_occupied = dists <= atom_radii[idx]

    # Calculate accessible points
    is_accessible = ~is_occupied
    accessible_points = grid_coords[is_accessible]
    accessible_count = len(accessible_points)

    # Calculate volume
    total_points = len(grid_coords)
    sphere_volume = (4 / 3) * np.pi * radius ** 3
    cavity_volume = (accessible_count / total_points) * sphere_volume

    # Only do additional analyses if we have accessible points
    if accessible_count > 0:
        # Calculate pocket openness - shoot rays in different directions
        n_rays = 26  # Fewer rays for speed (26 = ~20° resolution)
        ray_length = radius
        ray_step = radius / 5.0  # Fewer steps along ray

        # Generate fibonacci spiral for uniform ray directions (faster to reuse the same approach)
        indices = np.arange(n_rays)
        z = 1 - (2 * indices) / (n_rays - 1)
        radius_xy = np.sqrt(1 - z * z)
        theta = 2 * np.pi * indices / phi

        x = radius_xy * np.cos(theta)
        y = radius_xy * np.sin(theta)

        ray_directions = np.column_stack([x, y, z])

        # Normalize directions
        norms = np.sqrt(np.sum(ray_directions ** 2, axis=1)).reshape(-1, 1)
        ray_directions = ray_directions / norms

        # Shoot rays, count exits
        exits = 0

        # Create all ray points at once for vectorization
        steps = np.arange(1, 6) * ray_step  # 5 steps per ray

        for direction in ray_directions:
            # Points along this ray
            ray_points = hydroxyl_coord + direction.reshape(1, 3) * steps.reshape(-1, 1)

            # Find nearest atoms for all points at once
            dists, idx = atom_tree.query(ray_points)

            # Check if points are blocked
            is_blocked = dists <= atom_radii[idx]

            # If last point is not blocked, ray found an exit
            if not is_blocked[-1]:
                exits += 1

        # Calculate openness score
        pocket_openness = exits / n_rays

        # Calculate pocket depth - maximum distance from hydroxyl to accessible points
        depths = np.sqrt(np.sum((accessible_points - hydroxyl_coord) ** 2, axis=1))
        pocket_depth = np.max(depths) if len(depths) > 0 else 0.0

        # Calculate cavity polarity - count polar atoms (O, N) in pocket
        polar_atoms = np.sum((atom_elements == 'O') | (atom_elements == 'N'))
        cavity_polarity = polar_atoms / len(atom_elements) if len(atom_elements) > 0 else 0.0

        # Calculate pocket sphericity - how close to a sphere is the cavity
        if len(accessible_points) > 3:
            # Calculate center of mass of accessible points
            com = np.mean(accessible_points, axis=0)

            # Calculate average and std deviation of distances from COM
            distances = np.sqrt(np.sum((accessible_points - com) ** 2, axis=1))
            avg_dist = np.mean(distances)
            std_dist = np.std(distances)

            # Sphericity: lower std deviation = more spherical
            sphericity = 1.0 - min(1.0, std_dist / avg_dist)
        else:
            sphericity = 0.0
    else:
        pocket_openness = 0.0
        pocket_depth = 0.0
        cavity_polarity = 0.0
        sphericity = 0.0

    # Return cavity metrics
    return {
        'cavity_volume': cavity_volume,
        'pocket_openness': pocket_openness,
        'cavity_depth': pocket_depth,
        'cavity_polarity': cavity_polarity,
        'cavity_sphericity': sphericity
    }


def calculate_core_rim_surface(residue, neighbors, model, total_sasa=None, hydroxyl_sasa=None):
    """Classify residue as core, rim, or surface based on exposure and packing

    Args:
        residue: The target residue
        neighbors: List of neighboring residues
        model: BioPython Model object
        total_sasa: Pre-calculated total SASA (optional)
        hydroxyl_sasa: Pre-calculated hydroxyl SASA (optional)

    Returns:
        dict: Dictionary of core-rim-surface metrics
    """
    # Use provided SASA values if available, otherwise calculate them
    if total_sasa is None or hydroxyl_sasa is None:
        total_sasa, hydroxyl_sasa = calculate_sasa_vectorized(residue, neighbors)

    # Get total possible SASA for isolated residue (approximate)
    residue_type = three_to_one(residue.resname)
    max_sasa = {
        'A': 129.0, 'R': 274.0, 'N': 195.0, 'D': 193.0, 'C': 167.0,
        'Q': 225.0, 'E': 223.0, 'G': 104.0, 'H': 224.0, 'I': 197.0,
        'L': 201.0, 'K': 236.0, 'M': 224.0, 'F': 240.0, 'P': 159.0,
        'S': 155.0, 'T': 172.0, 'W': 285.0, 'Y': 263.0, 'V': 174.0
    }

    # Calculate relative SASA
    rel_sasa = total_sasa / max_sasa.get(residue_type, 200.0) if residue_type in max_sasa else None

    # Calculate packing density (vectorized version)
    try:
        # Collect atom coordinates
        residue_atoms = np.array([atom.coord for atom in residue])

        # Count atoms only once from all neighbors
        neighbor_atoms = []
        for neighbor in neighbors:
            for atom in neighbor:
                neighbor_atoms.append(atom.coord)

        if not residue_atoms.size or not neighbor_atoms:
            packing_density = 0
        else:
            # Convert to numpy array
            neighbor_atoms = np.array(neighbor_atoms)

            # Calculate all pairwise distances at once
            # For each residue atom, calculate distance to all neighbor atoms
            packing_density = 0
            cutoff_squared = 5.0 * 5.0  # Square the cutoff for faster comparison

            for res_atom in residue_atoms:
                # Calculate squared distances
                diffs = neighbor_atoms - res_atom
                sq_dists = np.sum(diffs * diffs, axis=1)

                # Count neighbors within cutoff
                packing_density += np.sum(sq_dists <= cutoff_squared)

    except Exception:
        packing_density = 0

    # Normalize packing density by residue size
    atom_count = len(list(residue.get_atoms()))
    normalized_packing = packing_density / atom_count if atom_count > 0 else 0

    # Classify the residue
    classification = None
    if rel_sasa is not None:
        if rel_sasa < 0.1:  # Buried
            classification = 'core'
        elif rel_sasa < 0.4:  # Partially exposed
            classification = 'rim'
        else:  # Highly exposed
            classification = 'surface'

    # Protein centroid calculation (can be optimized if pre-calculated)
    # This part is expensive if model has many chains and residues
    try:
        ca_atom = residue['CA']
        ca_coord = np.array(ca_atom.coord)

        # For efficiency, we can use a model attribute to cache the centroid
        if hasattr(model, '_cached_centroid') and hasattr(model, '_cached_max_dist'):
            # Use cached values
            centroid = model._cached_centroid
            max_dist = model._cached_max_dist
        else:
            # Collect CA coordinates
            all_ca_coords = []
            for chain in model:
                for res in chain:
                    if is_aa(res) and 'CA' in res:
                        all_ca_coords.append(res['CA'].coord)

            if all_ca_coords:
                all_ca_coords = np.array(all_ca_coords)
                centroid = np.mean(all_ca_coords, axis=0)

                # Vectorized calculation of max distance
                diffs = all_ca_coords - centroid
                sq_dists = np.sum(diffs * diffs, axis=1)
                max_dist = np.sqrt(np.max(sq_dists))

                # Cache for future calls
                model._cached_centroid = centroid
                model._cached_max_dist = max_dist
            else:
                centroid = None
                max_dist = None

        # Calculate distance to centroid
        if centroid is not None and max_dist is not None:
            dist_to_centroid = np.sqrt(np.sum((ca_coord - centroid) ** 2))
            normalized_centrality = dist_to_centroid / max_dist if max_dist > 0 else None
        else:
            dist_to_centroid = None
            normalized_centrality = None

    except KeyError:
        dist_to_centroid = None
        normalized_centrality = None

    return {
        'core_rim_surface': classification,
        'relative_sasa': rel_sasa,
        'normalized_packing': normalized_packing,
        'dist_to_centroid': dist_to_centroid,
        'normalized_centrality': normalized_centrality
    }


def calculate_disorder_features(residue, chain_data, plddt_dict=None, window=7):
    """Calculate disorder-related features using a simplified predictor

    Args:
        residue: The target residue
        chain_data: Pre-computed chain data
        plddt_dict: Optional dictionary of pLDDT values (AlphaFold)
        window: Size of sequence window for calculations

    Returns:
        dict: Dictionary of disorder-related features
    """
    # Get residue ID and sequence position
    res_id = residue.id[1]
    if res_id not in chain_data['seq_ids']:
        return {
            'disorder_score': None,
            'disorder_gradient': None,
            'order_disorder_boundary': None,
            'binding_induced_folding': None
        }

    seq_pos = chain_data['seq_ids'].index(res_id)
    chain = residue.get_parent()
    sequence = chain_data['sequence']

    # 1. Simple sequence-based disorder propensity
    # Disorder propensities from various scales (higher = more disordered)
    disorder_propensity = {
        'A': 0.06, 'R': 0.18, 'N': 0.23, 'D': 0.29, 'C': -0.26,
        'Q': 0.32, 'E': 0.27, 'G': 0.16, 'H': 0.23, 'I': -0.66,
        'L': -0.57, 'K': 0.26, 'M': -0.40, 'F': -0.61, 'P': 0.17,
        'S': 0.13, 'T': 0.05, 'W': -0.42, 'Y': -0.20, 'V': -0.50
    }

    # Calculate window sequence disorder propensity
    window_start = max(0, seq_pos - window)
    window_end = min(len(sequence), seq_pos + window + 1)
    window_seq = sequence[window_start:window_end]

    window_disorder = sum(disorder_propensity.get(aa, 0) for aa in window_seq) / len(window_seq) if window_seq else None

    # 2. pLDDT-based disorder (if available)
    plddt_disorder = None
    plddt_gradient = None

    if plddt_dict:
        # Get pLDDT scores in window
        window_plddt = []
        for i in range(window_start, window_end):
            if i < len(chain_data['seq_ids']):
                win_res_id = chain_data['seq_ids'][i]
                if win_res_id in plddt_dict:
                    window_plddt.append((i - seq_pos, plddt_dict[win_res_id]))

        # Calculate pLDDT-based disorder score
        # Lower pLDDT = higher disorder
        if window_plddt:
            center_plddt = plddt_dict.get(res_id, 100)
            plddt_disorder = (100 - center_plddt) / 100  # 0-1 scale, higher = more disordered

            # Calculate pLDDT gradient
            if len(window_plddt) > 1:
                # Linear regression for gradient
                x = np.array([pos for pos, _ in window_plddt])
                y = np.array([score for _, score in window_plddt])

                # Calculate slope
                if np.var(x) > 0:
                    slope = np.cov(x, y)[0, 1] / np.var(x)
                    plddt_gradient = slope

    # 3. Structure-based disorder metrics
    # B-factor/flexibility
    b_factor_disorder = None
    try:
        b_factors = [atom.bfactor for atom in residue]
        avg_b = np.mean(b_factors) if b_factors else None

        # Normalize B-factor (approximate)
        if avg_b is not None:
            b_factor_disorder = min(1.0, avg_b / 100.0)  # 0-1 scale, higher = more disordered
    except Exception:
        pass

    # 4. Combined disorder score (weighted average)
    disorder_score = None
    components = []

    if window_disorder is not None:
        components.append((window_disorder, 0.3))  # 30% weight
    if plddt_disorder is not None:
        components.append((plddt_disorder, 0.5))  # 50% weight
    if b_factor_disorder is not None:
        components.append((b_factor_disorder, 0.2))  # 20% weight

    if components:
        weighted_sum = sum(score * weight for score, weight in components)
        total_weight = sum(weight for _, weight in components)
        disorder_score = weighted_sum / total_weight if total_weight > 0 else None

    # 5. Order-disorder boundary detection
    order_disorder_boundary = None

    if plddt_dict and len(sequence) > 2 * window:
        # Calculate disorder scores along the sequence
        seq_disorder = []
        for i in range(len(chain_data['seq_ids'])):
            curr_id = chain_data['seq_ids'][i]
            curr_plddt = plddt_dict.get(curr_id, 100)
            curr_disorder = (100 - curr_plddt) / 100
            seq_disorder.append(curr_disorder)

        # Look for sharp transitions (order-disorder boundaries)
        if seq_pos > window and seq_pos < len(seq_disorder) - window:
            # Calculate average disorder before and after current position
            before_avg = np.mean(seq_disorder[seq_pos - window:seq_pos])
            after_avg = np.mean(seq_disorder[seq_pos:seq_pos + window])

            # Calculate boundary score (higher = sharper transition)
            boundary_score = abs(after_avg - before_avg)

            # Threshold for boundary detection
            order_disorder_boundary = boundary_score if boundary_score > 0.2 else 0.0

    # 6. Binding-induced folding potential
    # Simple estimate based on disorder score and presence of linear motifs
    binding_induced_folding = None

    if disorder_score is not None:
        # Patterns for known binding motifs (simplified)
        motif_patterns = [
            r'[ILMV]..[ILMV]',  # Hydrophobic anchors
            r'P..P',  # Proline-rich regions
            r'[ST]P..[KR]',  # CDK phosphorylation sites
            r'[RK][RK][RK][RK]',  # NLS-like
            r'L..[LI]'  # Leucine-rich
        ]

        # Calculate motif matches
        motif_match = False
        for pattern in motif_patterns:
            if re.search(pattern, window_seq):
                motif_match = True
                break

        # Higher disorder + motif presence = higher binding potential
        if motif_match:
            binding_induced_folding = disorder_score
        else:
            binding_induced_folding = disorder_score * 0.5

    return {
        'disorder_score': disorder_score,
        'disorder_gradient': plddt_gradient,
        'order_disorder_boundary': order_disorder_boundary,
        'binding_induced_folding': binding_induced_folding
    }


def calculate_plddt_gradient_features(residue, chain_data, plddt_dict, window=7):
    """Calculate detailed pLDDT gradient features

    Args:
        residue: The target residue
        chain_data: Pre-computed chain data
        plddt_dict: Dictionary of pLDDT values (AlphaFold)
        window: Size of sequence window for calculations

    Returns:
        dict: Dictionary of pLDDT gradient features
    """
    # Get residue ID and sequence position
    res_id = residue.id[1]
    if res_id not in chain_data['seq_ids'] or not plddt_dict:
        return {
            'plddt_gradient_magnitude': None,
            'plddt_gradient_direction': None,
            'plddt_local_variance': None,
            'plddt_relative_to_domain': None,
            'plddt_transition_strength': None
        }

    seq_pos = chain_data['seq_ids'].index(res_id)

    # Get local pLDDT values in sequence window
    window_start = max(0, seq_pos - window)
    window_end = min(len(chain_data['seq_ids']), seq_pos + window + 1)

    window_plddt = []
    for i in range(window_start, window_end):
        win_res_id = chain_data['seq_ids'][i]
        if win_res_id in plddt_dict:
            window_plddt.append((i - seq_pos, plddt_dict[win_res_id]))

    if not window_plddt:
        return {
            'plddt_gradient_magnitude': None,
            'plddt_gradient_direction': None,
            'plddt_local_variance': None,
            'plddt_relative_to_domain': None,
            'plddt_transition_strength': None
        }

    # Calculate gradient using linear regression
    x = np.array([pos for pos, _ in window_plddt])
    y = np.array([score for _, score in window_plddt])

    # Calculate gradient magnitude and direction
    if len(window_plddt) > 2 and np.var(x) > 0:
        slope = np.cov(x, y)[0, 1] / np.var(x)
        gradient_magnitude = abs(slope)
        gradient_direction = np.sign(slope)
    else:
        gradient_magnitude = None
        gradient_direction = None

    # Calculate local variance
    local_variance = np.var(y) if len(y) > 1 else None

    # Calculate relative to domain average
    # Estimate domain as ±30 residues
    domain_start = max(0, seq_pos - 30)
    domain_end = min(len(chain_data['seq_ids']), seq_pos + 31)

    domain_plddt = []
    for i in range(domain_start, domain_end):
        if i < len(chain_data['seq_ids']):
            dom_res_id = chain_data['seq_ids'][i]
            if dom_res_id in plddt_dict:
                domain_plddt.append(plddt_dict[dom_res_id])

    domain_avg = np.mean(domain_plddt) if domain_plddt else None
    residue_plddt = plddt_dict.get(res_id, None)

    relative_to_domain = None
    if domain_avg is not None and residue_plddt is not None:
        relative_to_domain = residue_plddt - domain_avg

    # Calculate transition strength
    # Look for abrupt changes in pLDDT values
    transition_strength = None
    if len(window_plddt) > 3:
        # Calculate maximum absolute difference between consecutive values
        consecutive_diffs = []
        positions = sorted(window_plddt, key=lambda x: x[0])

        for i in range(len(positions) - 1):
            diff = abs(positions[i + 1][1] - positions[i][1])
            consecutive_diffs.append(diff)

        transition_strength = max(consecutive_diffs) if consecutive_diffs else None

    return {
        'plddt_gradient_magnitude': gradient_magnitude,
        'plddt_gradient_direction': gradient_direction,
        'plddt_local_variance': local_variance,
        'plddt_relative_to_domain': relative_to_domain,
        'plddt_transition_strength': transition_strength
    }


def calculate_hydration_stability(residue, neighbors, grid_spacing=0.7, radius=5.0):
    """Calculate hydration stability metrics around STY hydroxyl

    Args:
        residue: The target residue
        neighbors: List of neighboring residues
        grid_spacing: Grid resolution for hydration analysis
        radius: Radius around hydroxyl to analyze

    Returns:
        dict: Dictionary of hydration stability metrics
    """
    # Get hydroxyl oxygen atom
    hydroxyl_atom = None
    if residue.resname == 'SER' and 'OG' in residue:
        hydroxyl_atom = residue['OG']
    elif residue.resname == 'THR' and 'OG1' in residue:
        hydroxyl_atom = residue['OG1']
    elif residue.resname == 'TYR' and 'OH' in residue:
        hydroxyl_atom = residue['OH']

    if not hydroxyl_atom:
        return {
            'hydration_sites': 0,
            'hydration_stability': None,
            'hydration_residence_score': None,
            'bridged_hbond_network': None
        }

    # Get hydroxyl coordinates
    hydroxyl_coord = np.array(hydroxyl_atom.coord)

    # Extract all polar atoms (potential H-bonding partners)
    polar_atoms = []
    for neighbor in neighbors:
        for atom in neighbor:
            if atom.element in ['O', 'N']:
                polar_atoms.append(atom.coord)

    # Generate grid points in spherical shell around hydroxyl
    # (where water molecules are likely to be)
    water_shell_min = 2.5  # Minimum water distance
    water_shell_max = 3.5  # Maximum water distance

    # Create uniform grid around hydroxyl
    grid_margin = water_shell_max + 0.5  # Add small margin
    nx = int(2 * grid_margin / grid_spacing) + 1
    ny = int(2 * grid_margin / grid_spacing) + 1
    nz = int(2 * grid_margin / grid_spacing) + 1

    grid_points = []
    for i in range(nx):
        x = (i * grid_spacing) - grid_margin + hydroxyl_coord[0]
        for j in range(ny):
            y = (j * grid_spacing) - grid_margin + hydroxyl_coord[1]
            for k in range(nz):
                z = (k * grid_spacing) - grid_margin + hydroxyl_coord[2]

                point = np.array([x, y, z])

                # Calculate distance to hydroxyl
                dist_to_hydroxyl = np.sqrt(np.sum((point - hydroxyl_coord) ** 2))

                # Check if within water shell
                if water_shell_min <= dist_to_hydroxyl <= water_shell_max:
                    grid_points.append(point)

    # Score each potential hydration site
    hydration_sites = []

    for point in grid_points:
        # Check if this point is too close to any atom in the protein
        too_close = False
        for neighbor in neighbors:
            for atom in neighbor:
                dist = np.sqrt(np.sum((point - atom.coord) ** 2))
                if dist < 1.6:  # Minimum atom-water distance
                    too_close = True
                    break
            if too_close:
                break

        if too_close:
            continue

        # Calculate score based on hydrogen bonding potential
        h_bond_score = 0

        # Check hydroxyl as H-bond partner
        dist_to_hydroxyl = np.sqrt(np.sum((point - hydroxyl_coord) ** 2))
        angle_ideal = 3.14159 / 4  # ~45 degrees - water prefers tetrahedral geometry
        h_bond_score += 1.0 * np.exp(-(dist_to_hydroxyl - 2.8) ** 2)  # Distance component

        # Check other polar atoms as H-bond partners
        for polar_atom_coord in polar_atoms:
            dist = np.sqrt(np.sum((point - polar_atom_coord) ** 2))

            if 2.5 <= dist <= 3.5:  # Typical H-bond distance range
                # Calculate angle between hydroxyl-water-polar_atom
                v1 = hydroxyl_coord - point
                v2 = polar_atom_coord - point

                # Normalize vectors
                v1_norm = np.sqrt(np.sum(v1 ** 2))
                v2_norm = np.sqrt(np.sum(v2 ** 2))

                if v1_norm > 0 and v2_norm > 0:
                    v1 = v1 / v1_norm
                    v2 = v2 / v2_norm

                    # Calculate cosine of angle
                    cos_angle = np.dot(v1, v2)

                    # Ensure numerical stability
                    cos_angle = max(-1.0, min(1.0, cos_angle))

                    # Calculate angle in radians
                    angle = np.arccos(cos_angle)

                    # Score based on tetrahedral angles
                    # Perfect tetrahedral angle is ~109.5 degrees (1.91 radians)
                    angle_score = np.exp(-(angle - 1.91) ** 2 / 0.3)

                    h_bond_score += 0.5 * angle_score  # Add to score (with lower weight)

        # If good hydration site, add to list
        if h_bond_score > 0.5:
            hydration_sites.append((point, h_bond_score))

    # Count hydration sites and calculate stability
    num_sites = len(hydration_sites)

    # Calculate average stability score
    hydration_stability = None
    if num_sites > 0:
        hydration_stability = sum(score for _, score in hydration_sites) / num_sites

    # Calculate residence score (estimate of water residence time)
    # Higher score = longer residence time
    hydration_residence_score = None
    if num_sites > 0:
        # Get maximum score (most stable site)
        max_score = max(score for _, score in hydration_sites)
        hydration_residence_score = max_score

    # Calculate bridged H-bond network score
    # Look for cases where water can bridge between hydroxyl and other polar groups
    bridged_hbond_network = None
    if num_sites > 0:
        bridge_count = 0

        for _, site_score in hydration_sites:
            if site_score > 1.0:  # Only consider stable sites
                bridge_count += 1

        bridged_hbond_network = bridge_count / 3.0 if bridge_count < 3 else 1.0  # Normalize

    return {
        'hydration_sites': num_sites,
        'hydration_stability': hydration_stability,
        'hydration_residence_score': hydration_residence_score,
        'bridged_hbond_network': bridged_hbond_network
    }


def calculate_motif_plddt(chain_data, residue, plddt_dict, window=7):
    """Calculate mean pLDDT for the sequence motif around a residue

    Args:
        chain_data: Pre-computed chain data
        residue: The target residue
        plddt_dict: Dictionary mapping residue IDs to pLDDT values
        window: Size of the sequence window on each side (default: 7)

    Returns:
        float: Mean pLDDT of residues in the motif
    """
    res_id = residue.id[1]

    # Get the sequence IDs from chain data
    seq_ids = chain_data['seq_ids']

    # Find the position in the sequence
    if res_id not in seq_ids:
        return None

    seq_pos = seq_ids.index(res_id)

    # Get residue IDs in the motif window
    start = max(0, seq_pos - window)
    end = min(len(seq_ids), seq_pos + window + 1)

    motif_res_ids = seq_ids[start:end]

    # Get pLDDT values for each residue in the motif
    motif_plddt = [plddt_dict.get(res_id, 0) for res_id in motif_res_ids]

    # Remove missing values (0)
    motif_plddt = [p for p in motif_plddt if p > 0]

    if not motif_plddt:
        return None

    return np.mean(motif_plddt)


def distance_to_nearest_vectorized(residue, chain_data):
    """Calculate distance to the nearest residue of each type and property category

    Args:
        residue: The target residue
        chain_data: Pre-computed chain data with KD-tree

    Returns:
        dict: Distances to nearest residue of each type and property category
    """
    # Get CA atom coordinates for the central residue
    try:
        ca_atom = residue['CA']
    except KeyError:
        # Return default values if CA atom is missing
        return {f"dist_to_{aa}": None for aa in "ACDEFGHIKLMNPQRSTVWY"}, {
            "dist_to_hydrophobic": None,
            "dist_to_polar": None,
            "dist_to_charged": None,
            "dist_to_acidic": None,
            "dist_to_basic": None,
            "dist_to_small": None,
            "dist_to_medium": None,
            "dist_to_large": None
        }

    ca_coord = np.array(ca_atom.coord)

    # Get distances to all residues in the chain
    distances = {}
    property_distances = {
        "hydrophobic": [],
        "polar": [],
        "charged": [],
        "acidic": [],
        "basic": [],
        "small": [],
        "medium": [],
        "large": []
    }

    # Collect all residue coordinates by type
    residue_coords = {aa: [] for aa in "ACDEFGHIKLMNPQRSTVWY"}

    # Iterate through all residues in the chain
    for res in chain_data.get('all_residues', []):
        # Skip the target residue
        if res.id == residue.id and res.get_parent() == residue.get_parent():
            continue

        # Skip non-amino acids
        if not is_aa(res):
            continue

        # Get residue type
        res_type = three_to_one(res.resname)
        if res_type == 'X':
            continue

        # Get CA coordinates
        try:
            res_ca = res['CA']
            coords = res_ca.coord

            # Add to appropriate lists
            residue_coords[res_type].append(coords)

            # Add to property-based lists
            props = aa_properties.get(res_type, {})
            if props.get('hydrophobic', False):
                property_distances["hydrophobic"].append(coords)
            if props.get('polar', False):
                property_distances["polar"].append(coords)
            if props.get('charged', False):
                property_distances["charged"].append(coords)

                # Acidic (negative) or basic (positive)
                charge = props.get('charge', 0)
                if charge < 0:
                    property_distances["acidic"].append(coords)
                elif charge > 0:
                    property_distances["basic"].append(coords)

            # Size properties
            size = props.get('size', '')
            if size == 'small':
                property_distances["small"].append(coords)
            elif size == 'medium':
                property_distances["medium"].append(coords)
            elif size == 'large':
                property_distances["large"].append(coords)

        except KeyError:
            # Skip residues without CA atoms
            continue

    # Calculate distances for each amino acid type
    aa_distances = {}
    for aa, coords_list in residue_coords.items():
        if coords_list:
            # Convert to numpy array
            coords_array = np.array(coords_list)

            # Calculate Euclidean distances
            vectors = coords_array - ca_coord
            dist_sq = np.sum(vectors * vectors, axis=1)
            min_dist = np.sqrt(np.min(dist_sq)) if len(dist_sq) > 0 else None

            aa_distances[f"dist_to_{aa}"] = min_dist
        else:
            aa_distances[f"dist_to_{aa}"] = None

    # Calculate distances for each property category
    property_results = {}
    for prop, coords_list in property_distances.items():
        if coords_list:
            # Convert to numpy array
            coords_array = np.array(coords_list)

            # Calculate Euclidean distances
            vectors = coords_array - ca_coord
            dist_sq = np.sum(vectors * vectors, axis=1)
            min_dist = np.sqrt(np.min(dist_sq)) if len(dist_sq) > 0 else None

            property_results[f"dist_to_{prop}"] = min_dist
        else:
            property_results[f"dist_to_{prop}"] = None

    return aa_distances, property_results


def num_nearest_type_vectorized(residue, neighbors):
    """Count the number of each amino acid type and property within 10Å

    Args:
        residue: The target residue
        neighbors: List of neighboring residues within 10Å

    Returns:
        dict: Counts of each amino acid type and property category
    """
    # Initialize counts for each amino acid
    aa_counts = {f"num_{aa}_10A": 0 for aa in "ACDEFGHIKLMNPQRSTVWY"}

    # Initialize counts for property categories
    property_counts = {
        "num_hydrophobic_10A": 0,
        "num_polar_10A": 0,
        "num_charged_10A": 0,
        "num_acidic_10A": 0,
        "num_basic_10A": 0,
        "num_small_10A": 0,
        "num_medium_10A": 0,
        "num_large_10A": 0
    }

    # Count residues by type and property
    for res in neighbors:
        # Skip non-amino acids
        if not is_aa(res):
            continue

        # Get residue type
        res_type = three_to_one(res.resname)
        if res_type == 'X':
            continue

        # Increment type count
        aa_counts[f"num_{res_type}_10A"] += 1

        # Increment property counts
        props = aa_properties.get(res_type, {})
        if props.get('hydrophobic', False):
            property_counts["num_hydrophobic_10A"] += 1
        if props.get('polar', False):
            property_counts["num_polar_10A"] += 1
        if props.get('charged', False):
            property_counts["num_charged_10A"] += 1

            # Acidic (negative) or basic (positive)
            charge = props.get('charge', 0)
            if charge < 0:
                property_counts["num_acidic_10A"] += 1
            elif charge > 0:
                property_counts["num_basic_10A"] += 1

        # Size properties
        size = props.get('size', '')
        if size == 'small':
            property_counts["num_small_10A"] += 1
        elif size == 'medium':
            property_counts["num_medium_10A"] += 1
        elif size == 'large':
            property_counts["num_large_10A"] += 1

    return aa_counts, property_counts


def process_structure_file(structure_file, total_files, file_index):
    """Process a single structure file with timing for each major function"""
    import time
    from collections import defaultdict

    results = []
    # Dictionary to store timing data for individual functions
    function_timings = {}
    # Dictionary to store timing data for function categories
    category_timings = defaultdict(float)

    total_start_time = time.time()

    # Extract UniProt ID from filename
    filename = os.path.basename(structure_file)
    if 'AF-' in filename and '-F' in filename:
        uniprot_id = filename.split('AF-')[1].split('-F')[0]
    else:
        uniprot_id = filename.split('.')[0]  # Fallback

    print(f"[{file_index}/{total_files}] Processing {filename}...")

    # Determine file type and select appropriate parser
    if filename.endswith('.cif.gz') or filename.endswith('.cif'):
        parser = MMCIFParser(QUIET=True)
        is_cif = True
    else:
        parser = PDBParser(QUIET=True)
        is_cif = False

    total_sty_found = 0

    try:
        # ------ FILE PARSING PHASE ------
        parsing_start = time.time()

        # Extract secondary structure information
        ss_dict = {}
        if is_cif:
            func_start = time.time()
            ss_dict = extract_secondary_structure(structure_file)
            function_timings['extract_secondary_structure'] = time.time() - func_start

        # Handle gzipped files and parse structure
        func_start = time.time()
        if structure_file.endswith('.gz'):
            with gzip.open(structure_file, 'rt') as f:
                temp_file = f"{structure_file}.temp"
                with open(temp_file, 'w') as temp:
                    temp.write(f.read())
                structure = parser.get_structure(uniprot_id, temp_file)
                os.remove(temp_file)  # Clean up
        else:
            structure = parser.get_structure(uniprot_id, structure_file)
        function_timings['parse_structure'] = time.time() - func_start

        model = structure[0]

        # Get pLDDT scores
        func_start = time.time()
        plddt_dict = parse_plddt_from_structure_file(structure_file)
        function_timings['parse_plddt'] = time.time() - func_start

        # Cache secondary structures
        func_start = time.time()
        ss_data = cache_all_secondary_structures(model, ss_dict)
        function_timings['cache_secondary_structures'] = time.time() - func_start

        category_timings['file_parsing'] = time.time() - parsing_start

        # Process each chain
        for chain_idx, chain in enumerate(model):
            chain_id = chain.id
            print(f"  Processing chain {chain_id} ({chain_idx + 1}/{len(model)} chains)...")

            # ------ CHAIN SETUP PHASE ------
            chain_setup_start = time.time()

            # Build chain data
            func_start = time.time()
            chain_data = build_chain_data(chain)
            function_timings['build_chain_data'] = time.time() - func_start

            # Add residues to chain_data
            chain_data['all_residues'] = list(chain.get_residues())

            # Get STY residues
            sty_residues = chain_data['sty_residues']

            if not sty_residues:
                print(f"  No S/T/Y residues found in chain {chain_id}")
                continue

            print(f"  Found {len(sty_residues)} S/T/Y residues in chain {chain_id}")

            # Calculate HSE data
            func_start = time.time()
            hse_data = calculate_hse_fast(model, sty_residues, chain_data)
            function_timings['calculate_hse'] = time.time() - func_start

            category_timings['chain_setup'] = time.time() - chain_setup_start

            # Process each STY residue
            for i, residue in enumerate(sty_residues):
                if (i + 1) % 100 == 0:
                    print(
                        f"    Processed {i + 1}/{len(sty_residues)} STY residues in chain {chain_id} (current: {residue.id[1]})")

                res_id = residue.id[1]
                res_code = three_to_one(residue.resname)

                # Get pLDDT score
                plddt = plddt_dict.get(res_id, None)

                # ------ NEIGHBOR IDENTIFICATION PHASE ------
                neighbor_start = time.time()

                # Find neighbors using KD-tree
                func_start = time.time()
                neighbors = find_neighbors_fast(chain_data, residue)
                function_timings['find_neighbors'] = function_timings.get('find_neighbors', 0) + (
                            time.time() - func_start)

                category_timings['neighbor_identification'] += time.time() - neighbor_start

                # ------ SEQUENCE FEATURES PHASE ------
                sequence_start = time.time()

                # Get sequence motif
                func_start = time.time()
                motif = get_sequence_motif(chain, residue)
                function_timings['get_sequence_motif'] = function_timings.get('get_sequence_motif', 0) + (
                            time.time() - func_start)

                # Calculate motif pLDDT
                func_start = time.time()
                motif_plddt = calculate_motif_plddt(chain_data, residue, plddt_dict)
                function_timings['calculate_motif_plddt'] = function_timings.get('calculate_motif_plddt', 0) + (
                            time.time() - func_start)

                category_timings['sequence_features'] += time.time() - sequence_start

                # ------ HYDROXYL IDENTIFICATION PHASE ------
                hydroxyl_start = time.time()

                # Get hydroxyl oxygen atom
                hydroxyl_atom = None
                if residue.resname == 'SER' and 'OG' in residue:
                    hydroxyl_atom = residue['OG']
                elif residue.resname == 'THR' and 'OG1' in residue:
                    hydroxyl_atom = residue['OG1']
                elif residue.resname == 'TYR' and 'OH' in residue:
                    hydroxyl_atom = residue['OH']

                category_timings['hydroxyl_identification'] += time.time() - hydroxyl_start

                # ------ EXPOSURE CALCULATION PHASE ------
                exposure_start = time.time()

                # Calculate hydroxyl exposure
                func_start = time.time()
                hydroxyl_exposure, backbone_contacts = calculate_hydroxyl_exposure_vectorized(hydroxyl_atom, neighbors)
                function_timings['calculate_hydroxyl_exposure'] = function_timings.get('calculate_hydroxyl_exposure',
                                                                                       0) + (time.time() - func_start)

                # Calculate SASA
                func_start = time.time()
                total_sasa, hydroxyl_sasa = calculate_sasa_vectorized(residue, neighbors)
                function_timings['calculate_sasa'] = function_timings.get('calculate_sasa', 0) + (
                            time.time() - func_start)

                # Calculate residue depth
                func_start = time.time()
                depth_results, min_depth, avg_depth, ca_depth, hydroxyl_depth = calculate_residue_depth_multi(residue,
                                                                                                              model)
                function_timings['calculate_residue_depth'] = function_timings.get('calculate_residue_depth', 0) + (
                            time.time() - func_start)

                category_timings['exposure_calculations'] += time.time() - exposure_start

                # ------ POCKET ANALYSIS PHASE ------
                pocket_start = time.time()

                # Analyze pocket properties
                func_start = time.time()
                pocket_props = analyze_pocket_properties(neighbors)
                function_timings['analyze_pocket_properties'] = function_timings.get('analyze_pocket_properties', 0) + (
                            time.time() - func_start)

                # Analyze pocket shape
                func_start = time.time()
                pocket_shape = analyze_pocket_shape(residue, neighbors)
                function_timings['analyze_pocket_shape'] = function_timings.get('analyze_pocket_shape', 0) + (
                            time.time() - func_start)

                # Calculate cavity volume
                func_start = time.time()
                cavity_features = calculate_cavity_features(hydroxyl_atom, neighbors)
                function_timings['calculate_cavity_features'] = function_timings.get('calculate_cavity_features', 0) + (
                            time.time() - func_start)

                # Calculate pocket solvent features
                func_start = time.time()
                solvent_features = calculate_pocket_solvent_features(residue, neighbors)
                function_timings['calculate_pocket_solvent'] = function_timings.get('calculate_pocket_solvent', 0) + (
                            time.time() - func_start)

                category_timings['pocket_analysis'] += time.time() - pocket_start

                # ------ HYDROGEN BOND CALCULATION PHASE ------
                hbond_start = time.time()

                # Find hydrogen bonds
                func_start = time.time()
                h_bonds = find_hydrogen_bonds_vectorized(hydroxyl_atom, neighbors)
                function_timings['find_hydrogen_bonds'] = function_timings.get('find_hydrogen_bonds', 0) + (
                            time.time() - func_start)

                # Calculate hydrogen bond energy
                func_start = time.time()
                hbond_energy_data = calculate_hydrogen_bond_energy(residue, neighbors)
                function_timings['calculate_hydrogen_bond_energy'] = function_timings.get(
                    'calculate_hydrogen_bond_energy', 0) + (time.time() - func_start)

                category_timings['hydrogen_bond_calculations'] += time.time() - hbond_start

                # ------ SECONDARY STRUCTURE PHASE ------
                ss_start = time.time()

                # Get secondary structure
                ss_key = (chain_id, res_id)
                secondary_structure = ss_dict.get(ss_key, 'Unknown')

                # Map secondary structure codes
                ss_descriptions = {
                    'HELX_LH_PP_P': 'Alpha Helix',
                    'HELX_RH_PP_P': 'Right-handed Helix',
                    'HELX_RH_3T_P': '3/10 Helix',
                    'HELX_RH_PI_P': 'Pi Helix',
                    'HELX_LH_P': 'Left-handed Helix',
                    'STRN': 'Beta Strand',
                    'TURN_TY1_P': 'Type I Turn',
                    'TURN_TY2_P': 'Type II Turn',
                    'TURN_TY1_PM': 'Type I Prime Turn',
                    'TURN_TY2_PM': 'Type II Prime Turn',
                    'TURN_TY3_P': 'Type III Turn',
                    'BEND': 'Bend',
                    'SHEET': 'Beta Sheet',
                    'COIL': 'Coil'
                }
                ss_description = ss_descriptions.get(secondary_structure, secondary_structure)

                # Calculate secondary structure distances
                func_start = time.time()
                ss_distances = calculate_secondary_structure_distances(residue, chain_data['seq_ids'], ss_data)
                function_timings['calculate_ss_distances'] = function_timings.get('calculate_ss_distances', 0) + (
                            time.time() - func_start)

                category_timings['secondary_structure_analysis'] += time.time() - ss_start

                # ------ SIDE CHAIN ANALYSIS PHASE ------
                side_chain_start = time.time()

                # Calculate side chain angles
                func_start = time.time()
                side_chain_angles = calculate_side_chain_angles(residue)
                function_timings['calculate_side_chain_angles'] = function_timings.get('calculate_side_chain_angles',
                                                                                       0) + (time.time() - func_start)

                category_timings['side_chain_analysis'] += time.time() - side_chain_start

                # ------ DENSITY AND PACKING PHASE ------
                packing_start = time.time()

                # Calculate packing density
                func_start = time.time()
                packing_density = calculate_packing_density(residue, neighbors)
                function_timings['calculate_packing_density'] = function_timings.get('calculate_packing_density', 0) + (
                            time.time() - func_start)

                category_timings['density_packing_calculations'] += time.time() - packing_start

                # ------ RESIDUE CONTACTS PHASE ------
                contacts_start = time.time()

                # Calculate contact order
                func_start = time.time()
                contact_order_data = calculate_contact_order(residue, chain_data, neighbors)
                function_timings['calculate_contact_order'] = function_timings.get('calculate_contact_order', 0) + (
                            time.time() - func_start)

                # Calculate distance to nearest residue types
                func_start = time.time()
                nearest_aa_distances, nearest_property_distances = distance_to_nearest_vectorized(residue, chain_data)
                function_timings['calculate_nearest_distances'] = function_timings.get('calculate_nearest_distances',
                                                                                       0) + (time.time() - func_start)

                # Count nearby residues by type
                func_start = time.time()
                aa_counts, property_counts = num_nearest_type_vectorized(residue, neighbors)
                function_timings['count_nearest_residues'] = function_timings.get('count_nearest_residues', 0) + (
                            time.time() - func_start)

                category_timings['residue_contact_calculations'] += time.time() - contacts_start

                # ------ FLEXIBILITY ANALYSIS PHASE ------
                flexibility_start = time.time()

                # Calculate B-factor statistics
                func_start = time.time()
                bfactor_stats = calculate_b_factor_statistics(residue, neighbors, plddt_dict)
                function_timings['calculate_bfactor_stats'] = function_timings.get('calculate_bfactor_stats', 0) + (
                            time.time() - func_start)

                # Calculate local structural features
                func_start = time.time()
                local_structure = calculate_local_structural_features(residue, chain_data)
                function_timings['calculate_local_structure'] = function_timings.get('calculate_local_structure', 0) + (
                            time.time() - func_start)

                # Calculate domain boundary features
                func_start = time.time()
                domain_features = calculate_domain_boundary_features(residue, chain_data)
                function_timings['calculate_domain_boundary'] = function_timings.get('calculate_domain_boundary', 0) + (
                            time.time() - func_start)

                category_timings['flexibility_analysis'] += time.time() - flexibility_start

                # Calculate disorder features
                func_start = time.time()
                disorder_features = calculate_disorder_features(residue, chain_data, plddt_dict)
                function_timings['calculate_disorder_features'] = function_timings.get('calculate_disorder_features',
                                                                                       0) + (time.time() - func_start)

                # Calculate pLDDT gradient features
                func_start = time.time()
                plddt_gradient_features = calculate_plddt_gradient_features(residue, chain_data, plddt_dict)
                function_timings['calculate_plddt_gradient'] = function_timings.get('calculate_plddt_gradient', 0) + (
                            time.time() - func_start)

                # Calculate core-rim-surface features
                func_start = time.time()
                core_rim_features = calculate_core_rim_surface(residue, neighbors, model, total_sasa, hydroxyl_sasa)
                function_timings['calculate_core_rim_surface'] = function_timings.get('calculate_core_rim_surface',
                                                                                      0) + (time.time() - func_start)

                # Calculate bottleneck metrics
                # func_start = time.time()
                # bottleneck_metrics = calculate_bottleneck_metrics(residue, chain_data, neighbors)
                # function_timings['calculate_bottleneck'] = function_timings.get('calculate_bottleneck', 0) + (time.time() - func_start)

                # Calculate hydration stability
                # func_start = time.time()
                # hydration_features = calculate_hydration_stability(residue, neighbors)
                # function_timings['calculate_hydration'] = function_timings.get('calculate_hydration', 0) + (time.time() - func_start)

                # ------ CHEMICAL ENVIRONMENT PHASE ------
                chemical_start = time.time()

                # Calculate hydrophobic features
                func_start = time.time()
                hydrophobic_features = calculate_hydrophobic_feature(residue, neighbors)
                function_timings['calculate_hydrophobic'] = function_timings.get('calculate_hydrophobic', 0) + (
                            time.time() - func_start)

                # Calculate electrostatic features
                func_start = time.time()
                electrostatic_features = calculate_electrostatic_features(residue, neighbors)
                function_timings['calculate_electrostatic'] = function_timings.get('calculate_electrostatic', 0) + (
                            time.time() - func_start)

                # Calculate aromatic features
                func_start = time.time()
                aromatic_features = calculate_aromatic_features(residue, neighbors)
                function_timings['calculate_aromatic'] = function_timings.get('calculate_aromatic', 0) + (
                            time.time() - func_start)

                category_timings['chemical_environment_analysis'] += time.time() - chemical_start

                # Get HSE data
                residue_hse = hse_data.get(res_id, {})

                # ------ RESULT COMPILATION PHASE ------
                # Create result dictionary
                result = {
                    'UniProtID': uniprot_id,
                    'ResidueNumber': res_id,
                    'Site': f"{uniprot_id}_{res_id}",
                    'ResidueType': res_code,
                    'Motif': motif,
                    'pLDDT': plddt,
                    'NeighborCount': len(neighbors),
                    'ChainID': chain_id,
                    'SecondaryStructure': secondary_structure,
                    'SecondaryStructureDesc': ss_description,
                    'HydroxylExposure': hydroxyl_exposure,
                    'BackboneContacts': backbone_contacts,
                    'HydrogenBonds': h_bonds,
                    'MotifPLDDT': motif_plddt,

                    # SASA features
                    'SASA_Total': total_sasa,
                    'SASA_Hydroxyl': hydroxyl_sasa,
                    'SASA_Ratio': hydroxyl_sasa / total_sasa if total_sasa > 0 else None,

                    # Residue depth features
                    'MinDepth': min_depth,
                    'AvgDepth': avg_depth,
                    'CADepth': ca_depth,
                    'HydroxylDepth': hydroxyl_depth,

                    # Add secondary structure distances
                    'SeqDistToHelix': ss_distances.get('helix', (None, None))[0],
                    'SpatialDistToHelix': ss_distances.get('helix', (None, None))[1],
                    'SeqDistToStrand': ss_distances.get('strand', (None, None))[0],
                    'SpatialDistToStrand': ss_distances.get('strand', (None, None))[1],
                    'SeqDistToTurn': ss_distances.get('turn', (None, None))[0],
                    'SpatialDistToTurn': ss_distances.get('turn', (None, None))[1],
                    'SeqDistToCoil': ss_distances.get('coil', (None, None))[0],
                    'SpatialDistToCoil': ss_distances.get('coil', (None, None))[1],
                }

                # Add side chain angles
                result.update({
                    'Chi1': side_chain_angles['chi1'],
                    'Chi2': side_chain_angles['chi2'],
                    'HydroxylAngle': side_chain_angles['hydroxyl_angle'],
                    'CA_Hydroxyl_Distance': side_chain_angles['ca_cb_og_dist'],
                })

                # Add all other feature dictionaries
                result.update(packing_density)
                result.update(contact_order_data)
                result.update(hbond_energy_data)
                result.update(bfactor_stats)
                result.update(pocket_shape)
                result.update(local_structure)
                result.update(hydrophobic_features)
                result.update(electrostatic_features)
                result.update(aromatic_features)
                result.update(solvent_features)
                result.update(domain_features)
                result.update(pocket_props)
                result.update(nearest_aa_distances)
                result.update(nearest_property_distances)
                result.update(aa_counts)
                result.update(property_counts)

                # Update the result dictionary with new features
                result.update(disorder_features)
                result.update(plddt_gradient_features)
                result.update(core_rim_features)
                result.update({
                    'MinDepth': min_depth,
                    'AvgDepth': avg_depth,
                    'CADepth': ca_depth,
                    'HydroxylDepth': hydroxyl_depth
                })
                result.update(cavity_features)
                # result.update(bottleneck_metrics)
                # result.update(hydration_features)

                # Add HSE data
                result.update({
                    'HSE_CA_U': residue_hse.get('HSE_CA_U'),
                    'HSE_CA_D': residue_hse.get('HSE_CA_D'),
                    'HSE_CA_RATIO': residue_hse.get('HSE_CA_RATIO'),
                    'HSE_CB_U': residue_hse.get('HSE_CB_U'),
                    'HSE_CB_D': residue_hse.get('HSE_CB_D'),
                    'HSE_CB_RATIO': residue_hse.get('HSE_CB_RATIO'),
                })

                results.append(result)
                total_sty_found += 1

    except Exception as e:
        print(f"Error processing {structure_file}: {e}")
        traceback.print_exc()

    # Calculate total execution time
    total_time = time.time() - total_start_time

    # Print timing summary
    print("\n----- Timing Summary -----")

    # Print category timings
    print("\nFunction Category Timings:")
    print(f"{'Category':<30} {'Time (s)':<10} {'% of Total':<10}")
    print("-" * 50)
    for category, duration in sorted(category_timings.items(), key=lambda x: x[1], reverse=True):
        percent = (duration / total_time) * 100
        print(f"{category:<30} {duration:<10.3f} {percent:<10.1f}%")

    # Print individual function timings
    print("\nIndividual Function Timings:")
    print(f"{'Function':<40} {'Time (s)':<10} {'% of Total':<10}")
    print("-" * 60)
    for func, duration in sorted(function_timings.items(), key=lambda x: x[1], reverse=True):
        percent = (duration / total_time) * 100
        print(f"{func:<40} {duration:<10.3f} {percent:<10.1f}%")

    print(f"\nTotal processing time: {total_time:.3f} seconds")
    print(f"Found {total_sty_found} STY sites")

    return results


def process_structure_file_parallel(args):
    """Wrapper function for multiprocessing that unpacks arguments"""
    structure_file, file_idx, total_files, temp_dir = args
    try:
        # Call your existing process_structure_file function
        results = process_structure_file(structure_file, total_files, file_idx)

        # Save results to a temp file
        if results:
            temp_file = os.path.join(temp_dir, f"temp_results_{file_idx}.csv")
            temp_df = pd.DataFrame(results)
            temp_df.to_csv(temp_file, index=False)
            return len(results)
        return 0
    except Exception as e:
        print(f"Error processing {structure_file}: {e}")
        traceback.print_exc()
        return 0


# Main function implementation for parallel processing
def main():
    parser = argparse.ArgumentParser(
        description='Analyze STY residues in protein structure files with enhanced features')
    parser.add_argument('structure_dir', nargs='?', default='Complete_AF_Proteome/',
                        help='Directory containing structure files (.pdb, .cif, or gzipped versions)')
    parser.add_argument('--output', '-o', default='sty_analysis_enhanced.csv', help='Output CSV file')
    parser.add_argument('--temp-dir', '-d', default='./temp_results', help='Directory to store temporary result files')
    parser.add_argument('--file-types', '-t', default='cif',
                        help='Comma-separated list of file types to process (default: cif)')
    parser.add_argument('--max-files', '-m', type=int, default=0, help='Maximum number of files to process (0 for all)')
    parser.add_argument('--summary', '-s', action='store_true', help='Generate summary statistics after analysis')
    parser.add_argument('--fallback-pdb', '-f', action='store_true', help='Fallback to PDB files if CIF not available')
    parser.add_argument('--num-processes', '-p', type=int, default=0,
                        help='Number of processes to use (0 for automatic detection)')
    parser.add_argument('--batch-size', '-b', type=int, default=10,
                        help='Number of files to process in each batch')
    args = parser.parse_args()

    # Create temp directory if it doesn't exist
    if not os.path.exists(args.temp_dir):
        print(f"Creating temporary directory: {args.temp_dir}")
        os.makedirs(args.temp_dir)

    # Parse file types to look for
    file_extensions = args.file_types.split(',')

    # Get list of structure files, prioritizing CIF over PDB for the same UniProt ID
    structure_files = []
    uniprot_files = {}  # Maps UniProt ID to file path

    # First pass: collect all files by UniProt ID
    print("Scanning directory for structure files...")
    for f in os.listdir(args.structure_dir):
        # Extract UniProt ID from filename
        if 'AF-' in f and '-F' in f:
            uniprot_id = f.split('AF-')[1].split('-F')[0]
            file_path = os.path.join(args.structure_dir, f)

            # Check if this is a supported file type
            is_supported = any(f.endswith(ext) or f.endswith(f"{ext}.gz") for ext in file_extensions)

            # If fallback is enabled, also check for PDB files
            if args.fallback_pdb and not is_supported:
                is_supported = any(f.endswith(ext) or f.endswith(f"{ext}.gz") for ext in ['pdb'])

            if is_supported:
                if uniprot_id not in uniprot_files:
                    uniprot_files[uniprot_id] = []
                uniprot_files[uniprot_id].append(file_path)

    # Second pass: prioritize CIF files over PDB
    for uniprot_id, files in uniprot_files.items():
        # Sort files by extension preference
        sorted_files = sorted(files,
                              key=lambda x: 0 if x.endswith('.cif') or x.endswith('.cif.gz') else 1)

        # Add the highest priority file for this UniProt ID
        if sorted_files:
            structure_files.append(sorted_files[0])

    # Limit number of files if requested
    if args.max_files > 0 and args.max_files < len(structure_files):
        print(f"Limiting analysis to first {args.max_files} files out of {len(structure_files)} total")
        structure_files = structure_files[:args.max_files]

    total_files = len(structure_files)
    print(f"Found {total_files} unique structure files.")

    # Determine optimal number of processes
    if args.num_processes <= 0:
        # Get number of available CPU cores
        available_cores = mp.cpu_count()
        # Use slightly fewer cores than available to leave some headroom
        args.num_processes = max(1, min(120, int(available_cores * 0.9)))

    print(f"Running with {args.num_processes} parallel processes")

    # Track overall progress
    start_time = time.time()
    total_sty_sites = 0
    processed_file_count = 0

    try:
        # Import tqdm for progress bar
        from tqdm import tqdm
    except ImportError:
        # Define a simple fallback if tqdm is not installed
        class tqdm:
            def __init__(self, total, desc):
                self.total = total
                self.n = 0
                self.desc = desc
                print(f"{desc}: 0/{total}")

            def update(self, n):
                self.n += n
                print(f"{self.desc}: {self.n}/{self.total}")

    # Prepare arguments for parallel processing
    process_args = []
    for i, structure_file in enumerate(structure_files):
        process_args.append((structure_file, i + 1, total_files, args.temp_dir))

    # Create a pool of worker processes
    pool = mp.Pool(processes=args.num_processes)

    # Process files in parallel with progress tracking
    with tqdm(total=len(process_args), desc="Processing files") as pbar:
        for result in pool.imap_unordered(process_structure_file_parallel, process_args):
            total_sty_sites += result
            if result > 0:
                processed_file_count += 1
            pbar.update(1)

    # Close the pool
    pool.close()
    pool.join()

    # Calculate total execution time
    total_time = time.time() - start_time

    # Combine all temp files into final result
    print(f"\nCombining {processed_file_count} temporary files into final output...")

    temp_files = [os.path.join(args.temp_dir, f) for f in os.listdir(args.temp_dir)
                  if f.startswith("temp_results_") and f.endswith(".csv")]

    if temp_files:
        # Read and combine all temporary CSV files in smaller batches to manage memory
        batch_size = args.batch_size
        num_batches = (len(temp_files) + batch_size - 1) // batch_size

        # Process in batches
        combined_df = None

        for batch_idx in range(num_batches):
            start_idx = batch_idx * batch_size
            end_idx = min((batch_idx + 1) * batch_size, len(temp_files))

            batch_files = temp_files[start_idx:end_idx]
            batch_dfs = []

            print(f"Processing batch {batch_idx + 1}/{num_batches} ({len(batch_files)} files)")

            for temp_file in batch_files:
                try:
                    df = pd.read_csv(temp_file)
                    batch_dfs.append(df)
                except Exception as e:
                    print(f"Error reading {temp_file}: {e}")

            if batch_dfs:
                # Combine dataframes from this batch
                batch_combined = pd.concat(batch_dfs, ignore_index=True)

                # Add to the full dataset
                if combined_df is None:
                    combined_df = batch_combined
                else:
                    combined_df = pd.concat([combined_df, batch_combined], ignore_index=True)

                # Clear memory
                batch_dfs = None
                batch_combined = None

        if combined_df is not None and not combined_df.empty:
            combined_df.to_csv(args.output, index=False)

            print(f"\nAnalysis complete. Combined {len(combined_df)} STY sites from {processed_file_count} structures.")
            print(f"Results saved to {args.output}")
            print(f"Total execution time: {total_time / 60:.1f} minutes ({total_time / 3600:.2f} hours)")

            # Generate summary statistics if requested
            if args.summary:
                print("\nGenerating summary statistics...")

                # Count of each residue type
                residue_counts = combined_df['ResidueType'].value_counts()
                print("\nResidue Type Counts:")
                for res, count in residue_counts.items():
                    print(f"  {res}: {count} ({count / len(combined_df) * 100:.1f}%)")

                # Secondary structure distribution
                if 'SecondaryStructure' in combined_df.columns:
                    ss_counts = combined_df['SecondaryStructure'].value_counts()
                    print("\nSecondary Structure Distribution:")
                    for ss, count in ss_counts.items():
                        print(f"  {ss}: {count} ({count / len(combined_df) * 100:.1f}%)")

                # Write summary to file
                summary_file = f"{os.path.splitext(args.output)[0]}_summary.txt"
                with open(summary_file, 'w') as f:
                    f.write(f"Enhanced STY Analysis Summary\n")
                    f.write(f"===================\n\n")
                    f.write(f"Total STY sites analyzed: {len(combined_df)}\n")
                    f.write(f"Structures processed: {processed_file_count}\n\n")

                    f.write("Residue Type Counts:\n")
                    for res, count in residue_counts.items():
                        f.write(f"  {res}: {count} ({count / len(combined_df) * 100:.1f}%)\n")

                    if 'SecondaryStructure' in combined_df.columns:
                        f.write("\nSecondary Structure Distribution:\n")
                        for ss, count in ss_counts.items():
                            f.write(f"  {ss}: {count} ({count / len(combined_df) * 100:.1f}%)\n")

                print(f"Summary statistics saved to {summary_file}")
        else:
            print("\nNo results were found after combining temporary files.")
    else:
        print("\nNo temporary result files were found.")

    print("\nAnalysis complete!")


if __name__ == "__main__":
    # Set up multiprocessing to use 'spawn' method for better compatibility
    try:
        mp.set_start_method('spawn')
    except RuntimeError:
        # Method already set
        pass


    # Register signal handlers
    def signal_handler(sig, frame):
        print('Received signal to terminate. Cleaning up...')
        sys.exit(0)


    signal.signal(signal.SIGINT, signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)

    main()