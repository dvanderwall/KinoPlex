#!/usr/bin/env python3
"""
KinoPlex Database Builder - OPTIMIZED VERSION

This script creates the KinoPlex database from phosphosite and protein metadata feather files.
Builds a normalized schema optimized for the next-generation phosphoproteomics platform.
Features memory-efficient kinase score loading and robust data validation.

Usage:
    python build_kinoplex_db.py

Requirements:
    pandas, sqlite3 (built-in), tqdm (optional for progress bars)
"""

import pandas as pd
import sqlite3
import os
import re
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import sys
import gc

# Optional: Progress bars (install with: pip install tqdm)
try:
    from tqdm import tqdm

    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False


    def tqdm(iterable, **kwargs):
        return iterable

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class KinoPlexDatabaseBuilder:
    """Build and populate the KinoPlex database from feather files."""

    def __init__(self,
                 phosphosite_file: str,
                 metadata_file: str,
                 db_path: str = "kinoplex.db",
                 overwrite: bool = True):
        """
        Initialize the database builder.

        Args:
            phosphosite_file: Path to phosphosite data feather file
            metadata_file: Path to protein metadata feather file
            db_path: Path for output database file
            overwrite: Whether to overwrite existing database
        """
        self.phosphosite_file = Path(phosphosite_file)
        self.metadata_file = Path(metadata_file)
        self.db_path = Path(db_path)
        self.overwrite = overwrite

        # Validate input files
        if not self.phosphosite_file.exists():
            raise FileNotFoundError(f"Phosphosite file not found: {phosphosite_file}")
        if not self.metadata_file.exists():
            raise FileNotFoundError(f"Metadata file not found: {metadata_file}")

        # Handle existing database
        if self.db_path.exists() and not overwrite:
            raise FileExistsError(f"Database exists and overwrite=False: {db_path}")
        elif self.db_path.exists():
            self.db_path.unlink()
            logger.info(f"Removed existing database: {db_path}")

        self.conn = None

    def __enter__(self):
        """Context manager entry."""
        self.conn = sqlite3.connect(str(self.db_path))
        self.conn.execute("PRAGMA journal_mode=WAL")  # Better performance
        self.conn.execute("PRAGMA synchronous=NORMAL")
        self.conn.execute("PRAGMA temp_store=MEMORY")
        self.conn.execute("PRAGMA mmap_size=268435456")  # 256MB mmap
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        if self.conn:
            self.conn.close()

    def create_schema(self) -> None:
        """Create the database schema with all tables and indexes."""
        logger.info("Creating database schema...")

        schema_sql = """
        -- Proteins table (using gene_symbol as primary identifier)
        CREATE TABLE proteins (
            gene_symbol VARCHAR(50) PRIMARY KEY,
            ncbi_gene_id INTEGER,
            ensembl_gene_id VARCHAR(50),
            protein_count INTEGER DEFAULT 1,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        );

        -- Phosphosites table (main data) - WIDE FORMAT with all kinase scores
        CREATE TABLE phosphosites (
            site_id VARCHAR(50) PRIMARY KEY,
            uniprot_id VARCHAR(20),
            gene_symbol VARCHAR(50) NOT NULL,
            position INTEGER NOT NULL,
            residue_type VARCHAR(1) CHECK (residue_type IN ('S', 'T')),
            motif VARCHAR(15),

            -- Prediction scores
            known_positive BOOLEAN DEFAULT FALSE,
            predicted_prob_raw FLOAT,
            predicted_prob_calibrated FLOAT,
            predicted_raw INTEGER,
            predicted_calibrated INTEGER,
            qvalue FLOAT,
            significant_fdr05 BOOLEAN DEFAULT FALSE,
            significant_fdr10 BOOLEAN DEFAULT FALSE,
            significant_fdr20 BOOLEAN DEFAULT FALSE,

            -- All 303 Kinase Scores (wide format)
            AAK1_MotifScore FLOAT,
            ACVR2A_MotifScore FLOAT,
            ACVR2B_MotifScore FLOAT,
            AKT1_MotifScore FLOAT,
            AKT2_MotifScore FLOAT,
            AKT3_MotifScore FLOAT,
            ALK2_MotifScore FLOAT,
            ALK4_MotifScore FLOAT,
            ALPHAK3_MotifScore FLOAT,
            AMPKA1_MotifScore FLOAT,
            AMPKA2_MotifScore FLOAT,
            ANKRD3_MotifScore FLOAT,
            ASK1_MotifScore FLOAT,
            ATM_MotifScore FLOAT,
            ATR_MotifScore FLOAT,
            AURA_MotifScore FLOAT,
            AURB_MotifScore FLOAT,
            AURC_MotifScore FLOAT,
            BCKDK_MotifScore FLOAT,
            BIKE_MotifScore FLOAT,
            BMPR1A_MotifScore FLOAT,
            BMPR1B_MotifScore FLOAT,
            BMPR2_MotifScore FLOAT,
            BRAF_MotifScore FLOAT,
            BRSK1_MotifScore FLOAT,
            BRSK2_MotifScore FLOAT,
            BUB1_MotifScore FLOAT,
            CAMK1A_MotifScore FLOAT,
            CAMK1B_MotifScore FLOAT,
            CAMK1D_MotifScore FLOAT,
            CAMK1G_MotifScore FLOAT,
            CAMK2A_MotifScore FLOAT,
            CAMK2B_MotifScore FLOAT,
            CAMK2D_MotifScore FLOAT,
            CAMK2G_MotifScore FLOAT,
            CAMK4_MotifScore FLOAT,
            CAMKK1_MotifScore FLOAT,
            CAMKK2_MotifScore FLOAT,
            CAMLCK_MotifScore FLOAT,
            CDC7_MotifScore FLOAT,
            CDK1_MotifScore FLOAT,
            CDK10_MotifScore FLOAT,
            CDK12_MotifScore FLOAT,
            CDK13_MotifScore FLOAT,
            CDK14_MotifScore FLOAT,
            CDK16_MotifScore FLOAT,
            CDK17_MotifScore FLOAT,
            CDK18_MotifScore FLOAT,
            CDK19_MotifScore FLOAT,
            CDK2_MotifScore FLOAT,
            CDK3_MotifScore FLOAT,
            CDK4_MotifScore FLOAT,
            CDK5_MotifScore FLOAT,
            CDK6_MotifScore FLOAT,
            CDK7_MotifScore FLOAT,
            CDK8_MotifScore FLOAT,
            CDK9_MotifScore FLOAT,
            CDKL1_MotifScore FLOAT,
            CDKL5_MotifScore FLOAT,
            CHAK1_MotifScore FLOAT,
            CHAK2_MotifScore FLOAT,
            CHK1_MotifScore FLOAT,
            CHK2_MotifScore FLOAT,
            CK1A_MotifScore FLOAT,
            CK1A2_MotifScore FLOAT,
            CK1D_MotifScore FLOAT,
            CK1E_MotifScore FLOAT,
            CK1G1_MotifScore FLOAT,
            CK1G2_MotifScore FLOAT,
            CK1G3_MotifScore FLOAT,
            CK2A1_MotifScore FLOAT,
            CK2A2_MotifScore FLOAT,
            CLK1_MotifScore FLOAT,
            CLK2_MotifScore FLOAT,
            CLK3_MotifScore FLOAT,
            CLK4_MotifScore FLOAT,
            COT_MotifScore FLOAT,
            CRIK_MotifScore FLOAT,
            DAPK1_MotifScore FLOAT,
            DAPK2_MotifScore FLOAT,
            DAPK3_MotifScore FLOAT,
            DCAMKL1_MotifScore FLOAT,
            DCAMKL2_MotifScore FLOAT,
            DLK_MotifScore FLOAT,
            DMPK1_MotifScore FLOAT,
            DNAPK_MotifScore FLOAT,
            DRAK1_MotifScore FLOAT,
            DSTYK_MotifScore FLOAT,
            DYRK1A_MotifScore FLOAT,
            DYRK1B_MotifScore FLOAT,
            DYRK2_MotifScore FLOAT,
            DYRK3_MotifScore FLOAT,
            DYRK4_MotifScore FLOAT,
            EEF2K_MotifScore FLOAT,
            ERK1_MotifScore FLOAT,
            ERK2_MotifScore FLOAT,
            ERK5_MotifScore FLOAT,
            ERK7_MotifScore FLOAT,
            FAM20C_MotifScore FLOAT,
            GAK_MotifScore FLOAT,
            GCK_MotifScore FLOAT,
            GCN2_MotifScore FLOAT,
            GRK1_MotifScore FLOAT,
            GRK2_MotifScore FLOAT,
            GRK3_MotifScore FLOAT,
            GRK4_MotifScore FLOAT,
            GRK5_MotifScore FLOAT,
            GRK6_MotifScore FLOAT,
            GRK7_MotifScore FLOAT,
            GSK3A_MotifScore FLOAT,
            GSK3B_MotifScore FLOAT,
            HASPIN_MotifScore FLOAT,
            HGK_MotifScore FLOAT,
            HIPK1_MotifScore FLOAT,
            HIPK2_MotifScore FLOAT,
            HIPK3_MotifScore FLOAT,
            HIPK4_MotifScore FLOAT,
            HPK1_MotifScore FLOAT,
            HRI_MotifScore FLOAT,
            HUNK_MotifScore FLOAT,
            ICK_MotifScore FLOAT,
            IKKA_MotifScore FLOAT,
            IKKB_MotifScore FLOAT,
            IKKE_MotifScore FLOAT,
            IRAK1_MotifScore FLOAT,
            IRAK4_MotifScore FLOAT,
            IRE1_MotifScore FLOAT,
            IRE2_MotifScore FLOAT,
            JNK1_MotifScore FLOAT,
            JNK2_MotifScore FLOAT,
            JNK3_MotifScore FLOAT,
            KHS1_MotifScore FLOAT,
            KHS2_MotifScore FLOAT,
            KIS_MotifScore FLOAT,
            LATS1_MotifScore FLOAT,
            LATS2_MotifScore FLOAT,
            LKB1_MotifScore FLOAT,
            LOK_MotifScore FLOAT,
            LRRK2_MotifScore FLOAT,
            MAK_MotifScore FLOAT,
            MAP3K15_MotifScore FLOAT,
            MAPKAPK2_MotifScore FLOAT,
            MAPKAPK3_MotifScore FLOAT,
            MAPKAPK5_MotifScore FLOAT,
            MARK1_MotifScore FLOAT,
            MARK2_MotifScore FLOAT,
            MARK3_MotifScore FLOAT,
            MARK4_MotifScore FLOAT,
            MASTL_MotifScore FLOAT,
            MEK1_MotifScore FLOAT,
            MEK2_MotifScore FLOAT,
            MEK5_MotifScore FLOAT,
            MEKK1_MotifScore FLOAT,
            MEKK2_MotifScore FLOAT,
            MEKK3_MotifScore FLOAT,
            MEKK6_MotifScore FLOAT,
            MELK_MotifScore FLOAT,
            MINK_MotifScore FLOAT,
            MLK1_MotifScore FLOAT,
            MLK2_MotifScore FLOAT,
            MLK3_MotifScore FLOAT,
            MLK4_MotifScore FLOAT,
            MNK1_MotifScore FLOAT,
            MNK2_MotifScore FLOAT,
            MOK_MotifScore FLOAT,
            MOS_MotifScore FLOAT,
            MPSK1_MotifScore FLOAT,
            MRCKA_MotifScore FLOAT,
            MRCKB_MotifScore FLOAT,
            MSK1_MotifScore FLOAT,
            MSK2_MotifScore FLOAT,
            MST1_MotifScore FLOAT,
            MST2_MotifScore FLOAT,
            MST3_MotifScore FLOAT,
            MST4_MotifScore FLOAT,
            MTOR_MotifScore FLOAT,
            MYLK4_MotifScore FLOAT,
            MYO3A_MotifScore FLOAT,
            MYO3B_MotifScore FLOAT,
            NDR1_MotifScore FLOAT,
            NDR2_MotifScore FLOAT,
            NEK1_MotifScore FLOAT,
            NEK11_MotifScore FLOAT,
            NEK2_MotifScore FLOAT,
            NEK3_MotifScore FLOAT,
            NEK4_MotifScore FLOAT,
            NEK5_MotifScore FLOAT,
            NEK6_MotifScore FLOAT,
            NEK7_MotifScore FLOAT,
            NEK8_MotifScore FLOAT,
            NEK9_MotifScore FLOAT,
            NIK_MotifScore FLOAT,
            NIM1_MotifScore FLOAT,
            NLK_MotifScore FLOAT,
            NUAK1_MotifScore FLOAT,
            NUAK2_MotifScore FLOAT,
            OSR1_MotifScore FLOAT,
            P38A_MotifScore FLOAT,
            P38B_MotifScore FLOAT,
            P38D_MotifScore FLOAT,
            P38G_MotifScore FLOAT,
            P70S6K_MotifScore FLOAT,
            P70S6KB_MotifScore FLOAT,
            P90RSK_MotifScore FLOAT,
            PAK1_MotifScore FLOAT,
            PAK2_MotifScore FLOAT,
            PAK3_MotifScore FLOAT,
            PAK4_MotifScore FLOAT,
            PAK5_MotifScore FLOAT,
            PAK6_MotifScore FLOAT,
            PASK_MotifScore FLOAT,
            PBK_MotifScore FLOAT,
            PDHK1_MotifScore FLOAT,
            PDHK4_MotifScore FLOAT,
            PDK1_MotifScore FLOAT,
            PERK_MotifScore FLOAT,
            PHKG1_MotifScore FLOAT,
            PHKG2_MotifScore FLOAT,
            PIM1_MotifScore FLOAT,
            PIM2_MotifScore FLOAT,
            PIM3_MotifScore FLOAT,
            PINK1_MotifScore FLOAT,
            PKACA_MotifScore FLOAT,
            PKACB_MotifScore FLOAT,
            PKACG_MotifScore FLOAT,
            PKCA_MotifScore FLOAT,
            PKCB_MotifScore FLOAT,
            PKCD_MotifScore FLOAT,
            PKCE_MotifScore FLOAT,
            PKCG_MotifScore FLOAT,
            PKCH_MotifScore FLOAT,
            PKCI_MotifScore FLOAT,
            PKCT_MotifScore FLOAT,
            PKCZ_MotifScore FLOAT,
            PKG1_MotifScore FLOAT,
            PKG2_MotifScore FLOAT,
            PKN1_MotifScore FLOAT,
            PKN2_MotifScore FLOAT,
            PKN3_MotifScore FLOAT,
            PKR_MotifScore FLOAT,
            PLK1_MotifScore FLOAT,
            PLK2_MotifScore FLOAT,
            PLK3_MotifScore FLOAT,
            PLK4_MotifScore FLOAT,
            PRKD1_MotifScore FLOAT,
            PRKD2_MotifScore FLOAT,
            PRKD3_MotifScore FLOAT,
            PRKX_MotifScore FLOAT,
            PRP4_MotifScore FLOAT,
            PRPK_MotifScore FLOAT,
            QIK_MotifScore FLOAT,
            QSK_MotifScore FLOAT,
            RAF1_MotifScore FLOAT,
            RIPK1_MotifScore FLOAT,
            RIPK2_MotifScore FLOAT,
            RIPK3_MotifScore FLOAT,
            ROCK1_MotifScore FLOAT,
            ROCK2_MotifScore FLOAT,
            RSK2_MotifScore FLOAT,
            RSK3_MotifScore FLOAT,
            RSK4_MotifScore FLOAT,
            SBK_MotifScore FLOAT,
            SGK1_MotifScore FLOAT,
            SGK3_MotifScore FLOAT,
            SIK_MotifScore FLOAT,
            SKMLCK_MotifScore FLOAT,
            SLK_MotifScore FLOAT,
            SMG1_MotifScore FLOAT,
            SMMLCK_MotifScore FLOAT,
            SNRK_MotifScore FLOAT,
            SRPK1_MotifScore FLOAT,
            SRPK2_MotifScore FLOAT,
            SRPK3_MotifScore FLOAT,
            SSTK_MotifScore FLOAT,
            STK33_MotifScore FLOAT,
            STLK3_MotifScore FLOAT,
            TAK1_MotifScore FLOAT,
            TAO1_MotifScore FLOAT,
            TAO2_MotifScore FLOAT,
            TAO3_MotifScore FLOAT,
            TBK1_MotifScore FLOAT,
            TGFBR1_MotifScore FLOAT,
            TGFBR2_MotifScore FLOAT,
            TLK1_MotifScore FLOAT,
            TLK2_MotifScore FLOAT,
            TNIK_MotifScore FLOAT,
            TSSK1_MotifScore FLOAT,
            TSSK2_MotifScore FLOAT,
            TTBK1_MotifScore FLOAT,
            TTBK2_MotifScore FLOAT,
            TTK_MotifScore FLOAT,
            ULK1_MotifScore FLOAT,
            ULK2_MotifScore FLOAT,
            VRK1_MotifScore FLOAT,
            VRK2_MotifScore FLOAT,
            WNK1_MotifScore FLOAT,
            WNK3_MotifScore FLOAT,
            WNK4_MotifScore FLOAT,
            YANK2_MotifScore FLOAT,
            YANK3_MotifScore FLOAT,
            YSK1_MotifScore FLOAT,
            YSK4_MotifScore FLOAT,
            ZAK_MotifScore FLOAT,

            -- Structural features
            secondary_structure VARCHAR(10),
            secondary_structure_desc TEXT,
            seq_dist_to_helix INTEGER,
            spatial_dist_to_helix FLOAT,
            seq_dist_to_strand FLOAT,
            spatial_dist_to_strand FLOAT,
            seq_dist_to_turn INTEGER,
            spatial_dist_to_turn FLOAT,
            plddt FLOAT,
            neighbor_count INTEGER,
            hydroxyl_exposure FLOAT,
            hydrogen_bonds INTEGER,
            motif_plddt FLOAT,
            sasa_hydroxyl FLOAT,
            sasa_ratio FLOAT,

            -- Packing density features
            packing_density_5a FLOAT,
            packing_density_8a FLOAT,
            packing_density_10a FLOAT,
            packing_density_12a FLOAT,
            packing_density_15a FLOAT,
            hydroxyl_packing_5a FLOAT,
            hydroxyl_packing_8a FLOAT,
            hydroxyl_packing_10a FLOAT,
            hydroxyl_packing_12a FLOAT,
            hydroxyl_packing_15a FLOAT,

            -- Additional structural features  
            contact_order FLOAT,
            long_range_contacts INTEGER,
            short_range_contacts INTEGER,
            backbone_hbonds INTEGER,
            sidechain_hbonds INTEGER,
            avg_bfactor FLOAT,
            neighbor_avg_bfactor FLOAT,
            pocket_width_y FLOAT,
            pocket_width_z FLOAT,
            pocket_volume FLOAT,
            pocket_asymmetry FLOAT,
            pocket_atom_count INTEGER,
            local_rmsd FLOAT,
            hydrophobic_patch_size INTEGER,
            aromatic_count INTEGER,
            water_channel_score FLOAT,
            polar_atom_count INTEGER,
            structural_discontinuity FLOAT,
            hydrophobic_count INTEGER,
            polar_count INTEGER,
            small_count INTEGER,
            large_count INTEGER,
            aliphatic_count INTEGER,

            -- Distance features
            dist_to_c FLOAT,
            dist_to_f FLOAT,
            dist_to_w FLOAT,
            dist_to_y FLOAT,
            num_hydrophobic_10a INTEGER,
            num_polar_10a INTEGER,
            num_charged_10a INTEGER,
            num_acidic_10a INTEGER,
            num_basic_10a INTEGER,

            -- Additional features
            disorder_score FLOAT,
            cavity_volume FLOAT,
            pocket_openness FLOAT,
            cavity_depth FLOAT,
            cavity_polarity FLOAT,
            cavity_sphericity FLOAT,
            hse_ca_u INTEGER,
            hse_ca_d INTEGER,
            hse_ca_ratio FLOAT,
            hse_cb_u INTEGER,
            hse_cb_d INTEGER,
            hse_cb_ratio FLOAT,

            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,

            FOREIGN KEY (gene_symbol) REFERENCES proteins(gene_symbol)
        );

        -- Remove kinase_scores table - not needed in wide format

        -- Gene sets table (unchanged)
        CREATE TABLE gene_sets (
            gs_id INTEGER PRIMARY KEY AUTOINCREMENT,
            gs_name VARCHAR(200) NOT NULL,
            gs_collection VARCHAR(100),
            gs_description TEXT,
            gs_pmid INTEGER,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,

            UNIQUE(gs_name, gs_collection)
        );

        -- Gene set memberships (unchanged)
        CREATE TABLE gene_set_memberships (
            membership_id INTEGER PRIMARY KEY AUTOINCREMENT,
            gene_symbol VARCHAR(50) NOT NULL,
            ncbi_gene_id INTEGER,
            ensembl_gene_id VARCHAR(50),
            gs_id INTEGER NOT NULL,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,

            FOREIGN KEY (gene_symbol) REFERENCES proteins(gene_symbol),
            FOREIGN KEY (gs_id) REFERENCES gene_sets(gs_id),
            UNIQUE(gene_symbol, gs_id)
        );

        -- Create indexes for optimal query performance
        CREATE INDEX idx_phosphosites_gene_symbol ON phosphosites(gene_symbol);
        CREATE INDEX idx_phosphosites_qvalue ON phosphosites(qvalue);
        CREATE INDEX idx_phosphosites_significance ON phosphosites(significant_fdr05, significant_fdr10);
        CREATE INDEX idx_phosphosites_position ON phosphosites(gene_symbol, position);
        CREATE INDEX idx_phosphosites_residue ON phosphosites(residue_type);
        CREATE INDEX idx_phosphosites_calibrated_prob ON phosphosites(predicted_prob_calibrated DESC);

        -- No kinase score indexes needed - they're now columns

        CREATE INDEX idx_memberships_gene ON gene_set_memberships(gene_symbol);
        CREATE INDEX idx_memberships_geneset ON gene_set_memberships(gs_id);
        CREATE INDEX idx_memberships_ncbi ON gene_set_memberships(ncbi_gene_id);

        CREATE INDEX idx_gene_sets_collection ON gene_sets(gs_collection);
        CREATE INDEX idx_gene_sets_pmid ON gene_sets(gs_pmid);
        CREATE INDEX idx_gene_sets_name ON gene_sets(gs_name);
        """

        # Execute schema creation
        self.conn.executescript(schema_sql)
        self.conn.commit()
        logger.info("Schema created successfully")

    def load_data(self) -> None:
        """Load all data from feather files into the database."""
        logger.info("Loading data from feather files...")

        # Load the datasets
        logger.info("Reading phosphosite data...")
        phospho_df = pd.read_feather(self.phosphosite_file)
        logger.info(f"Loaded {len(phospho_df):,} phosphosites")

        logger.info("Reading protein metadata...")
        metadata_df = pd.read_feather(self.metadata_file)
        logger.info(f"Loaded {len(metadata_df):,} metadata records")

        # Load each table - NO KINASE TRANSFORMATION NEEDED!
        self._load_proteins(phospho_df, metadata_df)
        self._load_gene_sets(metadata_df)
        self._load_gene_set_memberships(metadata_df)
        self._load_phosphosites(phospho_df)  # Includes all kinase scores as columns

        logger.info("All data loaded successfully - wide format preserved!")

    def _load_proteins(self, phospho_df: pd.DataFrame, metadata_df: pd.DataFrame) -> None:
        """Load proteins table from both datasets."""
        logger.info("Loading proteins table...")

        # Get unique genes from phosphosite data
        phospho_genes = phospho_df[['gene_symbol', 'uniprotID']].drop_duplicates()
        phospho_genes = phospho_genes.rename(columns={'uniprotID': 'uniprot_id'})

        # Get unique genes from metadata
        metadata_genes = metadata_df[['gene_symbol', 'ncbi_gene', 'ensembl_gene']].drop_duplicates()
        metadata_genes = metadata_genes.rename(columns={
            'ncbi_gene': 'ncbi_gene_id',
            'ensembl_gene': 'ensembl_gene_id'
        })

        # Merge and get unique proteins
        proteins = pd.merge(
            phospho_genes, metadata_genes,
            on='gene_symbol', how='outer'
        ).drop_duplicates(subset=['gene_symbol'])

        # Count proteins per gene symbol (in case of duplicates)
        protein_counts = phospho_genes.groupby('gene_symbol').size().reset_index(name='protein_count')
        proteins = pd.merge(proteins, protein_counts, on='gene_symbol', how='left')
        proteins['protein_count'] = proteins['protein_count'].fillna(1).astype(int)

        # Select columns for database
        protein_cols = ['gene_symbol', 'ncbi_gene_id', 'ensembl_gene_id', 'protein_count']
        proteins = proteins[protein_cols]

        # Load to database
        proteins.to_sql('proteins', self.conn, if_exists='append', index=False)
        logger.info(f"Loaded {len(proteins):,} proteins")

    def _load_gene_sets(self, metadata_df: pd.DataFrame) -> None:
        """Load gene sets table."""
        logger.info("Loading gene sets...")

        # Get unique gene sets
        gene_sets = metadata_df[[
            'gs_name', 'gs_collection_name', 'gs_description', 'gs_pmid'
        ]].drop_duplicates()

        gene_sets = gene_sets.rename(columns={
            'gs_collection_name': 'gs_collection'
        })

        # Load to database
        gene_sets.to_sql('gene_sets', self.conn, if_exists='append', index=False)
        logger.info(f"Loaded {len(gene_sets):,} gene sets")

    def _load_gene_set_memberships(self, metadata_df: pd.DataFrame) -> None:
        """Load gene set memberships."""
        logger.info("Loading gene set memberships...")

        # Prepare memberships data
        memberships = metadata_df[[
            'gene_symbol', 'ncbi_gene', 'ensembl_gene', 'gs_name', 'gs_collection_name'
        ]].copy()

        # Get gene set IDs
        gene_sets_map = pd.read_sql_query(
            "SELECT gs_id, gs_name, gs_collection FROM gene_sets",
            self.conn
        )
        gene_sets_map['gs_collection_name'] = gene_sets_map['gs_collection']

        # Merge to get gs_id
        memberships = pd.merge(
            memberships,
            gene_sets_map[['gs_id', 'gs_name', 'gs_collection_name']],
            on=['gs_name', 'gs_collection_name']
        )

        # Select final columns
        membership_cols = ['gene_symbol', 'ncbi_gene', 'ensembl_gene', 'gs_id']
        memberships = memberships[membership_cols].rename(columns={
            'ncbi_gene': 'ncbi_gene_id',
            'ensembl_gene': 'ensembl_gene_id'
        })

        # Debug: Check for duplicates before removal
        duplicates_before = len(memberships)
        unique_constraint_dupes = memberships.duplicated(subset=['gene_symbol', 'gs_id']).sum()
        logger.info(f"Found {unique_constraint_dupes:,} duplicate (gene_symbol, gs_id) combinations")

        # Handle duplicates by consolidating metadata for the same gene-geneset pair
        # Keep the first non-null value for each ID type
        memberships_consolidated = memberships.groupby(['gene_symbol', 'gs_id']).agg({
            'ncbi_gene_id': 'first',  # Take first non-null value
            'ensembl_gene_id': 'first'  # Take first non-null value
        }).reset_index()

        duplicates_after = len(memberships_consolidated)
        logger.info(f"Consolidated {duplicates_before:,} rows to {duplicates_after:,} unique memberships")

        # Load to database
        memberships_consolidated.to_sql('gene_set_memberships', self.conn, if_exists='append', index=False)
        logger.info(f"Loaded {len(memberships_consolidated):,} gene set memberships")

    def _load_phosphosites(self, phospho_df: pd.DataFrame) -> None:
        """Load phosphosites with all structural features."""
        logger.info("Loading phosphosites...")

        # Prepare phosphosite data
        phospho_data = phospho_df.copy()

        # Check for missing critical values
        missing_gene_symbols = phospho_data['gene_symbol'].isnull().sum()
        missing_sites = phospho_data['Site'].isnull().sum()

        if missing_gene_symbols > 0:
            logger.warning(
                f"Found {missing_gene_symbols:,} phosphosites with missing gene_symbol - these will be excluded")
            phospho_data = phospho_data.dropna(subset=['gene_symbol'])

        if missing_sites > 0:
            logger.warning(f"Found {missing_sites:,} phosphosites with missing Site ID - these will be excluded")
            phospho_data = phospho_data.dropna(subset=['Site'])

        logger.info(f"Processing {len(phospho_data):,} valid phosphosites")

        # Check for duplicate site_ids (primary key constraint)
        duplicate_sites = phospho_data['Site'].duplicated().sum()
        if duplicate_sites > 0:
            logger.warning(f"Found {duplicate_sites:,} duplicate site_ids - consolidating...")

            # For duplicate sites, keep the first occurrence and warn about data loss
            original_count = len(phospho_data)
            phospho_data = phospho_data.drop_duplicates(subset=['Site'], keep='first')
            final_count = len(phospho_data)
            logger.info(f"Consolidated {original_count:,} rows to {final_count:,} unique sites")

        logger.info(f"Loading {len(phospho_data):,} unique phosphosites to database")

        # Rename columns to match database schema (convert to lowercase with underscores)
        column_mapping = {
            'Site': 'site_id',
            'uniprotID': 'uniprot_id',
            'gene_symbol': 'gene_symbol',
            'position': 'position',
            'ResidueType': 'residue_type',
            'Motif': 'motif',
            'KnownPositive': 'known_positive',
            'PredictedProb_Raw': 'predicted_prob_raw',
            'PredictedProb_Calibrated': 'predicted_prob_calibrated',
            'Predicted_Raw': 'predicted_raw',
            'Predicted_Calibrated': 'predicted_calibrated',
            'qvalue': 'qvalue',
            'significant_fdr05': 'significant_fdr05',
            'significant_fdr10': 'significant_fdr10',
            'significant_fdr20': 'significant_fdr20',
            'SecondaryStructure': 'secondary_structure',
            'SecondaryStructureDesc': 'secondary_structure_desc',
            'SeqDistToHelix': 'seq_dist_to_helix',
            'SpatialDistToHelix': 'spatial_dist_to_helix',
            'SeqDistToStrand': 'seq_dist_to_strand',
            'SpatialDistToStrand': 'spatial_dist_to_strand',
            'SeqDistToTurn': 'seq_dist_to_turn',
            'SpatialDistToTurn': 'spatial_dist_to_turn',
            'pLDDT': 'plddt',
            'NeighborCount': 'neighbor_count',
            'HydroxylExposure': 'hydroxyl_exposure',
            'HydrogenBonds': 'hydrogen_bonds',
            'MotifPLDDT': 'motif_plddt',
            'SASA_Hydroxyl': 'sasa_hydroxyl',
            'SASA_Ratio': 'sasa_ratio',
            'packing_density_5.0A': 'packing_density_5a',
            'packing_density_8.0A': 'packing_density_8a',
            'packing_density_10.0A': 'packing_density_10a',
            'packing_density_12.0A': 'packing_density_12a',
            'packing_density_15.0A': 'packing_density_15a',
            'hydroxyl_packing_5.0A': 'hydroxyl_packing_5a',
            'hydroxyl_packing_8.0A': 'hydroxyl_packing_8a',
            'hydroxyl_packing_10.0A': 'hydroxyl_packing_10a',
            'hydroxyl_packing_12.0A': 'hydroxyl_packing_12a',
            'hydroxyl_packing_15.0A': 'hydroxyl_packing_15a',
            'contact_order': 'contact_order',
            'long_range_contacts': 'long_range_contacts',
            'short_range_contacts': 'short_range_contacts',
            'backbone_hbonds': 'backbone_hbonds',
            'sidechain_hbonds': 'sidechain_hbonds',
            'avg_bfactor': 'avg_bfactor',
            'neighbor_avg_bfactor': 'neighbor_avg_bfactor',
            'pocket_width_y': 'pocket_width_y',
            'pocket_width_z': 'pocket_width_z',
            'pocket_volume': 'pocket_volume',
            'pocket_asymmetry': 'pocket_asymmetry',
            'pocket_atom_count': 'pocket_atom_count',
            'local_rmsd': 'local_rmsd',
            'hydrophobic_patch_size': 'hydrophobic_patch_size',
            'aromatic_count': 'aromatic_count',
            'water_channel_score': 'water_channel_score',
            'polar_atom_count': 'polar_atom_count',
            'structural_discontinuity': 'structural_discontinuity',
            'hydrophobic_count': 'hydrophobic_count',
            'polar_count': 'polar_count',
            'small_count': 'small_count',
            'large_count': 'large_count',
            'aliphatic_count': 'aliphatic_count',
            'dist_to_C': 'dist_to_c',
            'dist_to_F': 'dist_to_f',
            'dist_to_W': 'dist_to_w',
            'dist_to_Y': 'dist_to_y',
            'num_hydrophobic_10A': 'num_hydrophobic_10a',
            'num_polar_10A': 'num_polar_10a',
            'num_charged_10A': 'num_charged_10a',
            'num_acidic_10A': 'num_acidic_10a',
            'num_basic_10A': 'num_basic_10a',
            'disorder_score': 'disorder_score',
            'cavity_volume': 'cavity_volume',
            'pocket_openness': 'pocket_openness',
            'cavity_depth': 'cavity_depth',
            'cavity_polarity': 'cavity_polarity',
            'cavity_sphericity': 'cavity_sphericity',
            'HSE_CA_U': 'hse_ca_u',
            'HSE_CA_D': 'hse_ca_d',
            'HSE_CA_RATIO': 'hse_ca_ratio',
            'HSE_CB_U': 'hse_cb_u',
            'HSE_CB_D': 'hse_cb_d',
            'HSE_CB_RATIO': 'hse_cb_ratio'
        }

        # Apply column mapping for columns that exist - includes ALL kinase scores
        existing_columns = {k: v for k, v in column_mapping.items() if k in phospho_data.columns}

        # Add kinase score columns (they keep their original names in wide format)
        kinase_columns = [col for col in phospho_data.columns if col.endswith('_MotifScore')]
        logger.info(f"Found {len(kinase_columns)} kinase score columns - keeping in wide format")

        # Add kinase columns to the column mapping (they keep their names)
        for kinase_col in kinase_columns:
            existing_columns[kinase_col] = kinase_col

        phospho_data = phospho_data.rename(columns=existing_columns)

        # Get columns that match our database schema
        db_columns = [col for col in existing_columns.values() if col in phospho_data.columns]

        # Select only the columns we need
        phospho_final = phospho_data[db_columns]

        # Load to database in chunks for better performance
        chunk_size = 1000
        total_chunks = len(phospho_final) // chunk_size + 1

        for i in tqdm(range(total_chunks), desc="Loading phosphosites"):
            start_idx = i * chunk_size
            end_idx = min((i + 1) * chunk_size, len(phospho_final))
            chunk = phospho_final.iloc[start_idx:end_idx]

            if not chunk.empty:
                chunk.to_sql('phosphosites', self.conn, if_exists='append', index=False)

        self.conn.commit()
        logger.info(f"Loaded {len(phospho_final):,} phosphosites")

    def _load_kinase_scores(self, phospho_df: pd.DataFrame) -> None:
        """Extract and load kinase scores - ULTRA-FAST VERSION using SQLite native import."""
        logger.info("Loading kinase scores...")

        # Filter data same as phosphosites to maintain consistency
        phospho_data = phospho_df.copy()

        # Apply same filters as phosphosites
        missing_gene_symbols = phospho_data['gene_symbol'].isnull().sum()
        missing_sites = phospho_data['Site'].isnull().sum()

        if missing_gene_symbols > 0:
            phospho_data = phospho_data.dropna(subset=['gene_symbol'])

        if missing_sites > 0:
            phospho_data = phospho_data.dropna(subset=['Site'])

        # Apply same duplicate site consolidation as phosphosites
        duplicate_sites = phospho_data['Site'].duplicated().sum()
        if duplicate_sites > 0:
            phospho_data = phospho_data.drop_duplicates(subset=['Site'], keep='first')

        # Find kinase score columns
        kinase_columns = [col for col in phospho_data.columns if col.endswith('_MotifScore')]
        logger.info(f"Found {len(kinase_columns)} kinase score columns")

        # ULTRA-FAST APPROACH: Write to CSV and use SQLite's native import
        import tempfile
        import csv

        # Create temporary CSV file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False, newline='') as temp_file:
            csv_path = temp_file.name
            writer = csv.writer(temp_file)

            # Write header
            writer.writerow(['site_id', 'kinase_name', 'motif_score'])

            logger.info("Writing kinase scores to temporary CSV...")

            # SMART FILTERING: Only include scores above threshold to reduce data size
            score_threshold = 0.1  # Adjust this threshold as needed
            total_rows_written = 0

            # Get site IDs and convert to numpy for speed
            site_ids = phospho_data['Site'].values

            # Process each kinase column efficiently using numpy
            for i, kinase_col in enumerate(tqdm(kinase_columns, desc="Processing kinases")):
                kinase_name = kinase_col.replace('_MotifScore', '')
                scores = phospho_data[kinase_col].values

                # Use numpy boolean indexing for speed - only keep scores above threshold
                valid_mask = (~pd.isna(scores)) & (scores >= score_threshold)
                valid_sites = site_ids[valid_mask]
                valid_scores = scores[valid_mask]

                # Write valid scores for this kinase
                for site_id, score in zip(valid_sites, valid_scores):
                    writer.writerow([site_id, kinase_name, float(score)])
                    total_rows_written += 1

        logger.info(f"Wrote {total_rows_written:,} kinase scores to CSV (threshold >= {score_threshold})")

        # Drop indexes before bulk import for maximum speed
        logger.info("Temporarily dropping indexes for faster import...")
        try:
            self.conn.execute("DROP INDEX IF EXISTS idx_kinase_scores_site")
            self.conn.execute("DROP INDEX IF EXISTS idx_kinase_scores_kinase")
            self.conn.execute("DROP INDEX IF EXISTS idx_kinase_scores_score")
            self.conn.execute("DROP INDEX IF EXISTS idx_kinase_scores_top_scores")
        except:
            pass  # Indexes might not exist yet

        # Use SQLite's native CSV import - MUCH faster than pandas
        logger.info("Importing CSV into SQLite using native import...")

        # Configure SQLite for MAXIMUM import speed
        self.conn.execute("PRAGMA synchronous = OFF")  # Disable disk sync
        self.conn.execute("PRAGMA journal_mode = OFF")  # Disable journal completely
        self.conn.execute("PRAGMA temp_store = MEMORY")  # Keep temp data in memory
        self.conn.execute("PRAGMA cache_size = 2000000")  # 8GB cache for 1M batches
        self.conn.execute("PRAGMA locking_mode = EXCLUSIVE")  # Exclusive access
        self.conn.execute("PRAGMA count_changes = OFF")  # Disable change counting

        # Begin single large transaction for entire import
        self.conn.execute("BEGIN IMMEDIATE TRANSACTION")

        # ULTRA-FAST: Read CSV and insert in massive batches
        with open(csv_path, 'r') as f:
            reader = csv.DictReader(f)

            # MAXIMUM BATCH SIZE for ultimate speed
            batch_size = 1000000  # 1M rows per batch
            batch = []
            rows_imported = 0

            # Prepare statement once for reuse
            insert_sql = "INSERT INTO kinase_scores (site_id, kinase_name, motif_score) VALUES (?, ?, ?)"

            for row in reader:
                batch.append((row['site_id'], row['kinase_name'], float(row['motif_score'])))

                if len(batch) >= batch_size:
                    self.conn.executemany(insert_sql, batch)
                    rows_imported += len(batch)
                    batch = []

                    # Less frequent progress updates for speed
                    if rows_imported % 2000000 == 0:  # Every 2M rows
                        logger.info(f"Imported {rows_imported:,} rows...")

            # Insert remaining batch
            if batch:
                self.conn.executemany(insert_sql, batch)
                rows_imported += len(batch)

        # Commit the entire transaction at once
        self.conn.execute("COMMIT")
        logger.info(f"Imported {rows_imported:,} total rows")

        # Restore SQLite settings to safe defaults
        self.conn.execute("PRAGMA synchronous = NORMAL")
        self.conn.execute("PRAGMA journal_mode = WAL")
        self.conn.execute("PRAGMA locking_mode = NORMAL")
        self.conn.execute("PRAGMA count_changes = ON")

        # Rebuild indexes after import
        logger.info("Rebuilding indexes...")
        self.conn.execute("CREATE INDEX idx_kinase_scores_site ON kinase_scores(site_id)")
        self.conn.execute("CREATE INDEX idx_kinase_scores_kinase ON kinase_scores(kinase_name)")
        self.conn.execute("CREATE INDEX idx_kinase_scores_score ON kinase_scores(kinase_name, motif_score DESC)")
        self.conn.execute("CREATE INDEX idx_kinase_scores_top_scores ON kinase_scores(motif_score DESC)")

        self.conn.commit()

        # Clean up temporary file
        os.unlink(csv_path)

        logger.info(f"Loaded {total_rows_written:,} kinase scores (filtered by threshold >= {score_threshold})")
        logger.info("Note: Low-scoring predictions excluded to optimize database size and query performance")

    def create_views(self) -> None:
        """Create useful views for common queries - optimized for wide format."""
        logger.info("Creating database views...")

        views_sql = """
        -- Phosphosite pathway view (unchanged)
        CREATE VIEW phosphosite_pathways AS
        SELECT 
            p.site_id,
            p.gene_symbol,
            p.position,
            p.residue_type,
            p.qvalue,
            p.predicted_prob_calibrated,
            p.significant_fdr05,
            p.significant_fdr10,
            gs.gs_name,
            gs.gs_collection,
            gs.gs_description,
            gs.gs_pmid
        FROM phosphosites p
        JOIN gene_set_memberships gsm ON p.gene_symbol = gsm.gene_symbol
        JOIN gene_sets gs ON gsm.gs_id = gs.gs_id;

        -- Significant phosphosites with pathway info
        CREATE VIEW significant_phosphosites AS
        SELECT 
            p.*,
            COUNT(gsm.gs_id) as pathway_count,
            GROUP_CONCAT(DISTINCT gs.gs_collection) as collections
        FROM phosphosites p
        LEFT JOIN gene_set_memberships gsm ON p.gene_symbol = gsm.gene_symbol
        LEFT JOIN gene_sets gs ON gsm.gs_id = gs.gs_id
        WHERE p.significant_fdr05 = TRUE
        GROUP BY p.site_id;

        -- High confidence predictions view (structural + kinase scores)
        CREATE VIEW high_confidence_sites AS
        SELECT 
            site_id,
            gene_symbol,
            position,
            residue_type,
            predicted_prob_calibrated,
            qvalue,
            -- Top kinase scores (example - can be customized)
            AKT1_MotifScore,
            CDK1_MotifScore,
            GSK3B_MotifScore,
            ERK1_MotifScore,
            ERK2_MotifScore
        FROM phosphosites
        WHERE significant_fdr05 = TRUE 
          AND predicted_prob_calibrated > 0.7
        ORDER BY predicted_prob_calibrated DESC;
        """

        self.conn.executescript(views_sql)
        self.conn.commit()
        logger.info("Views created successfully")

    def validate_database(self) -> None:
        """Validate the database integrity and provide summary statistics."""
        logger.info("Validating database...")

        # Get table counts
        tables = ['proteins', 'phosphosites', 'gene_sets', 'gene_set_memberships', 'kinase_scores']

        print("\n" + "=" * 60)
        print("KinoPlex Database Summary")
        print("=" * 60)

        for table in tables:
            count = self.conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
            print(f"{table:25}: {count:,} records")

        # Additional statistics
        significant_sites = self.conn.execute(
            "SELECT COUNT(*) FROM phosphosites WHERE significant_fdr05 = TRUE"
        ).fetchone()[0]

        unique_genes = self.conn.execute(
            "SELECT COUNT(DISTINCT gene_symbol) FROM phosphosites"
        ).fetchone()[0]

        unique_kinases = self.conn.execute(
            "SELECT COUNT(DISTINCT kinase_name) FROM kinase_scores"
        ).fetchone()[0]

        avg_kinases_per_site = self.conn.execute(
            "SELECT AVG(kinase_count) FROM (SELECT COUNT(*) as kinase_count FROM kinase_scores GROUP BY site_id)"
        ).fetchone()[0]

        top_kinase = self.conn.execute(
            "SELECT kinase_name, COUNT(*) as target_count FROM kinase_scores GROUP BY kinase_name ORDER BY target_count DESC LIMIT 1"
        ).fetchone()

        print("\nKey Statistics:")
        print(f"{'Significant sites (FDR<0.05)':35}: {significant_sites:,}")
        print(f"{'Unique genes with phosphosites':35}: {unique_genes:,}")
        print(f"{'Unique kinases':35}: {unique_kinases:,}")
        print(f"{'Avg kinases per site':35}: {avg_kinases_per_site:.1f}")
        print(f"{'Top kinase by targets':35}: {top_kinase[0]} ({top_kinase[1]:,} sites)")

        # Check for potential issues
        orphaned_phosphosites = self.conn.execute("""
            SELECT COUNT(*) FROM phosphosites p 
            WHERE NOT EXISTS (SELECT 1 FROM proteins pr WHERE pr.gene_symbol = p.gene_symbol)
        """).fetchone()[0]

        orphaned_memberships = self.conn.execute("""
            SELECT COUNT(*) FROM gene_set_memberships gsm 
            WHERE NOT EXISTS (SELECT 1 FROM proteins p WHERE p.gene_symbol = gsm.gene_symbol)
        """).fetchone()[0]

        orphaned_kinase_scores = self.conn.execute("""
            SELECT COUNT(*) FROM kinase_scores ks 
            WHERE NOT EXISTS (SELECT 1 FROM phosphosites p WHERE p.site_id = ks.site_id)
        """).fetchone()[0]

        if orphaned_phosphosites > 0 or orphaned_memberships > 0 or orphaned_kinase_scores > 0:
            print(f"\n⚠️  Data Issues:")
            if orphaned_phosphosites > 0:
                print(f"   Phosphosites without proteins: {orphaned_phosphosites}")
            if orphaned_memberships > 0:
                print(f"   Memberships without proteins: {orphaned_memberships}")
            if orphaned_kinase_scores > 0:
                print(f"   Kinase scores without phosphosites: {orphaned_kinase_scores}")
        else:
            print(f"\n✅ Data integrity checks passed")

        # Database size
        db_size = self.db_path.stat().st_size / (1024 ** 3)  # GB
        print(f"\nDatabase file size: {db_size:.2f} GB")
        print("=" * 60)


def main():
    """Main function to build the database."""
    # File paths - update these to match your system
    phosphosite_file = "/Users/davidvanderwall/Desktop/Combined_ST_Data_Predictions_Kinases_Structures.feather"
    metadata_file = "/Users/davidvanderwall/Desktop/Protein_Level_MetaData.feather"
    db_path = "kinoplex.db"

    try:
        # Build the database
        with KinoPlexDatabaseBuilder(phosphosite_file, metadata_file, db_path) as builder:
            builder.create_schema()
            builder.load_data()
            builder.create_views()
            builder.validate_database()

        logger.info(f"KinoPlex database successfully created: {db_path}")
        print(f"\n🎉 Database ready! You can now connect to: {db_path}")

        # Sample queries to test
        print(f"\n📋 Sample queries to test your database:")
        print(f"   - Find TP53 phosphosites: SELECT * FROM phosphosites WHERE gene_symbol = 'TP53';")
        print(f"   - Get site with all kinase scores: SELECT * FROM phosphosites WHERE site_id = 'P04637_15';")
        print(f"   - High confidence sites: SELECT * FROM high_confidence_sites LIMIT 10;")
        print(f"   - Pathway analysis: SELECT * FROM phosphosite_pathways WHERE gs_collection = 'Hallmark';")
        print(
            f"   - Top AKT1 targets: SELECT site_id, gene_symbol, AKT1_MotifScore FROM phosphosites WHERE AKT1_MotifScore > 0.8 ORDER BY AKT1_MotifScore DESC;")

    except Exception as e:
        logger.error(f"Failed to build database: {e}")
        raise


if __name__ == "__main__":
    main()