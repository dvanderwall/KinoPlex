#!/usr/bin/env python3
"""
KinoPlex Database Enhancer - Add Kinase Mapping Tables from Feather Files

This script adds protein-kinase and pathway-kinase mapping tables to the existing
KinoPlex database, enabling direct kinase enrichment lookups at both protein and
pathway levels.

Usage:
    python add_kinase_mappings.py --protein-feather protein_level_kinase_enrichment.feather \
                                   --pathway-feather pathway_kinase_integrated_scoremapping.feather

Requirements:
    pandas, sqlite3 (built-in), tqdm (optional for progress bars)
"""

import pandas as pd
import sqlite3
import argparse
import logging
from pathlib import Path
from typing import Optional
import sys

# Optional: Progress bars
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


class KinaseMappingAdder:
    """Add kinase mapping tables to existing KinoPlex database."""

    def __init__(self,
                 db_path: str = "kinoplex.db",
                 protein_feather: Optional[str] = None,
                 pathway_feather: Optional[str] = None):
        """
        Initialize the kinase mapping adder.

        Args:
            db_path: Path to existing KinoPlex database
            protein_feather: Path to protein-kinase enrichment feather file
            pathway_feather: Path to pathway-kinase mapping feather file
        """
        self.db_path = Path(db_path)
        self.protein_feather = Path(protein_feather) if protein_feather else None
        self.pathway_feather = Path(pathway_feather) if pathway_feather else None

        # Validate database exists
        if not self.db_path.exists():
            raise FileNotFoundError(f"Database not found: {db_path}")

        # Validate feather files
        if self.protein_feather and not self.protein_feather.exists():
            raise FileNotFoundError(f"Protein feather file not found: {protein_feather}")
        if self.pathway_feather and not self.pathway_feather.exists():
            raise FileNotFoundError(f"Pathway feather file not found: {pathway_feather}")

        self.conn = None

    def __enter__(self):
        """Context manager entry."""
        self.conn = sqlite3.connect(str(self.db_path))
        self.conn.execute("PRAGMA journal_mode=WAL")
        self.conn.execute("PRAGMA synchronous=NORMAL")
        self.conn.execute("PRAGMA temp_store=MEMORY")
        self.conn.execute("PRAGMA mmap_size=268435456")  # 256MB mmap
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        if self.conn:
            self.conn.close()

    def create_protein_kinase_table(self) -> None:
        """Create the protein-kinase enrichment table."""
        logger.info("Creating protein_kinase_enrichment table...")

        # Drop existing table if it exists
        self.conn.execute("DROP TABLE IF EXISTS protein_kinase_enrichment")
        self.conn.execute("DROP VIEW IF EXISTS top_protein_kinases")
        self.conn.execute("DROP VIEW IF EXISTS kinase_summary")
        self.conn.execute("DROP VIEW IF EXISTS integrated_kinase_predictions")

        schema_sql = """
        CREATE TABLE protein_kinase_enrichment (
            enrichment_id INTEGER PRIMARY KEY AUTOINCREMENT,
            gene_symbol VARCHAR(50) NOT NULL,
            kinase VARCHAR(50) NOT NULL,
            n_sites INTEGER,
            mean_score FLOAT,
            max_score FLOAT,
            z_score_kinase FLOAT,
            fold_over_median FLOAT,
            site_enrichment FLOAT,
            site_fraction_protein FLOAT,
            n_above_q75 INTEGER,
            n_above_q90 INTEGER,
            mean_rank FLOAT,
            robust_score FLOAT,
            enrichment_score FLOAT,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,

            FOREIGN KEY (gene_symbol) REFERENCES proteins(gene_symbol),
            UNIQUE(gene_symbol, kinase)
        );

        -- Create indexes for optimal query performance
        CREATE INDEX idx_protein_kinase_gene ON protein_kinase_enrichment(gene_symbol);
        CREATE INDEX idx_protein_kinase_kinase ON protein_kinase_enrichment(kinase);
        CREATE INDEX idx_protein_kinase_enrichment ON protein_kinase_enrichment(enrichment_score DESC);
        CREATE INDEX idx_protein_kinase_zscore ON protein_kinase_enrichment(z_score_kinase DESC);
        CREATE INDEX idx_protein_kinase_combined ON protein_kinase_enrichment(gene_symbol, enrichment_score DESC);
        """

        self.conn.executescript(schema_sql)
        self.conn.commit()
        logger.info("protein_kinase_enrichment table created successfully")

    def create_pathway_kinase_table(self) -> None:
        """Create the pathway-kinase mapping table."""
        logger.info("Creating pathway_kinase_mapping table...")

        # Drop existing table and views if they exist
        self.conn.execute("DROP TABLE IF EXISTS pathway_kinase_mapping")
        self.conn.execute("DROP VIEW IF EXISTS top_pathway_kinases")

        schema_sql = """
        CREATE TABLE pathway_kinase_mapping (
            mapping_id INTEGER PRIMARY KEY AUTOINCREMENT,
            gs_name VARCHAR(200) NOT NULL,
            kinase VARCHAR(50) NOT NULL,

            -- Core enrichment metrics
            n_genes_with_kinase INTEGER,
            n_total_sites INTEGER,
            mean_enrichment FLOAT,
            median_enrichment FLOAT,
            max_enrichment FLOAT,
            mean_z_score FLOAT,
            prop_significant FLOAT,
            kinase_mean_bg FLOAT,

            -- Gene coverage metrics
            n_genes INTEGER,
            gene_coverage FLOAT,
            total_genes_x INTEGER,
            total_genes_y INTEGER,

            -- Statistical scores
            fisher_p FLOAT,
            enrichment_ratio FLOAT,
            score_zscore_only FLOAT,
            score_integrated FLOAT,
            score_ratio_based FLOAT,

            -- Rankings
            rank_integrated INTEGER,
            rank_zscore INTEGER,
            rank_ratio INTEGER,

            -- Site distribution metrics
            mean_sites_per_gene FLOAT,
            median_sites_per_gene INTEGER,
            sd_sites_per_gene FLOAT,
            mean_enrichment_all FLOAT,
            sd_enrichment_all FLOAT,

            -- Quantile metrics
            q75_enrichment FLOAT,
            q90_enrichment FLOAT,
            n_genes_above_kinase_q75 INTEGER,
            n_genes_above_kinase_q90 INTEGER,
            n_sites_above_kinase_q75 INTEGER,
            n_sites_above_kinase_q90 INTEGER,

            -- Additional enrichment metrics
            total_sites_pathway_kinase INTEGER,
            sites_vs_background_fold FLOAT,
            prop_genes_above_q75 FLOAT,
            prop_genes_above_q90 FLOAT,
            prop_sites_above_q75 FLOAT,
            prop_sites_above_q90 FLOAT,
            fold_enrichment_genes_q90 FLOAT,
            fold_enrichment_sites_q90 FLOAT,
            cohens_d FLOAT,

            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,

            UNIQUE(gs_name, kinase)
        );

        -- Create indexes for optimal query performance
        CREATE INDEX idx_pathway_kinase_pathway ON pathway_kinase_mapping(gs_name);
        CREATE INDEX idx_pathway_kinase_kinase ON pathway_kinase_mapping(kinase);
        CREATE INDEX idx_pathway_kinase_integrated ON pathway_kinase_mapping(score_integrated DESC);
        CREATE INDEX idx_pathway_kinase_zscore ON pathway_kinase_mapping(score_zscore_only DESC);
        CREATE INDEX idx_pathway_kinase_ratio ON pathway_kinase_mapping(score_ratio_based DESC);
        CREATE INDEX idx_pathway_kinase_rank ON pathway_kinase_mapping(rank_integrated);
        CREATE INDEX idx_pathway_kinase_combined ON pathway_kinase_mapping(gs_name, score_integrated DESC);
        CREATE INDEX idx_pathway_kinase_fisher ON pathway_kinase_mapping(fisher_p);
        """

        self.conn.executescript(schema_sql)
        self.conn.commit()
        logger.info("pathway_kinase_mapping table created successfully")

    def load_protein_kinase_data(self) -> None:
        """Load protein-kinase enrichment data from feather file."""
        if not self.protein_feather:
            logger.warning("No protein feather file provided, skipping protein-kinase data")
            return

        logger.info(f"Loading protein-kinase data from {self.protein_feather}")

        # Read feather file
        df = pd.read_feather(self.protein_feather)
        logger.info(f"Read {len(df)} records from feather file")

        # Drop unnamed index column if present
        for col in df.columns:
            if col.startswith('Unnamed:') or col == '':
                df = df.drop(col, axis=1)
                logger.info(f"Dropped index column: {col}")

        # Validate required columns
        required_cols = ['gene_symbol', 'kinase']
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            raise ValueError(f"Missing required columns: {missing_cols}")

        # Check for duplicates
        duplicates = df.duplicated(subset=['gene_symbol', 'kinase']).sum()
        if duplicates > 0:
            logger.warning(f"Found {duplicates} duplicate (gene_symbol, kinase) pairs, keeping first occurrence")
            df = df.drop_duplicates(subset=['gene_symbol', 'kinase'], keep='first')

        # Validate gene symbols exist in database
        existing_genes = pd.read_sql_query(
            "SELECT DISTINCT gene_symbol FROM proteins",
            self.conn
        )['gene_symbol'].tolist()

        invalid_genes = df[~df['gene_symbol'].isin(existing_genes)]['gene_symbol'].unique()
        if len(invalid_genes) > 0:
            logger.warning(f"Found {len(invalid_genes)} gene symbols not in database, these will be skipped")
            if len(invalid_genes) <= 10:
                logger.info(f"Invalid genes: {list(invalid_genes)}")
            else:
                logger.info(f"First 10 invalid genes: {list(invalid_genes[:10])}...")
            df = df[df['gene_symbol'].isin(existing_genes)]

        logger.info(f"Loading {len(df)} protein-kinase enrichment records to database")

        # Load to database in chunks for better performance
        chunk_size = 1000
        total_chunks = len(df) // chunk_size + 1

        for i in tqdm(range(total_chunks), desc="Loading protein-kinase data"):
            start_idx = i * chunk_size
            end_idx = min((i + 1) * chunk_size, len(df))
            chunk = df.iloc[start_idx:end_idx]

            if not chunk.empty:
                chunk.to_sql('protein_kinase_enrichment', self.conn, if_exists='append', index=False)

        self.conn.commit()
        logger.info(f"Successfully loaded {len(df)} protein-kinase enrichment records")

    def load_pathway_kinase_data(self) -> None:
        """Load pathway-kinase mapping data from feather file."""
        if not self.pathway_feather:
            logger.warning("No pathway feather file provided, skipping pathway-kinase data")
            return

        logger.info(f"Loading pathway-kinase data from {self.pathway_feather}")

        # Read feather file
        df = pd.read_feather(self.pathway_feather)
        logger.info(f"Read {len(df)} records from feather file")

        # Drop unnamed index column if present
        for col in df.columns:
            if col.startswith('Unnamed:') or col == '':
                df = df.drop(col, axis=1)
                logger.info(f"Dropped index column: {col}")

        # Handle column name variations (dots in column names)
        column_mapping = {
            'total_genes.x': 'total_genes_x',
            'total_genes.y': 'total_genes_y'
        }
        df = df.rename(columns=column_mapping)

        # Validate required columns
        required_cols = ['gs_name', 'kinase']
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            raise ValueError(f"Missing required columns: {missing_cols}")

        # Check for duplicates
        duplicates = df.duplicated(subset=['gs_name', 'kinase']).sum()
        if duplicates > 0:
            logger.warning(f"Found {duplicates} duplicate (gs_name, kinase) pairs, keeping first occurrence")
            df = df.drop_duplicates(subset=['gs_name', 'kinase'], keep='first')

        logger.info(f"Loading {len(df)} pathway-kinase mapping records to database")

        # Select only columns that exist in our table schema
        table_columns = [
            'gs_name', 'kinase', 'n_genes_with_kinase', 'n_total_sites',
            'mean_enrichment', 'median_enrichment', 'max_enrichment', 'mean_z_score',
            'prop_significant', 'kinase_mean_bg', 'n_genes', 'gene_coverage',
            'total_genes_x', 'total_genes_y', 'fisher_p', 'enrichment_ratio',
            'score_zscore_only', 'score_integrated', 'score_ratio_based',
            'rank_integrated', 'rank_zscore', 'rank_ratio', 'mean_sites_per_gene',
            'median_sites_per_gene', 'sd_sites_per_gene', 'mean_enrichment_all',
            'sd_enrichment_all', 'q75_enrichment', 'q90_enrichment',
            'n_genes_above_kinase_q75', 'n_genes_above_kinase_q90',
            'n_sites_above_kinase_q75', 'n_sites_above_kinase_q90',
            'total_sites_pathway_kinase', 'sites_vs_background_fold',
            'prop_genes_above_q75', 'prop_genes_above_q90', 'prop_sites_above_q75',
            'prop_sites_above_q90', 'fold_enrichment_genes_q90',
            'fold_enrichment_sites_q90', 'cohens_d'
        ]

        # Filter to only columns that exist in the dataframe
        available_columns = [col for col in table_columns if col in df.columns]
        missing_from_df = [col for col in table_columns if col not in df.columns]

        if missing_from_df:
            logger.warning(f"Some columns not found in feather file: {missing_from_df}")

        df_filtered = df[available_columns]

        # Load to database in chunks for better performance
        chunk_size = 1000
        total_chunks = len(df_filtered) // chunk_size + 1

        for i in tqdm(range(total_chunks), desc="Loading pathway-kinase data"):
            start_idx = i * chunk_size
            end_idx = min((i + 1) * chunk_size, len(df_filtered))
            chunk = df_filtered.iloc[start_idx:end_idx]

            if not chunk.empty:
                chunk.to_sql('pathway_kinase_mapping', self.conn, if_exists='append', index=False)

        self.conn.commit()
        logger.info(f"Successfully loaded {len(df_filtered)} pathway-kinase mapping records")

    def create_views(self) -> None:
        """Create useful views for the new tables."""
        logger.info("Creating kinase mapping views...")

        views_sql = """
        -- Top kinases per protein
        CREATE VIEW IF NOT EXISTS top_protein_kinases AS
        SELECT 
            pk.gene_symbol,
            pk.kinase,
            pk.enrichment_score,
            pk.z_score_kinase,
            pk.n_sites,
            pk.n_above_q90,
            p.protein_count
        FROM protein_kinase_enrichment pk
        JOIN proteins p ON pk.gene_symbol = p.gene_symbol
        WHERE pk.enrichment_score IS NOT NULL
        ORDER BY pk.gene_symbol, pk.enrichment_score DESC;

        -- Top kinases per pathway
        CREATE VIEW IF NOT EXISTS top_pathway_kinases AS
        SELECT 
            pkm.gs_name,
            pkm.kinase,
            pkm.score_integrated,
            pkm.rank_integrated,
            pkm.n_genes_with_kinase,
            pkm.n_total_sites,
            pkm.fisher_p,
            gs.gs_collection,
            gs.gs_description
        FROM pathway_kinase_mapping pkm
        LEFT JOIN gene_sets gs ON pkm.gs_name = gs.gs_name
        WHERE pkm.score_integrated IS NOT NULL
        ORDER BY pkm.gs_name, pkm.score_integrated DESC;

        -- Kinase summary statistics
        CREATE VIEW IF NOT EXISTS kinase_summary AS
        SELECT 
            kinase,
            COUNT(DISTINCT pk.gene_symbol) as n_target_proteins,
            COUNT(DISTINCT pkm.gs_name) as n_target_pathways,
            AVG(pk.enrichment_score) as avg_protein_enrichment,
            AVG(pkm.score_integrated) as avg_pathway_score,
            MAX(pk.enrichment_score) as max_protein_enrichment,
            MAX(pkm.score_integrated) as max_pathway_score
        FROM (
            SELECT DISTINCT kinase FROM protein_kinase_enrichment
            UNION 
            SELECT DISTINCT kinase FROM pathway_kinase_mapping
        ) k
        LEFT JOIN protein_kinase_enrichment pk ON k.kinase = pk.kinase
        LEFT JOIN pathway_kinase_mapping pkm ON k.kinase = pkm.kinase
        GROUP BY k.kinase;

        -- Integrated kinase view linking all three levels
        CREATE VIEW IF NOT EXISTS integrated_kinase_predictions AS
        SELECT 
            p.gene_symbol,
            p.site_id,
            p.position,
            p.residue_type,
            p.predicted_prob_calibrated,
            pk.kinase,
            pk.enrichment_score as protein_enrichment,
            pk.z_score_kinase as protein_z_score,
            pkm.score_integrated as pathway_score,
            pkm.rank_integrated as pathway_rank
        FROM phosphosites p
        LEFT JOIN protein_kinase_enrichment pk ON p.gene_symbol = pk.gene_symbol
        LEFT JOIN gene_set_memberships gsm ON p.gene_symbol = gsm.gene_symbol
        LEFT JOIN gene_sets gs ON gsm.gs_id = gs.gs_id
        LEFT JOIN pathway_kinase_mapping pkm ON gs.gs_name = pkm.gs_name AND pk.kinase = pkm.kinase
        WHERE p.significant_fdr05 = TRUE;
        """

        self.conn.executescript(views_sql)
        self.conn.commit()
        logger.info("Views created successfully")

    def validate_additions(self) -> None:
        """Validate the newly added tables and provide statistics."""
        logger.info("Validating new tables...")

        print("\n" + "=" * 60)
        print("KinoPlex Database Enhancement Summary")
        print("=" * 60)

        # Check protein-kinase table
        if self.table_exists('protein_kinase_enrichment'):
            protein_count = self.conn.execute(
                "SELECT COUNT(*) FROM protein_kinase_enrichment"
            ).fetchone()[0]

            unique_proteins = self.conn.execute(
                "SELECT COUNT(DISTINCT gene_symbol) FROM protein_kinase_enrichment"
            ).fetchone()[0]

            unique_kinases_protein = self.conn.execute(
                "SELECT COUNT(DISTINCT kinase) FROM protein_kinase_enrichment"
            ).fetchone()[0]

            top_enriched = self.conn.execute("""
                SELECT gene_symbol, kinase, enrichment_score 
                FROM protein_kinase_enrichment 
                WHERE enrichment_score IS NOT NULL
                ORDER BY enrichment_score DESC 
                LIMIT 3
            """).fetchall()

            print(f"\nProtein-Kinase Enrichment Table:")
            print(f"  Total records: {protein_count:,}")
            print(f"  Unique proteins: {unique_proteins:,}")
            print(f"  Unique kinases: {unique_kinases_protein:,}")
            if top_enriched:
                print(f"  Top enrichments:")
                for gene, kinase, score in top_enriched:
                    if score is not None:
                        print(f"    {gene} - {kinase}: {score:.3f}")

        # Check pathway-kinase table
        if self.table_exists('pathway_kinase_mapping'):
            pathway_count = self.conn.execute(
                "SELECT COUNT(*) FROM pathway_kinase_mapping"
            ).fetchone()[0]

            unique_pathways = self.conn.execute(
                "SELECT COUNT(DISTINCT gs_name) FROM pathway_kinase_mapping"
            ).fetchone()[0]

            unique_kinases_pathway = self.conn.execute(
                "SELECT COUNT(DISTINCT kinase) FROM pathway_kinase_mapping"
            ).fetchone()[0]

            top_integrated = self.conn.execute("""
                SELECT gs_name, kinase, score_integrated 
                FROM pathway_kinase_mapping 
                WHERE score_integrated IS NOT NULL
                ORDER BY score_integrated DESC 
                LIMIT 3
            """).fetchall()

            print(f"\nPathway-Kinase Mapping Table:")
            print(f"  Total records: {pathway_count:,}")
            print(f"  Unique pathways: {unique_pathways:,}")
            print(f"  Unique kinases: {unique_kinases_pathway:,}")
            if top_integrated:
                print(f"  Top integrated scores:")
                for pathway, kinase, score in top_integrated:
                    pathway_short = pathway[:50] + "..." if len(pathway) > 50 else pathway
                    print(f"    {pathway_short} - {kinase}: {score:.3f}")

        # Test integration with existing tables
        if self.table_exists('protein_kinase_enrichment'):
            integrated_test = self.conn.execute("""
                SELECT COUNT(*) 
                FROM protein_kinase_enrichment pk
                JOIN proteins p ON pk.gene_symbol = p.gene_symbol
            """).fetchone()[0]

            print(f"\nIntegration Tests:")
            print(f"  Protein-kinase records with valid proteins: {integrated_test:,}")

        print("=" * 60)

    def table_exists(self, table_name: str) -> bool:
        """Check if a table exists in the database."""
        result = self.conn.execute(
            "SELECT name FROM sqlite_master WHERE type='table' AND name=?",
            (table_name,)
        ).fetchone()
        return result is not None


def main():
    """Main function to add kinase mappings to database."""
    parser = argparse.ArgumentParser(
        description="Add kinase mapping tables to KinoPlex database from feather files"
    )
    parser.add_argument(
        "--db",
        default="kinoplex.db",
        help="Path to KinoPlex database (default: kinoplex.db)"
    )
    parser.add_argument(
        "--protein-feather",
        help="Path to protein-kinase enrichment feather file"
    )
    parser.add_argument(
        "--pathway-feather",
        help="Path to pathway-kinase mapping feather file"
    )
    parser.add_argument(
        "--skip-views",
        action="store_true",
        help="Skip creating database views"
    )

    args = parser.parse_args()

    if not args.protein_feather and not args.pathway_feather:
        logger.error("At least one feather file must be provided")
        sys.exit(1)

    try:
        with KinaseMappingAdder(
                db_path=args.db,
                protein_feather=args.protein_feather,
                pathway_feather=args.pathway_feather
        ) as adder:

            # Create tables
            if args.protein_feather:
                adder.create_protein_kinase_table()
                adder.load_protein_kinase_data()

            if args.pathway_feather:
                adder.create_pathway_kinase_table()
                adder.load_pathway_kinase_data()

            # Create views
            if not args.skip_views:
                adder.create_views()

            # Validate
            adder.validate_additions()

        logger.info("✅ Kinase mapping tables successfully added to database")

        print("\n📊 Sample queries to test the new tables:")
        print("  - Top kinases for TP53:")
        print(
            "    SELECT * FROM protein_kinase_enrichment WHERE gene_symbol = 'TP53' ORDER BY enrichment_score DESC LIMIT 5;")
        print("\n  - Top pathways for AKT1:")
        print(
            "    SELECT gs_name, score_integrated FROM pathway_kinase_mapping WHERE kinase = 'AKT1' ORDER BY score_integrated DESC LIMIT 5;")
        print("\n  - Kinase summary:")
        print("    SELECT * FROM kinase_summary ORDER BY n_target_proteins DESC LIMIT 10;")
        print("\n  - Integrated view for EGFR:")
        print("    SELECT * FROM top_protein_kinases WHERE gene_symbol = 'EGFR' LIMIT 5;")

    except Exception as e:
        logger.error(f"Failed to add kinase mappings: {e}")
        raise


if __name__ == "__main__":
    main()