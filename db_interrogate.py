#!/usr/bin/env python3
"""
KinoPlex Database Interrogation Script

This script provides comprehensive analysis of the KinoPlex database structure,
contents, and data quality to help diagnose API issues and optimize queries.
"""

import sqlite3
import pandas as pd
from pathlib import Path
import json


class DatabaseInterrogator:
    """Comprehensive database analysis tool."""

    def __init__(self, db_path: str = "kinoplex.db"):
        self.db_path = Path(db_path)
        if not self.db_path.exists():
            raise FileNotFoundError(f"Database not found: {db_path}")

        self.conn = sqlite3.connect(str(self.db_path))
        self.conn.row_factory = sqlite3.Row

    def analyze_database(self):
        """Perform comprehensive database analysis."""
        print("=" * 80)
        print("KINOPLEX DATABASE COMPREHENSIVE ANALYSIS")
        print("=" * 80)

        # Basic database info
        self.get_database_info()

        # Table analysis
        self.analyze_all_tables()

        # Enrichment tables deep dive
        self.analyze_enrichment_tables()

        # Data quality checks
        self.check_data_quality()

        # Query performance analysis
        self.analyze_query_patterns()

        # Sample data extraction
        self.extract_sample_data()

    def get_database_info(self):
        """Get basic database information."""
        print("\n" + "=" * 60)
        print("DATABASE OVERVIEW")
        print("=" * 60)

        # Database size
        db_size = self.db_path.stat().st_size / (1024 ** 3)
        print(f"Database file size: {db_size:.2f} GB")

        # SQLite version and settings
        version = self.conn.execute("SELECT sqlite_version()").fetchone()[0]
        print(f"SQLite version: {version}")

        # Journal mode
        journal_mode = self.conn.execute("PRAGMA journal_mode").fetchone()[0]
        print(f"Journal mode: {journal_mode}")

        # Page size
        page_size = self.conn.execute("PRAGMA page_size").fetchone()[0]
        print(f"Page size: {page_size} bytes")

        # Total pages
        page_count = self.conn.execute("PRAGMA page_count").fetchone()[0]
        print(f"Total pages: {page_count:,}")

    def analyze_all_tables(self):
        """Analyze all tables in the database."""
        print("\n" + "=" * 60)
        print("TABLE STRUCTURE ANALYSIS")
        print("=" * 60)

        # Get all tables
        tables_query = """
        SELECT name, type FROM sqlite_master 
        WHERE type IN ('table', 'view') 
        ORDER BY type, name
        """
        tables = self.conn.execute(tables_query).fetchall()

        print(f"Found {len(tables)} tables and views:")

        for table in tables:
            table_name = table['name']
            table_type = table['type']

            print(f"\n{table_type.upper()}: {table_name}")
            print("-" * 40)

            # Get column information
            columns_info = self.conn.execute(f"PRAGMA table_info({table_name})").fetchall()

            if table_type == 'table':
                # Row count for tables only
                try:
                    row_count = self.conn.execute(f"SELECT COUNT(*) as count FROM {table_name}").fetchone()['count']
                    print(f"Rows: {row_count:,}")
                except:
                    print("Rows: Unable to count")

            print(f"Columns ({len(columns_info)}):")
            for col in columns_info[:10]:  # Show first 10 columns
                col_info = f"  {col['name']} ({col['type']})"
                if col['pk']:
                    col_info += " [PRIMARY KEY]"
                if col['notnull']:
                    col_info += " [NOT NULL]"
                print(col_info)

            if len(columns_info) > 10:
                print(f"  ... and {len(columns_info) - 10} more columns")

            # Special analysis for key tables
            if table_name == 'phosphosites':
                self.analyze_phosphosites_table()
            elif table_name == 'protein_kinase_enrichment':
                self.analyze_protein_kinase_table()
            elif table_name == 'pathway_kinase_mapping':
                self.analyze_pathway_kinase_table()

    def analyze_phosphosites_table(self):
        """Detailed analysis of phosphosites table."""
        print("\n  PHOSPHOSITES DETAILED ANALYSIS:")

        # Count kinase score columns
        columns = self.conn.execute("PRAGMA table_info(phosphosites)").fetchall()
        kinase_columns = [col['name'] for col in columns if col['name'].endswith('_MotifScore')]

        print(f"  - Kinase score columns: {len(kinase_columns)}")
        if len(kinase_columns) > 0:
            print(f"  - First 5 kinases: {kinase_columns[:5]}")
            print(f"  - Last 5 kinases: {kinase_columns[-5:]}")

        # Check for data in kinase columns
        if kinase_columns:
            sample_kinase = kinase_columns[0]
            non_null_scores = self.conn.execute(
                f"SELECT COUNT(*) as count FROM phosphosites WHERE {sample_kinase} IS NOT NULL AND {sample_kinase} > 0"
            ).fetchone()['count']
            print(f"  - Sites with {sample_kinase} scores > 0: {non_null_scores:,}")

        # Significance analysis
        try:
            sig_fdr05 = self.conn.execute(
                "SELECT COUNT(*) as count FROM phosphosites WHERE significant_fdr05 = 1"
            ).fetchone()['count']
            print(f"  - Significant sites (FDR < 0.05): {sig_fdr05:,}")
        except:
            print("  - No significance column found")

        # Confidence distribution
        try:
            conf_stats = self.conn.execute("""
                SELECT 
                    AVG(predicted_prob_calibrated) as avg_conf,
                    MIN(predicted_prob_calibrated) as min_conf,
                    MAX(predicted_prob_calibrated) as max_conf,
                    COUNT(*) as total
                FROM phosphosites 
                WHERE predicted_prob_calibrated IS NOT NULL
            """).fetchone()

            if conf_stats['total'] > 0:
                print(f"  - Confidence range: {conf_stats['min_conf']:.3f} - {conf_stats['max_conf']:.3f}")
                print(f"  - Average confidence: {conf_stats['avg_conf']:.3f}")
        except:
            print("  - No confidence scores found")

    def analyze_protein_kinase_table(self):
        """Detailed analysis of protein_kinase_enrichment table."""
        print("\n  PROTEIN-KINASE ENRICHMENT DETAILED ANALYSIS:")

        # Unique counts
        unique_proteins = self.conn.execute(
            "SELECT COUNT(DISTINCT gene_symbol) as count FROM protein_kinase_enrichment"
        ).fetchone()['count']

        unique_kinases = self.conn.execute(
            "SELECT COUNT(DISTINCT kinase) as count FROM protein_kinase_enrichment"
        ).fetchone()['count']

        print(f"  - Unique proteins: {unique_proteins:,}")
        print(f"  - Unique kinases: {unique_kinases:,}")

        # Score distribution
        score_stats = self.conn.execute("""
            SELECT 
                AVG(enrichment_score) as avg_score,
                MIN(enrichment_score) as min_score,
                MAX(enrichment_score) as max_score,
                COUNT(*) as total_records
            FROM protein_kinase_enrichment 
            WHERE enrichment_score IS NOT NULL
        """).fetchone()

        if score_stats['total_records'] > 0:
            print(f"  - Enrichment score range: {score_stats['min_score']:.3f} - {score_stats['max_score']:.3f}")
            print(f"  - Average enrichment score: {score_stats['avg_score']:.3f}")

        # Top enriched protein-kinase pairs
        top_enriched = self.conn.execute("""
            SELECT gene_symbol, kinase, enrichment_score 
            FROM protein_kinase_enrichment 
            WHERE enrichment_score IS NOT NULL
            ORDER BY enrichment_score DESC 
            LIMIT 3
        """).fetchall()

        print("  - Top 3 enriched pairs:")
        for row in top_enriched:
            print(f"    {row['gene_symbol']} - {row['kinase']}: {row['enrichment_score']:.3f}")

    def analyze_pathway_kinase_table(self):
        """Detailed analysis of pathway_kinase_mapping table."""
        print("\n  PATHWAY-KINASE MAPPING DETAILED ANALYSIS:")

        # Unique counts
        unique_pathways = self.conn.execute(
            "SELECT COUNT(DISTINCT gs_name) as count FROM pathway_kinase_mapping"
        ).fetchone()['count']

        unique_kinases = self.conn.execute(
            "SELECT COUNT(DISTINCT kinase) as count FROM pathway_kinase_mapping"
        ).fetchone()['count']

        print(f"  - Unique pathways: {unique_pathways:,}")
        print(f"  - Unique kinases: {unique_kinases:,}")

        # Score distribution
        score_stats = self.conn.execute("""
            SELECT 
                AVG(score_integrated) as avg_score,
                MIN(score_integrated) as min_score,
                MAX(score_integrated) as max_score,
                COUNT(*) as total_records
            FROM pathway_kinase_mapping 
            WHERE score_integrated IS NOT NULL
        """).fetchone()

        if score_stats['total_records'] > 0:
            print(f"  - Integrated score range: {score_stats['min_score']:.3f} - {score_stats['max_score']:.3f}")
            print(f"  - Average integrated score: {score_stats['avg_score']:.3f}")

        # Top enriched pathway-kinase pairs
        top_enriched = self.conn.execute("""
            SELECT gs_name, kinase, score_integrated 
            FROM pathway_kinase_mapping 
            WHERE score_integrated IS NOT NULL
            ORDER BY score_integrated DESC 
            LIMIT 3
        """).fetchall()

        print("  - Top 3 enriched pairs:")
        for row in top_enriched:
            pathway_short = row['gs_name'][:50] + "..." if len(row['gs_name']) > 50 else row['gs_name']
            print(f"    {pathway_short} - {row['kinase']}: {row['score_integrated']:.3f}")

    def analyze_enrichment_tables(self):
        """Deep dive into enrichment table relationships."""
        print("\n" + "=" * 60)
        print("ENRICHMENT TABLES INTEGRATION ANALYSIS")
        print("=" * 60)

        # Check if enrichment tables exist
        enrichment_tables = ['protein_kinase_enrichment', 'pathway_kinase_mapping']
        existing_enrichment = []

        for table in enrichment_tables:
            exists = self.conn.execute(
                "SELECT name FROM sqlite_master WHERE type='table' AND name=?",
                (table,)
            ).fetchone()

            if exists:
                existing_enrichment.append(table)
                print(f"✓ {table} exists")
            else:
                print(f"✗ {table} MISSING")

        if not existing_enrichment:
            print("\nWARNING: No enrichment tables found!")
            print("The database may not have been updated with enrichment data.")
            return

        # Cross-table consistency checks
        if 'protein_kinase_enrichment' in existing_enrichment:
            print("\nProtein-Kinase Enrichment Integration:")

            # Check protein overlap with main tables
            protein_overlap = self.conn.execute("""
                SELECT COUNT(DISTINCT pk.gene_symbol) as enrichment_proteins,
                       COUNT(DISTINCT p.gene_symbol) as main_proteins,
                       COUNT(DISTINCT CASE WHEN pr.gene_symbol IS NOT NULL THEN pk.gene_symbol END) as overlap
                FROM protein_kinase_enrichment pk
                LEFT JOIN proteins pr ON pk.gene_symbol = pr.gene_symbol
                CROSS JOIN (SELECT COUNT(DISTINCT gene_symbol) as gene_symbol FROM proteins) p
            """).fetchone()

            print(f"  - Proteins in enrichment table: {protein_overlap['enrichment_proteins']:,}")
            print(f"  - Overlap with proteins table: {protein_overlap['overlap']:,}")

        if 'pathway_kinase_mapping' in existing_enrichment:
            print("\nPathway-Kinase Mapping Integration:")

            # Check pathway name consistency
            pathway_overlap = self.conn.execute("""
                SELECT COUNT(DISTINCT pkm.gs_name) as mapping_pathways,
                       COUNT(DISTINCT gs.gs_name) as main_pathways,
                       COUNT(DISTINCT CASE WHEN gs.gs_name IS NOT NULL THEN pkm.gs_name END) as overlap
                FROM pathway_kinase_mapping pkm
                LEFT JOIN gene_sets gs ON pkm.gs_name = gs.gs_name
                CROSS JOIN (SELECT COUNT(DISTINCT gs_name) as gs_name FROM gene_sets) g
            """).fetchone()

            print(f"  - Pathways in mapping table: {pathway_overlap['mapping_pathways']:,}")
            print(f"  - Overlap with gene_sets table: {pathway_overlap['overlap']:,}")

    def check_data_quality(self):
        """Perform data quality checks."""
        print("\n" + "=" * 60)
        print("DATA QUALITY ANALYSIS")
        print("=" * 60)

        # Check for orphaned records
        print("Checking for orphaned records:")

        # Orphaned phosphosites
        orphaned_phosphosites = self.conn.execute("""
            SELECT COUNT(*) as count FROM phosphosites p 
            WHERE NOT EXISTS (SELECT 1 FROM proteins pr WHERE pr.gene_symbol = p.gene_symbol)
        """).fetchone()['count']

        print(f"  - Phosphosites without proteins: {orphaned_phosphosites:,}")

        # Orphaned memberships
        orphaned_memberships = self.conn.execute("""
            SELECT COUNT(*) as count FROM gene_set_memberships gsm 
            WHERE NOT EXISTS (SELECT 1 FROM proteins p WHERE p.gene_symbol = gsm.gene_symbol)
        """).fetchone()['count']

        print(f"  - Memberships without proteins: {orphaned_memberships:,}")

        # Check enrichment table consistency
        try:
            orphaned_protein_enrichment = self.conn.execute("""
                SELECT COUNT(*) as count FROM protein_kinase_enrichment pk 
                WHERE NOT EXISTS (SELECT 1 FROM proteins p WHERE p.gene_symbol = pk.gene_symbol)
            """).fetchone()['count']
            print(f"  - Protein enrichments without proteins: {orphaned_protein_enrichment:,}")
        except:
            print("  - Protein enrichment table not found for consistency check")

        # Check for null values in critical columns
        print("\nNull value analysis:")
        critical_columns = [
            ('phosphosites', 'gene_symbol'),
            ('phosphosites', 'position'),
            ('phosphosites', 'predicted_prob_calibrated'),
            ('proteins', 'gene_symbol'),
            ('gene_sets', 'gs_name')
        ]

        for table, column in critical_columns:
            try:
                null_count = self.conn.execute(
                    f"SELECT COUNT(*) as count FROM {table} WHERE {column} IS NULL"
                ).fetchone()['count']
                print(f"  - {table}.{column}: {null_count:,} nulls")
            except:
                print(f"  - {table}.{column}: Column not found")

    def analyze_query_patterns(self):
        """Analyze common query patterns and performance."""
        print("\n" + "=" * 60)
        print("QUERY PERFORMANCE ANALYSIS")
        print("=" * 60)

        # Check indexes
        print("Existing indexes:")
        indexes = self.conn.execute("""
            SELECT name, tbl_name, sql FROM sqlite_master 
            WHERE type = 'index' AND name NOT LIKE 'sqlite_%'
            ORDER BY tbl_name, name
        """).fetchall()

        for idx in indexes:
            print(f"  - {idx['name']} on {idx['tbl_name']}")

        # Test key query patterns
        print("\nQuery pattern testing:")

        import time

        test_queries = [
            ("Protein lookup", "SELECT * FROM phosphosites WHERE gene_symbol = 'TP53'"),
            ("Pathway search", "SELECT * FROM gene_sets WHERE gs_name LIKE '%APOPTOSIS%' LIMIT 10"),
            ("Significance filter", "SELECT COUNT(*) FROM phosphosites WHERE significant_fdr05 = 1")
        ]

        for query_name, query in test_queries:
            try:
                start_time = time.time()
                result = self.conn.execute(query).fetchall()
                end_time = time.time()
                duration = (end_time - start_time) * 1000
                print(f"  - {query_name}: {duration:.1f}ms ({len(result)} rows)")
            except Exception as e:
                print(f"  - {query_name}: ERROR - {str(e)}")

    def extract_sample_data(self):
        """Extract sample data for API development."""
        print("\n" + "=" * 60)
        print("SAMPLE DATA EXTRACTION")
        print("=" * 60)

        # Sample pathway for testing
        sample_pathway = self.conn.execute("""
            SELECT gs_id, gs_name, gs_collection, gs_description 
            FROM gene_sets 
            WHERE gs_name LIKE '%COMPLEMENT%' 
            LIMIT 1
        """).fetchone()

        if sample_pathway:
            print(f"Sample pathway found: {sample_pathway['gs_name']} (ID: {sample_pathway['gs_id']})")

            # Get proteins in this pathway
            pathway_proteins = self.conn.execute("""
                SELECT COUNT(DISTINCT gene_symbol) as count
                FROM gene_set_memberships
                WHERE gs_id = ?
            """, (sample_pathway['gs_id'],)).fetchone()['count']

            print(f"  - Proteins in pathway: {pathway_proteins:,}")

            # Check if we have enrichment data for this pathway
            try:
                enrichment_data = self.conn.execute("""
                    SELECT COUNT(*) as count
                    FROM pathway_kinase_mapping
                    WHERE gs_name = ?
                """, (sample_pathway['gs_name'],)).fetchone()['count']

                print(f"  - Kinase mappings for pathway: {enrichment_data:,}")

                if enrichment_data > 0:
                    # Get sample kinase data
                    sample_kinases = self.conn.execute("""
                        SELECT kinase, score_integrated, rank_integrated
                        FROM pathway_kinase_mapping
                        WHERE gs_name = ? AND score_integrated IS NOT NULL
                        ORDER BY score_integrated DESC
                        LIMIT 5
                    """, (sample_pathway['gs_name'],)).fetchall()

                    print("  - Top 5 kinases for this pathway:")
                    for kinase in sample_kinases:
                        print(
                            f"    {kinase['kinase']}: {kinase['score_integrated']:.3f} (rank {kinase['rank_integrated']})")

            except Exception as e:
                print(f"  - Could not check enrichment data: {e}")

        else:
            print("No COMPLEMENT pathways found - testing with any pathway")
            any_pathway = self.conn.execute("""
                SELECT gs_id, gs_name, gs_collection 
                FROM gene_sets 
                LIMIT 1
            """).fetchone()

            if any_pathway:
                print(f"Using pathway: {any_pathway['gs_name']} (ID: {any_pathway['gs_id']})")

    def close(self):
        """Close database connection."""
        self.conn.close()


def main():
    """Run the database interrogation."""
    try:
        interrogator = DatabaseInterrogator("kinoplex.db")
        interrogator.analyze_database()
        interrogator.close()

        print("\n" + "=" * 80)
        print("DATABASE ANALYSIS COMPLETE")
        print("=" * 80)
        print("\nThis analysis will help identify:")
        print("1. Whether enrichment tables exist and are populated")
        print("2. Data consistency between tables")
        print("3. Query performance issues")
        print("4. Sample data for API testing")

    except FileNotFoundError:
        print("ERROR: kinoplex.db not found in current directory")
        print("Make sure you're running this script in the same directory as your database file")
    except Exception as e:
        print(f"ERROR: {e}")


if __name__ == "__main__":
    main()