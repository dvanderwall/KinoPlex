#!/usr/bin/env python3
"""
KinoPlex Pathway Analysis Debugger - FIXED VERSION
Comprehensive debugging tool for pathway analysis issues
"""

import sqlite3
import json
import sys
from pathlib import Path
from typing import Dict, List, Any
import pandas as pd


class PathwayAnalysisDebugger:
    def __init__(self, db_path: str):
        """Initialize the debugger with database connection."""
        self.db_path = Path(db_path)
        if not self.db_path.exists():
            print(f"❌ Database not found at: {db_path}")
            sys.exit(1)

        self.conn = sqlite3.connect(str(self.db_path))
        self.conn.row_factory = sqlite3.Row  # Enable column access by name
        print(f"✅ Connected to database: {db_path}")
        print("=" * 80)

    def debug_database_structure(self):
        """Check the actual database structure and column names."""
        print("\n🔍 DATABASE STRUCTURE ANALYSIS")
        print("-" * 80)

        # Get all tables
        tables = self.conn.execute(
            "SELECT name FROM sqlite_master WHERE type='table' ORDER BY name"
        ).fetchall()

        print("📊 Tables found in database:")
        for table in tables:
            count = self.conn.execute(f"SELECT COUNT(*) as cnt FROM {table['name']}").fetchone()['cnt']
            print(f"   - {table['name']}: {count:,} records")

        # Check phosphosites table structure in detail
        print("\n📋 Phosphosites Table Structure:")
        columns = self.conn.execute("PRAGMA table_info(phosphosites)").fetchall()

        # Group columns by type
        basic_cols = []
        kinase_cols = []
        structural_cols = []

        for col in columns:
            col_name = col['name']
            col_type = col['type']

            if '_MotifScore' in col_name:
                kinase_cols.append(col_name)
            elif col_name in ['site_id', 'gene_symbol', 'position', 'residue_type',
                              'predicted_prob_calibrated', 'qvalue', 'significant_fdr05']:
                basic_cols.append((col_name, col_type))
            else:
                structural_cols.append(col_name)

        print(f"\n   Basic columns ({len(basic_cols)}):")
        for col_name, col_type in basic_cols:
            print(f"      - {col_name} ({col_type})")

        print(f"\n   Kinase score columns found: {len(kinase_cols)}")
        if kinase_cols:
            print("   First 10 kinase columns:")
            for kinase in kinase_cols[:10]:
                print(f"      - {kinase}")
            if len(kinase_cols) > 10:
                print(f"      ... and {len(kinase_cols) - 10} more")
        else:
            print("   ⚠️  NO KINASE SCORE COLUMNS FOUND!")

        print(f"\n   Structural feature columns: {len(structural_cols)}")

        return kinase_cols

    def check_sample_data(self):
        """Check sample data from phosphosites table."""
        print("\n🔬 SAMPLE DATA CHECK")
        print("-" * 80)

        # Get a sample phosphosite record
        sample = self.conn.execute(
            "SELECT * FROM phosphosites LIMIT 1"
        ).fetchone()

        if not sample:
            print("❌ No phosphosite records found!")
            return

        print("Sample phosphosite record:")
        print(f"   site_id: {sample['site_id']}")
        print(f"   gene_symbol: {sample['gene_symbol']}")
        print(f"   position: {sample['position']}")
        print(f"   residue_type: {sample['residue_type']}")
        print(f"   predicted_prob_calibrated: {sample['predicted_prob_calibrated']}")
        print(f"   qvalue: {sample['qvalue']}")

        # Check if any kinase scores exist
        kinase_scores_found = []
        for key in sample.keys():
            if '_MotifScore' in key and sample[key] is not None:
                kinase_scores_found.append((key, sample[key]))

        if kinase_scores_found:
            print(f"\n   Kinase scores found in this record: {len(kinase_scores_found)}")
            for kinase, score in kinase_scores_found[:5]:
                print(f"      - {kinase}: {score:.4f}")
        else:
            print("\n   ⚠️  No kinase scores found in this record!")

    def check_pathway_data(self):
        """Check gene_sets and gene_set_memberships tables."""
        print("\n🗂️ PATHWAY DATA CHECK")
        print("-" * 80)

        # Check gene_sets
        gene_sets_count = self.conn.execute(
            "SELECT COUNT(*) as cnt FROM gene_sets"
        ).fetchone()['cnt']

        print(f"Total gene sets: {gene_sets_count:,}")

        # FIXED: Specify which table's gs_id we're using
        sample_sets = self.conn.execute(
            """SELECT gs.gs_id, gs.gs_name, gs.gs_collection, 
               COUNT(DISTINCT gsm.gene_symbol) as gene_count
               FROM gene_sets gs
               LEFT JOIN gene_set_memberships gsm ON gs.gs_id = gsm.gs_id
               GROUP BY gs.gs_id, gs.gs_name, gs.gs_collection
               ORDER BY gene_count DESC
               LIMIT 5"""
        ).fetchall()

        print("\nTop 5 gene sets by gene count:")
        for gs in sample_sets:
            print(f"   [{gs['gs_id']}] {gs['gs_name'][:50]}")
            print(f"        Collection: {gs['gs_collection']}, Genes: {gs['gene_count']}")

        return sample_sets

    def test_specific_pathway(self, pathway_id: int = None):
        """Test with a specific pathway ID."""
        print(f"\n🧪 TESTING SPECIFIC PATHWAY")
        print("-" * 80)

        # If no pathway_id provided, use the first one with good gene count
        if pathway_id is None:
            result = self.conn.execute(
                """SELECT gs.gs_id 
                   FROM gene_sets gs
                   JOIN gene_set_memberships gsm ON gs.gs_id = gsm.gs_id
                   GROUP BY gs.gs_id
                   HAVING COUNT(DISTINCT gsm.gene_symbol) > 10
                   LIMIT 1"""
            ).fetchone()

            if result:
                pathway_id = result['gs_id']
            else:
                print("❌ No pathways found with sufficient genes")
                return

        print(f"Testing pathway ID: {pathway_id}")

        # Get pathway info
        pathway = self.conn.execute(
            "SELECT * FROM gene_sets WHERE gs_id = ?", (pathway_id,)
        ).fetchone()

        if not pathway:
            print(f"❌ Pathway {pathway_id} not found!")
            return

        print(f"Pathway: {pathway['gs_name']}")
        print(f"Collection: {pathway['gs_collection']}")

        # Get genes in this pathway
        genes = self.conn.execute(
            """SELECT DISTINCT gene_symbol 
               FROM gene_set_memberships 
               WHERE gs_id = ?
               LIMIT 10""",
            (pathway_id,)
        ).fetchall()

        print(f"\nFirst 10 genes in pathway:")
        gene_list = [g['gene_symbol'] for g in genes]
        for gene in gene_list:
            print(f"   - {gene}")

        # Now test the full analysis query
        print(f"\n📊 Testing pathway analysis query...")

        # First, check if we can get the basic phosphosite data
        basic_query = """
        SELECT 
            p.site_id,
            p.gene_symbol,
            p.position,
            p.residue_type,
            p.predicted_prob_calibrated,
            p.qvalue,
            p.significant_fdr05,
            p.secondary_structure,
            p.plddt,
            p.sasa_ratio,
            p.disorder_score,
            p.neighbor_count,
            p.motif
        FROM gene_set_memberships gsm
        JOIN phosphosites p ON gsm.gene_symbol = p.gene_symbol
        WHERE gsm.gs_id = ?
        LIMIT 5
        """

        basic_results = self.conn.execute(basic_query, (pathway_id,)).fetchall()

        if not basic_results:
            print("❌ No phosphosites found for genes in this pathway!")

            # Debug: Check if genes exist in phosphosites table
            for gene in gene_list[:3]:
                count = self.conn.execute(
                    "SELECT COUNT(*) as cnt FROM phosphosites WHERE gene_symbol = ?",
                    (gene,)
                ).fetchone()['cnt']
                print(f"   Sites for {gene}: {count}")
        else:
            print(f"✅ Found phosphosites for pathway genes")
            print(f"   Sample sites:")
            for site in basic_results[:3]:
                print(f"      - {site['gene_symbol']} {site['residue_type']}{site['position']} "
                      f"(conf: {site['predicted_prob_calibrated']:.3f if site['predicted_prob_calibrated'] else 0:.3f})")

        # Now test with kinase columns
        print(f"\n🔍 Testing kinase score retrieval...")

        # Get column names dynamically
        cursor = self.conn.execute("SELECT * FROM phosphosites LIMIT 0")
        all_columns = [description[0] for description in cursor.description]
        kinase_columns = [col for col in all_columns if '_MotifScore' in col]

        if not kinase_columns:
            print("❌ No kinase score columns found in phosphosites table!")
            return

        print(f"✅ Found {len(kinase_columns)} kinase score columns")

        # Build a query with first 3 kinase columns as a test
        test_kinases = kinase_columns[:3] if kinase_columns else []

        if test_kinases:
            kinase_cols_str = ', '.join([f'p.{col}' for col in test_kinases])
            test_query = f"""
            SELECT 
                p.site_id,
                p.gene_symbol,
                {kinase_cols_str}
            FROM gene_set_memberships gsm
            JOIN phosphosites p ON gsm.gene_symbol = p.gene_symbol
            WHERE gsm.gs_id = ?
            AND p.significant_fdr05 = 1
            LIMIT 3
            """

            try:
                kinase_results = self.conn.execute(test_query, (pathway_id,)).fetchall()

                if kinase_results:
                    print(f"✅ Successfully retrieved kinase scores")
                    for site in kinase_results:
                        print(f"\n   {site['site_id']} ({site['gene_symbol']}):")
                        for kinase_col in test_kinases:
                            score = site[kinase_col]
                            if score is not None and score > 0:
                                print(f"      {kinase_col}: {score:.4f}")
                else:
                    print("⚠️  Query executed but no results with significant_fdr05")

                    # Try without the significant filter
                    test_query_no_filter = f"""
                    SELECT 
                        p.site_id,
                        p.gene_symbol,
                        p.significant_fdr05,
                        {kinase_cols_str}
                    FROM gene_set_memberships gsm
                    JOIN phosphosites p ON gsm.gene_symbol = p.gene_symbol
                    WHERE gsm.gs_id = ?
                    LIMIT 3
                    """

                    kinase_results = self.conn.execute(test_query_no_filter, (pathway_id,)).fetchall()
                    if kinase_results:
                        print(f"\n✅ Found sites without significant filter:")
                        for site in kinase_results:
                            print(f"   {site['site_id']} (significant: {site['significant_fdr05']})")

            except Exception as e:
                print(f"❌ Error retrieving kinase scores: {e}")

        return pathway_id

    def test_api_endpoint_simulation(self, pathway_id: int):
        """Simulate what the API endpoint should return."""
        print(f"\n🌐 API ENDPOINT SIMULATION")
        print("-" * 80)

        # This simulates the exact query from the routes.py file
        cursor = self.conn.execute("SELECT * FROM phosphosites LIMIT 0")
        all_columns = [description[0] for description in cursor.description]
        kinase_columns = [col for col in all_columns if '_MotifScore' in col]

        print(f"Simulating API call for pathway {pathway_id}")
        print(f"Expected to retrieve {len(kinase_columns)} kinase scores per site")

        # Count total sites for this pathway
        count_query = """
        SELECT COUNT(DISTINCT p.site_id) as total
        FROM gene_set_memberships gsm
        JOIN phosphosites p ON gsm.gene_symbol = p.gene_symbol
        WHERE gsm.gs_id = ?
        """

        total_sites = self.conn.execute(count_query, (pathway_id,)).fetchone()['total']
        print(f"Total phosphosites in pathway: {total_sites:,}")

        # Get phosphocompetent sites count (matching the frontend filter logic)
        competent_query = """
        SELECT COUNT(DISTINCT p.site_id) as total
        FROM gene_set_memberships gsm
        JOIN phosphosites p ON gsm.gene_symbol = p.gene_symbol
        WHERE gsm.gs_id = ?
        AND p.predicted_prob_calibrated >= 0.5
        AND p.qvalue <= 0.05
        """

        competent_sites = self.conn.execute(competent_query, (pathway_id,)).fetchone()['total']
        print(f"Phosphocompetent sites (prob>=0.5, q<=0.05): {competent_sites:,}")

        # Memory estimate
        memory_estimate = (total_sites * len(kinase_columns) * 8) / (1024 * 1024)  # MB
        print(f"Estimated memory for kinase scores: {memory_estimate:.2f} MB")

        if memory_estimate > 100:
            print("⚠️  Large dataset - may cause performance issues")

        # Now test the actual API query structure
        print("\n📋 Testing full API query (limited to 2 sites for debugging)...")

        # Build the full kinase columns string
        kinase_cols_str = ', '.join([f'p.{col}' for col in kinase_columns])

        full_query = f"""
        SELECT 
            p.site_id,
            p.gene_symbol,
            p.position,
            p.residue_type,
            p.motif,
            p.predicted_prob_calibrated,
            p.qvalue,
            p.significant_fdr05,
            p.significant_fdr10,
            p.significant_fdr20,
            p.plddt,
            p.secondary_structure,
            p.sasa_ratio,
            p.disorder_score,
            p.neighbor_count,
            p.hydrogen_bonds,
            p.hydroxyl_exposure,
            {kinase_cols_str}
        FROM gene_set_memberships gsm
        JOIN phosphosites p ON gsm.gene_symbol = p.gene_symbol
        WHERE gsm.gs_id = ?
        ORDER BY p.gene_symbol, p.position
        LIMIT 2
        """

        try:
            results = self.conn.execute(full_query, (pathway_id,)).fetchall()

            if results:
                print(f"✅ Full query successful! Retrieved {len(results)} sites")

                # Analyze the first site
                site = results[0]
                print(f"\nAnalyzing site: {site['site_id']}")

                # Count non-null kinase scores
                non_null_kinases = 0
                top_kinases = []

                for kinase_col in kinase_columns:
                    score = site[kinase_col]
                    if score is not None and score > 0:
                        non_null_kinases += 1
                        kinase_name = kinase_col.replace('_MotifScore', '')
                        top_kinases.append((kinase_name, score))

                top_kinases.sort(key=lambda x: x[1], reverse=True)

                print(f"   Non-null kinase scores: {non_null_kinases}/{len(kinase_columns)}")
                print(f"   Top 3 kinases:")
                for kinase, score in top_kinases[:3]:
                    print(f"      - {kinase}: {score:.4f}")

            else:
                print("❌ Full query returned no results")

        except Exception as e:
            print(f"❌ Error running full API query: {e}")
            print(f"   This is likely the cause of the pathway analysis failure")

    def check_data_consistency(self):
        """Check for data consistency issues."""
        print(f"\n🔧 DATA CONSISTENCY CHECKS")
        print("-" * 80)

        # Check for orphaned phosphosites
        orphan_query = """
        SELECT COUNT(*) as cnt
        FROM phosphosites p
        WHERE NOT EXISTS (
            SELECT 1 FROM proteins pr 
            WHERE pr.gene_symbol = p.gene_symbol
        )
        """
        orphans = self.conn.execute(orphan_query).fetchone()['cnt']

        if orphans > 0:
            print(f"⚠️  Found {orphans:,} phosphosites without matching protein records")
        else:
            print("✅ All phosphosites have matching protein records")

        # Check for genes without phosphosites
        no_sites_query = """
        SELECT COUNT(*) as cnt
        FROM proteins pr
        WHERE NOT EXISTS (
            SELECT 1 FROM phosphosites p 
            WHERE p.gene_symbol = pr.gene_symbol
        )
        """
        no_sites = self.conn.execute(no_sites_query).fetchone()['cnt']

        if no_sites > 0:
            print(f"ℹ️  {no_sites:,} proteins have no phosphosites")

        # Check site_id format
        sample_sites = self.conn.execute(
            "SELECT DISTINCT site_id FROM phosphosites LIMIT 10"
        ).fetchall()

        print("\nSample site_id formats:")
        for site in sample_sites:
            print(f"   - {site['site_id']}")

        # Check for NULL values in critical columns
        null_checks = [
            ('gene_symbol', 'gene_symbol'),
            ('position', 'position'),
            ('residue_type', 'residue_type'),
            ('predicted_prob_calibrated', 'confidence score')
        ]

        print("\nNULL value checks:")
        for col, desc in null_checks:
            null_count = self.conn.execute(
                f"SELECT COUNT(*) as cnt FROM phosphosites WHERE {col} IS NULL"
            ).fetchone()['cnt']

            if null_count > 0:
                print(f"   ⚠️  {null_count:,} sites have NULL {desc}")
            else:
                print(f"   ✅ No NULL values in {desc}")

    def test_pathway_endpoint_direct(self, pathway_id: int):
        """Test the exact endpoint that's failing."""
        print(f"\n🔍 TESTING EXACT API ENDPOINT FOR PATHWAY {pathway_id}")
        print("-" * 80)

        # Get pathway info first
        pathway = self.conn.execute(
            "SELECT * FROM gene_sets WHERE gs_id = ?",
            (pathway_id,)
        ).fetchone()

        if pathway:
            print(f"Testing: {pathway['gs_name']}")
            print(f"Collection: {pathway['gs_collection']}")

        # Test if we're getting the pathway statistics correctly
        stats_query = """
        SELECT 
            COUNT(DISTINCT gsm.gene_symbol) as total_proteins,
            COUNT(DISTINCT p.site_id) as total_sites,
            COUNT(DISTINCT CASE WHEN p.significant_fdr05 = 1 THEN p.site_id END) as significant_sites
        FROM gene_set_memberships gsm
        LEFT JOIN phosphosites p ON gsm.gene_symbol = p.gene_symbol
        WHERE gsm.gs_id = ?
        """

        stats = self.conn.execute(stats_query, (pathway_id,)).fetchone()
        print(f"\nPathway statistics:")
        print(f"   Total proteins: {stats['total_proteins']}")
        print(f"   Total sites: {stats['total_sites']}")
        print(f"   Significant sites: {stats['significant_sites']}")

        if stats['total_sites'] == 0:
            print("\n❌ PROBLEM: No phosphosites found for this pathway!")
            print("   This would cause the analysis page to show nothing")

    def generate_fix_suggestions(self):
        """Generate suggestions for fixing issues."""
        print(f"\n💡 RECOMMENDATIONS")
        print("-" * 80)

        # Check if kinase scores exist
        cursor = self.conn.execute("SELECT * FROM phosphosites LIMIT 0")
        all_columns = [description[0] for description in cursor.description]
        kinase_columns = [col for col in all_columns if '_MotifScore' in col]

        if len(kinase_columns) == 303:
            print(f"✅ All {len(kinase_columns)} kinase score columns present")
            print("\nThe database structure looks correct!")
            print("\nPotential issues to check:")
            print("1. The pathway analysis page might be timing out on large pathways")
            print("2. There might be a JavaScript error in the frontend")
            print("3. The API endpoint might be returning too much data at once")
            print("\nRecommended debugging steps:")
            print("   1. Open browser Developer Tools (F12)")
            print("   2. Go to Network tab")
            print("   3. Navigate to a pathway analysis page")
            print("   4. Look for the API call to /api/pathway/[id]/analysis")
            print("   5. Check if it completes or times out")
            print("   6. Check Console tab for JavaScript errors")

    def close(self):
        """Close database connection."""
        self.conn.close()
        print("\n✅ Database connection closed")


def main():
    """Run the complete debugging suite."""
    print("=" * 80)
    print("🔬 KINOPLEX PATHWAY ANALYSIS DEBUGGER")
    print("=" * 80)

    # Find the database
    db_paths = [
        'kinoplex.db',
        '/Users/davidvanderwall/Desktop/KinoPlex_App/kinoplex.db',
        'app/kinoplex.db'
    ]

    db_path = None
    for path in db_paths:
        if Path(path).exists():
            db_path = path
            break

    if not db_path:
        print("❌ Could not find kinoplex.db")
        print("Please specify the path to your database file")
        sys.exit(1)

    # Run debugging
    debugger = PathwayAnalysisDebugger(db_path)

    try:
        # Run all checks
        kinase_cols = debugger.debug_database_structure()
        debugger.check_sample_data()
        debugger.check_pathway_data()
        pathway_id = debugger.test_specific_pathway()

        if pathway_id:
            debugger.test_api_endpoint_simulation(pathway_id)
            debugger.test_pathway_endpoint_direct(pathway_id)

        debugger.check_data_consistency()
        debugger.generate_fix_suggestions()

    finally:
        debugger.close()

    print("\n" + "=" * 80)
    print("Debugging complete!")
    print("=" * 80)


if __name__ == "__main__":
    main()