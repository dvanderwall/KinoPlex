"""Fixed API endpoints for pathway and protein search."""

from flask import request, jsonify, current_app, g
from app.api import bp
import sqlite3
import json
import time


def get_db():
    """Get database connection."""
    if 'db' not in g:
        g.db = sqlite3.connect(current_app.config['DATABASE_PATH'])
        g.db.row_factory = sqlite3.Row
    return g.db


def execute_query(query, params=(), fetch_one=False):
    """Execute database query with error handling."""
    try:
        db = get_db()
        cursor = db.execute(query, params)
        if fetch_one:
            return cursor.fetchone()
        return cursor.fetchall()
    except Exception as e:
        current_app.logger.error(f"Database query error: {e}")
        return None if fetch_one else []


@bp.route('/stats/overview')
def get_overview_stats():
    """Get comprehensive database statistics."""
    try:
        # Core counts
        total_sites_result = execute_query("SELECT COUNT(*) as count FROM phosphosites", fetch_one=True)
        total_sites = total_sites_result['count'] if total_sites_result else 0

        # Check if tables exist before querying
        tables_query = "SELECT name FROM sqlite_master WHERE type='table'"
        existing_tables = [row['name'] for row in execute_query(tables_query)]

        total_proteins = 0
        if 'proteins' in existing_tables:
            result = execute_query("SELECT COUNT(*) as count FROM proteins", fetch_one=True)
            total_proteins = result['count'] if result else 0

        total_pathways = 0
        if 'gene_sets' in existing_tables:
            result = execute_query("SELECT COUNT(*) as count FROM gene_sets", fetch_one=True)
            total_pathways = result['count'] if result else 0

        # Try to get significant sites
        significant_sites = 0
        try:
            result = execute_query("SELECT COUNT(*) as count FROM phosphosites WHERE significant_fdr05 = 1",
                                   fetch_one=True)
            significant_sites = result['count'] if result else 0
        except:
            # Column might not exist
            pass

        # Quality distribution (simplified)
        quality_dist = {'very_high': 0, 'high': 0, 'moderate': 0, 'low': 0}
        try:
            quality_query = """
            SELECT 
                COUNT(CASE WHEN predicted_prob_calibrated > 0.9 THEN 1 END) as very_high,
                COUNT(CASE WHEN predicted_prob_calibrated BETWEEN 0.7 AND 0.9 THEN 1 END) as high,
                COUNT(CASE WHEN predicted_prob_calibrated BETWEEN 0.5 AND 0.7 THEN 1 END) as moderate,
                COUNT(CASE WHEN predicted_prob_calibrated < 0.5 THEN 1 END) as low
            FROM phosphosites
            """
            result = execute_query(quality_query, fetch_one=True)
            if result:
                quality_dist = dict(result)
        except:
            pass

        # Pathway collections (simplified)
        pathway_collections = []
        if 'gene_sets' in existing_tables:
            try:
                collections_query = """
                SELECT gs_collection, COUNT(*) as count
                FROM gene_sets 
                GROUP BY gs_collection 
                ORDER BY count DESC
                LIMIT 10
                """
                results = execute_query(collections_query)
                pathway_collections = [dict(row) for row in results] if results else []
            except:
                pass

        # Residue distribution
        residue_distribution = []
        try:
            residue_query = """
            SELECT residue_type, COUNT(*) as count
            FROM phosphosites 
            GROUP BY residue_type
            ORDER BY count DESC
            """
            results = execute_query(residue_query)
            residue_distribution = [dict(row) for row in results] if results else []
        except:
            pass

        return jsonify({
            'total_sites': total_sites,
            'total_proteins': total_proteins,
            'total_pathways': total_pathways,
            'significant_sites': significant_sites,
            'quality_distribution': quality_dist,
            'pathway_collections': pathway_collections,
            'residue_distribution': residue_distribution,
            'generated_at': time.time()
        })

    except Exception as e:
        current_app.logger.error(f"Error getting overview stats: {e}")
        return jsonify({'error': 'Failed to retrieve statistics'}), 500


@bp.route('/search/pathways')
def search_pathways():
    """Advanced pathway search with filters - FIXED VERSION."""
    query = request.args.get('q', '').strip()
    collection = request.args.get('collection', '')
    min_genes = int(request.args.get('min_genes', 0))
    limit = min(int(request.args.get('limit', 50)), 100)

    if not query or len(query) < 2:
        return jsonify({'pathways': []})

    try:
        # Check if required tables exist
        tables_query = "SELECT name FROM sqlite_master WHERE type='table'"
        existing_tables = [row['name'] for row in execute_query(tables_query)]

        if 'gene_sets' not in existing_tables:
            return jsonify({'pathways': [], 'message': 'Gene sets table not found'})

        # Build conditions for the WHERE clause
        conditions = ["(gs.gs_name LIKE ? OR gs.gs_description LIKE ?)"]
        params = [f'%{query}%', f'%{query}%']

        if collection:
            conditions.append("gs.gs_collection = ?")
            params.append(collection)

        # FIXED: Properly join tables to get actual gene counts
        sql = f"""
        SELECT 
            gs.gs_id, 
            gs.gs_name, 
            gs.gs_collection, 
            gs.gs_description, 
            gs.gs_pmid,
            COUNT(DISTINCT gsm.gene_symbol) as gene_count,
            COUNT(DISTINCT CASE WHEN p.significant_fdr05 = 1 THEN p.site_id END) as significant_sites
        FROM gene_sets gs
        LEFT JOIN gene_set_memberships gsm ON gs.gs_id = gsm.gs_id
        LEFT JOIN phosphosites p ON gsm.gene_symbol = p.gene_symbol
        WHERE {' AND '.join(conditions)}
        GROUP BY gs.gs_id, gs.gs_name, gs.gs_collection, gs.gs_description, gs.gs_pmid
        HAVING COUNT(DISTINCT gsm.gene_symbol) >= ?
        ORDER BY significant_sites DESC, gene_count DESC
        LIMIT ?
        """

        params.extend([min_genes, limit])

        results = execute_query(sql, tuple(params))
        pathways = [dict(row) for row in results] if results else []

        return jsonify({
            'pathways': pathways,
            'query': query,
            'filters': {
                'collection': collection,
                'min_genes': min_genes
            }
        })

    except Exception as e:
        current_app.logger.error(f"Error searching pathways: {e}")
        return jsonify({'error': f'Search failed: {str(e)}'}), 500


@bp.route('/search/proteins')
def search_proteins():
    """Advanced protein search with site statistics - FIXED VERSION."""
    query = request.args.get('q', '').strip().upper()
    limit = min(int(request.args.get('limit', 50)), 100)

    if not query or len(query) < 1:
        return jsonify({'proteins': []})

    try:
        # FIXED: Better protein search with proper statistics
        sql = """
        SELECT 
            gene_symbol,
            COUNT(site_id) as site_count,
            COUNT(CASE WHEN significant_fdr05 = 1 THEN 1 END) as significant_sites,
            AVG(CASE WHEN predicted_prob_calibrated IS NOT NULL THEN predicted_prob_calibrated END) as avg_confidence,
            MAX(CASE WHEN predicted_prob_calibrated IS NOT NULL THEN predicted_prob_calibrated END) as max_confidence
        FROM phosphosites 
        WHERE gene_symbol LIKE ?
        GROUP BY gene_symbol
        ORDER BY significant_sites DESC, site_count DESC
        LIMIT ?
        """

        results = execute_query(sql, (f'{query}%', limit))
        proteins = []

        for row in results:
            protein_dict = dict(row)

            # Get pathway count for this protein
            pathway_count_sql = """
            SELECT COUNT(DISTINCT gs_id) as pathway_count
            FROM gene_set_memberships
            WHERE gene_symbol = ?
            """
            pathway_result = execute_query(pathway_count_sql, (protein_dict['gene_symbol'],), fetch_one=True)
            protein_dict['pathway_count'] = pathway_result['pathway_count'] if pathway_result else 0

            proteins.append(protein_dict)

        return jsonify({'proteins': proteins})

    except Exception as e:
        current_app.logger.error(f"Error searching proteins: {e}")
        return jsonify({'error': f'Search failed: {str(e)}'}), 500


@bp.route('/pathway/<int:pathway_id>')
def get_pathway_details(pathway_id):
    """Get comprehensive pathway analysis - FIXED VERSION."""
    try:
        # Get pathway basic info
        pathway_query = """
        SELECT gs_id, gs_name, gs_collection, gs_description, gs_pmid
        FROM gene_sets 
        WHERE gs_id = ?
        """
        pathway = execute_query(pathway_query, (pathway_id,), fetch_one=True)

        if not pathway:
            return jsonify({'error': 'Pathway not found'}), 404

        # Get pathway statistics
        stats_query = """
        SELECT 
            COUNT(DISTINCT gsm.gene_symbol) as total_genes,
            COUNT(DISTINCT p.site_id) as total_sites,
            COUNT(DISTINCT CASE WHEN p.significant_fdr05 = 1 THEN p.site_id END) as significant_sites
        FROM gene_set_memberships gsm
        LEFT JOIN phosphosites p ON gsm.gene_symbol = p.gene_symbol
        WHERE gsm.gs_id = ?
        """
        stats_result = execute_query(stats_query, (pathway_id,), fetch_one=True)
        stats = dict(stats_result) if stats_result else {'total_genes': 0, 'total_sites': 0, 'significant_sites': 0}

        # Get associated proteins with their stats
        proteins_query = """
        SELECT 
            gsm.gene_symbol,
            COUNT(p.site_id) as total_sites,
            COUNT(CASE WHEN p.significant_fdr05 = 1 THEN 1 END) as significant_sites,
            AVG(CASE WHEN p.predicted_prob_calibrated IS NOT NULL THEN p.predicted_prob_calibrated END) as avg_confidence
        FROM gene_set_memberships gsm
        LEFT JOIN phosphosites p ON gsm.gene_symbol = p.gene_symbol
        WHERE gsm.gs_id = ?
        GROUP BY gsm.gene_symbol
        ORDER BY significant_sites DESC, total_sites DESC
        LIMIT 50
        """
        proteins_results = execute_query(proteins_query, (pathway_id,))
        proteins = [dict(row) for row in proteins_results] if proteins_results else []

        return jsonify({
            'pathway': dict(pathway),
            'proteins': proteins,
            'stats': stats
        })

    except Exception as e:
        current_app.logger.error(f"Error getting pathway details: {e}")
        return jsonify({'error': 'Failed to retrieve pathway details'}), 500


@bp.route('/protein/<gene_symbol>')
def get_protein_details(gene_symbol):
    """Get comprehensive protein analysis - FIXED VERSION."""
    gene_symbol = gene_symbol.upper()

    try:
        # Get all phosphosites for this protein
        sites_query = """
        SELECT 
            site_id, position, residue_type, motif,
            predicted_prob_calibrated, qvalue,
            significant_fdr05, significant_fdr10, significant_fdr20,
            plddt, secondary_structure,
            sasa_ratio, disorder_score, neighbor_count
        FROM phosphosites 
        WHERE gene_symbol = ?
        ORDER BY position
        """

        sites_results = execute_query(sites_query, (gene_symbol,))
        if not sites_results:
            return jsonify({'error': 'Protein not found'}), 404

        sites = [dict(row) for row in sites_results]

        # Get protein metadata if available
        protein_query = """
        SELECT gene_symbol, ncbi_gene_id, ensembl_gene_id
        FROM proteins
        WHERE gene_symbol = ?
        """
        protein_result = execute_query(protein_query, (gene_symbol,), fetch_one=True)

        protein = dict(protein_result) if protein_result else {
            'gene_symbol': gene_symbol,
            'ncbi_gene_id': None,
            'ensembl_gene_id': None
        }

        # Get pathway associations
        pathways_query = """
        SELECT gs.gs_id, gs.gs_name, gs.gs_collection, gs.gs_description
        FROM gene_set_memberships gsm
        JOIN gene_sets gs ON gsm.gs_id = gs.gs_id
        WHERE gsm.gene_symbol = ?
        ORDER BY gs.gs_collection, gs.gs_name
        LIMIT 20
        """
        pathways_results = execute_query(pathways_query, (gene_symbol,))
        pathways = [dict(row) for row in pathways_results] if pathways_results else []

        # Calculate statistics
        stats = {
            'total_sites': len(sites),
            'significant_sites': sum(1 for s in sites if s.get('significant_fdr05')),
            'high_confidence': sum(1 for s in sites if s.get('predicted_prob_calibrated', 0) > 0.8),
            'avg_confidence': sum(s.get('predicted_prob_calibrated', 0) for s in sites) / len(sites) if sites else 0,
            'pathway_count': len(pathways)
        }

        return jsonify({
            'protein': protein,
            'sites': sites,
            'pathways': pathways,
            'stats': stats
        })

    except Exception as e:
        current_app.logger.error(f"Error getting protein details: {e}")
        return jsonify({'error': 'Failed to retrieve protein details'}), 500


# app/api/routes.py - Updated pathway analysis endpoint
# app/api/routes.py - Complete pathway analysis endpoint with all kinase scores

@bp.route('/pathway/<int:pathway_id>/analysis')
def get_pathway_analysis(pathway_id):
    """Comprehensive pathway-level phosphoproteomics analysis with full kinase extraction."""
    try:
        # Get pathway basic info
        pathway_query = """
        SELECT gs_id, gs_name, gs_collection, gs_description, gs_pmid
        FROM gene_sets 
        WHERE gs_id = ?
        """
        pathway = execute_query(pathway_query, (pathway_id,), fetch_one=True)

        if not pathway:
            return jsonify({'error': 'Pathway not found'}), 404

        # Get list of all kinase columns from the database schema
        # This is the complete list of 303 kinases
        kinase_columns = [
            'AAK1_MotifScore', 'ACVR2A_MotifScore', 'ACVR2B_MotifScore', 'AKT1_MotifScore', 'AKT2_MotifScore',
            'AKT3_MotifScore', 'ALK2_MotifScore', 'ALK4_MotifScore', 'ALPHAK3_MotifScore', 'AMPKA1_MotifScore',
            'AMPKA2_MotifScore', 'ANKRD3_MotifScore', 'ASK1_MotifScore', 'ATM_MotifScore', 'ATR_MotifScore',
            'AURA_MotifScore', 'AURB_MotifScore', 'AURC_MotifScore', 'BCKDK_MotifScore', 'BIKE_MotifScore',
            'BMPR1A_MotifScore', 'BMPR1B_MotifScore', 'BMPR2_MotifScore', 'BRAF_MotifScore', 'BRSK1_MotifScore',
            'BRSK2_MotifScore', 'BUB1_MotifScore', 'CAMK1A_MotifScore', 'CAMK1B_MotifScore', 'CAMK1D_MotifScore',
            'CAMK1G_MotifScore', 'CAMK2A_MotifScore', 'CAMK2B_MotifScore', 'CAMK2D_MotifScore', 'CAMK2G_MotifScore',
            'CAMK4_MotifScore', 'CAMKK1_MotifScore', 'CAMKK2_MotifScore', 'CAMLCK_MotifScore', 'CDC7_MotifScore',
            'CDK1_MotifScore', 'CDK10_MotifScore', 'CDK12_MotifScore', 'CDK13_MotifScore', 'CDK14_MotifScore',
            'CDK16_MotifScore', 'CDK17_MotifScore', 'CDK18_MotifScore', 'CDK19_MotifScore', 'CDK2_MotifScore',
            'CDK3_MotifScore', 'CDK4_MotifScore', 'CDK5_MotifScore', 'CDK6_MotifScore', 'CDK7_MotifScore',
            'CDK8_MotifScore', 'CDK9_MotifScore', 'CDKL1_MotifScore', 'CDKL5_MotifScore', 'CHAK1_MotifScore',
            'CHAK2_MotifScore', 'CHK1_MotifScore', 'CHK2_MotifScore', 'CK1A_MotifScore', 'CK1A2_MotifScore',
            'CK1D_MotifScore', 'CK1E_MotifScore', 'CK1G1_MotifScore', 'CK1G2_MotifScore', 'CK1G3_MotifScore',
            'CK2A1_MotifScore', 'CK2A2_MotifScore', 'CLK1_MotifScore', 'CLK2_MotifScore', 'CLK3_MotifScore',
            'CLK4_MotifScore', 'COT_MotifScore', 'CRIK_MotifScore', 'DAPK1_MotifScore', 'DAPK2_MotifScore',
            'DAPK3_MotifScore', 'DCAMKL1_MotifScore', 'DCAMKL2_MotifScore', 'DLK_MotifScore', 'DMPK1_MotifScore',
            'DNAPK_MotifScore', 'DRAK1_MotifScore', 'DSTYK_MotifScore', 'DYRK1A_MotifScore', 'DYRK1B_MotifScore',
            'DYRK2_MotifScore', 'DYRK3_MotifScore', 'DYRK4_MotifScore', 'EEF2K_MotifScore', 'ERK1_MotifScore',
            'ERK2_MotifScore', 'ERK5_MotifScore', 'ERK7_MotifScore', 'FAM20C_MotifScore', 'GAK_MotifScore',
            'GCK_MotifScore', 'GCN2_MotifScore', 'GRK1_MotifScore', 'GRK2_MotifScore', 'GRK3_MotifScore',
            'GRK4_MotifScore', 'GRK5_MotifScore', 'GRK6_MotifScore', 'GRK7_MotifScore', 'GSK3A_MotifScore',
            'GSK3B_MotifScore', 'HASPIN_MotifScore', 'HGK_MotifScore', 'HIPK1_MotifScore', 'HIPK2_MotifScore',
            'HIPK3_MotifScore', 'HIPK4_MotifScore', 'HPK1_MotifScore', 'HRI_MotifScore', 'HUNK_MotifScore',
            'ICK_MotifScore', 'IKKA_MotifScore', 'IKKB_MotifScore', 'IKKE_MotifScore', 'IRAK1_MotifScore',
            'IRAK4_MotifScore', 'IRE1_MotifScore', 'IRE2_MotifScore', 'JNK1_MotifScore', 'JNK2_MotifScore',
            'JNK3_MotifScore', 'KHS1_MotifScore', 'KHS2_MotifScore', 'KIS_MotifScore', 'LATS1_MotifScore',
            'LATS2_MotifScore', 'LKB1_MotifScore', 'LOK_MotifScore', 'LRRK2_MotifScore', 'MAK_MotifScore',
            'MAP3K15_MotifScore', 'MAPKAPK2_MotifScore', 'MAPKAPK3_MotifScore', 'MAPKAPK5_MotifScore',
            'MARK1_MotifScore', 'MARK2_MotifScore', 'MARK3_MotifScore', 'MARK4_MotifScore', 'MASTL_MotifScore',
            'MEK1_MotifScore', 'MEK2_MotifScore', 'MEK5_MotifScore', 'MEKK1_MotifScore', 'MEKK2_MotifScore',
            'MEKK3_MotifScore', 'MEKK6_MotifScore', 'MELK_MotifScore', 'MINK_MotifScore', 'MLK1_MotifScore',
            'MLK2_MotifScore', 'MLK3_MotifScore', 'MLK4_MotifScore', 'MNK1_MotifScore', 'MNK2_MotifScore',
            'MOK_MotifScore', 'MOS_MotifScore', 'MPSK1_MotifScore', 'MRCKA_MotifScore', 'MRCKB_MotifScore',
            'MSK1_MotifScore', 'MSK2_MotifScore', 'MST1_MotifScore', 'MST2_MotifScore', 'MST3_MotifScore',
            'MST4_MotifScore', 'MTOR_MotifScore', 'MYLK4_MotifScore', 'MYO3A_MotifScore', 'MYO3B_MotifScore',
            'NDR1_MotifScore', 'NDR2_MotifScore', 'NEK1_MotifScore', 'NEK11_MotifScore', 'NEK2_MotifScore',
            'NEK3_MotifScore', 'NEK4_MotifScore', 'NEK5_MotifScore', 'NEK6_MotifScore', 'NEK7_MotifScore',
            'NEK8_MotifScore', 'NEK9_MotifScore', 'NIK_MotifScore', 'NIM1_MotifScore', 'NLK_MotifScore',
            'NUAK1_MotifScore', 'NUAK2_MotifScore', 'OSR1_MotifScore', 'P38A_MotifScore', 'P38B_MotifScore',
            'P38D_MotifScore', 'P38G_MotifScore', 'P70S6K_MotifScore', 'P70S6KB_MotifScore', 'P90RSK_MotifScore',
            'PAK1_MotifScore', 'PAK2_MotifScore', 'PAK3_MotifScore', 'PAK4_MotifScore', 'PAK5_MotifScore',
            'PAK6_MotifScore', 'PASK_MotifScore', 'PBK_MotifScore', 'PDHK1_MotifScore', 'PDHK4_MotifScore',
            'PDK1_MotifScore', 'PERK_MotifScore', 'PHKG1_MotifScore', 'PHKG2_MotifScore', 'PIM1_MotifScore',
            'PIM2_MotifScore', 'PIM3_MotifScore', 'PINK1_MotifScore', 'PKACA_MotifScore', 'PKACB_MotifScore',
            'PKACG_MotifScore', 'PKCA_MotifScore', 'PKCB_MotifScore', 'PKCD_MotifScore', 'PKCE_MotifScore',
            'PKCG_MotifScore', 'PKCH_MotifScore', 'PKCI_MotifScore', 'PKCT_MotifScore', 'PKCZ_MotifScore',
            'PKG1_MotifScore', 'PKG2_MotifScore', 'PKN1_MotifScore', 'PKN2_MotifScore', 'PKN3_MotifScore',
            'PKR_MotifScore', 'PLK1_MotifScore', 'PLK2_MotifScore', 'PLK3_MotifScore', 'PLK4_MotifScore',
            'PRKD1_MotifScore', 'PRKD2_MotifScore', 'PRKD3_MotifScore', 'PRKX_MotifScore', 'PRP4_MotifScore',
            'PRPK_MotifScore', 'QIK_MotifScore', 'QSK_MotifScore', 'RAF1_MotifScore', 'RIPK1_MotifScore',
            'RIPK2_MotifScore', 'RIPK3_MotifScore', 'ROCK1_MotifScore', 'ROCK2_MotifScore', 'RSK2_MotifScore',
            'RSK3_MotifScore', 'RSK4_MotifScore', 'SBK_MotifScore', 'SGK1_MotifScore', 'SGK3_MotifScore',
            'SIK_MotifScore', 'SKMLCK_MotifScore', 'SLK_MotifScore', 'SMG1_MotifScore', 'SMMLCK_MotifScore',
            'SNRK_MotifScore', 'SRPK1_MotifScore', 'SRPK2_MotifScore', 'SRPK3_MotifScore', 'SSTK_MotifScore',
            'STK33_MotifScore', 'STLK3_MotifScore', 'TAK1_MotifScore', 'TAO1_MotifScore', 'TAO2_MotifScore',
            'TAO3_MotifScore', 'TBK1_MotifScore', 'TGFBR1_MotifScore', 'TGFBR2_MotifScore', 'TLK1_MotifScore',
            'TLK2_MotifScore', 'TNIK_MotifScore', 'TSSK1_MotifScore', 'TSSK2_MotifScore', 'TTBK1_MotifScore',
            'TTBK2_MotifScore', 'TTK_MotifScore', 'ULK1_MotifScore', 'ULK2_MotifScore', 'VRK1_MotifScore',
            'VRK2_MotifScore', 'WNK1_MotifScore', 'WNK3_MotifScore', 'WNK4_MotifScore', 'YANK2_MotifScore',
            'YANK3_MotifScore', 'YSK1_MotifScore', 'YSK4_MotifScore', 'ZAK_MotifScore'
        ]

        # Build SQL query to get ALL data including all kinase scores
        kinase_columns_str = ', '.join([f'p.{col}' for col in kinase_columns])

        proteins_sites_query = f"""
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
            {kinase_columns_str}
        FROM gene_set_memberships gsm
        JOIN phosphosites p ON gsm.gene_symbol = p.gene_symbol
        WHERE gsm.gs_id = ?
        ORDER BY p.gene_symbol, p.position
        """

        sites_results = execute_query(proteins_sites_query, (pathway_id,))
        if not sites_results:
            return jsonify({'error': 'No phosphorylation sites found for this pathway'}), 404

        # Convert to list of dicts for easier processing
        all_sites = [dict(row) for row in sites_results]

        # Process each site to extract top 3 kinases
        for site in all_sites:
            kinase_scores = []
            for kinase_col in kinase_columns:
                score = site.get(kinase_col)
                if score is not None and score > 0:
                    kinase_name = kinase_col.replace('_MotifScore', '')
                    kinase_scores.append({'kinase': kinase_name, 'score': float(score)})

            # Sort by score and get top 3 kinases
            kinase_scores.sort(key=lambda x: x['score'], reverse=True)
            site['top_kinases'] = kinase_scores[:3]

            # Remove individual kinase columns from site data to reduce payload
            for kinase_col in kinase_columns:
                site.pop(kinase_col, None)

        # Organize by protein
        proteins_data = {}
        for site in all_sites:
            gene = site['gene_symbol']
            if gene not in proteins_data:
                proteins_data[gene] = {
                    'gene_symbol': gene,
                    'sites': [],
                    'stats': {
                        'total_sites': 0,
                        'significant_sites': 0,
                        'phosphocompetent_sites': 0,
                        'high_confidence_sites': 0,
                        'avg_confidence': 0,
                        'max_confidence': 0,
                        'min_position': float('inf'),
                        'max_position': 0
                    }
                }
            proteins_data[gene]['sites'].append(site)

        # Calculate protein-level statistics
        for gene, data in proteins_data.items():
            sites = data['sites']
            data['stats']['total_sites'] = len(sites)
            data['stats']['significant_sites'] = sum(1 for s in sites if s['significant_fdr05'])

            # Count phosphocompetent sites (default threshold: prob > 0.5, qvalue < 0.05)
            data['stats']['phosphocompetent_sites'] = sum(
                1 for s in sites
                if (s['predicted_prob_calibrated'] or 0) > 0.5 and (s['qvalue'] or 1) < 0.05
            )

            data['stats']['high_confidence_sites'] = sum(
                1 for s in sites if (s['predicted_prob_calibrated'] or 0) > 0.8
            )

            confidences = [s['predicted_prob_calibrated'] or 0 for s in sites]
            data['stats']['avg_confidence'] = sum(confidences) / len(confidences) if confidences else 0
            data['stats']['max_confidence'] = max(confidences) if confidences else 0

            positions = [s['position'] for s in sites if s['position']]
            if positions:
                data['stats']['min_position'] = min(positions)
                data['stats']['max_position'] = max(positions)

        # Calculate pathway-wide kinase enrichment
        kinase_enrichment = {}
        for site in all_sites:
            for kinase_data in site.get('top_kinases', []):
                kinase = kinase_data['kinase']
                if kinase not in kinase_enrichment:
                    kinase_enrichment[kinase] = {
                        'kinase': kinase,
                        'site_count': 0,
                        'protein_count': 0,
                        'phosphocompetent_sites': 0,
                        'avg_score': 0,
                        'max_score': 0,
                        'scores': [],
                        'proteins': set()
                    }

                kinase_enrichment[kinase]['site_count'] += 1
                kinase_enrichment[kinase]['scores'].append(kinase_data['score'])
                kinase_enrichment[kinase]['proteins'].add(site['gene_symbol'])

                # Count if this is a phosphocompetent site
                if (site['predicted_prob_calibrated'] or 0) > 0.5 and (site['qvalue'] or 1) < 0.05:
                    kinase_enrichment[kinase]['phosphocompetent_sites'] += 1

        # Calculate final kinase statistics
        enriched_kinases = []
        for kinase_data in kinase_enrichment.values():
            scores = kinase_data['scores']
            if scores:
                kinase_data['avg_score'] = sum(scores) / len(scores)
                kinase_data['max_score'] = max(scores)
                kinase_data['protein_count'] = len(kinase_data['proteins'])
                del kinase_data['proteins']  # Remove set object for JSON serialization
                del kinase_data['scores']  # Clean up
                enriched_kinases.append(kinase_data)

        # Sort kinases by phosphocompetent site count
        enriched_kinases.sort(key=lambda x: x['phosphocompetent_sites'], reverse=True)

        # Calculate pathway-wide statistics
        all_confidences = [s['predicted_prob_calibrated'] or 0 for s in all_sites]
        all_qvalues = [s['qvalue'] or 1 for s in all_sites if s['qvalue'] is not None]

        pathway_stats = {
            'total_proteins': len(proteins_data),
            'total_sites': len(all_sites),
            'significant_sites': sum(1 for s in all_sites if s['significant_fdr05']),
            'phosphocompetent_sites': sum(
                1 for s in all_sites
                if (s['predicted_prob_calibrated'] or 0) > 0.5 and (s['qvalue'] or 1) < 0.05
            ),
            'high_confidence_sites': sum(1 for s in all_sites if (s['predicted_prob_calibrated'] or 0) > 0.8),
            'avg_confidence': sum(all_confidences) / len(all_confidences) if all_confidences else 0,
            'median_qvalue': sorted(all_qvalues)[len(all_qvalues) // 2] if all_qvalues else 0,
            'serine_sites': sum(1 for s in all_sites if s['residue_type'] == 'S'),
            'threonine_sites': sum(1 for s in all_sites if s['residue_type'] == 'T'),
            'tyrosine_sites': sum(1 for s in all_sites if s['residue_type'] == 'Y')
        }

        # Calculate structural feature statistics
        valid_plddt = [s['plddt'] for s in all_sites if s['plddt'] is not None]
        valid_sasa = [s['sasa_ratio'] for s in all_sites if s['sasa_ratio'] is not None]
        valid_disorder = [s['disorder_score'] for s in all_sites if s['disorder_score'] is not None]

        structural_stats = {
            'avg_plddt': sum(valid_plddt) / len(valid_plddt) if valid_plddt else 0,
            'avg_sasa_ratio': sum(valid_sasa) / len(valid_sasa) if valid_sasa else 0,
            'avg_disorder': sum(valid_disorder) / len(valid_disorder) if valid_disorder else 0,
            'min_plddt': min(valid_plddt) if valid_plddt else 0,
            'max_plddt': max(valid_plddt) if valid_plddt else 0,
            'surface_exposed_sites': sum(1 for s in valid_sasa if s > 0.4),
            'disordered_sites': sum(1 for s in valid_disorder if s > 0.5)
        }

        # Secondary structure distribution
        structure_distribution = {}
        for site in all_sites:
            struct = site['secondary_structure'] or 'Unknown'
            structure_distribution[struct] = structure_distribution.get(struct, 0) + 1

        # Confidence distribution bins for phosphocompetent sites
        confidence_bins = {
            '0.0-0.3': 0, '0.3-0.5': 0, '0.5-0.7': 0,
            '0.7-0.9': 0, '0.9-1.0': 0
        }
        for site in all_sites:
            conf = site['predicted_prob_calibrated'] or 0
            if conf < 0.3:
                confidence_bins['0.0-0.3'] += 1
            elif conf < 0.5:
                confidence_bins['0.3-0.5'] += 1
            elif conf < 0.7:
                confidence_bins['0.5-0.7'] += 1
            elif conf < 0.9:
                confidence_bins['0.7-0.9'] += 1
            else:
                confidence_bins['0.9-1.0'] += 1

        return jsonify({
            'pathway': dict(pathway),
            'proteins': list(proteins_data.values()),
            'all_sites': all_sites,
            'enriched_kinases': enriched_kinases[:100],  # Top 100 kinases
            'pathway_stats': pathway_stats,
            'structural_stats': structural_stats,
            'structure_distribution': structure_distribution,
            'confidence_distribution': confidence_bins,
            'analysis_metadata': {
                'generated_at': time.time(),
                'kinases_analyzed': len(kinase_columns),
                'total_kinases_found': len(enriched_kinases),
                'analysis_type': 'comprehensive_pathway_phosphoproteomics'
            }
        })

    except Exception as e:
        current_app.logger.error(f"Error in pathway analysis: {e}")
        import traceback
        current_app.logger.error(traceback.format_exc())
        return jsonify({'error': f'Failed to perform pathway analysis: {str(e)}'}), 500

@bp.teardown_app_request
def close_db(error):
    """Close database connection."""
    db = g.pop('db', None)
    if db is not None:
        db.close()