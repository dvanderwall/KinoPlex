"""
FINAL CORRECTED KinoPlex API Routes

This fixes the critical pathway_id -> gs_name -> enrichment data chain
and ensures proper data structure for the frontend visualizations.
"""

from flask import request, jsonify, current_app, g
from app.api import bp
import sqlite3
import json
import time
from collections import defaultdict


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
            result = cursor.fetchone()
            return dict(result) if result else None
        return [dict(row) for row in cursor.fetchall()]
    except Exception as e:
        current_app.logger.error(f"Database query error: {e}")
        current_app.logger.error(f"Query: {query}")
        current_app.logger.error(f"Params: {params}")
        return None if fetch_one else []


@bp.route('/pathway/<int:pathway_id>/analysis')
def get_pathway_analysis(pathway_id):
    """FIXED: Pathway analysis with proper pathway_id -> gs_name -> enrichment chain."""
    try:
        current_app.logger.info(f"=== PATHWAY ANALYSIS START: pathway_id={pathway_id} ===")

        # STEP 1: Get pathway info using gs_id (pathway_id)
        pathway = execute_query("""
            SELECT gs_id, gs_name, 
                   COALESCE(gs_collection, 'Unknown') as gs_collection, 
                   COALESCE(gs_description, 'No description available') as gs_description,
                   gs_pmid
            FROM gene_sets 
            WHERE gs_id = ?
        """, (pathway_id,), fetch_one=True)

        if not pathway:
            current_app.logger.error(f"Pathway not found for gs_id: {pathway_id}")
            return jsonify({'error': 'Pathway not found'}), 404

        gs_name = pathway['gs_name']
        current_app.logger.info(f"Found pathway: '{gs_name}' (gs_id: {pathway_id})")

        # STEP 2: Get proteins in this pathway using gs_id
        proteins_in_pathway = execute_query("""
            SELECT DISTINCT gene_symbol
            FROM gene_set_memberships
            WHERE gs_id = ?
            LIMIT 100
        """, (pathway_id,))

        protein_list = [p['gene_symbol'] for p in proteins_in_pathway]
        current_app.logger.info(f"Found {len(protein_list)} proteins in pathway")

        # STEP 3: Get pathway-level kinase enrichment using gs_name
        # This is the critical query that needs to be fixed
        enriched_kinases = execute_query("""
            SELECT 
                kinase,
                score_integrated,
                rank_integrated,
                n_genes_with_kinase,
                gene_coverage,
                fisher_p,
                enrichment_ratio,
                score_zscore_only,
                score_ratio_based,
                mean_enrichment,
                max_enrichment,
                n_total_sites,
                prop_significant
            FROM pathway_kinase_mapping 
            WHERE gs_name = ? AND score_integrated IS NOT NULL
            ORDER BY score_integrated DESC
            LIMIT 50
        """, (gs_name,))

        current_app.logger.info(f"Found {len(enriched_kinases)} enriched kinases from pathway_kinase_mapping")

        # Add debug output for troubleshooting
        if len(enriched_kinases) == 0:
            current_app.logger.warning(f"No enriched kinases found for pathway: {gs_name}")
            # Try a broader query to debug
            test_query = execute_query("""
                SELECT COUNT(*) as count FROM pathway_kinase_mapping 
                WHERE gs_name LIKE ?
            """, (f"%{gs_name.split()[0]}%",), fetch_one=True)
            current_app.logger.info(f"Partial name match found {test_query['count']} rows")

        # STEP 4: Get protein-level enrichment for proteins in this pathway
        proteins_with_enrichment = []

        if protein_list:
            # Get protein-kinase enrichments for top proteins in this pathway
            placeholders = ','.join(['?' for _ in protein_list[:25]])  # Limit for performance
            protein_kinase_data = execute_query(f"""
                SELECT 
                    gene_symbol,
                    kinase,
                    enrichment_score,
                    z_score_kinase,
                    n_sites,
                    n_above_q75,
                    n_above_q90,
                    max_score,
                    mean_score
                FROM protein_kinase_enrichment
                WHERE gene_symbol IN ({placeholders}) AND enrichment_score IS NOT NULL
                ORDER BY gene_symbol, enrichment_score DESC
            """, tuple(protein_list[:25]))

            # Group by protein
            protein_kinase_map = {}
            for row in protein_kinase_data:
                if row['gene_symbol'] not in protein_kinase_map:
                    protein_kinase_map[row['gene_symbol']] = []
                protein_kinase_map[row['gene_symbol']].append(dict(row))

            # Create protein objects with phosphosite stats and kinase enrichments
            for gene_symbol in protein_list[:25]:  # Show top 25 proteins
                # Get basic phosphosite stats
                protein_stats = execute_query("""
                    SELECT 
                        COUNT(*) as total_sites,
                        SUM(CASE WHEN significant_fdr05 = 1 THEN 1 ELSE 0 END) as significant_sites,
                        MAX(COALESCE(predicted_prob_calibrated, 0)) as max_confidence
                    FROM phosphosites
                    WHERE gene_symbol = ?
                """, (gene_symbol,), fetch_one=True)

                protein_data = {
                    'gene_symbol': gene_symbol,
                    'total_sites': protein_stats['total_sites'] if protein_stats else 0,
                    'significant_sites': protein_stats['significant_sites'] if protein_stats else 0,
                    'max_confidence': protein_stats['max_confidence'] if protein_stats else 0,
                    'top_kinases': protein_kinase_map.get(gene_symbol, [])[:10]  # Top 10 kinases per protein
                }

                proteins_with_enrichment.append(protein_data)

        # STEP 5: Get phosphorylation sites from proteins in this pathway
        sample_sites = []
        if protein_list:
            sample_proteins = protein_list[:10]  # Sample from first 10 proteins
            placeholders = ','.join(['?' for _ in sample_proteins])

            sample_sites = execute_query(f"""
                SELECT 
                    site_id,
                    gene_symbol, 
                    position, 
                    residue_type,
                    predicted_prob_calibrated,
                    qvalue,
                    significant_fdr05,
                    secondary_structure,
                    plddt
                FROM phosphosites
                WHERE gene_symbol IN ({placeholders}) 
                  AND predicted_prob_calibrated >= 0.5
                ORDER BY predicted_prob_calibrated DESC
                LIMIT 200
            """, tuple(sample_proteins))

        # STEP 6: Calculate kinase family distribution from enriched kinases
        kinase_families = {}
        if enriched_kinases:
            family_scores = {}
            for kinase in enriched_kinases:
                family = get_kinase_family(kinase['kinase'])
                score = kinase['score_integrated'] or 0
                if family not in family_scores:
                    family_scores[family] = []
                family_scores[family].append(score)

            for family, scores in family_scores.items():
                kinase_families[family] = {
                    'avg_score': sum(scores) / len(scores),
                    'count': len(scores)
                }

        # Response data
        response_data = {
            'pathway': pathway,
            'proteins': proteins_with_enrichment,
            'sample_sites': sample_sites,
            'enriched_kinases': enriched_kinases,
            'kinase_families': kinase_families,
            'pathway_stats': {
                'total_proteins': len(protein_list),
                'total_sites': len(sample_sites),
                'significant_sites': sum(1 for s in sample_sites if s.get('significant_fdr05')),
                'avg_confidence': sum(s.get('predicted_prob_calibrated', 0) for s in sample_sites) / len(sample_sites) if sample_sites else 0,
                'residue_distribution': {
                    'serine': sum(1 for s in sample_sites if s.get('residue_type') == 'S'),
                    'threonine': sum(1 for s in sample_sites if s.get('residue_type') == 'T'),
                    'tyrosine': sum(1 for s in sample_sites if s.get('residue_type') == 'Y')
                }
            },
            'analysis_metadata': {
                'generated_at': time.time(),
                'sites_analyzed': len(sample_sites),
                'kinases_analyzed': len(enriched_kinases),
                'proteins_analyzed': len(proteins_with_enrichment)
            }
        }

        return jsonify(response_data)

    except Exception as e:
        current_app.logger.error(f"=== PATHWAY ANALYSIS ERROR: {e} ===")
        import traceback
        current_app.logger.error(traceback.format_exc())
        return jsonify({'error': f'Failed to perform pathway analysis: {str(e)}'}), 500

@bp.route('/search/pathways')
def search_pathways():
    """CORRECTED: Pathway search that properly queries the database."""
    query = request.args.get('q', '').strip()
    collection = request.args.get('collection', '')
    min_genes = int(request.args.get('min_genes', 0))
    limit = min(int(request.args.get('limit', 50)), 100)

    if not query or len(query) < 2:
        return jsonify({'pathways': [], 'message': 'Query too short'})

    try:
        current_app.logger.info(f"Searching pathways for: '{query}', collection: '{collection}', min_genes: {min_genes}")

        # Build search conditions
        conditions = ["(gs.gs_name LIKE ? OR COALESCE(gs.gs_description, '') LIKE ?)"]
        params = [f'%{query}%', f'%{query}%']

        if collection:
            conditions.append("gs.gs_collection = ?")
            params.append(collection)

        # Query with proper joins
        sql = f"""
        SELECT DISTINCT
            gs.gs_id, 
            gs.gs_name, 
            COALESCE(gs.gs_collection, 'Unknown') as gs_collection, 
            COALESCE(gs.gs_description, 'No description available') as gs_description, 
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

        current_app.logger.info(f"Found {len(results)} pathway results")
        return jsonify({
            'pathways': results,
            'query': query,
            'filters': {'collection': collection, 'min_genes': min_genes}
        })

    except Exception as e:
        current_app.logger.error(f"Error searching pathways: {e}")
        return jsonify({'error': f'Search failed: {str(e)}'}), 500


@bp.route('/search/proteins')
def search_proteins():
    """CORRECTED: Protein search with proper enrichment data."""
    query = request.args.get('q', '').strip().upper()
    limit = min(int(request.args.get('limit', 50)), 100)

    if not query or len(query) < 1:
        return jsonify({'proteins': [], 'message': 'Query required'})

    try:
        current_app.logger.info(f"Searching proteins for: '{query}'")

        # Get protein search results with phosphosite stats
        results = execute_query("""
            SELECT 
                gene_symbol,
                COUNT(site_id) as site_count,
                SUM(CASE WHEN significant_fdr05 = 1 THEN 1 ELSE 0 END) as significant_sites,
                AVG(predicted_prob_calibrated) as avg_confidence,
                MAX(predicted_prob_calibrated) as max_confidence
            FROM phosphosites 
            WHERE gene_symbol LIKE ?
            GROUP BY gene_symbol
            ORDER BY significant_sites DESC, site_count DESC
            LIMIT ?
        """, (f'{query}%', limit))

        proteins = []
        for row in results:
            protein_dict = dict(row)

            # Get pathway count
            pathway_result = execute_query("""
                SELECT COUNT(DISTINCT gs_id) as pathway_count
                FROM gene_set_memberships
                WHERE gene_symbol = ?
            """, (protein_dict['gene_symbol'],), fetch_one=True)

            protein_dict['pathway_count'] = pathway_result['pathway_count'] if pathway_result else 0
            proteins.append(protein_dict)

        current_app.logger.info(f"Found {len(proteins)} protein results")
        return jsonify({'proteins': proteins})

    except Exception as e:
        current_app.logger.error(f"Error searching proteins: {e}")
        return jsonify({'error': f'Search failed: {str(e)}'}), 500


@bp.route('/stats/overview')
def get_overview_stats():
    """CORRECTED: Overview stats from real database."""
    try:
        current_app.logger.info("Getting overview statistics")

        stats = {
            'total_sites': 0,
            'total_proteins': 0,
            'total_pathways': 0,
            'significant_sites': 0,
            'quality_distribution': {'very_high': 0, 'high': 0, 'moderate': 0, 'low': 0},
            'pathway_collections': [],
            'residue_distribution': [],
            'generated_at': time.time()
        }

        # Basic counts
        basic_stats = execute_query("""
            SELECT 
                (SELECT COUNT(*) FROM phosphosites) as total_sites,
                (SELECT COUNT(*) FROM proteins) as total_proteins,
                (SELECT COUNT(*) FROM gene_sets) as total_pathways,
                (SELECT COUNT(*) FROM phosphosites WHERE significant_fdr05 = 1) as significant_sites
        """, fetch_one=True)

        if basic_stats:
            stats.update(basic_stats)

        # Quality distribution
        quality_stats = execute_query("""
            SELECT 
                SUM(CASE WHEN predicted_prob_calibrated > 0.9 THEN 1 ELSE 0 END) as very_high,
                SUM(CASE WHEN predicted_prob_calibrated BETWEEN 0.7 AND 0.9 THEN 1 ELSE 0 END) as high,
                SUM(CASE WHEN predicted_prob_calibrated BETWEEN 0.5 AND 0.7 THEN 1 ELSE 0 END) as moderate,
                SUM(CASE WHEN predicted_prob_calibrated < 0.5 THEN 1 ELSE 0 END) as low
            FROM phosphosites
            WHERE predicted_prob_calibrated IS NOT NULL
        """, fetch_one=True)

        if quality_stats:
            stats['quality_distribution'] = dict(quality_stats)

        # Pathway collections
        collections = execute_query("""
            SELECT gs_collection, COUNT(*) as count
            FROM gene_sets 
            WHERE gs_collection IS NOT NULL
            GROUP BY gs_collection 
            ORDER BY count DESC
            LIMIT 10
        """)
        stats['pathway_collections'] = collections

        # Residue distribution
        residues = execute_query("""
            SELECT residue_type, COUNT(*) as count
            FROM phosphosites 
            WHERE residue_type IS NOT NULL
            GROUP BY residue_type
            ORDER BY count DESC
        """)
        stats['residue_distribution'] = residues

        current_app.logger.info("Overview statistics completed")
        return jsonify(stats)

    except Exception as e:
        current_app.logger.error(f"Error getting overview stats: {e}")
        return jsonify({'error': 'Failed to retrieve statistics'}), 500


@bp.route('/pathway/<int:pathway_id>')
def get_pathway_details(pathway_id):
    """Get basic pathway details."""
    try:
        # Get pathway info
        pathway = execute_query("""
            SELECT gs_id, gs_name, gs_collection, gs_description, gs_pmid
            FROM gene_sets WHERE gs_id = ?
        """, (pathway_id,), fetch_one=True)

        if not pathway:
            return jsonify({'error': 'Pathway not found'}), 404

        # Get pathway statistics
        stats = execute_query("""
            SELECT 
                COUNT(DISTINCT gsm.gene_symbol) as total_genes,
                COUNT(DISTINCT p.site_id) as total_sites,
                COUNT(DISTINCT CASE WHEN p.significant_fdr05 = 1 THEN p.site_id END) as significant_sites
            FROM gene_set_memberships gsm
            LEFT JOIN phosphosites p ON gsm.gene_symbol = p.gene_symbol
            WHERE gsm.gs_id = ?
        """, (pathway_id,), fetch_one=True) or {'total_genes': 0, 'total_sites': 0, 'significant_sites': 0}

        # Get top proteins
        proteins = execute_query("""
            SELECT 
                gsm.gene_symbol,
                COUNT(p.site_id) as total_sites,
                COUNT(CASE WHEN p.significant_fdr05 = 1 THEN 1 END) as significant_sites,
                AVG(p.predicted_prob_calibrated) as avg_confidence
            FROM gene_set_memberships gsm
            LEFT JOIN phosphosites p ON gsm.gene_symbol = p.gene_symbol
            WHERE gsm.gs_id = ?
            GROUP BY gsm.gene_symbol
            ORDER BY significant_sites DESC, total_sites DESC
            LIMIT 50
        """, (pathway_id,))

        return jsonify({
            'pathway': pathway,
            'proteins': proteins,
            'stats': stats
        })

    except Exception as e:
        current_app.logger.error(f"Error getting pathway details: {e}")
        return jsonify({'error': 'Failed to retrieve pathway details'}), 500


@bp.route('/protein/<gene_symbol>')
def get_protein_details(gene_symbol):
    """Get detailed protein information."""
    gene_symbol = gene_symbol.upper()

    try:
        # Get protein basic info
        protein = execute_query("""
            SELECT gene_symbol, ncbi_gene_id, ensembl_gene_id
            FROM proteins WHERE gene_symbol = ?
        """, (gene_symbol,), fetch_one=True) or {'gene_symbol': gene_symbol}

        # Get phosphosites
        sites = execute_query("""
            SELECT 
                site_id, position, residue_type, motif,
                predicted_prob_calibrated, qvalue,
                significant_fdr05, significant_fdr10, significant_fdr20,
                plddt, secondary_structure, sasa_ratio, disorder_score
            FROM phosphosites 
            WHERE gene_symbol = ?
            ORDER BY position
        """, (gene_symbol,))

        if not sites:
            return jsonify({'error': 'Protein not found'}), 404

        # Get pathways
        pathways = execute_query("""
            SELECT gs.gs_id, gs.gs_name, gs.gs_collection, gs.gs_description
            FROM gene_set_memberships gsm
            JOIN gene_sets gs ON gsm.gs_id = gs.gs_id
            WHERE gsm.gene_symbol = ?
            ORDER BY gs.gs_collection, gs.gs_name
            LIMIT 20
        """, (gene_symbol,))

        # Calculate statistics
        stats = {
            'total_sites': len(sites),
            'significant_sites': sum(1 for s in sites if s.get('significant_fdr05')),
            'high_confidence': sum(1 for s in sites if (s.get('predicted_prob_calibrated') or 0) > 0.8),
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


def get_kinase_family(kinase_name):
    """Get kinase family for classification."""
    families = {
        'AGC': ['AKT1', 'AKT2', 'AKT3', 'PKA', 'PKACA', 'PKACB', 'PKACG', 'PKG1', 'PKG2', 'PKN1', 'PKN2', 'PKN3', 'ROCK1', 'ROCK2', 'P70S6K', 'P70S6KB', 'RSK2', 'RSK3', 'RSK4', 'P90RSK', 'SGK1', 'SGK3'],
        'CAMK': ['CAMK1A', 'CAMK1B', 'CAMK1D', 'CAMK1G', 'CAMK2A', 'CAMK2B', 'CAMK2D', 'CAMK2G', 'CAMK4', 'CAMKK1', 'CAMKK2', 'CAMLCK', 'DAPK1', 'DAPK2', 'DAPK3', 'DCAMKL1', 'DCAMKL2', 'MARK1', 'MARK2', 'MARK3', 'MARK4'],
        'CK1': ['CK1A', 'CK1A2', 'CK1D', 'CK1E', 'CK1G1', 'CK1G2', 'CK1G3'],
        'CMGC': ['CDK1', 'CDK2', 'CDK3', 'CDK4', 'CDK5', 'CDK6', 'CDK7', 'CDK8', 'CDK9', 'CDK10', 'CDK12', 'CDK13', 'CDK14', 'CDK16', 'CDK17', 'CDK18', 'CDK19', 'CDKL1', 'CDKL5', 'CLK1', 'CLK2', 'CLK3', 'CLK4', 'DYRK1A', 'DYRK1B', 'DYRK2', 'DYRK3', 'DYRK4', 'ERK1', 'ERK2', 'ERK5', 'ERK7', 'GSK3A', 'GSK3B', 'JNK1', 'JNK2', 'JNK3', 'MAPKAPK2', 'MAPKAPK3', 'MAPKAPK5', 'P38A', 'P38B', 'P38D', 'P38G'],
        'STE': ['MEK1', 'MEK2', 'MEK5', 'MEKK1', 'MEKK2', 'MEKK3', 'MEKK6', 'PAK1', 'PAK2', 'PAK3', 'PAK4', 'PAK5', 'PAK6', 'ASK1', 'TAK1', 'TAO1', 'TAO2', 'TAO3'],
        'TK': ['ALK2', 'ALK4'],
        'TKL': ['ACVR2A', 'ACVR2B', 'BMPR1A', 'BMPR1B', 'BMPR2', 'TGFBR1', 'TGFBR2', 'BRAF', 'RAF1', 'MLK1', 'MLK2', 'MLK3', 'MLK4', 'LRRK2', 'RIPK1', 'RIPK2', 'RIPK3']
    }

    for family, members in families.items():
        if kinase_name in members:
            return family
    return 'Other'