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


@bp.route('/pathway/<int:pathway_id>/analysis')
def get_pathway_analysis(pathway_id):
    """Optimized pathway-level phosphoproteomics analysis."""
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
        cursor = get_db().execute("SELECT * FROM phosphosites LIMIT 0")
        all_columns = [description[0] for description in cursor.description]
        kinase_columns = [col for col in all_columns if '_MotifScore' in col]

        current_app.logger.info(f"Found {len(kinase_columns)} kinase columns for pathway {pathway_id}")

        # OPTIMIZATION 1: Only get phosphocompetent sites initially (prob > 0.5, qvalue < 0.05)
        # This dramatically reduces the initial data load
        sites_query = f"""
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
            p.hydroxyl_exposure
        FROM gene_set_memberships gsm
        JOIN phosphosites p ON gsm.gene_symbol = p.gene_symbol
        WHERE gsm.gs_id = ?
          AND p.predicted_prob_calibrated >= 0.5
          AND p.qvalue <= 0.05
        ORDER BY p.predicted_prob_calibrated DESC
        LIMIT 500
        """

        sites_results = execute_query(sites_query, (pathway_id,))

        if not sites_results:
            # If no phosphocompetent sites, get ANY sites (limited)
            sites_query_fallback = """
            SELECT 
                p.site_id,
                p.gene_symbol,
                p.position,
                p.residue_type,
                p.motif,
                p.predicted_prob_calibrated,
                p.qvalue,
                p.significant_fdr05,
                p.plddt,
                p.secondary_structure,
                p.sasa_ratio,
                p.disorder_score
            FROM gene_set_memberships gsm
            JOIN phosphosites p ON gsm.gene_symbol = p.gene_symbol
            WHERE gsm.gs_id = ?
            ORDER BY p.predicted_prob_calibrated DESC
            LIMIT 100
            """
            sites_results = execute_query(sites_query_fallback, (pathway_id,))

            if not sites_results:
                return jsonify({'error': 'No phosphorylation sites found for this pathway'}), 404

        # Convert to list of dicts
        all_sites = [dict(row) for row in sites_results]

        current_app.logger.info(f"Processing {len(all_sites)} phosphocompetent sites")

        # OPTIMIZATION 2: Only get top kinases for each site, not all 303
        # For each site, get only the top 3 kinase predictions
        for site in all_sites:
            site_id = site['site_id']

            # Build a query to get just this site's kinase scores
            kinase_cols_str = ', '.join(kinase_columns)
            kinase_query = f"""
            SELECT {kinase_cols_str}
            FROM phosphosites
            WHERE site_id = ?
            """

            kinase_result = execute_query(kinase_query, (site_id,), fetch_one=True)

            if kinase_result:
                # Extract and sort kinase scores
                kinase_scores = []
                for kinase_col in kinase_columns:
                    score = kinase_result[kinase_col]
                    if score is not None and score > 0:
                        kinase_name = kinase_col.replace('_MotifScore', '')
                        kinase_scores.append({'kinase': kinase_name, 'score': float(score)})

                # Sort and get top 3
                kinase_scores.sort(key=lambda x: x['score'], reverse=True)
                site['top_kinases'] = kinase_scores[:3]
            else:
                site['top_kinases'] = []

        # Organize by protein with statistics
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
                        'avg_confidence': 0
                    }
                }
            proteins_data[gene]['sites'].append(site)

        # Calculate protein-level statistics
        for gene, data in proteins_data.items():
            sites = data['sites']
            data['stats']['total_sites'] = len(sites)
            data['stats']['significant_sites'] = sum(1 for s in sites if s.get('significant_fdr05'))
            data['stats']['phosphocompetent_sites'] = len(sites)  # All are phosphocompetent due to filter
            data['stats']['high_confidence_sites'] = sum(
                1 for s in sites if (s.get('predicted_prob_calibrated', 0) > 0.8))

            confidences = [s.get('predicted_prob_calibrated', 0) for s in sites]
            data['stats']['avg_confidence'] = sum(confidences) / len(confidences) if confidences else 0

        # Calculate pathway-wide statistics
        pathway_stats = {
            'total_proteins': len(proteins_data),
            'total_sites': len(all_sites),
            'significant_sites': sum(1 for s in all_sites if s.get('significant_fdr05')),
            'phosphocompetent_sites': len(all_sites),
            'high_confidence_sites': sum(1 for s in all_sites if (s.get('predicted_prob_calibrated', 0) > 0.8)),
            'avg_confidence': sum(s.get('predicted_prob_calibrated', 0) for s in all_sites) / len(
                all_sites) if all_sites else 0,
            'serine_sites': sum(1 for s in all_sites if s.get('residue_type') == 'S'),
            'threonine_sites': sum(1 for s in all_sites if s.get('residue_type') == 'T'),
            'tyrosine_sites': sum(1 for s in all_sites if s.get('residue_type') == 'Y')
        }

        # Get enriched kinases across the pathway (limited to top 20)
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
                        'scores': [],
                        'proteins': set()
                    }

                kinase_enrichment[kinase]['site_count'] += 1
                kinase_enrichment[kinase]['scores'].append(kinase_data['score'])
                kinase_enrichment[kinase]['proteins'].add(site['gene_symbol'])
                kinase_enrichment[kinase]['phosphocompetent_sites'] += 1

        # Calculate final kinase statistics
        enriched_kinases = []
        for kinase_data in kinase_enrichment.values():
            scores = kinase_data['scores']
            if scores:
                kinase_data['avg_score'] = sum(scores) / len(scores)
                kinase_data['protein_count'] = len(kinase_data['proteins'])
                del kinase_data['proteins']  # Remove set for JSON serialization
                del kinase_data['scores']
                enriched_kinases.append(kinase_data)

        # Sort and limit to top 20 kinases
        enriched_kinases.sort(key=lambda x: x['phosphocompetent_sites'], reverse=True)
        enriched_kinases = enriched_kinases[:20]

        # Structure distribution
        structure_distribution = {}
        for site in all_sites:
            struct = site.get('secondary_structure', 'Unknown')
            structure_distribution[struct] = structure_distribution.get(struct, 0) + 1

        # Confidence distribution bins
        confidence_bins = {
            '0.5-0.7': 0,
            '0.7-0.9': 0,
            '0.9-1.0': 0
        }
        for site in all_sites:
            conf = site.get('predicted_prob_calibrated', 0)
            if conf < 0.7:
                confidence_bins['0.5-0.7'] += 1
            elif conf < 0.9:
                confidence_bins['0.7-0.9'] += 1
            else:
                confidence_bins['0.9-1.0'] += 1

        return jsonify({
            'pathway': dict(pathway),
            'proteins': list(proteins_data.values()),
            'all_sites': all_sites,
            'enriched_kinases': enriched_kinases,
            'pathway_stats': pathway_stats,
            'structure_distribution': structure_distribution,
            'confidence_distribution': confidence_bins,
            'analysis_metadata': {
                'generated_at': time.time(),
                'sites_analyzed': len(all_sites),
                'note': 'Showing top 500 phosphocompetent sites (prob≥0.5, q≤0.05)',
                'total_kinases_found': len(enriched_kinases),
                'analysis_type': 'optimized_phosphocompetent_focus'
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


# app/api/routes.py - ADD this enhanced endpoint to your existing routes.py

@bp.route('/pathway/<int:pathway_id>/analysis_v2')
def get_pathway_analysis_v2(pathway_id):
    """
    Optimized pathway analysis leveraging pre-computed kinase enrichment tables.
    This replaces the inefficient /pathway/<int:pathway_id>/analysis endpoint.
    """
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

        pathway_name = pathway['gs_name']

        # 1. Get pathway-level kinase enrichment (pre-computed!)
        pathway_kinases_query = """
        SELECT 
            kinase,
            score_integrated,
            score_zscore_only,
            enrichment_ratio,
            n_genes_with_kinase,
            n_total_sites,
            gene_coverage,
            fisher_p,
            cohens_d,
            q90_enrichment,
            n_genes_above_kinase_q90,
            n_sites_above_kinase_q90,
            rank_integrated
        FROM pathway_kinase_mapping
        WHERE gs_name = ?
        ORDER BY score_integrated DESC
        LIMIT 50
        """
        pathway_kinases = execute_query(pathway_kinases_query, (pathway_name,))

        # 2. Get pathway statistics
        stats_query = """
        SELECT 
            COUNT(DISTINCT gsm.gene_symbol) as total_genes,
            COUNT(DISTINCT p.site_id) as total_sites,
            COUNT(DISTINCT CASE WHEN p.significant_fdr05 = 1 THEN p.site_id END) as significant_sites,
            COUNT(DISTINCT CASE WHEN p.predicted_prob_calibrated > 0.8 THEN p.site_id END) as high_confidence_sites,
            AVG(p.predicted_prob_calibrated) as avg_confidence
        FROM gene_set_memberships gsm
        LEFT JOIN phosphosites p ON gsm.gene_symbol = p.gene_symbol
        WHERE gsm.gs_id = ?
        """
        stats_result = execute_query(stats_query, (pathway_id,), fetch_one=True)
        pathway_stats = dict(stats_result) if stats_result else {}

        # 3. Get proteins with their pre-computed kinase enrichments
        proteins_query = """
        SELECT DISTINCT
            gsm.gene_symbol,
            COUNT(DISTINCT p.site_id) as total_sites,
            COUNT(DISTINCT CASE WHEN p.significant_fdr05 = 1 THEN p.site_id END) as significant_sites,
            AVG(p.predicted_prob_calibrated) as avg_confidence,
            MAX(p.predicted_prob_calibrated) as max_confidence
        FROM gene_set_memberships gsm
        LEFT JOIN phosphosites p ON gsm.gene_symbol = p.gene_symbol
        WHERE gsm.gs_id = ?
        GROUP BY gsm.gene_symbol
        HAVING total_sites > 0
        ORDER BY significant_sites DESC, total_sites DESC
        LIMIT 100
        """
        proteins_results = execute_query(proteins_query, (pathway_id,))
        proteins_list = [dict(row) for row in proteins_results] if proteins_results else []

        # For each protein, get its top kinases from pre-computed table
        for protein in proteins_list[:20]:  # Limit to top 20 for performance
            kinase_query = """
            SELECT 
                kinase,
                enrichment_score,
                z_score_kinase,
                n_sites,
                n_above_q90,
                max_score
            FROM protein_kinase_enrichment
            WHERE gene_symbol = ?
            ORDER BY enrichment_score DESC
            LIMIT 5
            """
            protein['top_kinases'] = [dict(row) for row in execute_query(kinase_query, (protein['gene_symbol'],))]

        # 4. Get phosphocompetent sites (limited sample for performance)
        sites_query = """
        SELECT 
            p.site_id,
            p.gene_symbol,
            p.position,
            p.residue_type,
            p.motif,
            p.predicted_prob_calibrated,
            p.qvalue,
            p.significant_fdr05,
            p.secondary_structure,
            p.plddt,
            p.disorder_score,
            p.sasa_ratio
        FROM phosphosites p
        JOIN gene_set_memberships gsm ON p.gene_symbol = gsm.gene_symbol
        WHERE gsm.gs_id = ?
            AND p.predicted_prob_calibrated >= 0.5
            AND p.qvalue <= 0.05
        ORDER BY p.predicted_prob_calibrated DESC
        LIMIT 200
        """
        sites_results = execute_query(sites_query, (pathway_id,))
        sites = [dict(row) for row in sites_results] if sites_results else []

        # 5. Get kinase families distribution
        kinase_families = {
            'AGC': ['AKT1', 'AKT2', 'AKT3', 'PKA', 'PKG1', 'PKG2', 'ROCK1', 'ROCK2', 'SGK1', 'SGK3'],
            'CAMK': ['CAMK1A', 'CAMK2A', 'CAMK2B', 'CAMK2D', 'CAMK4', 'MARK1', 'MARK2'],
            'CMGC': ['CDK1', 'CDK2', 'CDK4', 'CDK5', 'GSK3A', 'GSK3B', 'ERK1', 'ERK2', 'JNK1', 'P38A'],
            'CK1': ['CK1A', 'CK1D', 'CK1E'],
            'STE': ['MEK1', 'MEK2', 'PAK1', 'PAK2'],
            'TKL': ['BRAF', 'RAF1', 'RIPK1', 'RIPK2']
        }

        family_distribution = {}
        for family, members in kinase_families.items():
            family_kinases = [k for k in pathway_kinases if dict(k)['kinase'] in members]
            if family_kinases:
                family_distribution[family] = {
                    'count': len(family_kinases),
                    'avg_score': sum(dict(k)['score_integrated'] for k in family_kinases) / len(family_kinases),
                    'top_kinase': dict(family_kinases[0])['kinase'] if family_kinases else None
                }

        # 6. Calculate kinase network connections (which kinases co-occur)
        if len(proteins_list) > 1:
            kinase_network = []
            top_kinases_list = list(set([dict(k)['kinase'] for k in pathway_kinases[:10]]))

            for i, k1 in enumerate(top_kinases_list):
                for k2 in top_kinases_list[i + 1:]:
                    # Count proteins where both kinases are enriched
                    cooccurrence_query = """
                    SELECT COUNT(DISTINCT pk1.gene_symbol) as shared_targets
                    FROM protein_kinase_enrichment pk1
                    JOIN protein_kinase_enrichment pk2 ON pk1.gene_symbol = pk2.gene_symbol
                    JOIN gene_set_memberships gsm ON pk1.gene_symbol = gsm.gene_symbol
                    WHERE gsm.gs_id = ?
                        AND pk1.kinase = ?
                        AND pk2.kinase = ?
                        AND pk1.enrichment_score > 1
                        AND pk2.enrichment_score > 1
                    """
                    result = execute_query(cooccurrence_query, (pathway_id, k1, k2), fetch_one=True)
                    if result and result['shared_targets'] > 0:
                        kinase_network.append({
                            'source': k1,
                            'target': k2,
                            'weight': result['shared_targets']
                        })
        else:
            kinase_network = []

        # 7. Time-series mock data for kinase activity (could be replaced with real experimental data)
        # This demonstrates how the pathway might respond to stimulation
        time_series_data = {
            'timepoints': [0, 5, 15, 30, 60, 120],  # minutes
            'kinases': {}
        }

        for kinase_data in pathway_kinases[:5]:
            kinase = dict(kinase_data)['kinase']
            base_activity = dict(kinase_data)['score_integrated']
            # Simulate activation curves
            if kinase in ['ERK1', 'ERK2', 'AKT1']:  # Early responders
                curve = [0, base_activity * 0.8, base_activity * 1.2, base_activity * 0.9, base_activity * 0.6,
                         base_activity * 0.4]
            elif kinase in ['GSK3B', 'CDK1']:  # Late responders
                curve = [0, base_activity * 0.2, base_activity * 0.5, base_activity * 0.8, base_activity * 1.0,
                         base_activity * 0.9]
            else:  # Sustained
                curve = [0, base_activity * 0.5, base_activity * 0.8, base_activity * 1.0, base_activity * 1.0,
                         base_activity * 0.95]
            time_series_data['kinases'][kinase] = curve

        return jsonify({
            'pathway': dict(pathway),
            'stats': pathway_stats,
            'enriched_kinases': [dict(row) for row in pathway_kinases],
            'proteins': proteins_list,
            'sample_sites': sites,
            'kinase_families': family_distribution,
            'kinase_network': kinase_network,
            'time_series': time_series_data,
            'analysis_metadata': {
                'generated_at': time.time(),
                'total_kinases_analyzed': len(pathway_kinases),
                'proteins_with_kinase_data': len([p for p in proteins_list if p.get('top_kinases')]),
                'version': '2.0'
            }
        })

    except Exception as e:
        current_app.logger.error(f"Error in pathway analysis v2: {e}")
        import traceback
        current_app.logger.error(traceback.format_exc())
        return jsonify({'error': f'Failed to perform pathway analysis: {str(e)}'}), 500