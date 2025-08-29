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


# Enhanced pathway analysis endpoint in routes.py
@bp.route('/pathway/<int:pathway_id>/analysis')
def get_pathway_analysis(pathway_id):
    """Comprehensive pathway-level phosphoproteomics analysis with enhanced metrics."""
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

        # Get list of all kinase columns for comprehensive analysis
        cursor = get_db().execute("SELECT * FROM phosphosites LIMIT 0")
        all_columns = [description[0] for description in cursor.description]
        kinase_columns = [col for col in all_columns if '_MotifScore' in col]

        current_app.logger.info(f"Analyzing pathway {pathway_id} with {len(kinase_columns)} kinases")

        # OPTIMIZED: Get high-quality phosphocompetent sites (prob ≥ 0.5, q ≤ 0.05)
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
            p.hydroxyl_exposure,
            {', '.join(kinase_columns)}
        FROM gene_set_memberships gsm
        JOIN phosphosites p ON gsm.gene_symbol = p.gene_symbol
        WHERE gsm.gs_id = ?
          AND p.predicted_prob_calibrated >= 0.5
          AND p.qvalue <= 0.05
        ORDER BY p.predicted_prob_calibrated DESC
        LIMIT 1000
        """

        sites_results = execute_query(sites_query, (pathway_id,))

        if not sites_results:
            return jsonify({'error': 'No high-quality phosphorylation sites found for this pathway'}), 404

        # Convert to list of dicts for easier processing
        all_sites = [dict(row) for row in sites_results]

        # Calculate kinase enrichment with 75th percentile analysis
        kinase_enrichment = {}
        kinase_percentiles = {}

        # First pass: collect all scores per kinase to calculate percentiles
        for kinase_col in kinase_columns:
            kinase_name = kinase_col.replace('_MotifScore', '')
            scores = []
            for site in all_sites:
                score = site.get(kinase_col)
                if score is not None and score > 0:
                    scores.append(float(score))

            if scores:
                scores.sort()
                # Calculate 75th percentile
                percentile_75_idx = int(len(scores) * 0.75)
                kinase_percentiles[kinase_name] = scores[percentile_75_idx] if percentile_75_idx < len(scores) else \
                scores[-1]

        # Second pass: analyze kinase enrichment with percentile thresholds
        for site in all_sites:
            for kinase_col in kinase_columns:
                kinase_name = kinase_col.replace('_MotifScore', '')
                score = site.get(kinase_col)

                if score is not None and score > 0:
                    if kinase_name not in kinase_enrichment:
                        kinase_enrichment[kinase_name] = {
                            'kinase': kinase_name,
                            'total_sites': 0,
                            'high_confidence_sites': 0,  # above 75th percentile
                            'proteins': set(),
                            'scores': [],
                            'avg_score': 0,
                            'percentile_75': kinase_percentiles.get(kinase_name, 0)
                        }

                    enrichment = kinase_enrichment[kinase_name]
                    enrichment['total_sites'] += 1
                    enrichment['scores'].append(float(score))
                    enrichment['proteins'].add(site['gene_symbol'])

                    # Check if above 75th percentile for this kinase
                    if score >= kinase_percentiles.get(kinase_name, 0):
                        enrichment['high_confidence_sites'] += 1

        # Finalize kinase enrichment statistics
        enriched_kinases = []
        for kinase_data in kinase_enrichment.values():
            if kinase_data['scores']:
                kinase_data['avg_score'] = sum(kinase_data['scores']) / len(kinase_data['scores'])
                kinase_data['protein_count'] = len(kinase_data['proteins'])
                # Remove the set for JSON serialization
                del kinase_data['proteins']
                del kinase_data['scores']
                enriched_kinases.append(kinase_data)

        # Sort kinases by high confidence sites (above 75th percentile)
        enriched_kinases.sort(key=lambda x: x['high_confidence_sites'], reverse=True)
        enriched_kinases = enriched_kinases[:30]  # Top 30 kinases

        # Organize data by protein with enhanced kinase analysis
        proteins_data = {}
        for site in all_sites:
            gene = site['gene_symbol']
            if gene not in proteins_data:
                proteins_data[gene] = {
                    'gene_symbol': gene,
                    'sites': [],
                    'kinase_associations': {},  # Store kinase-protein specific data
                    'stats': {
                        'total_sites': 0,
                        'significant_sites': 0,
                        'high_confidence_sites': 0,
                        'max_confidence': 0
                    }
                }

            # Add site to protein
            proteins_data[gene]['sites'].append(site)

            # Analyze kinase associations for this protein
            for kinase_col in kinase_columns:
                kinase_name = kinase_col.replace('_MotifScore', '')
                score = site.get(kinase_col)

                if score is not None and score > 0:
                    if kinase_name not in proteins_data[gene]['kinase_associations']:
                        proteins_data[gene]['kinase_associations'][kinase_name] = {
                            'kinase': kinase_name,
                            'sites': [],
                            'above_75th_count': 0,
                            'avg_score': 0,
                            'max_score': 0
                        }

                    assoc = proteins_data[gene]['kinase_associations'][kinase_name]
                    assoc['sites'].append({
                        'site_id': site['site_id'],
                        'position': site['position'],
                        'residue_type': site['residue_type'],
                        'score': float(score),
                        'confidence': site['predicted_prob_calibrated'],
                        'significant': site['significant_fdr05']
                    })

                    # Check if above 75th percentile
                    if score >= kinase_percentiles.get(kinase_name, 0):
                        assoc['above_75th_count'] += 1

                    assoc['max_score'] = max(assoc['max_score'], float(score))

        # Finalize protein statistics and kinase associations
        for gene, data in proteins_data.items():
            sites = data['sites']

            # Calculate protein stats
            data['stats']['total_sites'] = len(sites)
            data['stats']['significant_sites'] = sum(1 for s in sites if s.get('significant_fdr05'))
            data['stats']['high_confidence_sites'] = sum(
                1 for s in sites if (s.get('predicted_prob_calibrated', 0) > 0.8))

            confidences = [s.get('predicted_prob_calibrated', 0) for s in sites]
            data['stats']['max_confidence'] = max(confidences) if confidences else 0

            # Finalize kinase associations
            for kinase_name, assoc in data['kinase_associations'].items():
                if assoc['sites']:
                    assoc['avg_score'] = sum(s['score'] for s in assoc['sites']) / len(assoc['sites'])
                    assoc['site_count'] = len(assoc['sites'])

            # Get top 10 kinases by sites above 75th percentile
            top_kinases = sorted(
                data['kinase_associations'].values(),
                key=lambda x: x['above_75th_count'],
                reverse=True
            )[:10]

            data['top_kinases'] = top_kinases

        # Calculate comprehensive pathway statistics
        pathway_stats = {
            'total_proteins': len(proteins_data),
            'total_sites': len(all_sites),
            'significant_sites': sum(1 for s in all_sites if s.get('significant_fdr05')),
            'high_confidence_sites': sum(1 for s in all_sites if (s.get('predicted_prob_calibrated', 0) > 0.8)),
            'avg_confidence': sum(s.get('predicted_prob_calibrated', 0) for s in all_sites) / len(
                all_sites) if all_sites else 0,
            'residue_distribution': {
                'serine': sum(1 for s in all_sites if s.get('residue_type') == 'S'),
                'threonine': sum(1 for s in all_sites if s.get('residue_type') == 'T'),
                'tyrosine': sum(1 for s in all_sites if s.get('residue_type') == 'Y')
            }
        }

        # Secondary structure distribution
        structure_distribution = {}
        disorder_scores = []
        sasa_scores = []

        for site in all_sites:
            # Structure
            struct = site.get('secondary_structure', 'Unknown')
            structure_distribution[struct] = structure_distribution.get(struct, 0) + 1

            # Collect scores for statistics
            if site.get('disorder_score') is not None:
                disorder_scores.append(site['disorder_score'])
            if site.get('sasa_ratio') is not None:
                sasa_scores.append(site['sasa_ratio'])

        # Confidence distribution
        confidence_distribution = {'low': 0, 'medium': 0, 'high': 0}
        for site in all_sites:
            conf = site.get('predicted_prob_calibrated', 0)
            if conf < 0.7:
                confidence_distribution['low'] += 1
            elif conf < 0.9:
                confidence_distribution['medium'] += 1
            else:
                confidence_distribution['high'] += 1

        # Calculate kinase co-enrichment network data
        kinase_network = calculate_kinase_coenrichment(enriched_kinases[:20], all_sites, kinase_columns)

        return jsonify({
            'pathway': dict(pathway),
            'proteins': list(proteins_data.values()),
            'all_sites': all_sites,
            'enriched_kinases': enriched_kinases,
            'kinase_network': kinase_network,
            'pathway_stats': pathway_stats,
            'structure_distribution': structure_distribution,
            'confidence_distribution': confidence_distribution,
            'structural_stats': {
                'disorder': {
                    'mean': sum(disorder_scores) / len(disorder_scores) if disorder_scores else 0,
                    'median': sorted(disorder_scores)[len(disorder_scores) // 2] if disorder_scores else 0
                },
                'sasa': {
                    'mean': sum(sasa_scores) / len(sasa_scores) if sasa_scores else 0,
                    'median': sorted(sasa_scores)[len(sasa_scores) // 2] if sasa_scores else 0
                }
            },
            'analysis_metadata': {
                'generated_at': time.time(),
                'sites_analyzed': len(all_sites),
                'kinases_analyzed': len(enriched_kinases),
                'analysis_type': 'enhanced_phosphocompetent_analysis'
            }
        })

    except Exception as e:
        current_app.logger.error(f"Error in enhanced pathway analysis: {e}")
        import traceback
        current_app.logger.error(traceback.format_exc())
        return jsonify({'error': f'Failed to perform pathway analysis: {str(e)}'}), 500


def calculate_kinase_coenrichment(kinases, sites, kinase_columns):
    """Calculate kinase co-enrichment network for pathway visualization."""
    try:
        # Create protein-kinase association matrix
        protein_kinase_matrix = {}

        for site in sites:
            protein = site['gene_symbol']
            if protein not in protein_kinase_matrix:
                protein_kinase_matrix[protein] = {}

            for kinase_col in kinase_columns:
                kinase_name = kinase_col.replace('_MotifScore', '')
                if kinase_name in [k['kinase'] for k in kinases]:
                    score = site.get(kinase_col, 0)
                    if score > 0:
                        if kinase_name not in protein_kinase_matrix[protein]:
                            protein_kinase_matrix[protein][kinase_name] = []
                        protein_kinase_matrix[protein][kinase_name].append(score)

        # Calculate kinase-kinase co-occurrence
        kinase_names = [k['kinase'] for k in kinases]
        cooccurrence_matrix = {}

        for i, kinase1 in enumerate(kinase_names):
            cooccurrence_matrix[kinase1] = {}
            for j, kinase2 in enumerate(kinase_names):
                if i != j:
                    # Count proteins where both kinases are active
                    shared_proteins = 0
                    for protein, kinase_data in protein_kinase_matrix.items():
                        if kinase1 in kinase_data and kinase2 in kinase_data:
                            shared_proteins += 1

                    cooccurrence_matrix[kinase1][kinase2] = shared_proteins

        # Create network nodes and edges
        nodes = []
        for kinase in kinases:
            nodes.append({
                'id': kinase['kinase'],
                'label': kinase['kinase'],
                'size': kinase['high_confidence_sites'],
                'group': kinase['kinase'][:3],  # Group by first 3 letters for coloring
                'total_sites': kinase['total_sites'],
                'protein_count': kinase['protein_count']
            })

        edges = []
        for kinase1 in kinase_names:
            for kinase2 in kinase_names:
                if kinase1 != kinase2:
                    weight = cooccurrence_matrix[kinase1][kinase2]
                    if weight > 2:  # Only show meaningful connections
                        edges.append({
                            'source': kinase1,
                            'target': kinase2,
                            'weight': weight,
                            'id': f"{kinase1}-{kinase2}"
                        })

        return {
            'nodes': nodes,
            'edges': edges,
            'stats': {
                'node_count': len(nodes),
                'edge_count': len(edges),
                'max_connections': max(
                    [len([e for e in edges if e['source'] == n['id']]) for n in nodes]) if edges else 0
            }
        }

    except Exception as e:
        current_app.logger.error(f"Error calculating kinase co-enrichment: {e}")
        return {'nodes': [], 'edges': [], 'stats': {'node_count': 0, 'edge_count': 0, 'max_connections': 0}}


@bp.route('/pathway/<int:pathway_id>/protein/<gene_symbol>/kinases')
def get_protein_kinase_details(pathway_id, gene_symbol):
    """Get detailed kinase-protein associations for expanded view."""
    try:
        # Get kinase columns
        cursor = get_db().execute("SELECT * FROM phosphosites LIMIT 0")
        all_columns = [description[0] for description in cursor.description]
        kinase_columns = [col for col in all_columns if '_MotifScore' in col]

        # Get sites for this protein in this pathway
        kinase_cols_str = ', '.join(kinase_columns)
        query = f"""
        SELECT 
            p.site_id,
            p.position,
            p.residue_type,
            p.motif,
            p.predicted_prob_calibrated,
            p.qvalue,
            p.significant_fdr05,
            p.secondary_structure,
            p.plddt,
            {kinase_cols_str}
        FROM gene_set_memberships gsm
        JOIN phosphosites p ON gsm.gene_symbol = p.gene_symbol
        WHERE gsm.gs_id = ? AND p.gene_symbol = ?
          AND p.predicted_prob_calibrated >= 0.5
          AND p.qvalue <= 0.05
        ORDER BY p.predicted_prob_calibrated DESC
        """

        sites = execute_query(query, (pathway_id, gene_symbol.upper()))

        if not sites:
            return jsonify({'error': 'No sites found'}), 404

        # Process kinase associations
        kinase_data = {}
        for site in sites:
            site_dict = dict(site)
            for kinase_col in kinase_columns:
                kinase_name = kinase_col.replace('_MotifScore', '')
                score = site_dict.get(kinase_col)

                if score is not None and score > 0:
                    if kinase_name not in kinase_data:
                        kinase_data[kinase_name] = {
                            'kinase': kinase_name,
                            'sites': [],
                            'avg_score': 0,
                            'max_score': 0,
                            'site_count': 0
                        }

                    kinase_data[kinase_name]['sites'].append({
                        'site_id': site_dict['site_id'],
                        'position': site_dict['position'],
                        'residue_type': site_dict['residue_type'],
                        'motif': site_dict['motif'],
                        'score': float(score),
                        'confidence': site_dict['predicted_prob_calibrated'],
                        'qvalue': site_dict['qvalue'],
                        'significant': site_dict['significant_fdr05'],
                        'structure': site_dict['secondary_structure'],
                        'plddt': site_dict['plddt']
                    })

        # Calculate statistics
        for kinase_name, data in kinase_data.items():
            scores = [s['score'] for s in data['sites']]
            if scores:
                data['avg_score'] = sum(scores) / len(scores)
                data['max_score'] = max(scores)
                data['site_count'] = len(scores)

        # Sort by average score
        sorted_kinases = sorted(kinase_data.values(), key=lambda x: x['avg_score'], reverse=True)

        return jsonify({
            'gene_symbol': gene_symbol.upper(),
            'pathway_id': pathway_id,
            'kinase_associations': sorted_kinases,
            'total_kinases': len(sorted_kinases),
            'total_associations': sum(len(k['sites']) for k in sorted_kinases)
        })

    except Exception as e:
        current_app.logger.error(f"Error getting protein kinase details: {e}")
        return jsonify({'error': 'Failed to retrieve kinase details'}), 500