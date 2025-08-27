# app/api/analysis.py
"""Advanced analysis endpoints for discovery lab."""

from flask import request, jsonify, current_app
from app.api import bp
from app.database import DatabaseManager
from app import cache
import json
from collections import defaultdict, Counter
from typing import Dict, List, Any
import uuid


@bp.route('/analysis/compare-proteins', methods=['POST'])
def compare_proteins():
    """Compare multiple proteins side by side."""
    try:
        data = request.get_json()
        gene_symbols = [g.upper() for g in data.get('genes', [])]

        if not gene_symbols or len(gene_symbols) < 2:
            return jsonify({'error': 'At least 2 proteins required for comparison'}), 400

        if len(gene_symbols) > 10:
            return jsonify({'error': 'Maximum 10 proteins allowed for comparison'}), 400

        comparison_results = {}

        for gene in gene_symbols:
            protein_data = DatabaseManager.get_protein_summary(gene)
            if protein_data:
                # Get top kinases for this protein
                kinase_query = """
                SELECT kinase_scores
                FROM phosphosites
                WHERE gene_symbol = ? AND kinase_scores IS NOT NULL
                """

                kinase_results = DatabaseManager.execute_query(kinase_query, (gene,))

                # Aggregate kinase scores across all sites
                kinase_totals = defaultdict(list)
                for row in kinase_results:
                    try:
                        scores = json.loads(row['kinase_scores'])
                        for kinase, score in scores.items():
                            if score is not None:
                                kinase_totals[kinase].append(score)
                    except (json.JSONDecodeError, TypeError):
                        continue

                # Calculate average kinase scores
                avg_kinases = {}
                for kinase, scores in kinase_totals.items():
                    avg_kinases[kinase] = sum(scores) / len(scores)

                # Get top 10 kinases
                top_kinases = sorted(avg_kinases.items(), key=lambda x: x[1], reverse=True)[:10]

                comparison_results[gene] = {
                    **protein_data,
                    'top_kinases': [{'kinase': k, 'avg_score': round(v, 3)} for k, v in top_kinases]
                }

        # Find common pathways
        common_pathways = find_common_pathways(gene_symbols)

        # Calculate similarity metrics
        similarity_matrix = calculate_protein_similarity(gene_symbols, comparison_results)

        return jsonify({
            'comparison': comparison_results,
            'common_pathways': common_pathways,
            'similarity_matrix': similarity_matrix,
            'analysis_id': str(uuid.uuid4())
        })

    except Exception as e:
        current_app.logger.error(f"Error in protein comparison: {e}")
        return jsonify({'error': 'Failed to perform protein comparison'}), 500


@bp.route('/analysis/pathway-enrichment', methods=['POST'])
def pathway_enrichment_analysis():
    """Perform pathway enrichment analysis on a gene list."""
    try:
        data = request.get_json()
        gene_list = [g.upper() for g in data.get('genes', [])]
        min_overlap = data.get('min_overlap', 2)
        max_pathways = min(data.get('max_pathways', 50), 100)

        if not gene_list:
            return jsonify({'error': 'Gene list is required'}), 400

        enriched_pathways = DatabaseManager.get_pathway_enrichment(gene_list, min_overlap)

        # Calculate statistical significance (simplified Fisher's exact test approximation)
        total_genes_in_db = DatabaseManager.execute_query(
            "SELECT COUNT(DISTINCT gene_symbol) FROM proteins", fetch_one=True
        )['count']

        for pathway in enriched_pathways[:max_pathways]:
            # Simplified p-value calculation (hypergeometric distribution approximation)
            k = pathway['overlap_count']  # overlapping genes
            n = len(gene_list)  # genes in query
            K = pathway['total_genes']  # genes in pathway
            N = total_genes_in_db  # total genes in database

            # Simple enrichment score
            expected = (n * K) / N
            enrichment_score = k / expected if expected > 0 else 0
            pathway['enrichment_score'] = round(enrichment_score, 3)
            pathway['expected_overlap'] = round(expected, 2)

        return jsonify({
            'enriched_pathways': enriched_pathways[:max_pathways],
            'query_genes': gene_list,
            'parameters': {
                'min_overlap': min_overlap,
                'max_pathways': max_pathways
            },
            'analysis_id': str(uuid.uuid4())
        })

    except Exception as e:
        current_app.logger.error(f"Error in pathway enrichment analysis: {e}")
        return jsonify({'error': 'Failed to perform pathway enrichment analysis'}), 500


@bp.route('/analysis/structural-features', methods=['POST'])
def analyze_structural_features():
    """Analyze structural features across a set of sites."""
    try:
        data = request.get_json()
        gene_symbols = [g.upper() for g in data.get('genes', [])]
        confidence_threshold = float(data.get('confidence_threshold', 0.5))

        if not gene_symbols:
            return jsonify({'error': 'Gene list is required'}), 400

        # Get structural features for all sites
        placeholders = ','.join(['?' for _ in gene_symbols])
        query = f"""
        SELECT 
            gene_symbol,
            position,
            residue_type,
            predicted_prob_calibrated,
            secondary_structure,
            plddt,
            disorder_score,
            sasa_ratio,
            structural_features
        FROM phosphosites 
        WHERE gene_symbol IN ({placeholders})
        AND predicted_prob_calibrated >= ?
        ORDER BY predicted_prob_calibrated DESC
        """

        params = tuple(gene_symbols) + (confidence_threshold,)
        results = DatabaseManager.execute_query(query, params)

        # Analyze structural patterns
        structure_counts = Counter()
        confidence_by_structure = defaultdict(list)
        disorder_distribution = []
        sasa_distribution = []

        for row in results:
            site = dict(row)

            # Secondary structure analysis
            struct = site['secondary_structure'] or 'Unknown'
            structure_counts[struct] += 1
            confidence_by_structure[struct].append(site['predicted_prob_calibrated'] or 0)

            # Feature distributions
            if site['disorder_score'] is not None:
                disorder_distribution.append(site['disorder_score'])
            if site['sasa_ratio'] is not None:
                sasa_distribution.append(site['sasa_ratio'])

        # Calculate statistics
        structure_analysis = {}
        for struct, count in structure_counts.items():
            confidences = confidence_by_structure[struct]
            structure_analysis[struct] = {
                'count': count,
                'percentage': round((count / len(results)) * 100, 1),
                'avg_confidence': round(sum(confidences) / len(confidences), 3),
                'max_confidence': round(max(confidences), 3)
            }

        return jsonify({
            'total_sites': len(results),
            'structure_analysis': structure_analysis,
            'disorder_stats': {
                'mean': round(sum(disorder_distribution) / len(disorder_distribution),
                              3) if disorder_distribution else 0,
                'min': round(min(disorder_distribution), 3) if disorder_distribution else 0,
                'max': round(max(disorder_distribution), 3) if disorder_distribution else 0
            },
            'sasa_stats': {
                'mean': round(sum(sasa_distribution) / len(sasa_distribution), 3) if sasa_distribution else 0,
                'min': round(min(sasa_distribution), 3) if sasa_distribution else 0,
                'max': round(max(sasa_distribution), 3) if sasa_distribution else 0
            },
            'parameters': {
                'genes': gene_symbols,
                'confidence_threshold': confidence_threshold
            },
            'analysis_id': str(uuid.uuid4())
        })

    except Exception as e:
        current_app.logger.error(f"Error in structural analysis: {e}")
        return jsonify({'error': 'Failed to perform structural analysis'}), 500


def find_common_pathways(gene_symbols: List[str]) -> List[Dict]:
    """Find pathways common to multiple proteins."""
    if len(gene_symbols) < 2:
        return []

    placeholders = ','.join(['?' for _ in gene_symbols])
    query = f"""
    SELECT gs.gs_id, gs.gs_name, gs.gs_collection, gs.gs_description,
           COUNT(DISTINCT gsm.gene_symbol) as gene_overlap
    FROM gene_sets gs
    JOIN gene_set_memberships gsm ON gs.gs_id = gsm.gs_id
    WHERE gsm.gene_symbol IN ({placeholders})
    GROUP BY gs.gs_id
    HAVING gene_overlap >= 2
    ORDER BY gene_overlap DESC, gs.gs_name
    LIMIT 20
    """

    results = DatabaseManager.execute_query(query, tuple(gene_symbols))
    return [dict(row) for row in results]


def calculate_protein_similarity(gene_symbols: List[str], comparison_data: Dict) -> Dict:
    """Calculate similarity matrix between proteins."""
    similarity_matrix = {}

    for i, gene1 in enumerate(gene_symbols):
        similarity_matrix[gene1] = {}
        for j, gene2 in enumerate(gene_symbols):
            if i == j:
                similarity_matrix[gene1][gene2] = 1.0
            elif gene2 in similarity_matrix and gene1 in similarity_matrix[gene2]:
                similarity_matrix[gene1][gene2] = similarity_matrix[gene2][gene1]
            else:
                # Calculate similarity based on kinase profiles
                kinases1 = {k['kinase']: k['avg_score'] for k in comparison_data.get(gene1, {}).get('top_kinases', [])}
                kinases2 = {k['kinase']: k['avg_score'] for k in comparison_data.get(gene2, {}).get('top_kinases', [])}

                # Jaccard similarity coefficient for kinases
                all_kinases = set(kinases1.keys()) | set(kinases2.keys())
                if all_kinases:
                    intersection = sum(min(kinases1.get(k, 0), kinases2.get(k, 0)) for k in all_kinases)
                    union = sum(max(kinases1.get(k, 0), kinases2.get(k, 0)) for k in all_kinases)
                    similarity = intersection / union if union > 0 else 0
                else:
                    similarity = 0

                similarity_matrix[gene1][gene2] = round(similarity, 3)

    return similarity_matrix