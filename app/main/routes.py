"""Updated main web interface routes with protein detail view."""

from flask import render_template, request, current_app
from app.main import bp


@bp.route('/')
def index():
    """Homepage with dashboard overview."""
    return render_template('index.html')


@bp.route('/pathways')
def pathway_explorer():
    """Pathway-centric exploration interface."""
    return render_template('pathways.html')


@bp.route('/proteins')
def protein_explorer():
    """Protein-centric exploration interface."""
    return render_template('proteins.html')


@bp.route('/protein/<gene_symbol>')
def protein_detail(gene_symbol):
    """Detailed protein analysis view."""
    # Render the detailed protein template (site.html)
    return render_template('site.html', gene_symbol=gene_symbol.upper())


@bp.route('/pathway/<int:pathway_id>/analysis')
def pathway_analysis(pathway_id):
    """Comprehensive pathway phosphoproteomics analysis."""
    return render_template('pathway_analysis.html', pathway_id=pathway_id)


@bp.route('/kinases')
def kinase_explorer():
    """Kinase-centric exploration interface."""
    # Define kinase families for template
    kinase_families = {
        'AGC': ['AKT1', 'AKT2', 'AKT3', 'PKA', 'PKACA', 'PKACB', 'PKACG', 'PKG1', 'PKG2',
                'PKN1', 'PKN2', 'PKN3', 'ROCK1', 'ROCK2', 'P70S6K', 'P70S6KB', 'RSK2', 'RSK3', 'RSK4', 'P90RSK', 'SGK1',
                'SGK3'],
        'CAMK': ['CAMK1A', 'CAMK1B', 'CAMK1D', 'CAMK1G', 'CAMK2A', 'CAMK2B', 'CAMK2D', 'CAMK2G',
                 'CAMK4', 'CAMKK1', 'CAMKK2', 'CAMLCK', 'DAPK1', 'DAPK2', 'DAPK3', 'DCAMKL1', 'DCAMKL2',
                 'MARK1', 'MARK2', 'MARK3', 'MARK4'],
        'CK1': ['CK1A', 'CK1A2', 'CK1D', 'CK1E', 'CK1G1', 'CK1G2', 'CK1G3'],
        'CMGC': ['CDK1', 'CDK2', 'CDK3', 'CDK4', 'CDK5', 'CDK6', 'CDK7', 'CDK8', 'CDK9', 'CDK10',
                 'CDK12', 'CDK13', 'CDK14', 'CDK16', 'CDK17', 'CDK18', 'CDK19', 'CDKL1', 'CDKL5',
                 'CLK1', 'CLK2', 'CLK3', 'CLK4', 'DYRK1A', 'DYRK1B', 'DYRK2', 'DYRK3', 'DYRK4',
                 'ERK1', 'ERK2', 'ERK5', 'ERK7', 'GSK3A', 'GSK3B', 'JNK1', 'JNK2', 'JNK3',
                 'MAPKAPK2', 'MAPKAPK3', 'MAPKAPK5', 'P38A', 'P38B', 'P38D', 'P38G'],
        'STE': ['MEK1', 'MEK2', 'MEK5', 'MEKK1', 'MEKK2', 'MEKK3', 'MEKK6', 'PAK1', 'PAK2', 'PAK3',
                'PAK4', 'PAK5', 'PAK6', 'ASK1', 'TAK1', 'TAO1', 'TAO2', 'TAO3'],
        'TK': ['ALK2', 'ALK4'],
        'TKL': ['ACVR2A', 'ACVR2B', 'BMPR1A', 'BMPR1B', 'BMPR2', 'TGFBR1', 'TGFBR2', 'BRAF', 'RAF1',
                'MLK1', 'MLK2', 'MLK3', 'MLK4', 'LRRK2', 'RIPK1', 'RIPK2', 'RIPK3']
    }

    return render_template('kinases.html', kinase_families=kinase_families)


@bp.route('/discover')
def discovery_lab():
    """Discovery lab for advanced analysis."""
    return render_template('discover.html')


@bp.route('/docs')
def documentation():
    """API documentation and help."""
    return render_template('docs.html')