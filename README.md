# Structural Context Pipeline

A comprehensive pipeline for extracting structural context features from protein structures (AlphaFold models or experimental structures). This tool generates structural atlases and disorder/boundary panels for phosphorylation site analysis and other structural bioinformatics applications.

## Features

- **Structural Context Atlas**: Extract comprehensive structural features around residues of interest
- **Disorder/Boundary Panel**: Analyze disorder metrics, accessibility, and biophysical properties
- **Flexible Input**: Support for CIF and PDB format files (compressed or uncompressed)
- **Parallel Processing**: Multi-core support for efficient large-scale analysis
- **Merged Output**: Automatically combines features into a single analysis-ready CSV

## Installation

### Prerequisites

- Python 3.8 or higher
- pip package manager

### Quick Install

```bash
# Clone the repository
git clone https://github.com/yourusername/structural-context-pipeline.git
cd structural-context-pipeline

# Install with pip
pip install -e .
```

### Install from requirements.txt

```bash
# Install core dependencies
pip install -r requirements.txt

# For development (optional)
pip install -r requirements-dev.txt
```

## Repository Structure

```
structural-context-pipeline/
├── main/
│   ├── build_structure_context_library.py       # Structural atlas builder
│   ├── build_structure_context_library_update.py # Disorder/boundary panel
│   └── run_structural_context_pipeline.py       # Main pipeline script
├── Structure_Depot/                             # Place your structure files here
│   └── (your .cif/.pdb/.gz files)
├── setup.py                                      # Package setup
├── requirements.txt                              # Core dependencies
├── requirements-dev.txt                          # Development dependencies
└── README.md                                     # This file
```

## Usage

### Basic Usage

```bash
# Run the full pipeline (atlas + disorder panel + merge)
python main/run_structural_context_pipeline.py Structure_Depot -o results.csv

# Or if installed as package:
structural-context-pipeline Structure_Depot -o results.csv
```

### Common Options

```bash
# Specify file types and use PDB as fallback
python main/run_structural_context_pipeline.py Structure_Depot \
    -o combined_features.csv \
    -t cif \
    --fallback-pdb

# Use multiple CPU cores and custom batch size
python main/run_structural_context_pipeline.py Structure_Depot \
    -o results.csv \
    -p 8 \
    -b 16

# Process only a subset of files for testing
python main/run_structural_context_pipeline.py Structure_Depot \
    -o test_results.csv \
    -m 100

# Keep intermediate files for debugging
python main/run_structural_context_pipeline.py Structure_Depot \
    -o results.csv \
    --keep-intermediate
```

### Advanced Usage

```bash
# Run only the structural atlas
python main/run_structural_context_pipeline.py Structure_Depot \
    --only-atlas \
    -o atlas_only.csv

# Run only the disorder panel
python main/run_structural_context_pipeline.py Structure_Depot \
    --only-panel \
    -o panel_only.csv

# Custom window size and temp directory
python main/run_structural_context_pipeline.py Structure_Depot \
    -o results.csv \
    -w 10 \
    -d ./my_temp_results
```

## Command-Line Arguments

| Argument | Description | Default |
|----------|-------------|---------|
| `structure_dir` | Directory containing structure files (.cif/.pdb/.gz) | Required |
| `-o, --output` | Output CSV filename | `combined_structural_features.csv` |
| `-d, --temp-dir` | Directory for intermediate outputs | `./temp_results` |
| `-t, --file-types` | Comma-separated file types (cif,pdb) | `cif` |
| `-f, --fallback-pdb` | Include PDB when CIF not present | False |
| `-p, --num-processes` | Worker processes (0=auto) | 0 |
| `-b, --batch-size` | Batch size for file submissions | 16 |
| `-w, --window` | Half-window size (±window residues) | 7 |
| `-m, --max-files` | Max files to process (0=all) | 0 |
| `-q, --quiet` | Suppress per-file output | False |
| `--only-atlas` | Run only atlas builder | False |
| `--only-panel` | Run only disorder panel | False |
| `--keep-intermediate` | Keep intermediate CSVs | False |

## Output Format

The pipeline generates a CSV file with the following key columns:

- **UniProtID**: UniProt accession identifier
- **ChainID**: Protein chain identifier
- **ResidueNumber**: Residue position
- **ResidueType**: Single-letter amino acid code
- **Site**: Specific residue designation
- **[Structural Features]**: Various context-specific features from the atlas
- **[Disorder Features]**: Disorder metrics and boundary information from the panel

## Input Structure Files

Place your structure files in the `Structure_Depot/` directory (or specify a custom directory). Supported formats:

- CIF files (`.cif`)
- PDB files (`.pdb`)
- Compressed files (`.gz`)

The pipeline can process AlphaFold predicted structures or experimental structures from the PDB.

## Performance Tips

1. **Parallel Processing**: Use `-p` to specify CPU cores (e.g., `-p 8` for 8 cores)
2. **Batch Size**: Adjust `-b` based on available memory (larger = faster but more RAM)
3. **Test First**: Use `-m 10` to process only 10 files for testing
4. **File Types**: If you only have CIF files, use `-t cif` to skip PDB scanning

## Troubleshooting

### Missing Dependencies
```bash
pip install -r requirements.txt
```

### Memory Issues
- Reduce batch size: `-b 8`
- Reduce parallel processes: `-p 4`
- Process files in chunks: `-m 1000`

### Missing Structure Files
Ensure your structure files are in the correct directory and have proper naming conventions.

## Citation

If you use this pipeline in your research, please cite:

[Your paper citation here]

## License

[Specify your license here - e.g., MIT, GPL, etc.]

## Contact

For questions, issues, or contributions:
- GitHub Issues: https://github.com/yourusername/structural-context-pipeline/issues
- Email: your.email@example.com

## Acknowledgments

This pipeline was developed for phosphorylation site analysis and kinase-substrate prediction research using AlphaFold structural models.