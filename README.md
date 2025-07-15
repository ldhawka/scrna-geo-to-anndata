# scRNA GEO to AnnData Compiler

A Python tool for automated compilation and preprocessing of single-cell RNA-seq datasets downloaded from GEO databases into standardized AnnData objects.

## Overview

Processing single-cell RNA-seq data from GEO databases typically involves tedious manual steps: downloading raw files, matching metadata, handling different file formats, and running standard preprocessing pipelines. This tool has helped me automate the workflow, from raw GEO data to analysis-ready AnnData objects.

**Key Features:**
- **Automated GEO data processing** with flexible file format support
- **Metadata integration** from CSV/Excel files  
- **Memory-efficient handling** of large multi-sample datasets
- **Parameter optimization** for PCA and clustering analysis
- **Standardized scanpy workflows** with reproducible results
- **Configurable preprocessing** pipeline with sensible defaults

## Quick Start

### Installation

```bash
git clone https://github.com/ldhawka/scrna-geo-to-anndata.git
cd scrna-geo-to-anndata
pip install -r requirements.txt
```

### Basic Usage

```python
from anndata_compiler import GEOAnndataCompiler, create_config_template

# Create configuration
config = {
    'raw_data_dir': './GSE123456_RAW',           # Path to GEO RAW directory
    'metadata_file': './metadata.csv',           # Sample metadata file
    'output_file': './compiled_data.h5ad',       # Output AnnData file
    'sample_id_column': 'sample_id',             # Column for sample IDs
    'max_cells_per_sample': 750,                 # Subsample large samples
    'target_sum': 1e4,                           # Normalization target
    'n_top_genes': 3000,                         # Highly variable genes
    'optimize_params': True                      # Auto-optimize parameters
}

# Initialize and run
compiler = GEOAnndataCompiler(config)
adata = compiler.run_full_pipeline(plot_colors=['leiden', 'sample_id'])
```

## Metadata Format

The tool expects metadata in CSV or Excel format with sample information. Here's an example structure:

| Column | Description | Example |
|--------|-------------|---------|
| `Sample_name` | Full sample name from GEO | "GSM123_Sample1" |
| `Sample_geo_accession` | GEO accession number | "GSM123456" |
| `Source` | Sample source/tissue | "brain" |
| `Organism` | Species | "Homo sapiens" |
| `Disease_state` | Condition/treatment | "control", "disease" |
| `Sample_ID` | Unique sample identifier | "Sample1" |

**Note:** The `Sample_ID` column should match the sample IDs extracted from your filename format. You can customize the `sample_id_column` parameter to use a different column name.

**Getting metadata from GEO:**
1. Download sample metadata from your GEO series page
2. Convert to CSV format with the required columns
3. Ensure `Sample_ID` matches your filename pattern

## Configuration Options

| Parameter | Description | Default |
|-----------|-------------|---------|
| `raw_data_dir` | Path to GEO RAW directory | Required |
| `metadata_file` | CSV/Excel file with sample metadata | Required |
| `output_file` | Output h5ad file path | Required |
| `sample_id_column` | Column name for sample IDs in metadata | Required |
| `max_cells_per_sample` | Maximum cells to sample per file | 750 |
| `target_sum` | Normalization target sum | 10,000 |
| `n_top_genes` | Number of highly variable genes | 3,000 |
| `delimiter` | File delimiter ('whitespace' or character) | 'whitespace' |
| `random_state` | Random seed for reproducibility | 42 |
| `optimize_params` | Auto-optimize PCA/clustering parameters | False |

## Pipeline Steps

1. **Sample Processing**: Load and parse individual sample files
2. **Metadata Integration**: Match samples with metadata and annotations  
3. **Data Merging**: Combine samples into unified AnnData object
4. **Quality Control**: Calculate standard QC metrics
5. **Normalization**: Library size normalization and log transformation
6. **Feature Selection**: Identify highly variable genes
7. **Dimensionality Reduction**: PCA with optional parameter optimization
8. **Clustering**: Leiden clustering and UMAP visualization
9. **Output**: Save processed data and generate plots

## Advanced Features

### Custom Sample ID Extraction

Override the `extract_sample_id` method for custom filename parsing:

```python
class CustomCompiler(GEOAnndataCompiler):
    def extract_sample_id(self, filename):
        # Custom logic for your filename format
        return filename.split('_')[1]  # Example: use second part
```

### Parameter Optimization

When `optimize_params=True`, the tool automatically:
- Determines optimal number of PCs (targeting 80-90% variance explained)
- Optimizes neighbor count for clustering using elbow method
- Caps parameters to prevent overfitting with small datasets

### Memory Management

For large datasets, the tool includes:
- Automatic garbage collection between processing steps
- Data type optimization to reduce memory usage
- Progress tracking with informative logging

## Requirements

- Python ≥ 3.8
- pandas ≥ 2.3.0
- numpy ≥ 1.23.0  
- anndata ≥ 0.11.0
- scanpy ≥ 1.11.0
- scikit-learn ≥ 1.7.0
- matplotlib ≥ 3.10.0
- tqdm ≥ 4.66.0
- openpyxl ≥ 3.0.9

## Development

This tool was developed using Claude, combining domain expertise in computational biology with collaborative development to ensure robust, well-documented code following bioinformatics best practices.

## License

MIT License - see LICENSE file for details.

## Contributing

Issues and pull requests welcome! Feedback from the computational biology community is appreciated.
