# scRNA GEO to AnnData Compiler

A Python tool for compiling single-cell RNA-seq datasets from GEO into standardized, analysis-ready AnnData objects.

## Overview

Processing single-cell RNA-seq data from GEO typically involves tedious manual steps: downloading raw files, matching metadata, handling different file formats, and running standard preprocessing. This tool automates the workflow from raw GEO data to a QC-filtered, normalized h5ad file ready for interactive analysis.

**Key Features:**
- **Three data formats**: Simple CSV/TSV, paired count + metadata files, and 10x Cell Ranger MTX
- **Automatic format detection** from file patterns
- **QC filtering**: Mitochondrial %, min genes/cell, min cells/gene with configurable thresholds
- **Raw counts preservation** in `adata.layers['counts']`
- **Metadata integration** from CSV/Excel files
- **Memory-efficient handling** of large multi-sample datasets

## Quick Start

### Installation

```bash
git clone https://github.com/ldhawka/scrna-geo-to-anndata.git
cd scrna-geo-to-anndata
pip install -r requirements.txt
```

### Basic Usage

```python
import sys
sys.path.append('/path/to/scrna-geo-to-anndata')
from anndata_compiler import GEOAnndataCompiler

config = {
    'raw_data_dir': './GSE123456_RAW',
    'metadata_file': './metadata.csv',
    'output_file': './compiled_data.h5ad',
    'sample_id_column': 'sample_id',
}

compiler = GEOAnndataCompiler(config)
adata = compiler.run_full_pipeline()
```

Or install in development mode:
```bash
pip install -e .
```

## Pipeline Steps

1. **Sample Processing**: Load and parse individual sample files (CSV/TSV, paired metadata, or 10x MTX)
2. **Metadata Integration**: Match samples with metadata annotations
3. **Data Merging**: Combine samples into unified AnnData object
4. **QC Metrics**: Calculate genes/cell, counts/cell, mitochondrial %
5. **QC Filtering**: Remove low-quality cells and rarely-expressed genes
6. **Raw Counts Preservation**: Store in `adata.layers['counts']`
7. **Normalization**: CPM (target_sum=10k) + log1p
8. **Feature Selection**: Flag highly variable genes (all genes retained)
9. **Output**: Save analysis-ready h5ad file

## Configuration Options

| Parameter | Description | Default |
|-----------|-------------|---------|
| `raw_data_dir` | Path to GEO RAW directory | Required |
| `metadata_file` | CSV/Excel file with sample metadata | Required |
| `output_file` | Output h5ad file path | Required |
| `sample_id_column` | Column name for sample IDs in metadata | Required |
| `data_format` | `'auto'`, `'simple'`, `'with_cell_metadata'`, `'10x_mtx'` | `'auto'` |
| `max_cells_per_sample` | Maximum cells per sample (`None` = keep all) | `None` |
| `target_sum` | Normalization target sum | `10000` |
| `n_top_genes` | Number of highly variable genes | `3000` |
| `min_genes` | Min genes per cell (QC filter) | `200` |
| `min_cells` | Min cells per gene (QC filter) | `3` |
| `max_mito_pct` | Max mitochondrial % (QC filter) | `20.0` |
| `mito_prefix` | Mitochondrial gene prefix (`'MT-'` human, `'mt-'` mouse) | `'MT-'` |
| `delimiter` | File delimiter for CSV/TSV formats | `'whitespace'` |
| `random_state` | Random seed | `42` |
| `counts_pattern` | Glob pattern for count files | `None` (auto) |
| `cell_metadata_pattern` | Glob pattern for cell metadata files | `None` (auto) |
| `metadata_columns` | List of metadata columns to include (`None` = all) | `None` |

## Data Format Examples

### Simple Format (Traditional GEO)

Single expression file per sample:
```python
config = {
    'raw_data_dir': './GSE123456_RAW',
    'metadata_file': './metadata.csv',
    'output_file': './compiled_data.h5ad',
    'sample_id_column': 'Sample_ID',
    'data_format': 'simple',
}
```

### With Cell Metadata

Paired count and cell metadata files per sample:
```python
config = {
    'raw_data_dir': './GSE225948_RAW',
    'metadata_file': './sample_metadata.csv',
    'output_file': './compiled_with_cells.h5ad',
    'sample_id_column': 'Sample_ID',
    'data_format': 'with_cell_metadata',
    'delimiter': ',',
}
```

### 10x Cell Ranger MTX Format

Cell Ranger output with matrix.mtx + barcodes.tsv + features.tsv:
```python
config = {
    'raw_data_dir': './GSE288856_RAW',   # Subdirectories or prefixed files
    'metadata_file': './metadata.csv',
    'output_file': './compiled_10x.h5ad',
    'sample_id_column': 'sample_id',
    'data_format': '10x_mtx',
}
```

Supports two common GEO layouts:
- **Subdirectory**: Each sample in its own folder (e.g., `sample1/matrix.mtx.gz`)
- **Prefixed files**: Flat directory with sample-prefixed files (e.g., `GSM123_sample1.matrix.mtx.gz`)

## Metadata Format

The tool expects metadata in CSV or Excel format. Example:

| Sample_ID | Disease_state | Tissue | Age | Sex |
|-----------|--------------|--------|-----|-----|
| Sample1 | control | brain | 65 | M |
| Sample2 | disease | brain | 72 | F |

The `sample_id_column` should match sample IDs extracted from filenames. You can customize extraction by subclassing:

```python
class CustomCompiler(GEOAnndataCompiler):
    def extract_sample_id(self, filename):
        return filename.split('_')[1]  # Custom logic
```

## Output Format

The compiled h5ad contains:
- `adata.X` — Normalized, log-transformed expression
- `adata.layers['counts']` — Raw integer counts
- `adata.var['highly_variable']` — HVG boolean flag
- `adata.var['mt']` — Mitochondrial gene flag
- `adata.obs` — All metadata + QC metrics (`n_genes_by_counts`, `total_counts`, `pct_counts_mt`)

## Downstream Analysis

The compiled h5ad is ready for interactive analysis (PCA, UMAP, Leiden clustering). See **[Downstream Analysis Guide](docs/downstream_analysis_guide.md)** for a step-by-step walkthrough covering:

- QC visualization
- PCA and scree plot inspection
- Choosing neighbors and PCs
- Leiden clustering at multiple resolutions
- Evaluating clusters with marker genes

## Requirements

- Python >= 3.9
- anndata >= 0.11.0
- scanpy >= 1.11.0
- pandas >= 2.3.0
- numpy >= 1.23.0
- scipy >= 1.10.0
- matplotlib >= 3.10.0
- tqdm >= 4.66.0
- leidenalg >= 0.10.0 (for downstream clustering)
- openpyxl >= 3.0.9

## License

MIT License - see LICENSE file for details.

## Contributing

Issues, pull requests and feedback welcome!
