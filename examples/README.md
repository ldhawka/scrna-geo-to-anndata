# Examples

This directory contains example datasets and Jupyter notebooks demonstrating how to use the scRNA GEO to AnnData compiler with different data formats.

## Example Datasets

### 1. Simple Format - ALS PBMC Dataset
- **Location**: `data/simple_format_ALS/`
- **Source**: GSE244263 (subset of 3 samples)
- **Format**: Traditional GEO format with one expression file per sample
- **Files**: `GSM*.txt.gz` + `metadata.csv`
- **Samples**: 2 ALS patients, 1 control

### 2. With Cell Metadata - Stroke Dataset  
- **Location**: `data/with_cell_metadata_stroke/`
- **Source**: GSE225948 (subset of 3 samples)
- **Format**: Modern GEO format with paired count and cell metadata files
- **Files**: `*_counts.csv.gz` + `*_metadata.csv.gz` + `metadata.csv`
- **Samples**: 2 brain samples, 1 blood sample

## Jupyter Notebooks

### `01_simple_format_ALS_example.ipynb`
Demonstrates processing traditional GEO datasets:
- Single expression file per sample
- Sample-level metadata integration
- Basic scanpy preprocessing pipeline
- Visualization and quality control

### `02_stroke_with_cell_metadata.ipynb`
Demonstrates processing modern GEO datasets:
- Auto-detection of paired count/metadata files
- Integration of both sample-level and cell-level metadata
- Preservation of cell type annotations from GEO
- Selective metadata column inclusion

## Running the Examples

1. **Install requirements**:
   ```bash
   pip install -r requirements.txt
   ```

2. **Start Jupyter**:
   ```bash
   cd examples/notebooks
   jupyter notebook
   ```

3. **Run the notebooks** in order to see both data format examples.

## Data Sources

- **GSE244263**: ALS vs control peripheral blood mononuclear cells
- **GSE225948**: Stroke brain and blood samples with cell type annotations

Both datasets are subsets of the original GEO submissions, containing only a few samples for demonstration purposes.