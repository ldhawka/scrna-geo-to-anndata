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
Demonstrates compiling traditional GEO datasets:
- Single expression file per sample
- Sample-level metadata integration
- QC filtering (mito %, min/max genes, min cells)
- Per-sample Scrublet doublet detection
- Normalization with raw counts preservation
- Output: analysis-ready h5ad with optional downstream analysis preview

### `02_stroke_with_cell_metadata.ipynb`
Demonstrates compiling modern GEO datasets:
- Auto-detection of paired count/metadata files
- Integration of sample-level and cell-level metadata
- QC filtering and preprocessing
- Output: analysis-ready h5ad with cell annotations preserved

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

## A note on runtime

The two example datasets are real GEO data, not toy matrices — the simple-format ALS
files are ~77 MB gzipped dense text, and parsing them takes several minutes per sample.
That parse is the slow part, not the QC or the doublet detection. Be patient the first
time, and reuse the compiled h5ad afterwards.

## Downstream Analysis

After compilation:

- [docs/downstream_analysis_guide.md](../docs/downstream_analysis_guide.md) — PCA,
  Harmony batch integration, UMAP, and Leiden clustering
- [docs/annotation_guide.md](../docs/annotation_guide.md) — turning clusters into cell types
- [pipeline/README.md](../pipeline/README.md) — the same pipeline as numbered CLI steps

Note that the notebooks predate v0.3.0 and cover compilation only; they do not yet
show the Harmony or annotation steps.
