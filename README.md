# scRNA GEO to AnnData

A Python tool for taking single-cell RNA-seq datasets from GEO all the way to
annotated, analysis-ready AnnData objects.

## Overview

Processing single-cell RNA-seq data from GEO involves a lot of tedious,
error-prone, and rarely-written-down work: downloading raw files, matching
metadata, handling three incompatible file layouts, filtering, detecting
doublets, correcting the batch effect between samples, and — the part that
takes longest and is easiest to get wrong — deciding what each cluster is.

This package covers that whole path, with the defaults and the judgement calls
written down rather than left implicit.

**Key features:**
- **Three data formats** — simple CSV/TSV, paired count + cell metadata, 10x Cell Ranger MTX
- **Automatic format detection** from file patterns
- **Per-sample doublet detection** (Scrublet), flagged rather than silently dropped
- **QC filtering** — min/max genes per cell, mitochondrial %, min cells per gene
- **Raw counts preserved** in `adata.layers['counts']`
- **Harmony batch integration** — multi-sample GEO compilations almost always need it
- **Cluster annotation by triangulation** — four CellTypist runs, marker z-scores,
  unbiased DE, assembled into one table for a human to decide from
- **Memory-efficient** handling of large multi-sample datasets

## Quick Start

### Installation

```bash
git clone https://github.com/ldhawka/scrna-geo-to-anndata.git
cd scrna-geo-to-anndata
pip install -r requirements.txt
```

Or install in development mode:
```bash
pip install -e .
```

### Compile a dataset

```python
from anndata_compiler import GEOAnndataCompiler

config = {
    'raw_data_dir': './GSE123456_RAW',
    'metadata_file': './metadata.csv',
    'output_file': './compiled_data.h5ad',
    'sample_id_column': 'sample_id',
}

adata = GEOAnndataCompiler(config).run_full_pipeline()
```

### Integrate and cluster

```python
from anndata_compiler import scale_and_pca, harmony, cluster

adata = scale_and_pca(adata, n_comps=50)
adata = harmony(adata, batch_key='sample_id', n_pcs=30)
adata = cluster(adata, use_rep='X_pca_harmony', n_pcs=30, resolutions=(1.0,))
```

### Annotate

See the [Annotation Guide](docs/annotation_guide.md) — this is the step that
deserves the most care, and it ends with a human looking at a dotplot.

### Or run it as scripts

Every stage is also a numbered command-line step with disk checkpoints between
them. See [`pipeline/README.md`](pipeline/README.md).

## Pipeline Steps

| Step | What it does |
|---|---|
| **1a** | Load samples (CSV/TSV, paired metadata, or 10x MTX); match metadata; **Scrublet per sample**; merge; QC metrics; QC filtering; preserve raw counts; CP10K + log1p; flag HVGs |
| **1b** | Subset to HVGs, `scale(max_value=10)`, PCA (50 comps, randomized) |
| **1c** | **Harmony** batch integration → `obsm['X_pca_harmony']` |
| **1d** | Neighbours (15, 30 PCs), UMAP, Leiden (r=1.0) |
| **1e** | CellTypist: `Immune_All_Low` + `Immune_All_High`, on HVG **and** full gene set |
| **1g** | Marker panel z-scores, raw means, unbiased Wilcoxon DE, dotplot |
| **1h** | Assemble all evidence into an annotation template CSV, label columns blank |
| **1i** | Apply the filled-in labels; flag doublet/artifact clusters |

## Configuration Options

| Parameter | Description | Default |
|-----------|-------------|---------|
| `raw_data_dir` | Path to GEO RAW directory | Required |
| `metadata_file` | CSV/Excel file with sample metadata | Required |
| `output_file` | Output h5ad file path | Required |
| `sample_id_column` | Column name for sample IDs in metadata | Required |
| `data_format` | `'auto'`, `'simple'`, `'with_cell_metadata'`, `'10x_mtx'` | `'auto'` |
| `min_genes` | Min genes per cell | `500` |
| `max_genes` | Max genes per cell (`None` = no cap) | `5000` |
| `min_cells` | Min cells per gene | `3` |
| `max_mito_pct` | Max mitochondrial % | `15.0` |
| `mito_prefix` | Mitochondrial gene prefix (`'MT-'` human, `'mt-'` mouse) | `'MT-'` |
| `detect_doublets` | Run Scrublet per sample | `True` |
| `expected_doublet_rate` | Scrublet expected doublet rate | `0.06` |
| `filter_doublets` | Drop predicted doublets instead of flagging them | `False` |
| `target_sum` | Normalization target sum | `10000` |
| `n_top_genes` | Number of highly variable genes | `3000` |
| `max_cells_per_sample` | Maximum cells per sample (`None` = keep all) | `None` |
| `delimiter` | File delimiter for CSV/TSV formats | `'whitespace'` |
| `random_state` | Random seed | `42` |
| `counts_pattern` | Glob pattern for count files | `None` (auto) |
| `cell_metadata_pattern` | Glob pattern for cell metadata files | `None` (auto) |
| `metadata_columns` | List of metadata columns to include (`None` = all) | `None` |

### About these defaults

They come from PBMC 10x data (the Itou and Zhang/Gate cohorts) and are a
reasonable starting point for that. **Re-examine them for anything else** —
tissue data, single-nucleus data, and non-immune compartments all want
different thresholds. `max_mito_pct=15` in particular is too strict for many
tissues and far too loose for single-nucleus data.

The QC thresholds are deliberately more conservative than the scanpy tutorial
defaults (`min_genes=200`, no upper bound, `max_mito=20`). The `max_genes`
upper bound is there to catch doublets that Scrublet misses.

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

After step 1a:
- `adata.X` — Normalized, log-transformed expression
- `adata.layers['counts']` — Raw integer counts
- `adata.var['highly_variable']` — HVG boolean flag
- `adata.var['mt']` — Mitochondrial gene flag
- `adata.obs['doublet_score']`, `adata.obs['predicted_doublet']` — Scrublet, per sample
- `adata.obs` — All metadata + QC metrics (`n_genes_by_counts`, `total_counts`, `pct_counts_mt`)

After steps 1b–1d, additionally:
- `adata.raw` — Pre-scale log-normalized matrix, all genes (needed by 1e and 1g)
- `adata.obsm['X_pca']`, `adata.obsm['X_pca_harmony']`, `adata.obsm['X_umap']`
- `adata.obs['leiden_r1_0']` (one column per requested resolution)

After step 1i:
- `adata.obs['cell_type']`, `adata.obs['cell_type_flag']`

## Documentation

- **[Downstream Analysis Guide](docs/downstream_analysis_guide.md)** — QC review, PCA and
  the scree plot, Harmony, UMAP, Leiden, and the mistakes that fail silently
- **[Annotation Guide](docs/annotation_guide.md)** — the triangulation method, PBMC marker
  tables, the z-score trap, doublet and artifact clusters, annotation granularity across
  cohorts, and the compositional caveat on cell-type frequencies
- **[Pipeline scripts](pipeline/README.md)** — the same thing as numbered CLI steps

### Related

Microglia and CNS-myeloid markers are maintained separately as a hierarchical atlas —
9 umbrella families over 23 sub-states, human and mouse symbols, denoised gene sets,
DOIs: **https://ldhawka.github.io/microglia-subtype-markers/**. The pipeline can load
it directly (`--panel microglia`).

## Requirements

- Python >= 3.9
- anndata >= 0.11.0
- scanpy >= 1.11.0
- pandas >= 2.3.0
- numpy >= 1.23.0
- scipy >= 1.10.0
- scikit-image >= 0.19.0 (Scrublet)
- harmonypy >= 0.0.9 (batch integration)
- celltypist >= 1.6.0 (annotation)
- leidenalg >= 0.10.0, igraph >= 0.11.0 (clustering)
- matplotlib >= 3.10.0, seaborn >= 0.12.0
- tqdm >= 4.66.0
- openpyxl >= 3.0.9

## License

MIT License - see LICENSE file for details.

## Contributing

Issues, pull requests and feedback welcome!
