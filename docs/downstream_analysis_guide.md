# Downstream Analysis Guide

After compiling your data with `GEOAnndataCompiler` (step 1a), the output h5ad contains:

- **`adata.X`** — Normalized, log-transformed expression (`log1p(CP10K)`)
- **`adata.layers['counts']`** — Raw integer counts, preserved for methods that need them
- **`adata.var['highly_variable']`** — Boolean flag for highly variable genes
- **`adata.var['mt']`** — Boolean flag for mitochondrial genes
- **`adata.obs['doublet_score']`, `adata.obs['predicted_doublet']`** — Scrublet, run per sample
- **`adata.obs`** — Sample and cell metadata plus QC metrics (`n_genes_by_counts`,
  `total_counts`, `pct_counts_mt`)

This guide covers the interactive path from there: QC review, PCA, **Harmony batch
integration**, UMAP, and Leiden clustering. To run the same thing as scripts instead,
see [`pipeline/README.md`](../pipeline/README.md). To annotate the resulting clusters,
see [`annotation_guide.md`](annotation_guide.md).

---

## 1. Load and Inspect Compiled Data

```python
import scanpy as sc

adata = sc.read_h5ad('compiled_data.h5ad')
print(adata)
print(f"Layers: {list(adata.layers.keys())}")
print(f"Samples: {adata.obs['sample_id'].nunique()}")
```

## 2. QC Visualization

Verify the QC filtering was appropriate before proceeding:

```python
# Overall QC distributions (post-filtering)
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)

# Per-sample QC — check for batch effects or outlier samples
sc.pl.violin(adata, 'n_genes_by_counts', groupby='sample_id', rotation=45)
sc.pl.violin(adata, 'pct_counts_mt', groupby='sample_id', rotation=45)

# Scatter: genes vs counts, colored by mito %
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', color='pct_counts_mt')
```

**Look for a sample whose distributions sit well away from the rest.** A sample with
half the median genes/cell of the others will drive its own clusters no matter what
Harmony does, and is usually better excluded than corrected.

### Reviewing doublet calls

```python
sc.pl.violin(adata, 'doublet_score', groupby='sample_id', rotation=45)
print(adata.obs.groupby('sample_id', observed=True)['predicted_doublet'].mean())
```

The compiler **flags** doublets rather than dropping them, on purpose. Scrublet's
automatic threshold is unreliable on small or low-complexity samples and can flag an
implausible fraction of cells — if a sample comes back at 40% doublets, do not believe
it. Doublets that form their own *cluster* are also much easier to spot after
clustering, from co-expression of two lineages (see the annotation guide).

To drop them at compile time instead, pass `filter_doublets=True` (or
`--filter-doublets`).

## 3. Scale and PCA

```python
from anndata_compiler import scale_and_pca

# Stores the pre-scale log-normalized matrix in .raw, subsets to HVGs,
# scales with clipping, then runs PCA.
adata = scale_and_pca(adata, n_comps=50, max_value=10)

# Scree plot — look for the "elbow" where variance explained flattens
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)
```

Scaling to unit variance with `max_value=10` clipping stops a handful of very highly
expressed genes from dominating the components. The clip matters more than it looks:
without it a single outlier cell can define a PC on its own.

`.raw` is set before scaling so that CellTypist and the marker z-scores can still see
log-normalized expression for **all** genes later. Do not discard it.

**Choosing n_pcs:**
- Read the elbow from the scree plot
- Typical values: 20–40 for most scRNA-seq datasets; 30 is a reasonable default
- Too many adds noise; too few loses biological signal
- Whatever you pick, use the **same** number for Harmony and for `neighbors`

```python
# Inspect which metadata variables drive the top PCs
sc.pl.pca(adata, color=['sample_id', 'n_genes_by_counts'], ncols=2)
```

If `sample_id` visibly structures PC1/PC2, that is the batch effect you are about to
correct.

## 4. Harmony Batch Integration

Multi-sample GEO compilations almost always carry a per-sample technical effect —
different capture dates, chemistry versions, sequencing depths. **Without correction,
Leiden will happily return clusters that are one sample each.**

```python
from anndata_compiler import harmony

adata = harmony(adata, batch_key='sample_id', n_pcs=30)
# -> adata.obsm['X_pca_harmony']
```

**Choosing `batch_key`:**
- Normally the sample/library ID. This is the safe default.
- A coarser key (sequencing batch) only when samples within a batch are known to be
  technically identical.
- **Never** a key confounded with your biological contrast. Correcting on
  `disease_state` removes the signal you are looking for.

**Verify it worked.** Compare UMAPs before and after, coloured by sample. Harmony can
also over-merge genuinely distinct biology, so check that populations you expect to be
distinct have not collapsed.

## 5. Neighbors, UMAP, and Leiden

```python
from anndata_compiler import cluster

adata = cluster(
    adata,
    use_rep='X_pca_harmony',   # <- the corrected embedding
    n_neighbors=15,
    n_pcs=30,                  # same as Harmony
    resolutions=(0.5, 1.0, 1.5),
)
```

> **The most common silent mistake in this whole pipeline** is running `sc.pp.neighbors`
> without `use_rep='X_pca_harmony'` after running Harmony. It defaults to `X_pca`, the
> clusters come back uncorrected, the UMAP separates by sample, and nothing errors.

**Parameter guidance:**
- **`n_pcs`** — from the scree plot; match what Harmony used
- **`n_neighbors`** — controls local vs global structure
  - Lower (5–10): preserves fine local structure, may fragment clusters
  - Default (15): good balance for most datasets
  - Higher (20–30): smoother global structure, may merge small populations

**Choosing resolution:**
- **Lower (0.1–0.5)**: fewer, broader clusters — major cell types
- **Medium (0.5–1.0)**: balanced — typical starting point
- **Higher (1.0–2.0)**: finer clusters — subtypes and states

**`1.0` is the default here** for PBMC data at this granularity, and is the resolution
these defaults were validated at. Computing several resolutions is cheap and worth
keeping, but pick **one** for the annotation pass rather than annotating each — the
annotation effort scales with the number of clusters, and mixing resolutions across
cohorts makes them incomparable.

```python
sc.pl.umap(adata, color=['leiden_r1_0', 'sample_id'], ncols=2)

for res in [0.5, 1.0, 1.5]:
    key = f"leiden_r{str(res).replace('.', '_')}"
    print(f"Resolution {res}: {adata.obs[key].nunique()} clusters")
```

## 6. Evaluate Clustering

```python
# Cluster sizes
print(adata.obs['leiden_r1_0'].value_counts())

# Known markers
marker_genes = ['CD3D', 'CD14', 'MS4A1', 'FCGR3A', 'PPBP', 'LYZ']
sc.pl.dotplot(adata, marker_genes, groupby='leiden_r1_0')

# Unbiased marker genes per cluster
sc.tl.rank_genes_groups(adata, groupby='leiden_r1_0', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=10)
```

Sanity checks worth doing every time:
- Is any cluster dominated by a single sample? (batch effect survived, or a real
  sample-specific population — worth distinguishing)
- Does any cluster have markedly lower `n_genes_by_counts` than the rest? (possible
  technical split — see the annotation guide)
- Are top DE genes mostly ribosomal or mitochondrial? (possible artifact cluster)

## 7. Annotate

This is a large enough topic to have its own document. See
**[annotation_guide.md](annotation_guide.md)** for the triangulation method, the
marker tables, the interpretation traps, and the microglia atlas.

## 8. Save

```python
adata.write('analyzed_data.h5ad')
```

## Using Raw Counts

Some downstream methods require raw counts (differential expression, RNA velocity, etc.):

```python
# Access raw counts directly
raw_counts = adata.layers['counts']

# Example: create a raw-counts AnnData for pseudobulk DE
adata_raw = adata.copy()
adata_raw.X = adata_raw.layers['counts']
```

For differential expression across conditions, prefer **pseudobulk** (aggregate counts
per sample per cell type, then run a bulk method such as PyDESeq2) over per-cell tests.
Per-cell tests treat cells as independent replicates, which they are not — they are
pseudo-replicates within a donor, and the resulting p-values are anticonservative by
orders of magnitude.

If you go on to correlate **cell-type frequencies**, read §8 of the annotation guide
first — frequencies are compositional, and closure manufactures correlations that are
not biology.
