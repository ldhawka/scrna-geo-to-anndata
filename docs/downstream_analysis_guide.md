# Downstream Analysis Guide

After compiling your data with `GEOAnndataCompiler`, the output h5ad contains:

- **`adata.X`** — Normalized, log-transformed expression matrix
- **`adata.layers['counts']`** — Raw integer counts (preserved for methods that need them)
- **`adata.var['highly_variable']`** — Boolean flag for highly variable genes
- **`adata.var['mt']`** — Boolean flag for mitochondrial genes
- **`adata.obs`** — Sample and cell metadata, including QC metrics (`n_genes_by_counts`, `total_counts`, `pct_counts_mt`)

This guide walks through interactive downstream analysis: PCA, UMAP, and Leiden clustering.

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

## 3. PCA

Run PCA on highly variable genes and inspect the scree plot to choose the number of components:

```python
# Compute PCA (50 components as starting point)
sc.tl.pca(adata, svd_solver='arpack', n_comps=50, use_highly_variable=True)

# Scree plot — look for the "elbow" where variance explained flattens
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)
```

**How to choose n_comps:**
- Look for the elbow in the scree plot where additional PCs contribute little variance
- Typical values: 20-40 for most scRNA-seq datasets
- When in doubt, 30 is a reasonable default
- Using too many PCs adds noise; too few loses biological signal

```python
# Inspect which metadata variables drive the top PCs
sc.pl.pca(adata, color=['sample_id', 'n_genes_by_counts'], ncols=2)
```

## 4. Compute Neighbors

```python
n_pcs = 30  # Set based on your scree plot
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=n_pcs)
```

**Parameter guidance:**
- **`n_pcs`**: Set from the scree plot elbow (step 3)
- **`n_neighbors`**: Controls local vs global structure
  - Lower (5-10): Preserves fine local structure, may fragment clusters
  - Default (15): Good balance for most datasets
  - Higher (20-30): Smoother global structure, may merge small populations

## 5. UMAP

```python
sc.tl.umap(adata)

# Visualize by sample and metadata
sc.pl.umap(adata, color='sample_id')
```

## 6. Leiden Clustering at Multiple Resolutions

Rather than picking a single resolution, run several and compare:

```python
for res in [0.3, 0.5, 0.8, 1.0, 1.5]:
    sc.tl.leiden(adata, resolution=res, key_added=f'leiden_{res}')

# Compare resolutions side by side
sc.pl.umap(adata, color=['leiden_0.3', 'leiden_0.5', 'leiden_0.8', 'leiden_1.0'],
           ncols=2, wspace=0.3)

# Print cluster counts per resolution
for res in [0.3, 0.5, 0.8, 1.0, 1.5]:
    n_clusters = adata.obs[f'leiden_{res}'].nunique()
    print(f"Resolution {res}: {n_clusters} clusters")
```

**Choosing resolution:**
- **Lower (0.1-0.5)**: Fewer, broader clusters — good for major cell types
- **Medium (0.5-1.0)**: Balanced — typical starting point
- **Higher (1.0-2.0)**: Finer clusters — subtypes, states
- The "right" resolution depends on your biological question
- You can store multiple resolutions and use different ones for different analyses

## 7. Evaluate Clustering

```python
# Check cluster sizes
print(adata.obs['leiden_0.5'].value_counts())

# Visualize known marker genes
marker_genes = ['CD3D', 'CD14', 'MS4A1', 'FCGR3A']  # Adapt to your cell types
sc.pl.dotplot(adata, marker_genes, groupby='leiden_0.5')

# Find marker genes per cluster
sc.tl.rank_genes_groups(adata, groupby='leiden_0.5', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=10)
```

## 8. Save Analyzed Object

```python
adata.write('analyzed_data.h5ad')
```

## Using Raw Counts

Some downstream methods require raw counts (DEG testing, RNA velocity, etc.):

```python
# Access raw counts directly
raw_counts = adata.layers['counts']

# Example: create a raw-counts AnnData for DEG testing
adata_raw = adata.copy()
adata_raw.X = adata_raw.layers['counts']
```
