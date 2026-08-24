# Pipeline scripts

The same pipeline as the library API, as numbered command-line steps. Use these
when you want each stage checkpointed to disk and re-runnable in isolation —
which is what you want on a cluster, or any time a step takes long enough that
you do not want to redo it.

Run them in order. Every script takes `--help`.

```
1a  compile GEO files -> AnnData; Scrublet; QC; normalize; HVG
1b  subset to HVGs, scale, PCA          (writes .raw = pre-scale log-norm)
1c  Harmony batch integration           (writes obsm['X_pca_harmony'])
1d  neighbours, UMAP, Leiden
1e  CellTypist  x2 models x2 gene sets  (run twice: --mode hvg, --mode fullgene)
1g  marker z-scores, unbiased DE, dotplot
1h  build the annotation template CSV
--  you fill in manual_label_final      <- the human step
1i  apply labels back onto the cells
```

There is no `1f`. It is left free for the cohort-specific metadata work that
usually lands between clustering and annotation — parsing an author-supplied
`.rds`, reconciling GEO metadata against a clinical table, fixing a swapped
`disease_state` field.

## Worked example

```bash
REPO=/path/to/scrna-geo-to-anndata
cd $REPO/pipeline

# --- 1a: compile ---------------------------------------------------------
python step1a_compile_qc_normalize.py \
    --raw-dir  ./GSE123456_RAW \
    --metadata ./metadata.csv \
    --sample-id-column Sample_ID \
    --out      work/compiled_qc_hvg.h5ad

# --- 1b: scale + PCA -----------------------------------------------------
python step1b_scale_pca.py \
    --in  work/compiled_qc_hvg.h5ad \
    --out work/step1b_pca.h5ad \
    --n-comps 50 --fig-dir figures

# read the elbow off figures/pca_variance_ratio_scree.png, then:

# --- 1c: Harmony (updates the file in place) -----------------------------
python step1c_harmony.py --file work/step1b_pca.h5ad \
    --batch-key sample_id --n-pcs 30

# --- 1d: cluster ---------------------------------------------------------
python step1d_neighbors_umap_leiden.py --file work/step1b_pca.h5ad \
    --n-pcs 30 --resolution 1.0 --color sample_id disease_state \
    --fig-dir figures

# --- 1e: CellTypist, both modes ------------------------------------------
python step1e_celltypist.py --file work/step1b_pca.h5ad \
    --out-prefix work/data/celltypist --mode hvg
python step1e_celltypist.py --file work/step1b_pca.h5ad \
    --out-prefix work/data/celltypist --mode fullgene

# --- 1g: markers + unbiased DE + dotplot ---------------------------------
python step1g_marker_evidence.py --file work/step1b_pca.h5ad \
    --panel pbmc --out-dir figures

# --- 1h: assemble the template -------------------------------------------
python step1h_build_annotation_template.py --file work/step1b_pca.h5ad \
    --celltypist-dir work/data --evidence-dir figures \
    --out work/data/cluster_annotation_template.csv

# ---> open figures/dotplot__cluster_marker_dotplot.png next to the CSV
# ---> fill in manual_label_final for every cluster, notes for flags

# --- 1i: apply -----------------------------------------------------------
python step1i_apply_labels.py --file work/step1b_pca.h5ad \
    --template work/data/cluster_annotation_template.csv \
    --out work/annotated.h5ad
```

## Notes

**Steps 1c and 1d update the file in place.** They add to `obsm`/`obs` rather
than writing a new h5ad, so `work/step1b_pca.h5ad` accumulates the embedding
and the cluster columns. Keep the step 1a output as the checkpoint to fall back
to.

**Step 1b writes `.raw`.** Steps 1e and 1g restore the log-normalized matrix
from it — CellTypist needs `log1p(CP10K)` and will return confident nonsense if
handed scaled data, and marker z-scores need the full gene set so panel markers
are not missing merely because they were not highly variable. Do not strip
`.raw`.

**Defaults** are the ones used for the Itou and Zhang/Gate PBMC cohorts:
`min_genes=500`, `max_genes=5000`, `max_mito_pct=15`, Scrublet at
`expected_doublet_rate=0.06`, HVG 3000 (seurat flavor), `scale(max_value=10)`,
PCA 50 (randomized, seed 42), Harmony on sample over 30 PCs,
`neighbors(n_neighbors=15, n_pcs=30)`, Leiden **r=1.0**. They are a reasonable
starting point for 10x PBMC data and should be re-examined for anything else —
in particular tissue data, single-nucleus data, and non-immune compartments.

For mouse data, pass `--mito-prefix mt-` to step 1a.
