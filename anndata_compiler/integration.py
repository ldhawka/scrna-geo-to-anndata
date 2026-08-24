"""
Scaling, PCA, Harmony batch integration, and clustering.

This is the second half of the pipeline. `compiler.GEOAnndataCompiler` takes
raw GEO files to a QC'd, normalized, HVG-flagged h5ad; the functions here take
that object to a batch-corrected embedding with Leiden clusters, ready for
annotation (see `annotation.py`).

The defaults are the ones used for the Itou and Zhang/Gate PBMC cohorts. They
are a reasonable starting point for 10x PBMC data and should be re-examined for
anything else — in particular tissue data, nuclei, and non-immune compartments.
"""

import time

import numpy as np
import scanpy as sc

DEFAULT_SEED = 42
DEFAULT_N_COMPS = 50
DEFAULT_SCALE_MAX_VALUE = 10
DEFAULT_N_PCS = 30
DEFAULT_N_NEIGHBORS = 15
DEFAULT_RESOLUTION = 1.0


def scale_and_pca(
    adata,
    n_comps=DEFAULT_N_COMPS,
    max_value=DEFAULT_SCALE_MAX_VALUE,
    use_highly_variable=True,
    random_state=DEFAULT_SEED,
):
    """
    Subset to HVGs, scale, and run PCA.

    Scaling to unit variance with `max_value` clipping stops a handful of very
    highly expressed genes from dominating the principal components. The clip
    matters more than it looks: without it, a single outlier cell can define
    a PC on its own.

    `adata.X` is overwritten with the scaled matrix, so the log-normalized
    values must already be recoverable — `raw` is set here before scaling, and
    integer counts remain in `layers['counts']` from the compile step.

    Parameters
    ----------
    adata : AnnData
        Log-normalized, HVG-flagged object from the compile step.
    n_comps : int
        Number of principal components. 50 is a deliberate over-estimate; the
        number actually carried forward is `n_pcs` in `harmony` / `cluster`.
    max_value : float
        Clip scaled values above this. 10 is the scanpy convention.
    use_highly_variable : bool
        Subset to `var['highly_variable']` before scaling. Almost always True.
    random_state : int

    Returns
    -------
    AnnData with `obsm['X_pca']` and `uns['pca']`.
    """
    np.random.seed(random_state)

    if use_highly_variable:
        if 'highly_variable' not in adata.var.columns:
            raise ValueError(
                "var['highly_variable'] not found — run the compile step first, "
                "or pass use_highly_variable=False."
            )
        adata.raw = adata
        adata = adata[:, adata.var['highly_variable']].copy()
        print(f"  Subset to {adata.n_vars:,} highly variable genes")

    sc.pp.scale(adata, max_value=max_value)
    print(f"  Scaled (max_value={max_value})")

    sc.tl.pca(
        adata,
        n_comps=n_comps,
        svd_solver='randomized',
        random_state=random_state,
    )
    print(f"  PCA complete: obsm['X_pca'] {adata.obsm['X_pca'].shape}")

    return adata


def harmony(
    adata,
    batch_key,
    n_pcs=DEFAULT_N_PCS,
    max_iter_harmony=20,
    epsilon_harmony=1e-5,
    random_state=DEFAULT_SEED,
):
    """
    Batch-correct the PCA embedding with Harmony.

    Multi-sample GEO compilations almost always carry a per-sample technical
    effect — different capture dates, chemistry versions, sequencing depths.
    Without correction, Leiden will happily return clusters that are one sample
    each. Check `plot_post_harmony_diagnostics`-style UMAPs coloured by sample
    before and after to confirm the correction did something and did not
    over-merge genuinely distinct biology.

    `batch_key` should normally be the sample/library ID. Use a coarser key
    (e.g. sequencing batch) only when samples within a batch are known to be
    technically identical, and never use a key that is confounded with the
    biological contrast — correcting on `disease_state` removes the signal.

    Parameters
    ----------
    adata : AnnData
        Must have `obsm['X_pca']` from `scale_and_pca`.
    batch_key : str
        Column in `adata.obs` identifying the batch.
    n_pcs : int
        PCs to feed Harmony. This is also the number to pass to `cluster`.
    max_iter_harmony, epsilon_harmony : int, float
        Harmony convergence controls.
    random_state : int

    Returns
    -------
    AnnData with `obsm['X_pca_harmony']` (n_cells x n_pcs, float32).
    """
    import harmonypy as hm

    np.random.seed(random_state)

    if 'X_pca' not in adata.obsm:
        raise ValueError("obsm['X_pca'] not found — run scale_and_pca first.")
    if batch_key not in adata.obs.columns:
        raise ValueError(f"batch_key '{batch_key}' not in adata.obs")

    n_cells, n_avail = adata.obsm['X_pca'].shape
    if n_avail < n_pcs:
        raise ValueError(f"Only {n_avail} PCs available, need {n_pcs}")

    pca_input = np.asarray(adata.obsm['X_pca'][:, :n_pcs], dtype=np.float64)
    n_batches = adata.obs[batch_key].nunique()
    print(f"  Harmony input: {pca_input.shape}, batch_key={batch_key} ({n_batches} batches)")

    if n_batches < 2:
        raise ValueError(
            f"batch_key '{batch_key}' has {n_batches} unique value(s) — nothing to correct."
        )

    t0 = time.time()
    ho = hm.run_harmony(
        pca_input, adata.obs, [batch_key],
        max_iter_harmony=max_iter_harmony,
        epsilon_harmony=epsilon_harmony,
        random_state=random_state,
    )
    print(f"  Harmony complete in {(time.time() - t0) / 60:.1f} min")

    # harmonypy has returned Z_corr in both orientations across versions.
    z = np.asarray(ho.Z_corr)
    if z.shape == (n_pcs, n_cells):
        z = z.T
    if z.shape != (n_cells, n_pcs):
        raise ValueError(f"Unexpected Z_corr shape {z.shape}")

    adata.obsm['X_pca_harmony'] = np.asarray(z, dtype=np.float32)
    if not np.isfinite(adata.obsm['X_pca_harmony']).all():
        raise ValueError("Harmony returned non-finite values")

    print(f"  obsm['X_pca_harmony'] {adata.obsm['X_pca_harmony'].shape}")
    return adata


def cluster(
    adata,
    use_rep='X_pca_harmony',
    n_neighbors=DEFAULT_N_NEIGHBORS,
    n_pcs=DEFAULT_N_PCS,
    resolutions=(DEFAULT_RESOLUTION,),
    compute_umap=True,
    random_state=DEFAULT_SEED,
):
    """
    Neighbourhood graph, UMAP, and Leiden clustering.

    `use_rep` defaults to the Harmony embedding. Passing 'X_pca' here after
    running Harmony is a silent, common, and consequential mistake: the
    clusters come back uncorrected and the UMAP separates by sample.

    Several resolutions can be requested at once; each is stored as
    `leiden_r{res}` with dots replaced by underscores (e.g. `leiden_r1_0`).
    Multiple resolutions are cheap to compute and worth keeping, but pick one
    for the annotation pass rather than annotating each — see
    docs/annotation_guide.md.

    Parameters
    ----------
    adata : AnnData
    use_rep : str
        Key in `adata.obsm`. Use 'X_pca_harmony' after batch correction.
    n_neighbors : int
        Lower (5-10) preserves fine structure; higher (20-30) smooths.
    n_pcs : int
        Dimensions of `use_rep` to use. Match what was passed to `harmony`.
    resolutions : iterable of float
        Leiden resolutions. 1.0 is the PBMC default used here.
    compute_umap : bool
    random_state : int

    Returns
    -------
    AnnData with `obsm['X_umap']` and one `obs['leiden_r*']` column per resolution.
    """
    if use_rep not in adata.obsm:
        raise ValueError(
            f"obsm['{use_rep}'] not found. Available: {list(adata.obsm.keys())}"
        )

    sc.pp.neighbors(
        adata,
        n_neighbors=n_neighbors,
        n_pcs=n_pcs,
        use_rep=use_rep,
        random_state=random_state,
    )
    print(f"  Neighbours: n_neighbors={n_neighbors}, n_pcs={n_pcs}, use_rep={use_rep}")

    if compute_umap:
        sc.tl.umap(adata, random_state=random_state)
        print("  UMAP complete")

    for res in resolutions:
        key = f"leiden_r{str(res).replace('.', '_')}"
        sc.tl.leiden(
            adata,
            resolution=res,
            key_added=key,
            random_state=random_state,
            flavor='igraph',
            n_iterations=2,
            directed=False,
        )
        n_clusters = adata.obs[key].nunique()
        print(f"  Leiden r={res}: {n_clusters} clusters -> obs['{key}']")

    return adata
