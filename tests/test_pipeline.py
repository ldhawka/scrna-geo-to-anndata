"""
Tests for the integration and annotation steps.

Synthetic data throughout: three "cell types" with distinct marker blocks, and
three "samples" carrying a shared shift so there is a real batch effect for
Harmony to remove.

    pytest tests/ -v
"""

import numpy as np
import pandas as pd
import pytest
import scanpy as sc

from anndata_compiler import (
    PBMC_MARKER_PANEL,
    apply_labels,
    build_annotation_template,
    cluster,
    cluster_marker_zscores,
    harmony,
    scale_and_pca,
    unbiased_de,
)

N_CELLS = 600
N_TYPES = 3
N_SAMPLES = 3
MARKERS_PER_TYPE = 20
N_NOISE_GENES = 200


@pytest.fixture(scope='module')
def adata():
    """Log-normalized, HVG-flagged object — what step 1a produces."""
    import anndata as ad

    rng = np.random.default_rng(0)
    ct = rng.integers(0, N_TYPES, N_CELLS)
    smp = rng.integers(0, N_SAMPLES, N_CELLS)

    n_genes = N_TYPES * MARKERS_PER_TYPE + N_NOISE_GENES
    X = rng.poisson(0.5, (N_CELLS, n_genes)).astype(np.float32)

    for k in range(N_TYPES):
        lo, hi = k * MARKERS_PER_TYPE, (k + 1) * MARKERS_PER_TYPE
        X[ct == k, lo:hi] += rng.poisson(15.0, ((ct == k).sum(), MARKERS_PER_TYPE))

    # Batch effect: each sample lifts a shared block of otherwise-noise genes.
    batch_block = slice(N_TYPES * MARKERS_PER_TYPE, N_TYPES * MARKERS_PER_TYPE + 40)
    for s in range(N_SAMPLES):
        X[smp == s, batch_block] += rng.poisson(6.0, ((smp == s).sum(), 40))

    a = ad.AnnData(X)
    a.var_names = [f"G{i}" for i in range(n_genes)]
    a.obs_names = [f"C{i}" for i in range(N_CELLS)]
    a.obs['sample_id'] = pd.Categorical([f"S{s}" for s in smp])
    a.obs['true_ct'] = pd.Categorical([f"T{c}" for c in ct])
    a.layers['counts'] = a.X.copy()

    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    sc.pp.highly_variable_genes(a, flavor='seurat', n_top_genes=120, subset=False)
    return a


@pytest.fixture(scope='module')
def clustered(adata):
    a = scale_and_pca(adata.copy(), n_comps=20)
    a = harmony(a, batch_key='sample_id', n_pcs=10)
    a = cluster(a, use_rep='X_pca_harmony', n_pcs=10, resolutions=(1.0,))
    return a


# --- integration ----------------------------------------------------------

def test_scale_and_pca_preserves_raw(adata):
    a = scale_and_pca(adata.copy(), n_comps=20)
    assert 'X_pca' in a.obsm
    assert a.obsm['X_pca'].shape == (N_CELLS, 20)
    # .raw must hold the full, pre-scale gene set — steps 1e/1g depend on it.
    assert a.raw is not None
    assert a.raw.n_vars == adata.n_vars
    assert a.n_vars < adata.n_vars
    assert a.raw.X.min() >= 0, "raw should be log-normalized, not scaled"


def test_scale_and_pca_requires_hvg_flag(adata):
    a = adata.copy()
    del a.var['highly_variable']
    with pytest.raises(ValueError, match='highly_variable'):
        scale_and_pca(a)


def test_harmony_produces_embedding(adata):
    a = scale_and_pca(adata.copy(), n_comps=20)
    a = harmony(a, batch_key='sample_id', n_pcs=10)
    assert a.obsm['X_pca_harmony'].shape == (N_CELLS, 10)
    assert np.isfinite(a.obsm['X_pca_harmony']).all()


def test_harmony_reduces_batch_separation(adata):
    """The corrected embedding should mix samples better than raw PCA."""
    a = scale_and_pca(adata.copy(), n_comps=20)
    a = harmony(a, batch_key='sample_id', n_pcs=10)

    def batch_spread(emb):
        # Mean distance between per-sample centroids: lower = better mixed.
        cents = np.vstack([emb[a.obs['sample_id'] == s].mean(axis=0)
                           for s in a.obs['sample_id'].cat.categories])
        scale = emb.std()
        d = [np.linalg.norm(cents[i] - cents[j])
             for i in range(len(cents)) for j in range(i + 1, len(cents))]
        return np.mean(d) / scale

    assert batch_spread(a.obsm['X_pca_harmony']) < batch_spread(a.obsm['X_pca'][:, :10])


def test_harmony_rejects_single_batch(adata):
    a = scale_and_pca(adata.copy(), n_comps=20)
    a.obs['only_one'] = 'x'
    with pytest.raises(ValueError, match='nothing to correct'):
        harmony(a, batch_key='only_one', n_pcs=10)


def test_harmony_rejects_missing_pca(adata):
    with pytest.raises(ValueError, match='X_pca'):
        harmony(adata.copy(), batch_key='sample_id', n_pcs=10)


def test_cluster_recovers_cell_types(clustered):
    assert 'leiden_r1_0' in clustered.obs
    assert 'X_umap' in clustered.obsm
    # Each true type should be dominated by one Leiden cluster.
    xt = pd.crosstab(clustered.obs['true_ct'], clustered.obs['leiden_r1_0'])
    purity = (xt.max(axis=1) / xt.sum(axis=1)).min()
    assert purity > 0.85, f"lowest per-type purity {purity:.2f}"


def test_cluster_rejects_missing_rep(adata):
    a = scale_and_pca(adata.copy(), n_comps=20)
    with pytest.raises(ValueError, match='X_pca_harmony'):
        cluster(a, use_rep='X_pca_harmony', n_pcs=10)


def test_cluster_multiple_resolutions(adata):
    a = scale_and_pca(adata.copy(), n_comps=20)
    a = harmony(a, batch_key='sample_id', n_pcs=10)
    a = cluster(a, use_rep='X_pca_harmony', n_pcs=10,
                resolutions=(0.3, 1.0), compute_umap=False)
    assert 'leiden_r0_3' in a.obs and 'leiden_r1_0' in a.obs


# --- annotation -----------------------------------------------------------

def test_marker_zscores_shape_and_missing(clustered):
    full = clustered.raw.to_adata()
    full.obs = clustered.obs.copy()
    panel = {f"T{k}": [f"G{i}" for i in range(k * MARKERS_PER_TYPE,
                                              k * MARKERS_PER_TYPE + 5)]
             for k in range(N_TYPES)}
    panel['absent'] = ['NOT_A_REAL_GENE']

    z, raw, missing = cluster_marker_zscores(full, panel, 'leiden_r1_0')
    assert missing == ['NOT_A_REAL_GENE']
    assert z.shape == raw.shape
    assert z.shape[1] == N_TYPES * 5
    assert z.shape[0] == clustered.obs['leiden_r1_0'].nunique()
    # Raw means are returned so a high z on a near-zero gene can be caught.
    assert (raw.to_numpy() >= 0).all()


def test_marker_zscores_accepts_flat_list(clustered):
    full = clustered.raw.to_adata()
    full.obs = clustered.obs.copy()
    z, raw, missing = cluster_marker_zscores(full, ['G0', 'G1', 'G2'], 'leiden_r1_0')
    assert list(z.columns) == ['G0', 'G1', 'G2']
    assert missing == []


def test_unbiased_de_covers_all_clusters(clustered):
    full = clustered.raw.to_adata()
    full.obs = clustered.obs.copy()
    de = unbiased_de(full, 'leiden_r1_0', n_genes=5)
    assert len(de) == clustered.obs['leiden_r1_0'].nunique()
    assert 'top5_de_genes' in de.columns
    assert de['top5_de_genes'].str.count(';').eq(4).all()


def test_build_template_leaves_labels_blank(clustered, tmp_path):
    out = tmp_path / 'template.csv'
    t = build_annotation_template(clustered, 'leiden_r1_0', str(out))

    assert out.exists()
    for col in ('user_label', 'manual_label_final', 'notes'):
        assert col in t.columns
        assert (t[col] == '').all(), f"{col} should be blank for the human to fill"
    assert t['n_cells'].sum() == clustered.n_obs
    assert abs(t['pct_of_total'].sum() - 100) < 0.1


def test_build_template_merges_evidence(clustered, tmp_path):
    full = clustered.raw.to_adata()
    full.obs = clustered.obs.copy()
    panel = {f"T{k}": [f"G{i}" for i in range(k * MARKERS_PER_TYPE,
                                              k * MARKERS_PER_TYPE + 5)]
             for k in range(N_TYPES)}
    z, _, _ = cluster_marker_zscores(full, panel, 'leiden_r1_0')
    de = unbiased_de(full, 'leiden_r1_0', n_genes=5)

    out = tmp_path / 'template.csv'
    t = build_annotation_template(clustered, 'leiden_r1_0', str(out),
                                  zscores=z, de_genes=de)
    assert 'top_panel_markers_zscore' in t.columns
    assert 'top5_de_genes' in t.columns
    assert t['top_panel_markers_zscore'].notna().all()


def test_apply_labels_rejects_blank_template(clustered, tmp_path):
    out = tmp_path / 'template.csv'
    build_annotation_template(clustered, 'leiden_r1_0', str(out))
    with pytest.raises(ValueError, match='no manual_label_final'):
        apply_labels(clustered.copy(), str(out), 'leiden_r1_0')


def test_apply_labels_merges_sister_clusters(clustered, tmp_path):
    out = tmp_path / 'template.csv'
    t = build_annotation_template(clustered, 'leiden_r1_0', str(out))
    # Every cluster gets the same label — they must collapse into one category.
    t['manual_label_final'] = 'One Cell Type'
    t.to_csv(out, index=False)

    a = apply_labels(clustered.copy(), str(out), 'leiden_r1_0')
    assert a.obs['cell_type'].nunique() == 1
    assert a.obs['cell_type'].iloc[0] == 'One Cell Type'


def test_apply_labels_flags_and_drops(clustered, tmp_path):
    out = tmp_path / 'template.csv'
    t = build_annotation_template(clustered, 'leiden_r1_0', str(out))
    t['manual_label_final'] = [f"Type{i}" for i in range(len(t))]
    t['notes'] = ['doublet'] + [''] * (len(t) - 1)
    t.to_csv(out, index=False)

    flagged_cluster = str(t['cluster'].iloc[0])
    n_flagged = int((clustered.obs['leiden_r1_0'].astype(str) == flagged_cluster).sum())

    kept = apply_labels(clustered.copy(), str(out), 'leiden_r1_0', drop_flagged=False)
    assert int(kept.obs['cell_type_flag'].sum()) == n_flagged
    assert kept.n_obs == clustered.n_obs

    dropped = apply_labels(clustered.copy(), str(out), 'leiden_r1_0', drop_flagged=True)
    assert dropped.n_obs == clustered.n_obs - n_flagged


def test_apply_labels_rejects_unmapped_cluster(clustered, tmp_path):
    out = tmp_path / 'template.csv'
    t = build_annotation_template(clustered, 'leiden_r1_0', str(out))
    t['manual_label_final'] = 'X'
    t = t.iloc[:-1]  # drop a cluster the data still contains
    t.to_csv(out, index=False)
    with pytest.raises(ValueError, match='absent from template'):
        apply_labels(clustered.copy(), str(out), 'leiden_r1_0')


# --- panels ---------------------------------------------------------------

def test_pbmc_panel_is_well_formed():
    assert isinstance(PBMC_MARKER_PANEL, dict)
    all_genes = [g for v in PBMC_MARKER_PANEL.values() for g in v]
    assert len(all_genes) > 100
    assert all(isinstance(g, str) and g.strip() == g for g in all_genes)
    for essential in ('CD3D', 'CD14', 'MS4A1', 'FCGR3A', 'FOXP3', 'PPBP'):
        assert essential in all_genes, f"{essential} missing from the PBMC panel"


@pytest.mark.network
def test_microglia_panel_loads_from_atlas():
    """Requires network — the atlas is served from GitHub Pages."""
    from anndata_compiler import load_microglia_panel

    panel = load_microglia_panel()
    assert len(panel) == 23, "expected 23 sub-states"
    assert 'P2RY12' in panel['Homeostatic']

    umbrella = load_microglia_panel(level='umbrella')
    assert len(umbrella) == 9, "expected 9 umbrella families"

    mouse = load_microglia_panel(species='mouse', denoised=False)
    assert 'P2ry12' in mouse['Homeostatic'], "mouse symbols should be title-case"
