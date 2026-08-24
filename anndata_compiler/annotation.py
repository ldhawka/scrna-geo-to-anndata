"""
Cluster annotation by triangulation.

No single method annotates clusters reliably. Reference-mapping tools are
confident and systematically biased; marker panels are interpretable but
incomplete; unbiased differential expression is honest but unlabelled. This
module runs several of them, assembles the evidence into one table, and leaves
the label column blank for a human to fill in.

The evidence sources are:

  1. CellTypist `Immune_All_Low` and `Immune_All_High`, run twice — once on the
     HVG matrix and once on the full gene set. Four label columns. Where the two
     runs disagree, that disagreement is itself informative: a cluster whose
     HVG and full-gene calls flip between central-memory and effector-memory is
     usually a cluster containing both.
  2. Top panel markers by z-score, from a curated marker panel.
  3. Top unbiased Wilcoxon differential-expression genes, panel-free.
  4. Optionally, any labels the original authors deposited — carried as a
     column, not adopted. See docs/annotation_guide.md for why.

The human then reads a dotplot and decides. `build_annotation_template` writes
the table; `apply_labels` reads the filled-in version back.

See docs/annotation_guide.md for the marker tables, the interpretation traps,
and the granularity problem when comparing cohorts.
"""

import os

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse

# Curated PBMC panel. Human symbols. Grouped for dotplot readability.
PBMC_MARKER_PANEL = {
    'T (all)':            ['CD3D', 'CD3E', 'CD3G', 'TRAC'],
    'CD4 T':              ['CD4', 'IL7R', 'S100A4', 'IL32'],
    'Naive/CM':           ['CCR7', 'SELL', 'LEF1', 'TCF7'],
    'Th1':                ['CXCR3', 'TBX21', 'CCR5', 'IFNG'],
    'Th17':               ['RORC', 'CCR6', 'IL17A', 'IL17F', 'IL23R'],
    'Treg':               ['FOXP3', 'IL2RA', 'CTLA4', 'TIGIT'],
    'CD8 T':              ['CD8A', 'CD8B'],
    'Cytotoxic':          ['GZMA', 'GZMB', 'GZMH', 'GZMK', 'PRF1', 'NKG7', 'GNLY'],
    'Effector/TEMRA':     ['FGFBP2', 'KLRG1', 'CX3CR1', 'EOMES'],
    'Exhaustion':         ['PDCD1', 'LAG3', 'HAVCR2', 'TOX'],
    'MAIT':               ['SLC4A10', 'TRAV1-2', 'KLRB1', 'IL18R1'],
    'gdT':                ['TRDC', 'TRGC1', 'XCL1', 'XCL2'],
    'NK':                 ['NCAM1', 'FCGR3A', 'KLRF1', 'KLRD1'],
    'B (all)':            ['CD19', 'MS4A1', 'CD79A', 'CD79B'],
    'Naive B':            ['IGHD', 'IGHM', 'TCL1A', 'FCER2', 'CR2'],
    'Memory B':           ['CD27', 'AIM2', 'BANK1', 'TNFRSF13B'],
    'Plasma':             ['JCHAIN', 'XBP1', 'MZB1', 'PRDM1', 'SDC1', 'TNFRSF17'],
    'Classical mono':     ['CD14', 'VCAN', 'S100A8', 'S100A9', 'S100A12', 'LYZ', 'FCN1'],
    'Non-classical mono': ['CDKN1C', 'MS4A7', 'LST1', 'HES4'],
    'cDC1':               ['CLEC9A', 'XCR1', 'BATF3', 'IRF8'],
    'cDC2':               ['CD1C', 'CLEC10A', 'FCER1A'],
    'pDC':                ['TCF4', 'IRF7', 'LILRA4', 'CLEC4C', 'PLD4', 'BCL11A'],
    'HLA class II':       ['HLA-DRA', 'HLA-DPA1', 'HLA-DPB1', 'CD74'],
    'Platelet':           ['PPBP', 'PF4', 'GP1BB', 'TUBB1', 'NRGN', 'ITGA2B'],
    'Mast/Basophil':      ['TPSAB1', 'TPSB2', 'CPA3', 'MS4A2', 'HDC'],
    'HSPC':               ['CD34', 'SPINK2', 'PRSS57'],
    'Proliferation':      ['MKI67', 'TOP2A', 'STMN1', 'PCLAF'],
}

# Microglia / CNS-myeloid markers are maintained separately as a hierarchical
# atlas — 9 umbrella families over 23 sub-states, human and mouse symbols,
# up- and down-regulated sets, denoised of dissociation artifacts, with DOIs.
# Use it instead of writing a microglia panel by hand.
MICROGLIA_ATLAS_URL = "https://ldhawka.github.io/microglia-subtype-markers/"
MICROGLIA_ATLAS_JSON = (
    "https://ldhawka.github.io/microglia-subtype-markers/"
    "data/microglia_subtypes_hierarchical.json"
)


def load_microglia_panel(source=MICROGLIA_ATLAS_JSON, species='human',
                         denoised=True, level='sub_subtype'):
    """
    Build a marker panel dict from the microglia subtype atlas.

    Brain myeloid cells are not covered by `PBMC_MARKER_PANEL`, and the
    disease-associated states in particular are a continuum rather than
    discrete types — scoring per sub-state and rolling up to the umbrella
    family is more honest than forcing a one-of-23 call. See
    https://ldhawka.github.io/microglia-subtype-markers/ for the table, the
    validation, and the primary references.

    Parameters
    ----------
    source : str
        URL or local path to the hierarchical JSON.
    species : {'human', 'mouse'}
    denoised : bool
        Use the gene sets with dissociation-response genes removed. Recommended
        for anything that went through enzymatic dissociation.
    level : {'sub_subtype', 'umbrella'}
        Return one entry per sub-state, or pooled per umbrella family.

    Returns
    -------
    dict of {label: [genes]}, in the shape `cluster_marker_zscores` expects.
    """
    import json

    if source.startswith('http'):
        import urllib.request
        with urllib.request.urlopen(source) as fh:
            families = json.load(fh)
    else:
        with open(source) as fh:
            families = json.load(fh)

    # The hierarchical JSON is a list of 9 umbrella families, each carrying a
    # 'sub_subtypes' list. Denoised gene sets exist only for human.
    if denoised and species == 'human':
        key = 'upregulated_genes_denoised'
    else:
        key = f'upregulated_genes_{species}'

    panel = {}
    for family in families:
        umbrella = family.get('umbrella')
        for rec in family.get('sub_subtypes', []):
            raw = rec.get(key) or rec.get(f'upregulated_genes_{species}') or ''
            genes = [g.strip() for g in raw.replace(';', ',').split(',') if g.strip()]
            genes = [g for g in genes if g != '—']
            if not genes:
                continue
            label = umbrella if level == 'umbrella' else rec.get('sub_subtype')
            existing = panel.setdefault(label, [])
            existing.extend(g for g in genes if g not in existing)

    print(f"  Loaded {len(panel)} {level} panels from the microglia atlas "
          f"({species}, denoised={denoised})")
    return panel


def verify_lognorm(adata, max_expected=20):
    """
    Sanity-check that `adata.X` is log1p(CP10K), not raw counts or scaled data.

    CellTypist expects log-normalized input and will return confident nonsense
    if handed something else. This catches the two common mistakes: passing the
    post-`scale` object (negative values), or passing raw counts (large values).
    """
    sample = adata.X[:1000]
    if scipy.sparse.issparse(sample):
        sample = sample.toarray()
    smax = float(np.max(sample))
    smin = float(np.min(sample))

    if smin < 0:
        print(f"  WARNING: X min = {smin:.2f} < 0 — this looks like scaled data. "
              f"CellTypist needs log1p(CP10K); use .raw or re-load the pre-scale object.")
    elif smax > max_expected:
        print(f"  WARNING: X max = {smax:.2f} — expected <{max_expected} for log1p(CP10K). "
              f"These may be raw counts.")
    else:
        print(f"  log-norm check: X max = {smax:.2f} (<{max_expected}) OK")


def _summarize_per_cluster(obs, label_col, conf_col, cluster_col, prefix):
    """Majority-vote a per-cell label column over clusters, with top-3 breakdown."""
    rows = []
    clusters = obs[cluster_col].astype(str).unique()
    try:
        clusters = sorted(clusters, key=int)
    except ValueError:
        clusters = sorted(clusters)

    for cluster in clusters:
        mask = obs[cluster_col].astype(str) == cluster
        n = int(mask.sum())
        labels = obs.loc[mask, label_col].astype(str)
        top3 = labels.value_counts().head(3)
        row = {
            'cluster': cluster,
            'n_cells': n,
            f'{prefix}_majority_label': top3.index[0],
            f'{prefix}_majority_fraction': round(float(top3.iloc[0]) / n, 3),
            f'{prefix}_top3_labels': '; '.join(
                f"{l} ({c / n:.1%})" for l, c in top3.items()
            ),
        }
        if conf_col is not None and conf_col in obs.columns:
            confs = pd.to_numeric(obs.loc[mask, conf_col], errors='coerce')
            row[f'{prefix}_mean_confidence'] = round(float(confs.mean()), 3)
        rows.append(row)

    return pd.DataFrame(rows)


def run_celltypist(
    adata,
    cluster_key,
    models=('Immune_All_Low.pkl', 'Immune_All_High.pkl'),
    prefixes=('low', 'high'),
    majority_voting=True,
):
    """
    Run CellTypist and summarize per cluster.

    Run this twice — once on the HVG object and once on the full gene set —
    and merge with `prefixes=('fg_low', 'fg_high')` on the second pass. The
    full-gene run is usually the more confident of the two; where the two
    disagree, record the disagreement rather than picking a winner.

    Parameters
    ----------
    adata : AnnData
        Log-normalized (`log1p(CP10K)`). NOT scaled. Verified on entry.
    cluster_key : str
        e.g. 'leiden_r1_0'. Used as CellTypist's over-clustering for majority voting.
    models : tuple of str
        CellTypist model filenames. Downloaded on first use.
    prefixes : tuple of str
        Column prefixes in the output, one per model.
    majority_voting : bool

    Returns
    -------
    (per_cell_df, per_cluster_df)
    """
    import celltypist
    from celltypist import models as ct_models

    verify_lognorm(adata)
    ct_models.download_models(force_update=False, model=list(models))

    per_cell = pd.DataFrame(index=adata.obs_names)
    per_cluster = None

    for model, prefix in zip(models, prefixes):
        print(f"  CellTypist: {model} -> '{prefix}_*'")
        pred = celltypist.annotate(
            adata,
            model=model,
            majority_voting=majority_voting,
            over_clustering=cluster_key,
        )
        result = pred.predicted_labels
        label_col = 'majority_voting' if majority_voting else 'predicted_labels'

        obs = adata.obs.copy()
        obs[f'{prefix}_label'] = result[label_col].values
        per_cell[f'{prefix}_label'] = result[label_col].values

        # Confidence = the winning label's probability, per cell.
        conf_col = None
        if getattr(pred, 'probability_matrix', None) is not None:
            conf = pred.probability_matrix.max(axis=1)
            obs[f'{prefix}_conf'] = conf.values
            per_cell[f'{prefix}_conf'] = conf.values
            conf_col = f'{prefix}_conf'

        summary = _summarize_per_cluster(obs, f'{prefix}_label', conf_col,
                                         cluster_key, prefix)
        if per_cluster is None:
            per_cluster = summary
        else:
            drop = [c for c in summary.columns if c in ('n_cells',)]
            per_cluster = per_cluster.merge(
                summary.drop(columns=drop), on='cluster', how='outer'
            )

    return per_cell, per_cluster


def cluster_marker_zscores(adata, marker_panel, cluster_key, layer=None):
    """
    Mean expression of each panel marker per cluster, z-scored across clusters.

    The z-score is across clusters, so a value answers "is this marker higher in
    this cluster than in the others" — not "is this marker expressed". Those
    come apart badly for sparse genes: a marker at raw mean 0.01 everywhere and
    0.10 in one cluster scores z = +3, which reads as a strong positive and is
    noise. `raw_means` is returned alongside for exactly this reason; always
    check the raw value before acting on a high z. See docs/annotation_guide.md.

    Parameters
    ----------
    adata : AnnData
        Log-normalized. Use the full-gene object so panel markers are not
        missing merely because they were not highly variable.
    marker_panel : dict of {group: [genes]} or list of genes
    cluster_key : str
    layer : str or None

    Returns
    -------
    (zscores_df, raw_means_df, missing_genes) — clusters x markers.
    """
    if isinstance(marker_panel, dict):
        genes = [g for group in marker_panel.values() for g in group]
    else:
        genes = list(marker_panel)

    seen = set()
    genes = [g for g in genes if not (g in seen or seen.add(g))]

    present = [g for g in genes if g in adata.var_names]
    missing = [g for g in genes if g not in adata.var_names]
    if missing:
        print(f"  {len(missing)} panel genes not in var_names: {', '.join(missing[:10])}"
              + (' ...' if len(missing) > 10 else ''))

    X = adata[:, present].layers[layer] if layer else adata[:, present].X
    if scipy.sparse.issparse(X):
        X = X.toarray()

    df = pd.DataFrame(X, columns=present, index=adata.obs_names)
    df[cluster_key] = adata.obs[cluster_key].astype(str).values
    raw_means = df.groupby(cluster_key, observed=True).mean()

    try:
        raw_means = raw_means.loc[sorted(raw_means.index, key=int)]
    except ValueError:
        raw_means = raw_means.sort_index()

    sd = raw_means.std(axis=0).replace(0, np.nan)
    zscores = (raw_means - raw_means.mean(axis=0)) / sd

    return zscores, raw_means, missing


def unbiased_de(adata, cluster_key, n_genes=10, method='wilcoxon'):
    """
    Top differentially expressed genes per cluster, with no panel involved.

    This is the check on the panel: if a cluster's top DE genes are all
    ribosomal, or all mitochondrial, or all long intron-rich nuclear
    transcripts (MALAT1, ARHGAP15, CAMK4, INPP4B), the cluster may be split on
    a technical axis rather than a biological one, regardless of what the
    marker panel says about it.

    Returns
    -------
    DataFrame with columns: cluster, top{n}_de_genes.
    """
    sc.tl.rank_genes_groups(adata, groupby=cluster_key, method=method)

    rows = []
    for cluster in adata.obs[cluster_key].astype(str).unique():
        names = [adata.uns['rank_genes_groups']['names'][i][cluster]
                 for i in range(n_genes)]
        rows.append({'cluster': cluster,
                     f'top{n_genes}_de_genes': '; '.join(names)})

    df = pd.DataFrame(rows)
    try:
        df = df.iloc[np.argsort([int(c) for c in df['cluster']])].reset_index(drop=True)
    except ValueError:
        df = df.sort_values('cluster').reset_index(drop=True)
    return df


def _top_markers_per_cluster(zscores, n_top=8):
    rows = []
    for cluster, row in zscores.iterrows():
        top = row.sort_values(ascending=False).head(n_top)
        rows.append({
            'cluster': str(cluster),
            'top_panel_markers_zscore': '; '.join(f"{m} ({v:+.2f})" for m, v in top.items()),
        })
    return pd.DataFrame(rows)


def build_annotation_template(
    adata,
    cluster_key,
    output_path,
    celltypist_hvg=None,
    celltypist_fullgene=None,
    zscores=None,
    de_genes=None,
    author_labels=None,
    n_top_markers=8,
):
    """
    Assemble every evidence source into one CSV, with the label columns blank.

    The blank columns are the point. The table gathers what the automated
    methods think; it does not decide. Fill in `manual_label_final` from the
    dotplot and the evidence columns, put anything unresolved in `notes`
    (`doublet`, `artifact`, `low_confidence`, `merge_with_N`), then run
    `apply_labels`.

    Parameters
    ----------
    adata : AnnData
    cluster_key : str
    output_path : str
    celltypist_hvg, celltypist_fullgene : DataFrame or None
        Per-cluster summaries from `run_celltypist`.
    zscores : DataFrame or None
        From `cluster_marker_zscores`.
    de_genes : DataFrame or None
        From `unbiased_de`.
    author_labels : DataFrame or None
        Two columns: cluster, author_label. Carried for comparison, never adopted
        automatically — see docs/annotation_guide.md on annotation granularity.
    n_top_markers : int

    Returns
    -------
    The template DataFrame (also written to `output_path`).
    """
    counts = adata.obs[cluster_key].astype(str).value_counts()
    try:
        order = sorted(counts.index, key=int)
    except ValueError:
        order = sorted(counts.index)

    template = pd.DataFrame({
        'cluster': order,
        'n_cells': [int(counts[c]) for c in order],
    })
    template['pct_of_total'] = (
        template['n_cells'] / template['n_cells'].sum() * 100
    ).round(2)

    for df, name in ((celltypist_hvg, 'CellTypist (HVG)'),
                     (celltypist_fullgene, 'CellTypist (full-gene)')):
        if df is not None:
            df = df.copy()
            df['cluster'] = df['cluster'].astype(str)
            cols = [c for c in df.columns if c not in ('cluster', 'n_cells')]
            template = template.merge(df[['cluster'] + cols], on='cluster', how='left')
            print(f"  merged {name}: {len(cols)} columns")

    if zscores is not None:
        top = _top_markers_per_cluster(zscores, n_top=n_top_markers)
        template = template.merge(top, on='cluster', how='left')

    if de_genes is not None:
        de_genes = de_genes.copy()
        de_genes['cluster'] = de_genes['cluster'].astype(str)
        template = template.merge(de_genes, on='cluster', how='left')

    if author_labels is not None:
        author_labels = author_labels.copy()
        author_labels['cluster'] = author_labels['cluster'].astype(str)
        template = template.merge(author_labels, on='cluster', how='left')

    # Deliberately blank — the calls are the analyst's.
    template['user_label'] = ''
    template['manual_label_final'] = ''
    template['notes'] = ''

    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    template.to_csv(output_path, index=False)
    print(f"  Wrote annotation template: {output_path} ({template.shape[0]} clusters, "
          f"{template.shape[1]} columns)")
    print("  Fill in 'manual_label_final' and 'notes', then run apply_labels().")

    return template


def apply_labels(
    adata,
    template_path,
    cluster_key,
    label_column='manual_label_final',
    notes_column='notes',
    output_key='cell_type',
    drop_flagged=False,
    flag_values=('doublet', 'artifact'),
):
    """
    Map the filled-in template back onto cells.

    Clusters sharing a label collapse into one category, which is the intended
    behaviour — sister clusters of the same cell type should merge. Use exactly
    consistent label strings for that to work ('Naive CD4 T cells' and
    'Naive CD4 T' will not merge).

    Parameters
    ----------
    adata : AnnData
    template_path : str
    cluster_key : str
    label_column : str
    notes_column : str
        Notes matching `flag_values` mark clusters for exclusion.
    output_key : str
        New column in `adata.obs`.
    drop_flagged : bool
        If True, remove flagged cells; otherwise keep them with
        `obs[output_key + '_flag']` set.
    flag_values : tuple of str

    Returns
    -------
    AnnData with `obs[output_key]` (and `obs[output_key + '_flag']`).
    """
    template = pd.read_csv(template_path, dtype={'cluster': str})

    if label_column not in template.columns:
        raise ValueError(f"'{label_column}' not in template columns: {list(template.columns)}")

    unlabelled = template[template[label_column].isna() | (template[label_column] == '')]
    if len(unlabelled):
        raise ValueError(
            f"{len(unlabelled)} cluster(s) have no {label_column}: "
            f"{', '.join(unlabelled['cluster'])}. Fill these in before applying."
        )

    mapping = dict(zip(template['cluster'], template[label_column].str.strip()))
    clusters = adata.obs[cluster_key].astype(str)

    unmapped = set(clusters.unique()) - set(mapping)
    if unmapped:
        raise ValueError(f"Clusters present in data but absent from template: {sorted(unmapped)}")

    adata.obs[output_key] = clusters.map(mapping).astype('category')
    print(f"  Applied {len(set(mapping.values()))} unique labels from "
          f"{len(mapping)} clusters -> obs['{output_key}']")

    flag_key = f'{output_key}_flag'
    if notes_column in template.columns:
        notes = template[notes_column].fillna('').str.strip().str.lower()
        flagged = {
            c for c, note in zip(template['cluster'], notes)
            if any(f in note for f in flag_values)
        }
        adata.obs[flag_key] = clusters.isin(flagged).values
        n_flagged = int(adata.obs[flag_key].sum())
        if n_flagged:
            print(f"  {n_flagged:,} cells in {len(flagged)} flagged cluster(s): "
                  f"{sorted(flagged)}")
            if drop_flagged:
                adata = adata[~adata.obs[flag_key]].copy()
                print(f"  Dropped flagged cells -> {adata.n_obs:,} cells remain")

    return adata
