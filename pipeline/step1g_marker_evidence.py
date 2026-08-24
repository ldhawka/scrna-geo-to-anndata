#!/usr/bin/env python3
"""
Step 1g — marker z-scores, unbiased DE, and the dotplot you actually annotate from.

Three outputs, all per cluster:
  cluster_means_zscore.csv   panel markers, z-scored across clusters
  cluster_means_raw.csv      the same means, un-z-scored  <- check this
  cluster_top_de_genes.csv   top Wilcoxon DE genes, panel-free
plus a dotplot, a z-score heatmap, and an unbiased-DE dotplot.

The raw means matter. A z-score of +3 on a gene whose raw mean is 0.01 is
noise, and reads exactly like a real positive in the heatmap. See
docs/annotation_guide.md.

    # PBMC panel (built in)
    python step1g_marker_evidence.py --file work/step1b_pca.h5ad --panel pbmc

    # microglia / CNS myeloid, from the hierarchical atlas
    python step1g_marker_evidence.py --file work/step1b_pca.h5ad --panel microglia
"""
import argparse

import matplotlib
matplotlib.use("Agg")
import scanpy as sc

from _common import DEFAULT_LEIDEN_KEY, ensure_dir, get_logger, load
from anndata_compiler import (PBMC_MARKER_PANEL, cluster_marker_zscores,
                              load_microglia_panel, unbiased_de)

log = get_logger(__name__)


def resolve_panel(name, species, level):
    if name == "pbmc":
        return PBMC_MARKER_PANEL
    if name == "microglia":
        return load_microglia_panel(species=species, level=level)
    raise SystemExit(f"Unknown panel '{name}'. Use 'pbmc' or 'microglia', or "
                     f"import a dict of {{group: [genes]}} yourself.")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--file", required=True)
    p.add_argument("--panel", default="pbmc", help="'pbmc' or 'microglia'")
    p.add_argument("--species", default="human", choices=["human", "mouse"])
    p.add_argument("--level", default="sub_subtype",
                   choices=["sub_subtype", "umbrella"],
                   help="microglia panel granularity")
    p.add_argument("--cluster-key", default=DEFAULT_LEIDEN_KEY)
    p.add_argument("--out-dir", default="figures")
    p.add_argument("--n-de-genes", type=int, default=10)
    args = p.parse_args()

    adata = load(args.file, log)
    if adata.raw is None:
        raise SystemExit("No .raw — re-run step 1b.")

    # Full gene set, log-normalized: panel markers must not be missing merely
    # because they were not highly variable.
    full = adata.raw.to_adata()
    full.obs = adata.obs.copy()

    ensure_dir(args.out_dir)
    sc.settings.figdir = args.out_dir

    panel = resolve_panel(args.panel, args.species, args.level)
    z, raw, missing = cluster_marker_zscores(full, panel, args.cluster_key)
    z.to_csv(f"{args.out_dir}/cluster_means_zscore.csv")
    raw.to_csv(f"{args.out_dir}/cluster_means_raw.csv")
    log.info(f"Wrote cluster_means_zscore.csv and cluster_means_raw.csv "
             f"({z.shape[0]} clusters x {z.shape[1]} markers)")
    if missing:
        log.warning(f"{len(missing)} panel genes absent from the data: "
                    f"{', '.join(missing)}")

    sizes = (full.obs[args.cluster_key].astype(str).value_counts()
             .rename_axis("cluster").reset_index(name="n_cells"))
    sizes["pct_of_total"] = (sizes["n_cells"] / sizes["n_cells"].sum() * 100).round(2)
    sizes.to_csv(f"{args.out_dir}/cluster_sizes.csv", index=False)

    plot_panel = {g: [m for m in ms if m in full.var_names] for g, ms in panel.items()} \
        if isinstance(panel, dict) else panel
    plot_panel = {g: ms for g, ms in plot_panel.items() if ms}

    sc.pl.dotplot(full, plot_panel, groupby=args.cluster_key, standard_scale="var",
                  show=False, save="_cluster_marker_dotplot.png")
    sc.pl.matrixplot(full, plot_panel, groupby=args.cluster_key, standard_scale="var",
                     show=False, save="_cluster_marker_heatmap_zscore.png")
    log.info(f"Dotplot and heatmap -> {args.out_dir}/")

    de = unbiased_de(full, args.cluster_key, n_genes=args.n_de_genes)
    de.to_csv(f"{args.out_dir}/cluster_top_de_genes.csv", index=False)
    sc.pl.rank_genes_groups_dotplot(full, n_genes=5, groupby=args.cluster_key,
                                    show=False, save="_cluster_top_de_dotplot.png")
    log.info(f"Wrote cluster_top_de_genes.csv")
    log.info("Now open the dotplot. It is the primary evidence; the CSVs support it.")


if __name__ == "__main__":
    main()
