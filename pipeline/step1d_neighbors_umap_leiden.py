#!/usr/bin/env python3
"""
Step 1d — neighbourhood graph, UMAP, Leiden clustering.

Clusters on obsm['X_pca_harmony'] by default. Clustering on 'X_pca' after
running Harmony is the single most common silent mistake here — the clusters
come back uncorrected.

    python step1d_neighbors_umap_leiden.py --file work/step1b_pca.h5ad --resolution 1.0
"""
import argparse

import matplotlib
matplotlib.use("Agg")
import scanpy as sc

from _common import ensure_dir, get_logger, load, save
from anndata_compiler import cluster

log = get_logger(__name__)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--file", required=True, help="h5ad from step 1c; updated in place")
    p.add_argument("--use-rep", default="X_pca_harmony")
    p.add_argument("--n-neighbors", type=int, default=15)
    p.add_argument("--n-pcs", type=int, default=30)
    p.add_argument("--resolution", type=float, nargs="+", default=[1.0])
    p.add_argument("--color", nargs="*", default=["sample_id"],
                   help="Extra obs columns to colour the UMAP by")
    p.add_argument("--fig-dir", default="figures")
    p.add_argument("--seed", type=int, default=42)
    args = p.parse_args()

    adata = load(args.file, log)
    adata = cluster(adata, use_rep=args.use_rep, n_neighbors=args.n_neighbors,
                    n_pcs=args.n_pcs, resolutions=args.resolution,
                    random_state=args.seed)

    ensure_dir(args.fig_dir)
    sc.settings.figdir = args.fig_dir
    leiden_keys = [f"leiden_r{str(r).replace('.', '_')}" for r in args.resolution]
    colors = leiden_keys + [c for c in args.color if c in adata.obs.columns]
    sc.pl.umap(adata, color=colors, ncols=2, show=False, save="_clusters.png")
    log.info(f"UMAP -> {args.fig_dir}/umap_clusters.png")
    log.info("Check the sample-coloured panel: if clusters are one sample each, "
             "revisit step 1c.")

    save(adata, args.file, log)


if __name__ == "__main__":
    main()
