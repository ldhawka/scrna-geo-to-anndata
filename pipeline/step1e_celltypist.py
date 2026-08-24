#!/usr/bin/env python3
"""
Step 1e — CellTypist on HVG and full-gene matrices, majority-voted per cluster.

Two models (Immune_All_Low, Immune_All_High) x two gene sets (HVG, full-gene)
= four independent label columns. Run it twice:

    # HVG (the working object)
    python step1e_celltypist.py --file work/step1b_pca.h5ad \
        --out-prefix work/data/celltypist --mode hvg

    # full gene set (restored from .raw)
    python step1e_celltypist.py --file work/step1b_pca.h5ad \
        --out-prefix work/data/celltypist --mode fullgene

The full-gene run is usually more confident. Where the two disagree, record the
disagreement — it often means the cluster contains both states.

CellTypist needs log1p(CP10K), never scaled data. Both modes restore that: the
HVG mode subsets .raw to the HVGs rather than reusing the scaled .X.
"""
import argparse

from _common import DEFAULT_LEIDEN_KEY, ensure_dir, get_logger, load
from anndata_compiler import run_celltypist

log = get_logger(__name__)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--file", required=True)
    p.add_argument("--out-prefix", required=True,
                   help="e.g. work/data/celltypist -> *_cluster_summary.csv")
    p.add_argument("--mode", choices=["hvg", "fullgene"], default="hvg")
    p.add_argument("--cluster-key", default=DEFAULT_LEIDEN_KEY)
    p.add_argument("--models", nargs="+",
                   default=["Immune_All_Low.pkl", "Immune_All_High.pkl"])
    args = p.parse_args()

    adata = load(args.file, log)

    if adata.raw is None:
        raise SystemExit("No .raw on this object — re-run step 1b, which stores "
                         "the pre-scale log-normalized matrix in .raw.")

    lognorm = adata.raw.to_adata()
    lognorm.obs = adata.obs.copy()

    if args.mode == "hvg":
        hvg_genes = [g for g in adata.var_names if g in lognorm.var_names]
        lognorm = lognorm[:, hvg_genes].copy()
        prefixes = ["low", "high"]
        suffix = ""
    else:
        prefixes = ["fg_low", "fg_high"]
        suffix = "_fullgene"

    log.info(f"{args.mode}: {lognorm.n_obs:,} x {lognorm.n_vars:,}, "
             f"{lognorm.obs[args.cluster_key].nunique()} clusters")

    per_cell, per_cluster = run_celltypist(
        lognorm, args.cluster_key, models=tuple(args.models), prefixes=tuple(prefixes)
    )

    ensure_dir(args.out_prefix.rsplit("/", 1)[0] if "/" in args.out_prefix else ".")
    cells_path = f"{args.out_prefix}{suffix}_labels.csv"
    clusters_path = f"{args.out_prefix}{suffix}_cluster_summary.csv"
    per_cell.to_csv(cells_path)
    per_cluster.to_csv(clusters_path, index=False)
    log.info(f"Wrote {cells_path}")
    log.info(f"Wrote {clusters_path}")


if __name__ == "__main__":
    main()
