#!/usr/bin/env python3
"""
Step 1h — assemble every evidence source into one annotation template CSV.

Merges the four CellTypist columns (step 1e), the top panel markers by z-score
and the unbiased DE genes (step 1g), and optionally the original authors'
labels, into a single per-cluster table. The three label columns
(`user_label`, `manual_label_final`, `notes`) are left blank on purpose.

    python step1h_build_annotation_template.py \
        --file work/step1b_pca.h5ad \
        --celltypist-dir work/data --evidence-dir figures \
        --out work/data/cluster_annotation_template.csv

Then open the dotplot from step 1g next to this CSV, fill in
`manual_label_final` for every cluster, and run step 1i.
"""
import argparse
import os

import pandas as pd

from _common import DEFAULT_LEIDEN_KEY, get_logger, load
from anndata_compiler import build_annotation_template

log = get_logger(__name__)


def load_or_warn(path, name):
    if not os.path.exists(path):
        log.warning(f"{name} not found at {path} — skipping")
        return None
    df = pd.read_csv(path)
    log.info(f"Loaded {name}: {df.shape}")
    return df


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--file", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--celltypist-dir", default="work/data")
    p.add_argument("--celltypist-prefix", default="celltypist")
    p.add_argument("--evidence-dir", default="figures")
    p.add_argument("--cluster-key", default=DEFAULT_LEIDEN_KEY)
    p.add_argument("--author-labels", default=None,
                   help="CSV with columns cluster,author_label. Carried as a "
                        "comparison column, never adopted automatically.")
    p.add_argument("--n-top-markers", type=int, default=8)
    args = p.parse_args()

    adata = load(args.file, log)

    ct_hvg = load_or_warn(
        f"{args.celltypist_dir}/{args.celltypist_prefix}_cluster_summary.csv",
        "CellTypist summary (HVG)")
    ct_fg = load_or_warn(
        f"{args.celltypist_dir}/{args.celltypist_prefix}_fullgene_cluster_summary.csv",
        "CellTypist summary (full-gene)")

    zscores = load_or_warn(f"{args.evidence_dir}/cluster_means_zscore.csv",
                           "marker z-scores")
    if zscores is not None:
        zscores = zscores.set_index(zscores.columns[0])

    de = load_or_warn(f"{args.evidence_dir}/cluster_top_de_genes.csv",
                      "unbiased DE genes")

    authors = None
    if args.author_labels:
        authors = load_or_warn(args.author_labels, "authors' labels")

    build_annotation_template(
        adata, args.cluster_key, args.out,
        celltypist_hvg=ct_hvg, celltypist_fullgene=ct_fg,
        zscores=zscores, de_genes=de, author_labels=authors,
        n_top_markers=args.n_top_markers,
    )


if __name__ == "__main__":
    main()
