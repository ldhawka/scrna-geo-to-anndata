#!/usr/bin/env python3
"""
Step 1b — subset to HVGs, scale, run PCA.

Pre-scale values are kept in .raw so that CellTypist (step 1e) and the marker
z-scores (step 1g) can still see log-normalized expression for all genes.

    python step1b_scale_pca.py --in work/compiled_qc_hvg.h5ad --out work/step1b_pca.h5ad
"""
import argparse

import matplotlib
matplotlib.use("Agg")
import scanpy as sc

from _common import ensure_dir, get_logger, load, save
from anndata_compiler import scale_and_pca

log = get_logger(__name__)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--in", dest="inp", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--n-comps", type=int, default=50)
    p.add_argument("--max-value", type=float, default=10)
    p.add_argument("--fig-dir", default="figures")
    p.add_argument("--seed", type=int, default=42)
    args = p.parse_args()

    adata = load(args.inp, log)
    adata = scale_and_pca(adata, n_comps=args.n_comps,
                          max_value=args.max_value, random_state=args.seed)

    # Scree plot — the elbow sets --n-pcs for steps 1c and 1d.
    ensure_dir(args.fig_dir)
    sc.settings.figdir = args.fig_dir
    sc.pl.pca_variance_ratio(adata, n_pcs=args.n_comps, log=True,
                             show=False, save="_scree.png")
    log.info(f"Scree plot -> {args.fig_dir}/pca_variance_ratio_scree.png "
             f"(read the elbow, pass it as --n-pcs next)")

    save(adata, args.out, log)


if __name__ == "__main__":
    main()
