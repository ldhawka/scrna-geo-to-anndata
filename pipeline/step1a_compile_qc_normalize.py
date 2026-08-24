#!/usr/bin/env python3
"""
Step 1a — compile GEO raw files into one AnnData, detect doublets, QC, normalize.

Wraps `GEOAnndataCompiler`. Output is log1p(CP10K) in .X, integer counts in
layers['counts'], HVGs flagged in var['highly_variable'], and Scrublet scores
in obs.

    python step1a_compile_qc_normalize.py \
        --raw-dir ./GSE123456_RAW \
        --metadata ./metadata.csv \
        --sample-id-column sample_id \
        --out ./work/compiled_qc_hvg.h5ad
"""
import argparse

from _common import get_logger
from anndata_compiler import GEOAnndataCompiler

log = get_logger(__name__)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--raw-dir", required=True)
    p.add_argument("--metadata", required=True)
    p.add_argument("--sample-id-column", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--data-format", default="auto",
                   choices=["auto", "simple", "with_cell_metadata", "10x_mtx"])
    p.add_argument("--min-genes", type=int, default=500)
    p.add_argument("--max-genes", type=int, default=5000,
                   help="Upper bound on genes/cell; 0 disables")
    p.add_argument("--min-cells", type=int, default=3)
    p.add_argument("--max-mito-pct", type=float, default=15.0)
    p.add_argument("--mito-prefix", default="MT-", help="'MT-' human, 'mt-' mouse")
    p.add_argument("--n-top-genes", type=int, default=3000)
    p.add_argument("--expected-doublet-rate", type=float, default=0.06)
    p.add_argument("--no-doublets", action="store_true",
                   help="Skip Scrublet entirely")
    p.add_argument("--filter-doublets", action="store_true",
                   help="Drop predicted doublets instead of only flagging them")
    p.add_argument("--max-cells-per-sample", type=int, default=0,
                   help="Downsample each sample; 0 = keep all")
    p.add_argument("--seed", type=int, default=42)
    args = p.parse_args()

    config = {
        "raw_data_dir": args.raw_dir,
        "metadata_file": args.metadata,
        "output_file": args.out,
        "sample_id_column": args.sample_id_column,
        "data_format": args.data_format,
        "min_genes": args.min_genes,
        "max_genes": args.max_genes or None,
        "min_cells": args.min_cells,
        "max_mito_pct": args.max_mito_pct,
        "mito_prefix": args.mito_prefix,
        "n_top_genes": args.n_top_genes,
        "detect_doublets": not args.no_doublets,
        "expected_doublet_rate": args.expected_doublet_rate,
        "filter_doublets": args.filter_doublets,
        "max_cells_per_sample": args.max_cells_per_sample or None,
        "random_state": args.seed,
    }

    log.info(f"Compiling {args.raw_dir} -> {args.out}")
    adata = GEOAnndataCompiler(config).run_full_pipeline()
    log.info(f"Done: {adata.n_obs:,} cells x {adata.n_vars:,} genes")


if __name__ == "__main__":
    main()
