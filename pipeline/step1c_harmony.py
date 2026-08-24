#!/usr/bin/env python3
"""
Step 1c — Harmony batch integration on the PCA embedding.

Writes obsm['X_pca_harmony'] back into the same file. Skip this step only if
the data is genuinely a single library; multi-sample GEO compilations almost
always need it.

    python step1c_harmony.py --file work/step1b_pca.h5ad --batch-key sample_id
"""
import argparse

from _common import get_logger, load, save
from anndata_compiler import harmony

log = get_logger(__name__)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--file", required=True, help="h5ad from step 1b; updated in place")
    p.add_argument("--batch-key", default="sample_id")
    p.add_argument("--n-pcs", type=int, default=30)
    p.add_argument("--max-iter", type=int, default=20)
    p.add_argument("--epsilon", type=float, default=1e-5)
    p.add_argument("--seed", type=int, default=42)
    args = p.parse_args()

    adata = load(args.file, log)
    adata = harmony(adata, batch_key=args.batch_key, n_pcs=args.n_pcs,
                    max_iter_harmony=args.max_iter, epsilon_harmony=args.epsilon,
                    random_state=args.seed)
    save(adata, args.file, log)


if __name__ == "__main__":
    main()
