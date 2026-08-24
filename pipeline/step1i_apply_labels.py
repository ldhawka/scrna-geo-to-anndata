#!/usr/bin/env python3
"""
Step 1i — apply the filled-in template back onto the cells.

Fails loudly if any cluster is still unlabelled, rather than silently writing
NaN into obs. Clusters whose `notes` contain 'doublet' or 'artifact' are marked
in obs['cell_type_flag']; pass --drop-flagged to remove them.

    python step1i_apply_labels.py \
        --file work/step1b_pca.h5ad \
        --template work/data/cluster_annotation_template.csv \
        --out work/annotated.h5ad
"""
import argparse

from _common import DEFAULT_LEIDEN_KEY, get_logger, load, save
from anndata_compiler import apply_labels

log = get_logger(__name__)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--file", required=True)
    p.add_argument("--template", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--cluster-key", default=DEFAULT_LEIDEN_KEY)
    p.add_argument("--label-column", default="manual_label_final")
    p.add_argument("--output-key", default="cell_type")
    p.add_argument("--drop-flagged", action="store_true",
                   help="Remove cells in clusters flagged doublet/artifact")
    args = p.parse_args()

    adata = load(args.file, log)
    adata = apply_labels(adata, args.template, args.cluster_key,
                         label_column=args.label_column,
                         output_key=args.output_key,
                         drop_flagged=args.drop_flagged)

    counts = adata.obs[args.output_key].value_counts()
    log.info("Final cell-type counts:")
    for label, n in counts.items():
        log.info(f"  {label}: {n:,} ({n / adata.n_obs:.1%})")

    save(adata, args.out, log)


if __name__ == "__main__":
    main()
