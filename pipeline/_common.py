"""Shared helpers for the numbered pipeline steps."""

import logging
import os
import sys

# Make the repo importable when the scripts are run directly from a clone.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

DEFAULT_LEIDEN_KEY = "leiden_r1_0"


def get_logger(name):
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.StreamHandler(sys.stdout)],
        force=True,
    )
    return logging.getLogger(name)


def load(path, log):
    import scanpy as sc
    log.info(f"Reading {path}")
    adata = sc.read_h5ad(path)
    log.info(f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes; "
             f"obsm={list(adata.obsm.keys())}")
    return adata


def save(adata, path, log):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    adata.write(path)
    log.info(f"Saved -> {path}")


def ensure_dir(path):
    os.makedirs(os.path.abspath(path), exist_ok=True)
    return path
