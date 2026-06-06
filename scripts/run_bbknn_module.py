#!/usr/bin/env python3
"""Run BBKNN integration for one pair of .h5ad inputs."""

import argparse
import os
import random
import sys
import time

import anndata as ad
import bbknn
import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from feature_selection import load_features_file
from normalization import NORM_METHODS, apply_normalization


def parse_args():
    parser = argparse.ArgumentParser(description="BBKNN integration module")
    parser.add_argument("--input_a", required=True)
    parser.add_argument("--input_b", required=True)
    parser.add_argument("--sample_id", required=True)
    parser.add_argument("--species_a", required=True)
    parser.add_argument("--species_b", required=True)
    parser.add_argument(
        "--normalization", default="log_norm", choices=NORM_METHODS,
        help="Normalization applied before PCA/integration.",
    )
    parser.add_argument(
        "--features_file", default=None,
        help="Optional path to a gene list (one per line) from run_feature_selection.py. "
             "When provided, integration is restricted to those genes.",
    )
    return parser.parse_args()


def ensure_obs_column(adata, key, default_value):
    if key not in adata.obs.columns:
        adata.obs[key] = default_value


def main():
    args = parse_args()

    adata_a = sc.read_h5ad(args.input_a)
    adata_b = sc.read_h5ad(args.input_b)

    # Keep benchmark defaults aligned with existing scripts.
    sc.pp.filter_cells(adata_a, min_genes=200)
    sc.pp.filter_genes(adata_a, min_cells=10)
    sc.pp.filter_cells(adata_b, min_genes=200)
    sc.pp.filter_genes(adata_b, min_cells=10)

    common_genes = sorted(set(adata_a.var_names).intersection(adata_b.var_names))
    if args.features_file:
        selected = set(load_features_file(args.features_file))
        common_genes = [g for g in common_genes if g in selected]
    if len(common_genes) < 50:
        raise ValueError("Insufficient shared genes between inputs for BBKNN")

    adata_a = adata_a[:, adata_a.var_names.isin(common_genes)].copy()
    adata_b = adata_b[:, adata_b.var_names.isin(common_genes)].copy()

    ensure_obs_column(adata_a, "celltype", "unknown")
    ensure_obs_column(adata_b, "celltype", "unknown")

    adata_a.obs["batch"] = f"{args.sample_id}_{args.species_a}"
    adata_b.obs["batch"] = f"{args.sample_id}_{args.species_b}"

    adata_all = ad.concat([adata_a, adata_b], join="inner", merge="same")

    seed = random.randint(0, 2**32 - 1)
    np.random.seed(seed)

    start = time.time()
    apply_normalization(adata_all, args.normalization)
    sc.pp.scale(adata_all)
    sc.tl.pca(adata_all, svd_solver="arpack", n_comps=30)

    adata_bbknn = bbknn.bbknn(
        adata_all,
        batch_key="batch",
        copy=True,
        neighbors_within_batch=5,
        approx=False,
        trim=50,
    )
    sc.tl.pca(adata_bbknn, svd_solver="arpack", n_comps=30)
    adata_bbknn.obsm["X_pca"] *= -1
    elapsed = time.time() - start

    pca = pd.DataFrame(
        adata_bbknn.obsm["X_pca"],
        index=adata_bbknn.obs_names,
        columns=[f"PC{i+1}" for i in range(adata_bbknn.obsm["X_pca"].shape[1])],
    )
    pca.index.name = "cell"

    pca_file = f"{args.sample_id}_bbknn_embedding.tsv"
    h5ad_file = f"{args.sample_id}_bbknn_integration.h5ad"
    report_file = f"{args.sample_id}_bbknn_report.txt"

    pca.to_csv(pca_file, sep="\t")
    adata_bbknn.write_h5ad(h5ad_file)

    with open(report_file, "w", encoding="utf-8") as handle:
        handle.write(f"sample: {args.sample_id}\n")
        handle.write(f"species_a: {args.species_a}\n")
        handle.write(f"species_b: {args.species_b}\n")
        handle.write(f"cells: {adata_bbknn.n_obs}\n")
        handle.write(f"genes: {adata_bbknn.n_vars}\n")
        handle.write(f"seed: {seed}\n")
        handle.write(f"normalization: {args.normalization}\n")
        handle.write(f"features_file: {args.features_file or 'none'}\n")
        handle.write(f"n_genes_used: {adata_bbknn.n_vars}\n")
        handle.write(f"elapsed_seconds: {elapsed:.2f}\n")
        handle.write("status: ok\n")


if __name__ == "__main__":
    main()
