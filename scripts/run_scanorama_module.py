#!/usr/bin/env python3
"""Run Scanorama integration for one pair of .h5ad inputs."""

import argparse
import random
import time

import anndata as ad
import numpy as np
import pandas as pd
import scanorama
import scanpy as sc
import scipy.sparse as sp


def parse_args():
    parser = argparse.ArgumentParser(description="Scanorama integration module")
    parser.add_argument("--input_a", required=True)
    parser.add_argument("--input_b", required=True)
    parser.add_argument("--sample_id", required=True)
    parser.add_argument("--species_a", required=True)
    parser.add_argument("--species_b", required=True)
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
    if len(common_genes) < 50:
        raise ValueError("Insufficient shared genes between inputs for Scanorama")

    adata_a = adata_a[:, adata_a.var_names.isin(common_genes)].copy()
    adata_b = adata_b[:, adata_b.var_names.isin(common_genes)].copy()

    ensure_obs_column(adata_a, "celltype", "unknown")
    ensure_obs_column(adata_b, "celltype", "unknown")

    adata_a.obs["batch"] = f"{args.sample_id}_{args.species_a}"
    adata_b.obs["batch"] = f"{args.sample_id}_{args.species_b}"

    seed = random.randint(0, 2**32 - 1)
    np.random.seed(seed)

    adatas = [adata_a.copy(), adata_b.copy()]
    for item in adatas:
        sc.pp.normalize_total(item)
        sc.pp.log1p(item)
        if sp.issparse(item.X):
            item.X = item.X.tocsr()

    start = time.time()
    corrected = scanorama.correct_scanpy(adatas, return_dimred=True)
    elapsed = time.time() - start

    adata_scanorama = ad.concat(corrected, join="inner", merge="same")
    if "X_scanorama" not in adata_scanorama.obsm:
        raise ValueError("Scanorama did not produce X_scanorama embedding")

    embedding = adata_scanorama.obsm["X_scanorama"]
    adata_scanorama.obsm["X_pca"] = embedding

    pca = pd.DataFrame(
        embedding,
        index=adata_scanorama.obs_names,
        columns=[f"PC{i+1}" for i in range(embedding.shape[1])],
    )
    pca.index.name = "cell"

    pca_file = f"{args.sample_id}_scanorama_embedding.tsv"
    h5ad_file = f"{args.sample_id}_scanorama_integration.h5ad"
    report_file = f"{args.sample_id}_scanorama_report.txt"

    pca.to_csv(pca_file, sep="\t")
    adata_scanorama.write_h5ad(h5ad_file)

    with open(report_file, "w", encoding="utf-8") as handle:
        handle.write(f"sample: {args.sample_id}\n")
        handle.write(f"species_a: {args.species_a}\n")
        handle.write(f"species_b: {args.species_b}\n")
        handle.write(f"cells: {adata_scanorama.n_obs}\n")
        handle.write(f"genes: {adata_scanorama.n_vars}\n")
        handle.write(f"scanorama_dims: {embedding.shape[1]}\n")
        handle.write(f"seed: {seed}\n")
        handle.write(f"elapsed_seconds: {elapsed:.2f}\n")
        handle.write("status: ok\n")


if __name__ == "__main__":
    main()