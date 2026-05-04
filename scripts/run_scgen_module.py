#!/usr/bin/env python3
"""Run scGen integration for one pair of .h5ad inputs."""

import argparse
import random
import time

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import scgen
import scvi


def parse_args():
    parser = argparse.ArgumentParser(description="scGen integration module")
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

    adata_a.obs_names_make_unique()
    adata_b.obs_names_make_unique()

    # Keep benchmark defaults aligned with existing scripts.
    sc.pp.filter_cells(adata_a, min_genes=200)
    sc.pp.filter_genes(adata_a, min_cells=10)
    sc.pp.filter_cells(adata_b, min_genes=200)
    sc.pp.filter_genes(adata_b, min_cells=10)

    common_genes = sorted(set(adata_a.var_names).intersection(adata_b.var_names))
    if len(common_genes) < 50:
        raise ValueError("Insufficient shared genes between inputs for scGen")

    adata_a = adata_a[:, adata_a.var_names.isin(common_genes)].copy()
    adata_b = adata_b[:, adata_b.var_names.isin(common_genes)].copy()

    ensure_obs_column(adata_a, "celltype", "unknown")
    ensure_obs_column(adata_b, "celltype", "unknown")

    adata_a.obs["batch"] = f"{args.sample_id}_{args.species_a}"
    adata_b.obs["batch"] = f"{args.sample_id}_{args.species_b}"

    adata_all = ad.concat([adata_a, adata_b], join="inner", merge="same")
    adata_all.obs_names_make_unique()

    if sp.issparse(adata_all.X):
        adata_all.X = adata_all.X.tocsr()
    else:
        adata_all.X = np.asarray(adata_all.X)

    seed = random.randint(0, 2**32 - 1)
    random.seed(seed)
    np.random.seed(seed)
    scvi.settings.seed = seed

    # Normalise and log-transform before scGen training (matches benchmark).
    sc.pp.normalize_total(adata_all)
    sc.pp.log1p(adata_all)

    # scGen requires explicit celltype label; use "celltype" obs column.
    adata_all.obs["celltype"] = adata_all.obs["celltype"].astype(str)
    adata_all.obs["batch"] = adata_all.obs["batch"].astype(str)

    scgen.SCGEN.setup_anndata(adata_all, batch_key="batch", labels_key="celltype")

    start = time.time()
    model = scgen.SCGEN(adata_all)
    model.train(
        max_epochs=100,
        batch_size=32,
        early_stopping=True,
        early_stopping_patience=25,
    )

    corrected_adata = model.batch_removal()
    elapsed = time.time() - start

    latent = corrected_adata.obsm["corrected_latent"]

    pca = pd.DataFrame(
        latent,
        index=corrected_adata.obs_names,
        columns=[f"PC{i+1}" for i in range(latent.shape[1])],
    )
    pca.index.name = "cell"

    # Store embedding so downstream evaluation can find it.
    corrected_adata.obsm["X_pca"] = latent
    corrected_adata.obs["orig.ident"] = corrected_adata.obs["batch"]

    pca_file = f"{args.sample_id}_scgen_embedding.tsv"
    h5ad_file = f"{args.sample_id}_scgen_integration.h5ad"
    report_file = f"{args.sample_id}_scgen_report.txt"

    pca.to_csv(pca_file, sep="\t")
    corrected_adata.write_h5ad(h5ad_file)

    with open(report_file, "w", encoding="utf-8") as handle:
        handle.write(f"sample: {args.sample_id}\n")
        handle.write(f"species_a: {args.species_a}\n")
        handle.write(f"species_b: {args.species_b}\n")
        handle.write(f"cells: {corrected_adata.n_obs}\n")
        handle.write(f"genes: {corrected_adata.n_vars}\n")
        handle.write(f"latent_dims: {latent.shape[1]}\n")
        handle.write(f"seed: {seed}\n")
        handle.write(f"elapsed_seconds: {elapsed:.2f}\n")
        handle.write("status: ok\n")


if __name__ == "__main__":
    main()
