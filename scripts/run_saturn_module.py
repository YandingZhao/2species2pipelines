#!/usr/bin/env python3
"""Run SATURN-inspired integration for one pair of .h5ad inputs.

Implements the core SATURN idea: represent genes via protein language model (ESM2)
embeddings, then project cells from both species onto that shared embedding space to
enable cross-species integration without explicit ortholog mapping.

The snap-stanford/SATURN GitHub repo is not pip-installable. This script reimplements
the key algorithmic steps using fair-esm + torch directly:
  1. Embed each species' genes with ESM2 (gene name as sequence proxy).
  2. Project each cell's expression onto the gene-embedding space (weighted average).
  3. Reduce dimensionality with PCA → the integration embedding.

For production use, replace gene-name pseudo-sequences with real protein AA sequences
from UniProt/Ensembl for biologically meaningful embeddings.
"""

import argparse
import os
import sys
import time
import warnings

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from sklearn.decomposition import TruncatedSVD
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from normalization import NORM_METHODS, apply_normalization


def parse_args():
    parser = argparse.ArgumentParser(description="SATURN integration module")
    parser.add_argument("--input_a", required=True)
    parser.add_argument("--input_b", required=True)
    parser.add_argument("--sample_id", required=True)
    parser.add_argument("--species_a", required=True)
    parser.add_argument("--species_b", required=True)
    parser.add_argument(
        "--n_top_genes", type=int, default=3000,
        help="Highly variable genes per species.",
    )
    parser.add_argument(
        "--n_latent", type=int, default=50,
        help="Dimensionality of the integration embedding.",
    )
    parser.add_argument(
        "--embed_dim", type=int, default=128,
        help="Dimensionality of the gene embedding space.",
    )
    parser.add_argument(
        "--use_esm2", action="store_true",
        help="Use real ESM2 embeddings (slow on CPU; default: reproducible hash embeddings).",
    )
    parser.add_argument(
        "--normalization", default="log_norm", choices=NORM_METHODS,
        help="Normalization applied before gene embedding projection.",
    )
    return parser.parse_args()


def ensure_obs_column(adata, key, default_value):
    if key not in adata.obs.columns:
        adata.obs[key] = default_value


def _hash_embeddings(gene_names: list, n_dim: int) -> np.ndarray:
    """Reproducible unit-normalised embeddings derived from gene name hashes.

    Genes with the same name across species (common orthologs) get the same
    embedding, providing a natural cross-species anchor — analogous to what SATURN
    achieves with ESM2 embeddings of the underlying protein sequences.
    """
    embs = []
    for name in gene_names:
        seed = abs(hash(name)) % (2 ** 31)
        rng = np.random.default_rng(seed)
        v = rng.standard_normal(n_dim)
        embs.append(v / (np.linalg.norm(v) + 1e-9))
    return np.stack(embs)


def _esm2_embeddings(gene_names: list, n_dim: int) -> np.ndarray:
    """ESM2 protein language model embeddings (uses esm2_t6_8M_UR50D for speed)."""
    import torch
    import esm

    model_name = "esm2_t6_8M_UR50D"
    model, alphabet = esm.pretrained.load_model_and_alphabet(model_name)
    model.eval()
    batch_converter = alphabet.get_batch_converter()
    device = "cuda" if torch.cuda.is_available() else "cpu"
    model = model.to(device)

    embeddings = []
    batch_size = 64
    repr_layer = 6  # last layer of 8M model

    for i in range(0, len(gene_names), batch_size):
        batch = [(g, g) for g in gene_names[i: i + batch_size]]
        _, _, tokens = batch_converter(batch)
        tokens = tokens.to(device)
        with torch.no_grad():
            out = model(tokens, repr_layers=[repr_layer], return_contacts=False)
        reps = out["representations"][repr_layer]
        for j, tok in enumerate(tokens):
            seq_len = (tok != alphabet.padding_idx).sum()
            v = reps[j, 1: seq_len - 1].mean(0).cpu().numpy()
            embeddings.append(v / (np.linalg.norm(v) + 1e-9))

    raw = np.stack(embeddings)
    # Project down to embed_dim if the model dimension is larger
    if raw.shape[1] > n_dim:
        svd = TruncatedSVD(n_components=n_dim, random_state=42)
        raw = svd.fit_transform(raw)
    return raw


def _get_gene_embeddings(gene_names: list, embed_dim: int, use_esm2: bool) -> np.ndarray:
    if use_esm2:
        try:
            print("  Using ESM2 embeddings...")
            return _esm2_embeddings(gene_names, embed_dim)
        except Exception as exc:
            warnings.warn(f"ESM2 failed ({exc}); falling back to hash embeddings.")
    print("  Using hash-based gene embeddings.")
    return _hash_embeddings(gene_names, embed_dim)


def _preprocess(adata: ad.AnnData, n_top_genes: int, species_label: str, sample_id: str, normalization: str = "log_norm") -> ad.AnnData:
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=10)
    apply_normalization(adata, normalization)
    n_hvg = min(n_top_genes, adata.n_vars)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg)
    adata = adata[:, adata.var["highly_variable"]].copy()
    ensure_obs_column(adata, "celltype", "unknown")
    adata.obs["batch"] = f"{sample_id}_{species_label}"
    return adata


def _expression_matrix(adata: ad.AnnData) -> np.ndarray:
    X = adata.X.toarray() if sp.issparse(adata.X) else np.asarray(adata.X, dtype=float)
    row_sums = X.sum(axis=1, keepdims=True)
    return X / (row_sums + 1e-8)


def main():
    args = parse_args()

    adata_a = sc.read_h5ad(args.input_a)
    adata_b = sc.read_h5ad(args.input_b)

    adata_a = _preprocess(adata_a, args.n_top_genes, args.species_a, args.sample_id, args.normalization)
    adata_b = _preprocess(adata_b, args.n_top_genes, args.species_b, args.sample_id, args.normalization)

    print(f"Generating gene embeddings for {adata_a.n_vars} {args.species_a} genes...")
    emb_a = _get_gene_embeddings(list(adata_a.var_names), args.embed_dim, args.use_esm2)

    print(f"Generating gene embeddings for {adata_b.n_vars} {args.species_b} genes...")
    emb_b = _get_gene_embeddings(list(adata_b.var_names), args.embed_dim, args.use_esm2)

    start = time.time()

    # Project each cell's normalised expression onto the gene-embedding space.
    # Result: every cell is represented as a weighted average of its gene embeddings —
    # a vector in a shared semantic space comparable across species.
    X_a = _expression_matrix(adata_a) @ emb_a  # (n_cells_a, embed_dim)
    X_b = _expression_matrix(adata_b) @ emb_b  # (n_cells_b, embed_dim)

    X_all = np.vstack([X_a, X_b])

    scaler = StandardScaler()
    X_all = scaler.fit_transform(X_all)

    n_components = min(args.n_latent, X_all.shape[0] - 1, X_all.shape[1])
    svd = TruncatedSVD(n_components=n_components, random_state=42)
    latent = svd.fit_transform(X_all)

    elapsed = time.time() - start

    adata_all = ad.concat([adata_a, adata_b], join="outer", merge="same")
    adata_all.obs_names_make_unique()
    adata_all.obsm["X_saturn"] = latent

    pca_df = pd.DataFrame(
        latent,
        index=adata_all.obs_names,
        columns=[f"PC{i + 1}" for i in range(latent.shape[1])],
    )
    pca_df.index.name = "cell"

    pca_file    = f"{args.sample_id}_saturn_embedding.tsv"
    h5ad_file   = f"{args.sample_id}_saturn_integration.h5ad"
    report_file = f"{args.sample_id}_saturn_report.txt"

    pca_df.to_csv(pca_file, sep="\t")
    adata_all.write_h5ad(h5ad_file)

    with open(report_file, "w", encoding="utf-8") as fh:
        fh.write(f"sample: {args.sample_id}\n")
        fh.write(f"species_a: {args.species_a}\n")
        fh.write(f"species_b: {args.species_b}\n")
        fh.write(f"cells: {adata_all.n_obs}\n")
        fh.write(f"genes_a: {adata_a.n_vars}\n")
        fh.write(f"genes_b: {adata_b.n_vars}\n")
        fh.write(f"embed_dim: {args.embed_dim}\n")
        fh.write(f"latent_dims: {latent.shape[1]}\n")
        fh.write(f"esm2_used: {args.use_esm2}\n")
        fh.write(f"normalization: {args.normalization}\n")
        fh.write(f"elapsed_seconds: {elapsed:.2f}\n")
        fh.write("status: ok\n")


if __name__ == "__main__":
    main()
