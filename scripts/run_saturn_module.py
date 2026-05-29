#!/usr/bin/env python3
"""Run SATURN integration for one pair of .h5ad inputs.

SATURN (Single-cell Analysis across Taxa Using gene embeddings from pRotein laNguage models)
uses ESM2 protein language model embeddings to create a gene-level representation that is
comparable across species without explicit ortholog mapping. It showed the best or near-best
overall performance across most taxonomic levels in Zhong et al. 2025.

ESM2 protein embeddings are generated automatically from gene symbols using the ESM2
650M-parameter model (esm2_t33_650M_UR50D). This requires ~2 GB RAM and is slow on CPU;
GPU is strongly recommended.
"""

import argparse
import time
import warnings

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp


def parse_args():
    parser = argparse.ArgumentParser(description="SATURN integration module")
    parser.add_argument("--input_a", required=True)
    parser.add_argument("--input_b", required=True)
    parser.add_argument("--sample_id", required=True)
    parser.add_argument("--species_a", required=True)
    parser.add_argument("--species_b", required=True)
    parser.add_argument(
        "--n_top_genes", type=int, default=3000,
        help="Highly variable genes per species for SATURN.",
    )
    parser.add_argument(
        "--max_epochs", type=int, default=50,
        help="Maximum training epochs for SATURN.",
    )
    parser.add_argument(
        "--esm_model", default="esm2_t33_650M_UR50D",
        help="ESM2 model identifier for protein embeddings.",
    )
    return parser.parse_args()


def ensure_obs_column(adata, key, default_value):
    if key not in adata.obs.columns:
        adata.obs[key] = default_value


def _get_esm_embeddings(gene_names: list[str], esm_model: str) -> np.ndarray:
    """Generate ESM2 embeddings for a list of gene symbol strings.

    Uses gene symbols as pseudo-sequences: ESM2 is applied to the gene name string
    as a proxy for the protein sequence. For production use, replace with actual
    protein sequences retrieved from UniProt / Ensembl.
    """
    try:
        import esm
    except ImportError as exc:
        raise ImportError(
            "fair-esm is not installed. Install with: pip install fair-esm"
        ) from exc

    import torch

    model, alphabet = esm.pretrained.load_model_and_alphabet(esm_model)
    model.eval()
    batch_converter = alphabet.get_batch_converter()

    device = "cuda" if torch.cuda.is_available() else "cpu"
    model = model.to(device)

    embeddings = []
    batch_size = 32

    for i in range(0, len(gene_names), batch_size):
        batch_genes = gene_names[i : i + batch_size]
        # Use gene name as sequence proxy — replace with real AA sequence for production
        batch_data = [(g, g) for g in batch_genes]
        _, _, batch_tokens = batch_converter(batch_data)
        batch_tokens = batch_tokens.to(device)

        with torch.no_grad():
            results = model(batch_tokens, repr_layers=[33], return_contacts=False)

        # Mean-pool over sequence positions for each gene
        token_reps = results["representations"][33]
        for j, tokens in enumerate(batch_tokens):
            seq_len = (tokens != alphabet.padding_idx).sum()
            emb = token_reps[j, 1 : seq_len - 1].mean(0).cpu().numpy()
            embeddings.append(emb)

    return np.stack(embeddings)


def _preprocess(adata: ad.AnnData, n_top_genes: int, species_label: str, sample_id: str) -> ad.AnnData:
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=10)
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    n_hvg = min(n_top_genes, adata.n_vars)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg)
    ensure_obs_column(adata, "celltype", "unknown")
    adata.obs["batch"] = f"{sample_id}_{species_label}"
    return adata


def main():
    args = parse_args()

    try:
        import saturn
        from saturn.model import SATURN as SATURNModel
    except ImportError as exc:
        raise ImportError(
            "saturn-sc is not installed. Install with: pip install saturn-sc"
        ) from exc

    adata_a = sc.read_h5ad(args.input_a)
    adata_b = sc.read_h5ad(args.input_b)

    adata_a = _preprocess(adata_a, args.n_top_genes, args.species_a, args.sample_id)
    adata_b = _preprocess(adata_b, args.n_top_genes, args.species_b, args.sample_id)

    # Generate ESM2 gene embeddings per species
    print(f"Generating ESM2 embeddings for {adata_a.n_vars} genes ({args.species_a})...")
    emb_a = _get_esm_embeddings(list(adata_a.var_names), args.esm_model)

    print(f"Generating ESM2 embeddings for {adata_b.n_vars} genes ({args.species_b})...")
    emb_b = _get_esm_embeddings(list(adata_b.var_names), args.esm_model)

    # Attach embeddings to var DataFrames (SATURN expects them in adata.varm or as input matrix)
    adata_a.varm["gene_emb"] = emb_a
    adata_b.varm["gene_emb"] = emb_b

    start = time.time()

    # Prepare count matrices (SATURN expects raw-ish counts in X)
    X_a = adata_a.X.toarray() if sp.issparse(adata_a.X) else np.asarray(adata_a.X)
    X_b = adata_b.X.toarray() if sp.issparse(adata_b.X) else np.asarray(adata_b.X)

    model = SATURNModel(
        adatas=[adata_a, adata_b],
        gene_embeddings=[emb_a, emb_b],
    )
    model.train(max_epochs=args.max_epochs)

    latent = model.get_latent_representation()
    elapsed = time.time() - start

    # Reconstruct a combined AnnData with the latent space
    adata_all = ad.concat([adata_a, adata_b], join="outer", merge="same")
    adata_all.obs_names_make_unique()
    adata_all.obsm["X_saturn"] = latent

    pca = pd.DataFrame(
        latent,
        index=adata_all.obs_names,
        columns=[f"PC{i+1}" for i in range(latent.shape[1])],
    )
    pca.index.name = "cell"

    pca_file    = f"{args.sample_id}_saturn_embedding.tsv"
    h5ad_file   = f"{args.sample_id}_saturn_integration.h5ad"
    report_file = f"{args.sample_id}_saturn_report.txt"

    pca.to_csv(pca_file, sep="\t")
    adata_all.write_h5ad(h5ad_file)

    with open(report_file, "w", encoding="utf-8") as fh:
        fh.write(f"sample: {args.sample_id}\n")
        fh.write(f"species_a: {args.species_a}\n")
        fh.write(f"species_b: {args.species_b}\n")
        fh.write(f"cells: {adata_all.n_obs}\n")
        fh.write(f"genes_a: {adata_a.n_vars}\n")
        fh.write(f"genes_b: {adata_b.n_vars}\n")
        fh.write(f"latent_dims: {latent.shape[1]}\n")
        fh.write(f"esm_model: {args.esm_model}\n")
        fh.write(f"elapsed_seconds: {elapsed:.2f}\n")
        fh.write("status: ok\n")


if __name__ == "__main__":
    main()
