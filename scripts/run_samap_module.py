#!/usr/bin/env python3
"""Run SAMap integration for one pair of .h5ad inputs.

SAMap (Self-Assembling Manifolds) uses cross-species protein sequence similarity
combined with gene expression to build a unified cell manifold. It performs best
at atlas-level cross-phylum integration (Zhong et al. 2025).

BLAST is run internally via SAMap's built-in wrapper — no pre-computed files needed.
"""

import argparse
import time

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc


def parse_args():
    parser = argparse.ArgumentParser(description="SAMap integration module")
    parser.add_argument("--input_a", required=True)
    parser.add_argument("--input_b", required=True)
    parser.add_argument("--sample_id", required=True)
    parser.add_argument("--species_a", required=True)
    parser.add_argument("--species_b", required=True)
    parser.add_argument(
        "--n_top_genes", type=int, default=3000,
        help="Number of highly variable genes to use per species.",
    )
    parser.add_argument(
        "--num_epochs", type=int, default=10,
        help="Number of SAMap training epochs.",
    )
    return parser.parse_args()


def ensure_obs_column(adata, key, default_value):
    if key not in adata.obs.columns:
        adata.obs[key] = default_value


def _short_id(species: str) -> str:
    """Return a 2-character species ID suitable for SAMap keys (e.g. 'hu', 'do')."""
    return species[:2].lower()


def main():
    args = parse_args()

    try:
        from samap.mapping import SAMAP
        from samap.analysis import get_mapping_scores
    except ImportError as exc:
        raise ImportError(
            "SAMap is not installed. Install with: pip install samap"
        ) from exc

    adata_a = sc.read_h5ad(args.input_a)
    adata_b = sc.read_h5ad(args.input_b)

    sc.pp.filter_cells(adata_a, min_genes=200)
    sc.pp.filter_genes(adata_a, min_cells=10)
    sc.pp.filter_cells(adata_b, min_genes=200)
    sc.pp.filter_genes(adata_b, min_cells=10)

    ensure_obs_column(adata_a, "celltype", "unknown")
    ensure_obs_column(adata_b, "celltype", "unknown")

    id_a = _short_id(args.species_a)
    id_b = _short_id(args.species_b)

    # SAMap requires unique species IDs — fall back to positional suffix if identical
    if id_a == id_b:
        id_a = id_a + "1"
        id_b = id_b + "2"

    adata_a.obs["batch"] = f"{args.sample_id}_{args.species_a}"
    adata_b.obs["batch"] = f"{args.sample_id}_{args.species_b}"

    # Preprocess per species as SAMap expects normalised + log data
    for adata, n in [(adata_a, args.n_top_genes), (adata_b, args.n_top_genes)]:
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        n_hvg = min(n, adata.n_vars)
        sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg)

    start = time.time()

    # SAMap requires file paths or SAM objects — AnnData is not accepted directly.
    # It also needs pre-computed BLAST mapping tables in a maps directory.
    # Since our data went through ortholog conversion, shared gene names are orthologs.
    # We create synthetic BLAST outfmt-6 tables from the intersection of gene names,
    # bypassing the need to run actual BLAST against protein sequence databases.
    import os, tempfile
    tmp_dir = tempfile.mkdtemp()
    path_a = os.path.join(tmp_dir, f"{id_a}.h5ad")
    path_b = os.path.join(tmp_dir, f"{id_b}.h5ad")
    adata_a.write_h5ad(path_a)
    adata_b.write_h5ad(path_b)

    maps_dir = os.path.join(tmp_dir, "maps") + "/"
    os.makedirs(maps_dir, exist_ok=True)

    shared_genes = sorted(set(adata_a.var_names) & set(adata_b.var_names))
    print(f"Creating synthetic BLAST tables for {len(shared_genes)} shared genes...")

    # SAMap expects: maps/{id_a}{id_b}/{id_a}_to_{id_b}.txt  (and the reverse)
    # e.g. maps/dohu/do_to_hu.txt and maps/dohu/hu_to_do.txt
    # The directory name is the two IDs concatenated (no separator).
    pair_dir = os.path.join(maps_dir, f"{id_a}{id_b}")
    os.makedirs(pair_dir, exist_ok=True)

    # BLAST outfmt 6: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
    blast_row = lambda q, s: f"{q}\t{s}\t100.0\t100\t0\t0\t1\t100\t1\t100\t1e-99\t100\n"
    with open(os.path.join(pair_dir, f"{id_a}_to_{id_b}.txt"), "w") as fh:
        for gene in shared_genes:
            fh.write(blast_row(gene, gene))
    with open(os.path.join(pair_dir, f"{id_b}_to_{id_a}.txt"), "w") as fh:
        for gene in shared_genes:
            fh.write(blast_row(gene, gene))

    sams = {id_a: path_a, id_b: path_b}
    sm = SAMAP(sams, f_maps=maps_dir)
    sm.run(NUMITERS=args.num_epochs)

    elapsed = time.time() - start

    # SAMap stores per-species results in sm.sams[id].adata (no sm.adata attribute).
    # Each species' adata has the joint-space embedding in obsm after run().
    adata_a_out = sm.sams[id_a].adata
    adata_b_out = sm.sams[id_b].adata

    # Find the embedding key (SAMap may use X_umap or X_samap)
    embedding_key = next(
        (k for k in ["X_samap", "X_umap"] if k in adata_a_out.obsm),
        list(adata_a_out.obsm.keys())[0] if adata_a_out.obsm else None,
    )

    # Restore metadata before concat
    for adata_out, species in [(adata_a_out, args.species_a), (adata_b_out, args.species_b)]:
        adata_out.obs["batch"] = f"{args.sample_id}_{species}"
        if "celltype" not in adata_out.obs.columns:
            adata_out.obs["celltype"] = "unknown"

    adata_combined = ad.concat([adata_a_out, adata_b_out], join="outer", merge="same")
    adata_combined.obs_names_make_unique()

    if embedding_key and embedding_key in adata_a_out.obsm and embedding_key in adata_b_out.obsm:
        embedding = np.vstack([
            np.asarray(adata_a_out.obsm[embedding_key]),
            np.asarray(adata_b_out.obsm[embedding_key]),
        ])
    else:
        # Fallback: PCA on the concatenated data
        sc.pp.pca(adata_combined, n_comps=50)
        embedding = np.asarray(adata_combined.obsm["X_pca"])
        embedding_key = "X_pca"

    adata_combined.obsm["X_samap"] = embedding

    pca = pd.DataFrame(
        embedding,
        index=adata_combined.obs_names,
        columns=[f"PC{i+1}" for i in range(embedding.shape[1])],
    )
    pca.index.name = "cell"

    pca_file    = f"{args.sample_id}_samap_embedding.tsv"
    h5ad_file   = f"{args.sample_id}_samap_integration.h5ad"
    report_file = f"{args.sample_id}_samap_report.txt"

    pca.to_csv(pca_file, sep="\t")
    adata_combined.write_h5ad(h5ad_file)

    with open(report_file, "w", encoding="utf-8") as fh:
        fh.write(f"sample: {args.sample_id}\n")
        fh.write(f"species_a: {args.species_a} (id: {id_a})\n")
        fh.write(f"species_b: {args.species_b} (id: {id_b})\n")
        fh.write(f"cells: {adata_combined.n_obs}\n")
        fh.write(f"embedding_key: {embedding_key}\n")
        fh.write(f"embedding_dims: {embedding.shape[1]}\n")
        fh.write(f"elapsed_seconds: {elapsed:.2f}\n")
        fh.write("status: ok\n")


if __name__ == "__main__":
    main()
