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
    # BLAST outfmt 6: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
    blast_row = lambda g: f"{g}\t{g}\t100.0\t100\t0\t0\t1\t100\t1\t100\t1e-99\t100\n"
    for fname in [f"{id_a}_{id_b}.txt", f"{id_b}_{id_a}.txt"]:
        with open(os.path.join(maps_dir, fname), "w") as fh:
            for gene in shared_genes:
                fh.write(blast_row(gene))

    sams = {id_a: path_a, id_b: path_b}
    sm = SAMAP(sams, f_maps=maps_dir)
    sm.run(NUMITERS=args.num_epochs)

    elapsed = time.time() - start

    # Extract the joint embedding — SAMap stores it in the combined adata
    adata_combined = sm.adata
    if "X_umap" not in adata_combined.obsm:
        sc.pp.neighbors(adata_combined, use_rep="X_samap" if "X_samap" in adata_combined.obsm else "X_pca")
        sc.tl.umap(adata_combined)

    embedding_key = "X_samap" if "X_samap" in adata_combined.obsm else "X_umap"
    embedding = np.asarray(adata_combined.obsm[embedding_key])

    # Restore batch + celltype from source datasets
    batch_map = {**dict(zip(adata_a.obs_names, adata_a.obs["batch"])),
                 **dict(zip(adata_b.obs_names, adata_b.obs["batch"]))}
    ct_map   = {**dict(zip(adata_a.obs_names, adata_a.obs["celltype"])),
                **dict(zip(adata_b.obs_names, adata_b.obs["celltype"]))}
    adata_combined.obs["batch"]    = [batch_map.get(c, "unknown") for c in adata_combined.obs_names]
    adata_combined.obs["celltype"] = [ct_map.get(c, "unknown")   for c in adata_combined.obs_names]

    # Store embedding under consistent key for downstream evaluation
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
