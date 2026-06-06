#!/usr/bin/env python3
"""Standalone feature selection — outputs a gene list for downstream integration.

Reads two normalized h5ad files (one per species), merges them, selects features
using the requested method, and writes a plain-text gene list.

Usage
-----
python run_feature_selection.py \
    --input_a dog.h5ad --input_b human.h5ad \
    --sample_id SAMP --species_a dog --species_b human \
    --method seurat_v3 --n_features 2000

Outputs
-------
{sample_id}_{method}_{n_features}_features.txt   — one gene per line
{sample_id}_{method}_{n_features}_fs_report.txt  — run metadata
"""

import argparse
import os
import sys
import time
from pathlib import Path

import scanpy as sc

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from feature_selection import FS_METHODS, select_features


def parse_args():
    parser = argparse.ArgumentParser(description="Standalone feature selection")
    parser.add_argument("--input_a", required=True, help="h5ad for species A")
    parser.add_argument("--input_b", required=True, help="h5ad for species B")
    parser.add_argument("--sample_id", required=True)
    parser.add_argument("--species_a", required=True)
    parser.add_argument("--species_b", required=True)
    parser.add_argument(
        "--method", default="seurat_v3", choices=FS_METHODS,
        help="Feature selection method (default: seurat_v3)",
    )
    parser.add_argument(
        "--n_features", type=int, default=2000,
        help="Number of features to select; ignored for 'all' (default: 2000)",
    )
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    return parser.parse_args()


def main():
    args = parse_args()

    adata_a = sc.read_h5ad(args.input_a)
    adata_b = sc.read_h5ad(args.input_b)

    sc.pp.filter_cells(adata_a, min_genes=200)
    sc.pp.filter_genes(adata_a, min_cells=10)
    sc.pp.filter_cells(adata_b, min_genes=200)
    sc.pp.filter_genes(adata_b, min_cells=10)

    common_genes = sorted(set(adata_a.var_names).intersection(adata_b.var_names))
    if len(common_genes) < 50:
        raise ValueError(
            f"Only {len(common_genes)} shared genes between inputs — "
            "cannot run feature selection"
        )

    adata_a = adata_a[:, common_genes].copy()
    adata_b = adata_b[:, common_genes].copy()

    adata_a.obs["batch"] = f"{args.sample_id}_{args.species_a}"
    adata_b.obs["batch"] = f"{args.sample_id}_{args.species_b}"

    import anndata as ad
    adata_all = ad.concat([adata_a, adata_b], join="inner", merge="same")

    start = time.time()
    features = select_features(
        adata_all,
        method=args.method,
        n_features=args.n_features,
        batch_key="batch",
        seed=args.seed,
    )
    elapsed = time.time() - start

    tag = f"{args.method}_{args.n_features}"
    out_features = f"{args.sample_id}_{tag}_features.txt"
    out_report = f"{args.sample_id}_{tag}_fs_report.txt"

    Path(out_features).write_text("\n".join(features) + "\n")

    with open(out_report, "w", encoding="utf-8") as fh:
        fh.write(f"sample: {args.sample_id}\n")
        fh.write(f"species_a: {args.species_a}\n")
        fh.write(f"species_b: {args.species_b}\n")
        fh.write(f"method: {args.method}\n")
        fh.write(f"n_features_requested: {args.n_features}\n")
        fh.write(f"n_features_selected: {len(features)}\n")
        fh.write(f"n_common_genes: {len(common_genes)}\n")
        fh.write(f"seed: {args.seed}\n")
        fh.write(f"elapsed_seconds: {elapsed:.2f}\n")
        fh.write("status: ok\n")

    print(
        f"[feature_selection] {args.method}: {len(features)} features "
        f"from {len(common_genes)} common genes → {out_features}"
    )


if __name__ == "__main__":
    main()
