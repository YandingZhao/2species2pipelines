"""
SIMBA cross-species RNA integration module.

SIMBA embeds cells and genes jointly in a shared space via a cell-gene
bipartite graph trained with PyTorch BigGraph (PBG). After ortholog
conversion both species share the same gene vocabulary, so cross-species
integration is equivalent to batch correction in the shared embedding space.
"""

import argparse
import os
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--input_a",   required=True)
    p.add_argument("--input_b",   required=True)
    p.add_argument("--sample_id", required=True)
    p.add_argument("--species_a", required=True)
    p.add_argument("--species_b", required=True)
    p.add_argument("--normalization", default="log_norm")
    p.add_argument("--features_file", default=None)
    return p.parse_args()


def main():
    args = parse_args()
    import simba as si
    import torch

    t0 = time.time()

    # ── Load ─────────────────────────────────────────────────────────────────
    adata_a = ad.io.read_h5ad(args.input_a) if hasattr(ad.io, "read_h5ad") \
              else ad.read_h5ad(args.input_a)
    adata_b = ad.io.read_h5ad(args.input_b) if hasattr(ad.io, "read_h5ad") \
              else ad.read_h5ad(args.input_b)

    adata_a.obs["batch"]    = f"{args.sample_id}_{args.species_a}"
    adata_b.obs["batch"]    = f"{args.sample_id}_{args.species_b}"
    adata_a.obs["celltype"] = adata_a.obs.get("celltype", pd.Categorical(["unknown"] * adata_a.n_obs))
    adata_b.obs["celltype"] = adata_b.obs.get("celltype", pd.Categorical(["unknown"] * adata_b.n_obs))

    adata = ad.concat([adata_a, adata_b], join="inner", label="batch_src")

    # ── Preprocess ───────────────────────────────────────────────────────────
    si.pp.filter_genes(adata, min_n_cells=3)
    si.pp.cal_qc_rna(adata)

    if args.normalization not in ("raw_counts", "pre_normalized"):
        si.pp.normalize(adata, method="lib_size")
        si.pp.log_transform(adata)

    if args.features_file:
        genes = Path(args.features_file).read_text().splitlines()
        genes = [g for g in genes if g in adata.var_names]
        if len(genes) < 10:
            sys.exit("ERROR: features_file has fewer than 10 valid genes")
        adata.var["highly_variable"] = adata.var_names.isin(genes)
    else:
        si.pp.select_variable_genes(adata, n_top_genes=min(2000, adata.n_vars))

    # Discretize continuous expression into bins for graph edge weights
    si.tl.discretize(adata, n_bins=5)

    # ── Build cell-gene graph and train PBG ──────────────────────────────────
    graph_dir = f"simba_graph_{args.sample_id}"
    si.tl.gen_graph(
        list_CG=[adata],
        copy=False,
        use_highly_variable=True,
        dirname=graph_dir,
    )

    pbg_params = si.settings.pbg_params.copy()
    pbg_params["num_gpus"] = 1 if torch.cuda.is_available() else 0
    if not torch.cuda.is_available():
        pbg_params["workers"] = 2

    si.tl.pbg_train(pbg_params=pbg_params, auto_wd=True, save_wd=True, output="model")

    si.load_graph_stats()
    si.load_pbg_config()
    dict_adata = si.read_embedding()

    cell_emb = dict_adata["C"]  # cell embeddings (n_cells × d)

    # ── Write outputs ─────────────────────────────────────────────────────────
    adata.obsm["X_simba"] = np.array(cell_emb)

    emb_df = pd.DataFrame(
        cell_emb.values if hasattr(cell_emb, "values") else np.array(cell_emb),
        index=adata.obs_names,
        columns=[f"SIMBA_{i+1}" for i in range(np.array(cell_emb).shape[1])],
    )
    emb_df.index.name = "cell"
    emb_df.reset_index().to_csv(f"{args.sample_id}_simba_embedding.tsv", sep="\t", index=False)

    adata.write_h5ad(f"{args.sample_id}_simba_integration.h5ad")

    elapsed = time.time() - t0
    Path(f"{args.sample_id}_simba_report.txt").write_text("\n".join([
        f"sample: {args.sample_id}",
        f"species_a: {args.species_a}",
        f"species_b: {args.species_b}",
        f"normalization: {args.normalization}",
        f"features_file: {args.features_file or 'none'}",
        f"n_genes_used: {int(adata.var['highly_variable'].sum())}",
        f"cells: {adata.n_obs}",
        f"genes: {adata.n_vars}",
        f"elapsed_s: {elapsed:.1f}",
        "status: ok",
    ]) + "\n")


if __name__ == "__main__":
    main()
