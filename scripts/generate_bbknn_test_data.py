#!/usr/bin/env python3
"""Generate synthetic .h5ad inputs for BBKNN module tests."""

import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp


rng = np.random.default_rng(42)
shared_genes = [f"GENE{i:04d}" for i in range(300)]
unique_a = [f"A_UNIQ_{i:04d}" for i in range(120)]
unique_b = [f"B_UNIQ_{i:04d}" for i in range(120)]

genes_a = shared_genes + unique_a
genes_b = shared_genes + unique_b

n_cells = 90

X_a = rng.negative_binomial(5, 0.5, size=(n_cells, len(genes_a))).astype(np.float32)
X_b = rng.negative_binomial(5, 0.5, size=(n_cells, len(genes_b))).astype(np.float32)

obs_a = pd.DataFrame(
    {
        "batch": ["canine_batch1" if i < n_cells // 2 else "canine_batch2" for i in range(n_cells)],
        "celltype": [f"type{i % 4}" for i in range(n_cells)],
    },
    index=[f"cell_a_{i}" for i in range(n_cells)],
)

obs_b = pd.DataFrame(
    {
        "batch": ["human_batch1" if i < n_cells // 2 else "human_batch2" for i in range(n_cells)],
        "celltype": [f"type{i % 4}" for i in range(n_cells)],
    },
    index=[f"cell_b_{i}" for i in range(n_cells)],
)

adata_a = ad.AnnData(X=sp.csr_matrix(X_a), obs=obs_a, var=pd.DataFrame(index=genes_a))
adata_b = ad.AnnData(X=sp.csr_matrix(X_b), obs=obs_b, var=pd.DataFrame(index=genes_b))

adata_a.write_h5ad("tests/data/data/Dog_demo_canine.h5ad")
adata_b.write_h5ad("tests/data/data/Human_demo_human.h5ad")

print("Generated tests/data/data/Dog_demo_canine.h5ad")
print("Generated tests/data/data/Human_demo_human.h5ad")
