"""Shared feature selection utilities for the 2species2pipelines benchmark.

Supported methods
-----------------
seurat_v3          : Seurat v3 HVG (VST/NB-based) — top performer in Zappia et al. 2025
seurat_v3_batch    : Seurat v3 HVG, species-aware (per-batch selection, union)
seurat             : Seurat original HVG (mean/dispersion)
seurat_batch       : Seurat original HVG, species-aware
cell_ranger        : Cell Ranger HVG
cell_ranger_batch  : Cell Ranger HVG, species-aware
pearson            : Analytic Pearson residuals HVG
pearson_batch      : Pearson residuals HVG, species-aware
mean               : Top N genes by mean expression
variance           : Top N genes by variance
wilcoxon           : Wilcoxon rank-sum marker genes (requires >=2 cell type labels;
                     falls back to seurat_v3 when all labels are 'unknown')
random             : Random selection (negative control)
all                : All features (baseline — no selection)
"""

import random as _random
import warnings

import numpy as np
import scanpy as sc
import scipy.sparse as sp

FS_METHODS = (
    "seurat_v3",
    "seurat_v3_batch",
    "seurat",
    "seurat_batch",
    "cell_ranger",
    "cell_ranger_batch",
    "pearson",
    "pearson_batch",
    "mean",
    "variance",
    "wilcoxon",
    "random",
    "all",
)


def load_features_file(path: str) -> list:
    """Read a features.txt (one gene per line) and return a list."""
    with open(path) as fh:
        return [line.strip() for line in fh if line.strip()]


def _ensure_counts_layer(adata) -> None:
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()


def _lognorm_copy(adata):
    """Return a log-normalized copy without modifying the original."""
    tmp = adata.copy()
    sc.pp.normalize_total(tmp)
    sc.pp.log1p(tmp)
    return tmp


def select_features(
    adata,
    method: str,
    n_features: int = 2000,
    batch_key: str = "batch",
    seed: int = 42,
) -> list:
    """Select features from adata and return an ordered list of gene names.

    adata is expected to have raw counts in X (or in layers['counts'] if X is
    already normalized).  The function never modifies adata in-place.
    """
    if method not in FS_METHODS:
        raise ValueError(f"Unknown method {method!r}. Choose from {FS_METHODS}")

    n_genes = adata.n_vars
    _ensure_counts_layer(adata)

    # ------------------------------------------------------------------
    if method == "all":
        return list(adata.var_names)

    # ------------------------------------------------------------------
    if method == "random":
        rng = _random.Random(seed)
        k = min(n_features, n_genes)
        return sorted(rng.sample(list(adata.var_names), k))

    # ------------------------------------------------------------------
    # Seurat v3 / Pearson residuals — operate on raw counts
    if method in ("seurat_v3", "seurat_v3_batch", "pearson", "pearson_batch"):
        tmp = adata.copy()
        tmp.X = tmp.layers["counts"].copy()
        batch = batch_key if method.endswith("_batch") else None
        n = min(n_features, n_genes)

        if method.startswith("pearson"):
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    sc.experimental.pp.highly_variable_genes(
                        tmp,
                        flavor="pearson_residuals",
                        n_top_genes=n,
                        batch_key=batch,
                    )
            except Exception:
                # Older scanpy: fall back to seurat_v3
                sc.pp.highly_variable_genes(
                    tmp, flavor="seurat_v3", n_top_genes=n, batch_key=batch
                )
        else:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                sc.pp.highly_variable_genes(
                    tmp, flavor="seurat_v3", n_top_genes=n, batch_key=batch
                )

        return list(tmp.var_names[tmp.var["highly_variable"]])

    # ------------------------------------------------------------------
    # Seurat original / Cell Ranger — operate on log-normalized data
    if method in ("seurat", "seurat_batch", "cell_ranger", "cell_ranger_batch"):
        flavor = "seurat" if "seurat" in method else "cell_ranger"
        tmp = _lognorm_copy(adata)
        batch = batch_key if method.endswith("_batch") else None
        n = min(n_features, n_genes)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sc.pp.highly_variable_genes(
                tmp, flavor=flavor, n_top_genes=n, batch_key=batch
            )
        return list(tmp.var_names[tmp.var["highly_variable"]])

    # ------------------------------------------------------------------
    if method in ("mean", "variance"):
        tmp = _lognorm_copy(adata)
        X = tmp.X
        mat = X.toarray() if sp.issparse(X) else np.asarray(X)
        scores = mat.mean(axis=0) if method == "mean" else mat.var(axis=0)
        n = min(n_features, n_genes)
        idx = np.argsort(scores)[::-1][:n]
        return list(tmp.var_names[idx])

    # ------------------------------------------------------------------
    # Wilcoxon: supervised — requires >=2 distinct cell type labels.
    # Falls back to seurat_v3 when all labels are 'unknown' (typical in this
    # pipeline when no annotation is available).
    if method == "wilcoxon":
        labels = adata.obs.get("celltype", None)
        n_labels = labels.nunique() if labels is not None else 0
        if n_labels < 2:
            warnings.warn(
                "wilcoxon requires >=2 cell type labels; "
                "falling back to seurat_v3",
                stacklevel=2,
            )
            return select_features(adata, "seurat_v3", n_features, batch_key, seed)

        tmp = _lognorm_copy(adata)
        sc.tl.rank_genes_groups(tmp, groupby="celltype", method="wilcoxon")
        groups = tmp.uns["rank_genes_groups"]["names"].dtype.names
        per_group = max(1, n_features // len(groups))
        genes: list = []
        seen: set = set()
        for g in groups:
            for gene in tmp.uns["rank_genes_groups"]["names"][g][:per_group]:
                if gene not in seen:
                    seen.add(gene)
                    genes.append(gene)
        return genes[:n_features]

    raise ValueError(f"Unhandled method: {method!r}")
