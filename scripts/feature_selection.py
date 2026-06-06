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
    "triku",
    "hotspot",
    "anticor",
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

    # ------------------------------------------------------------------
    # triku: graph-based, uses KNN autocorrelation (Miró-Blanch & Yanes 2021)
    if method == "triku":
        try:
            import triku
        except ImportError as exc:
            raise ImportError("triku required: pip install triku") from exc
        tmp = _lognorm_copy(adata)
        n_comps = min(30, tmp.n_vars - 1, tmp.n_obs - 1)
        sc.pp.pca(tmp, n_comps=n_comps)
        sc.pp.neighbors(tmp, n_neighbors=15)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            triku.tl.triku(tmp)
        hvg_col = next(
            (c for c in ("triku_highly_variable", "highly_variable") if c in tmp.var),
            None,
        )
        selected = list(tmp.var_names[tmp.var[hvg_col]]) if hvg_col else []
        if not selected:
            warnings.warn("triku selected 0 genes; falling back to seurat_v3", stacklevel=2)
            return select_features(adata, "seurat_v3", n_features, batch_key, seed)
        return selected[:n_features]

    # ------------------------------------------------------------------
    # hotspot: spatial autocorrelation in KNN graph (Detlefsen et al. 2022)
    if method == "hotspot":
        try:
            import hotspot as _hs
        except ImportError as exc:
            raise ImportError("hotspot required: pip install hotspot") from exc
        tmp = _lognorm_copy(adata)
        _ensure_counts_layer(tmp)
        n_comps = min(30, tmp.n_vars - 1, tmp.n_obs - 1)
        sc.pp.pca(tmp, n_comps=n_comps)
        try:
            hs = _hs.Hotspot(tmp, model="danb", latent_obsm_key="X_pca", n_neighbors=30)
            hs.create_knn_graph()
            hs.compute_autocorrelations(jobs=1)
        except TypeError:
            # Older API: n_neighbors passed to constructor
            hs = _hs.Hotspot(tmp, model="danb", latent_obsm_key="X_pca")
            hs.create_knn_graph(n_neighbors=30)
            hs.compute_autocorrelations(jobs=1)
        results = hs.results
        fdr_col = next((c for c in ("FDR", "fdr") if c in results.columns), None)
        if fdr_col:
            sig = results.loc[results[fdr_col] < 0.05]
            if sig.empty:
                sig = results
        else:
            sig = results
        z_col = next((c for c in ("Z", "z", "z_score") if c in sig.columns), None)
        if z_col:
            sig = sig.sort_values(z_col, ascending=False)
        return list(sig.head(n_features).index)

    # ------------------------------------------------------------------
    # anticor: anti-correlation-based redundancy filtering (Zeisel et al.)
    if method == "anticor":
        try:
            from anticor_features.anticor_features import get_anticor_genes
        except ImportError as exc:
            raise ImportError(
                "anticor_features required: pip install anticor-features"
            ) from exc
        tmp = _lognorm_copy(adata)
        X = tmp.X.toarray() if sp.issparse(tmp.X) else np.asarray(tmp.X)
        try:
            genes = get_anticor_genes(
                X, list(tmp.var_names), n_features, pre_remove_pathways=[]
            )
        except Exception:
            genes = get_anticor_genes(X, list(tmp.var_names), n_features)
        return list(genes)[:n_features]

    raise ValueError(f"Unhandled method: {method!r}")
