"""Shared normalization utilities for all Python integration scripts.

Supported methods
-----------------
log_norm          : library-size normalization (CPM) + log1p.  Standard default.
pearson_residuals : analytic Pearson residuals (Lause et al. 2021).
                    Stabilizes variance without log-transform; unbiased HVG selection.
scran             : pooling-based size-factor normalization (Lun et al. 2016).
                    Must be pre-applied by run_normalize_scran.R; this function is a
                    no-op because X already contains scran log-normalized counts.
raw_counts        : no normalization — X is kept as raw integer counts.
                    Only valid for model-based methods with an explicit count
                    likelihood (scVI, scGen).
"""
import scanpy as sc

NORM_METHODS = ("log_norm", "pearson_residuals", "scran", "sctransform", "raw_counts")


def apply_normalization(adata, method: str) -> None:
    """Normalize adata.X in-place; always preserves raw counts in layers['counts'].

    For 'scran': expects X to already contain scran log-normalized counts written
    by run_normalize_scran.R; raw counts are expected in layers['counts'].
    """
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()

    if method == "log_norm":
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
    elif method == "pearson_residuals":
        sc.experimental.pp.normalize_pearson_residuals(adata)
    elif method == "scran":
        pass  # pre-normalized by run_normalize_scran.R; X already contains log-scran
    elif method == "sctransform":
        pass  # pre-normalized by run_normalize_sctransform.R; X already contains SCT residuals
    elif method == "raw_counts":
        pass  # keep X as integer counts
    else:
        raise ValueError(
            f"Unknown normalization {method!r}. Choose from {NORM_METHODS}"
        )
