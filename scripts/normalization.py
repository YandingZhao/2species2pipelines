"""Shared normalization utilities for all Python integration scripts."""
import scanpy as sc

NORM_METHODS = ("log_norm", "pearson_residuals", "raw_counts")


def apply_normalization(adata, method: str) -> None:
    """Normalize adata.X in-place; preserves raw counts in layers['counts'].

    Parameters
    ----------
    adata   : AnnData with raw (integer) counts in X.
    method  : One of 'log_norm', 'pearson_residuals', 'raw_counts'.
    """
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()

    if method == "log_norm":
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
    elif method == "pearson_residuals":
        sc.experimental.pp.normalize_pearson_residuals(adata)
    elif method == "raw_counts":
        pass  # keep X as integer counts
    else:
        raise ValueError(
            f"Unknown normalization {method!r}. Choose from {NORM_METHODS}"
        )
