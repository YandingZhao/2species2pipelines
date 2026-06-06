#!/usr/bin/env Rscript
# Apply SCTransform (regularized NB GLM) normalization to a pair of Seurat RDS files.
#
# SCTransform (Hafemeister & Satija 2019) fits a regularized negative-binomial
# regression per gene, regressing out sequencing depth, and returns Pearson
# residuals stored in the SCT assay's scale.data slot.  Unlike the analytic
# Pearson-residuals approximation in sc.experimental.pp.normalize_pearson_residuals
# (which uses a single global NB dispersion), SCTransform estimates gene-wise
# dispersion parameters using cross-gene regularization.
#
# Inputs : ortholog-converted Seurat RDS files (raw counts in the RNA assay)
# Outputs: h5ad files with
#            X             = SCTransform Pearson residuals (HVGs only)
#            layers/counts = raw integer counts (preserved for scVI / scGen)
#
# Only HVGs are written (typically 3000).  Python integration scripts that
# receive these files should pass --normalization sctransform (a no-op).

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, required = TRUE) {
  key <- paste0("--", name)
  idx <- which(args == key)
  if (length(idx) == 0 || idx == length(args)) {
    if (required) stop(paste("Missing required argument:", key), call. = FALSE)
    return(NULL)
  }
  args[idx + 1]
}

# ── Write h5ad via inline Python (no zellkonverter / basilisk needed) ─────────
write_h5ad <- function(resid_mtx_file, raw_mtx_file, obs_file, var_file, output_h5ad) {
  py_file <- tempfile(fileext = ".py")
  writeLines(c(
    "import sys, numpy as np, anndata as ad, pandas as pd, scipy.sparse as sp",
    "from scipy.io import mmread",
    "resid_mtx, raw_mtx, obs_f, var_f, out = sys.argv[1:6]",
    "X_resid = mmread(resid_mtx).tocsr().T",
    "X_raw   = mmread(raw_mtx ).tocsr().T",
    "obs = pd.read_csv(obs_f, sep='\\t', index_col=0)",
    "var = pd.read_csv(var_f, sep='\\t', index_col=0)",
    "obs.index = pd.Index(obs.index.astype(str), dtype='object')",
    "var.index = pd.Index(var.index.astype(str), dtype='object')",
    "for df in (obs, var):",
    "    for col in df.columns:",
    "        if pd.api.types.is_string_dtype(df[col].dtype):",
    "            df[col] = df[col].astype('object')",
    "adata = ad.AnnData(X=X_resid, obs=obs, var=var)",
    "adata.layers['counts'] = X_raw",
    "adata.write_h5ad(out)"
  ), con = py_file)

  py_exec <- Sys.getenv("RETICULATE_PYTHON")
  if (!nzchar(py_exec)) py_exec <- "python3"

  status <- system2(py_exec,
    c(py_file, resid_mtx_file, raw_mtx_file, obs_file, var_file, output_h5ad))
  unlink(py_file)
  if (!identical(status, 0L)) stop("Python h5ad write failed", call. = FALSE)
}

# ── Main normalization function ───────────────────────────────────────────────
normalize_one <- function(input_rds, output_h5ad) {
  message("Reading: ", input_rds)
  obj <- readRDS(input_rds)
  if (!inherits(obj, "Seurat"))
    stop(paste("Not a Seurat object:", input_rds), call. = FALSE)

  message("  Cells: ", ncol(obj), "  Genes: ", nrow(obj))
  message("  SCTransform...")
  obj <- SCTransform(obj, verbose = FALSE)

  DefaultAssay(obj) <- "SCT"
  # scale.data contains Pearson residuals for HVGs (genes × cells, dense)
  resid <- GetAssayData(obj, assay = "SCT", slot = "scale.data")
  hvg_genes <- rownames(resid)
  message("  HVGs: ", length(hvg_genes))

  # Raw counts restricted to the same HVG gene set for layers["counts"]
  raw <- GetAssayData(obj, assay = "RNA", slot = "counts")[hvg_genes, , drop = FALSE]
  if (!inherits(raw,   "dgCMatrix")) raw   <- as(raw,   "dgCMatrix")
  if (!inherits(resid, "dgCMatrix")) resid <- as(resid, "dgCMatrix")

  obs <- obj@meta.data[colnames(raw), , drop = FALSE]
  if (!("batch"    %in% colnames(obs))) obs$batch    <- "unknown"
  if (!("celltype" %in% colnames(obs))) obs$celltype <- "unknown"
  var <- data.frame(gene = hvg_genes, row.names = hvg_genes, stringsAsFactors = FALSE)

  resid_mtx_file <- tempfile(fileext = ".mtx")
  raw_mtx_file   <- tempfile(fileext = ".mtx")
  obs_file       <- tempfile(fileext = ".tsv")
  var_file       <- tempfile(fileext = ".tsv")

  writeMM(resid, file = resid_mtx_file)
  writeMM(raw,   file = raw_mtx_file)
  write.table(obs, file = obs_file, sep = "\t", quote = FALSE, col.names = NA)
  write.table(var, file = var_file, sep = "\t", quote = FALSE, col.names = NA)

  message("  Writing: ", output_h5ad)
  write_h5ad(resid_mtx_file, raw_mtx_file, obs_file, var_file, output_h5ad)
  unlink(c(resid_mtx_file, raw_mtx_file, obs_file, var_file))
}

# ── Entry point ───────────────────────────────────────────────────────────────
normalize_one(get_arg("input_a"), get_arg("output_a"))

tryCatch({
  input_b  <- get_arg("input_b",  required = FALSE)
  output_b <- get_arg("output_b", required = FALSE)
  if (!is.null(input_b) && !is.null(output_b))
    normalize_one(input_b, output_b)
  else
    message("No second input provided, skipping input_b.")
}, error = function(e) message("input_b skipped: ", conditionMessage(e)))
