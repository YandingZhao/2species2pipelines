#!/usr/bin/env Rscript
# Apply scran pooling-based normalization to a pair of Seurat RDS files.
#
# scran (Lun et al. 2016) estimates size factors via pooling-based
# deconvolution, handling composition bias between batches/species better
# than simple library-size normalization.
#
# Inputs : ortholog-converted Seurat RDS files (raw counts in the counts slot)
# Outputs: h5ad files with
#            X             = scran log-normalized counts
#            layers/counts = raw integer counts (preserved for scVI / scGen)
#
# Avoids zellkonverter::readH5AD() (which triggers basilisk → Python compile)
# by reading Seurat RDS directly and writing h5ad via an inline Python call —
# the same pattern used in run_seurat_to_anndata.R.

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(scran)
  library(scuttle)
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

# ── Extract raw counts from a Seurat object ──────────────────────────────────
extract_counts <- function(obj) {
  DefaultAssay(obj) <- "RNA"
  all_layers   <- tryCatch(Layers(obj[["RNA"]]), error = function(e) character(0))
  count_layers <- grep("^counts", all_layers, value = TRUE)

  if (length(count_layers) > 1) {
    matrices <- lapply(count_layers, function(lyr) {
      m <- tryCatch(LayerData(obj, assay = "RNA", layer = lyr), error = function(e) NULL)
      if (!is.null(m) && !inherits(m, "dgCMatrix")) m <- as(m, "dgCMatrix")
      m
    })
    matrices <- Filter(Negate(is.null), matrices)
    counts   <- do.call(cbind, matrices)[, colnames(obj), drop = FALSE]
  } else {
    counts <- tryCatch(
      LayerData(obj, assay = "RNA",
                layer = if (length(count_layers) == 1) count_layers else "counts"),
      error = function(e) GetAssayData(obj, assay = "RNA", slot = "counts")
    )
  }
  if (!inherits(counts, "dgCMatrix")) counts <- as(counts, "dgCMatrix")
  counts
}

# ── Write h5ad via inline Python (no zellkonverter / basilisk needed) ─────────
write_h5ad <- function(norm_mtx_file, raw_mtx_file, obs_file, var_file, output_h5ad) {
  py_file <- tempfile(fileext = ".py")
  writeLines(c(
    "import sys, numpy as np, anndata as ad, pandas as pd, scipy.sparse as sp",
    "from scipy.io import mmread",
    "norm_mtx, raw_mtx, obs_f, var_f, out = sys.argv[1:6]",
    "X_norm = mmread(norm_mtx).tocsr().T",
    "X_raw  = mmread(raw_mtx ).tocsr().T",
    "obs = pd.read_csv(obs_f, sep='\\t', index_col=0)",
    "var = pd.read_csv(var_f, sep='\\t', index_col=0)",
    "obs.index = pd.Index(obs.index.astype(str), dtype='object')",
    "var.index = pd.Index(var.index.astype(str), dtype='object')",
    "for df in (obs, var):",
    "    for col in df.columns:",
    "        if pd.api.types.is_string_dtype(df[col].dtype):",
    "            df[col] = df[col].astype('object')",
    "adata = ad.AnnData(X=X_norm, obs=obs, var=var)",
    "adata.layers['counts'] = X_raw",
    "adata.write_h5ad(out)"
  ), con = py_file)

  py_exec <- Sys.getenv("RETICULATE_PYTHON")
  if (!nzchar(py_exec)) py_exec <- "python3"

  status <- system2(py_exec,
    c(py_file, norm_mtx_file, raw_mtx_file, obs_file, var_file, output_h5ad))
  unlink(py_file)
  if (!identical(status, 0L)) stop("Python h5ad write failed", call. = FALSE)
}

# ── Main normalization function ───────────────────────────────────────────────
normalize_one <- function(input_rds, output_h5ad) {
  message("Reading: ", input_rds)
  obj <- readRDS(input_rds)
  if (!inherits(obj, "Seurat"))
    stop(paste("Not a Seurat object:", input_rds), call. = FALSE)

  raw <- extract_counts(obj)
  message("  Cells: ", ncol(raw), "  Genes: ", nrow(raw))

  # Build a minimal SingleCellExperiment for scran
  sce <- SingleCellExperiment(assays = list(counts = raw))

  n_cells  <- ncol(sce)
  min_size <- max(10L, as.integer(n_cells %/% 20L))
  message("  quickCluster (min.size=", min_size, ")...")
  clusters <- tryCatch(
    quickCluster(sce, assay.type = "counts", min.size = min_size),
    error = function(e) {
      message("  quickCluster failed (", conditionMessage(e),
              "); using single cluster")
      factor(rep("1", n_cells))
    }
  )

  message("  computeSumFactors...")
  sce <- computeSumFactors(sce, clusters = clusters, assay.type = "counts")
  sce <- logNormCounts(sce, assay.type = "counts", log = TRUE)

  norm_counts <- assay(sce, "logcounts")   # scran log-normalized
  raw_counts  <- assay(sce, "counts")      # original integers

  # Build obs / var data frames from the Seurat object metadata
  obs <- obj@meta.data[colnames(raw_counts), , drop = FALSE]
  if (!("batch"    %in% colnames(obs))) obs$batch    <- "unknown"
  if (!("celltype" %in% colnames(obs))) obs$celltype <- "unknown"
  var <- data.frame(gene = rownames(raw_counts),
                    row.names = rownames(raw_counts),
                    stringsAsFactors = FALSE)

  # Write temp files then h5ad
  norm_mtx_file <- tempfile(fileext = ".mtx")
  raw_mtx_file  <- tempfile(fileext = ".mtx")
  obs_file      <- tempfile(fileext = ".tsv")
  var_file      <- tempfile(fileext = ".tsv")

  writeMM(norm_counts, file = norm_mtx_file)
  writeMM(raw_counts,  file = raw_mtx_file)
  write.table(obs, file = obs_file, sep = "\t", quote = FALSE, col.names = NA)
  write.table(var, file = var_file, sep = "\t", quote = FALSE, col.names = NA)

  message("  Writing: ", output_h5ad)
  write_h5ad(norm_mtx_file, raw_mtx_file, obs_file, var_file, output_h5ad)
  unlink(c(norm_mtx_file, raw_mtx_file, obs_file, var_file))
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
