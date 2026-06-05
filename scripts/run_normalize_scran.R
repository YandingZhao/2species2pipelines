#!/usr/bin/env Rscript
# Apply scran pooling-based normalization to a pair of raw h5ad files.
#
# scran (Lun et al. 2016) estimates size factors via deconvolution from
# pooled cells, handling composition bias between batches/species better
# than simple library-size normalization.
#
# Output h5ad files have:
#   X             = scran log-normalized counts  (log1p of scran-scaled counts)
#   layers/counts = raw integer counts           (preserved for scVI / scGen)

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scran)
  library(scuttle)
  library(zellkonverter)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, required = TRUE) {
  key  <- paste0("--", name)
  idx  <- which(args == key)
  if (length(idx) == 0 || idx == length(args)) {
    if (required) stop(paste("Missing required argument:", key), call. = FALSE)
    return(NULL)
  }
  args[idx + 1]
}

normalize_one <- function(input_h5ad, output_h5ad) {
  message("  Reading: ", input_h5ad)
  sce <- readH5AD(input_h5ad, X_name = "X")

  # Coerce X to integer counts (run_seurat_to_anndata.R exports raw counts)
  raw <- assay(sce, "X")
  if (!is.integer(raw)) {
    raw <- round(raw)
    storage.mode(raw) <- "integer"
  }
  # Store raw under the canonical "counts" assay name
  assay(sce, "counts") <- raw

  # Quick clustering for robust size-factor estimation; fall back gracefully
  # if the dataset is too small for the default min.size.
  n_cells   <- ncol(sce)
  min_size  <- max(10L, as.integer(n_cells %/% 20L))
  message("  quickCluster (", n_cells, " cells, min.size=", min_size, ")...")
  clusters <- tryCatch(
    quickCluster(sce, assay.type = "counts", min.size = min_size),
    error = function(e) {
      message("  quickCluster failed (", conditionMessage(e),
              "); falling back to single cluster.")
      factor(rep("1", n_cells))
    }
  )

  message("  computeSumFactors...")
  sce <- computeSumFactors(sce, clusters = clusters, assay.type = "counts")

  # logNormCounts uses the computed sizeFactors() and stores result in "logcounts"
  sce <- logNormCounts(sce, assay.type = "counts", log = TRUE)

  # Rename "logcounts" → "X" so zellkonverter writes it as AnnData.X
  assay(sce, "X") <- assay(sce, "logcounts")
  assay(sce, "logcounts") <- NULL

  message("  Writing: ", output_h5ad)
  writeH5AD(sce, output_h5ad)
}

input_a  <- get_arg("input_a")
output_a <- get_arg("output_a")
normalize_one(input_a, output_a)

input_b  <- get_arg("input_b",  required = FALSE)
output_b <- get_arg("output_b", required = FALSE)
if (!is.null(input_b) && !is.null(output_b)) {
  normalize_one(input_b, output_b)
} else {
  message("No second input provided, skipping input_b.")
}
