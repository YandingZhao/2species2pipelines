#!/usr/bin/env Rscript
# Standalone normalization for the R branch of the pipeline.
#
# Normalizes each input Seurat RDS independently and writes a normalized RDS
# that R integration scripts can consume with --normalization pre_normalized
# (skipping their internal normalization step).
#
# Supported methods
# -----------------
# log_norm      : library-size normalization + log1p  (NormalizeData)
# scran         : pooling-based size-factor normalization (Lun et al. 2016)
# sctransform   : regularized NB GLM; writes SCT assay (Hafemeister & Satija 2019)
# raw_counts    : no normalization (identity)
#
# Usage
# -----
# Rscript run_normalize_r.R \
#   --input_a A.rds --input_b B.rds \
#   --method log_norm \
#   --output_a A_norm.rds --output_b B_norm.rds

suppressPackageStartupMessages(library(Seurat))

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

# ── scran helper (same logic as in integration modules) ─────────────────────
scran_lognorm <- function(obj) {
  suppressPackageStartupMessages({
    library(SingleCellExperiment); library(scran); library(scuttle)
  })
  obj        <- JoinLayers(obj, assay = "RNA")
  counts_mat <- LayerData(obj, assay = "RNA", layer = "counts")
  sce   <- SingleCellExperiment(assays = list(counts = counts_mat))
  n     <- ncol(sce)
  minsz <- max(10L, as.integer(n %/% 20L))
  clust <- tryCatch(quickCluster(sce, assay.type = "counts", min.size = minsz),
                    error = function(e) factor(rep("1", n)))
  sce   <- computeSumFactors(sce, clusters = clust, assay.type = "counts")
  sce   <- logNormCounts(sce, assay.type = "counts", log = TRUE)
  LayerData(obj, assay = "RNA", layer = "data") <- logcounts(sce)
  obj
}

# ── Per-object normalization ──────────────────────────────────────────────────
normalize_one <- function(input_rds, output_rds, method) {
  message("Normalizing: ", input_rds, "  method=", method)
  obj <- readRDS(input_rds)
  if (!inherits(obj, "Seurat"))
    stop(paste("Not a Seurat object:", input_rds), call. = FALSE)

  if (method == "log_norm") {
    obj <- NormalizeData(obj, verbose = FALSE)
  } else if (method == "scran") {
    obj <- scran_lognorm(obj)
  } else if (method == "sctransform") {
    obj <- SCTransform(obj, verbose = FALSE)
  } else if (method == "raw_counts") {
    # no-op: identity
  } else {
    stop(paste("Unknown method:", method), call. = FALSE)
  }

  message("  Writing: ", output_rds)
  saveRDS(obj, file = output_rds)
}

# ── Entry point ───────────────────────────────────────────────────────────────
method   <- get_arg("method")
input_a  <- get_arg("input_a")
output_a <- get_arg("output_a")
normalize_one(input_a, output_a, method)

tryCatch({
  input_b  <- get_arg("input_b",  required = FALSE)
  output_b <- get_arg("output_b", required = FALSE)
  if (!is.null(input_b) && !is.null(output_b))
    normalize_one(input_b, output_b, method)
  else
    message("No second input provided, skipping input_b.")
}, error = function(e) message("input_b skipped: ", conditionMessage(e)))
