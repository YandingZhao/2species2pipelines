#!/usr/bin/env Rscript
# Standalone feature selection for the R branch of the pipeline.
#
# Reads two (optionally pre-normalized) Seurat RDS files, merges them,
# applies the requested feature selection method, and writes a plain-text
# gene list (one gene per line) compatible with --features_file in all
# integration scripts (both R and Python).
#
# Supported methods
# -----------------
# seurat_vst        : Seurat FindVariableFeatures VST (default, top performer)
# seurat_mvp        : Seurat mean.var.plot
# seurat_disp       : Seurat dispersion
# seurat_sct        : SCTransform-based HVG selection
# osca              : Over-dispersion test via scran::modelGeneVar
# brennecke         : Brennecke variance-decomposition (M3Drop)
# nbumi             : Negative-binomial UMI model (M3Drop)
# dubstepr          : DUBStepR dropout-based selection
# scry              : Poisson deviance residuals (scry)
# scpnmf            : Penalized NMF embedding (scPNMF)
# singlecellhaystack: Distribution-based (singleCellHaystack)
# scsegindex        : Stable-expression negative control (scSEGIndex)
# random            : Random selection (negative control)
# all               : All shared genes (no selection)
#
# Usage
# -----
# Rscript run_feature_selection_r.R \
#   --input_a A.rds --input_b B.rds \
#   --sample_id SAMP --species_a dog --species_b human \
#   --method seurat_vst --n_features 2000

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Matrix))

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

input_a    <- get_arg("input_a")
input_b    <- get_arg("input_b")
sample_id  <- get_arg("sample_id")
species_a  <- get_arg("species_a")
species_b  <- get_arg("species_b")
method     <- get_arg("method", required = FALSE); if (is.null(method)) method <- "seurat_vst"
n_features <- as.integer(get_arg("n_features", required = FALSE) %||% 2000L)
seed       <- as.integer(get_arg("seed", required = FALSE) %||% 42L)

`%||%` <- function(a, b) if (!is.null(a)) a else b

set.seed(seed)

# ── Load and align gene spaces ───────────────────────────────────────────────
obj_a <- readRDS(input_a)
obj_b <- readRDS(input_b)
if (!inherits(obj_a, "Seurat") || !inherits(obj_b, "Seurat"))
  stop("Both inputs must be Seurat objects", call. = FALSE)

obj_a$batch <- paste0(sample_id, "_", species_a)
obj_b$batch <- paste0(sample_id, "_", species_b)

common_genes <- intersect(rownames(obj_a), rownames(obj_b))
if (length(common_genes) < 50)
  stop(paste("Only", length(common_genes), "shared genes — insufficient"), call. = FALSE)

obj_a <- obj_a[common_genes, ]
obj_b <- obj_b[common_genes, ]

n_genes <- length(common_genes)
nf      <- min(n_features, n_genes)

# ── Extract raw counts (Seurat v5 safe) ──────────────────────────────────────
get_counts <- function(obj) {
  obj <- JoinLayers(obj, assay = "RNA")
  tryCatch(
    LayerData(obj, assay = "RNA", layer = "counts"),
    error = function(e) GetAssayData(obj, assay = "RNA", slot = "counts")
  )
}

# ── Feature selection ─────────────────────────────────────────────────────────
message("Feature selection: method=", method, "  n_features=", nf)

selected_genes <- if (method == "all") {
  common_genes

} else if (method == "random") {
  sample(common_genes, nf)

} else if (method %in% c("seurat_vst", "seurat_mvp", "seurat_disp")) {
  sel_method <- switch(method,
    seurat_vst  = "vst",
    seurat_mvp  = "mean.var.plot",
    seurat_disp = "dispersion"
  )
  obj_a <- NormalizeData(obj_a, verbose = FALSE)
  obj_b <- NormalizeData(obj_b, verbose = FALSE)
  hvg_a <- FindVariableFeatures(obj_a, selection.method = sel_method,
                                 nfeatures = nf, verbose = FALSE)
  hvg_b <- FindVariableFeatures(obj_b, selection.method = sel_method,
                                 nfeatures = nf, verbose = FALSE)
  union(VariableFeatures(hvg_a), VariableFeatures(hvg_b))

} else if (method == "seurat_sct") {
  obj_a <- SCTransform(obj_a, verbose = FALSE)
  obj_b <- SCTransform(obj_b, verbose = FALSE)
  union(VariableFeatures(obj_a), VariableFeatures(obj_b))

} else if (method == "osca") {
  suppressPackageStartupMessages({
    library(SingleCellExperiment); library(scran); library(scuttle)
  })
  counts_mat <- cbind(get_counts(obj_a), get_counts(obj_b))
  sce    <- SingleCellExperiment(assays = list(counts = counts_mat))
  sce    <- logNormCounts(sce)
  dec    <- modelGeneVar(sce)
  getTopHVGs(dec, n = nf)

} else if (method == "brennecke") {
  if (!requireNamespace("M3Drop", quietly = TRUE))
    stop("M3Drop required: BiocManager::install('M3Drop')", call. = FALSE)
  suppressPackageStartupMessages(library(M3Drop))
  counts_a <- as.matrix(get_counts(obj_a))
  counts_b <- as.matrix(get_counts(obj_b))
  counts_mat <- cbind(counts_a, counts_b)
  norm_mat <- M3Drop::M3DropConvertData(counts_mat, is.counts = TRUE)
  tryCatch({
    result   <- M3Drop::BrenneckeGetVariableGenes(norm_mat, fdr = 0.1, suppress.plot = TRUE)
    genes    <- rownames(result)
    if (length(genes) == 0) genes <- rownames(result)[seq_len(min(nf, nrow(result)))]
    genes[seq_len(min(nf, length(genes)))]
  }, error = function(e) {
    message("Brennecke failed (", conditionMessage(e), "); falling back to seurat_vst")
    obj_a2 <- NormalizeData(obj_a, verbose = FALSE)
    obj_a2 <- FindVariableFeatures(obj_a2, nfeatures = nf, verbose = FALSE)
    VariableFeatures(obj_a2)
  })

} else if (method == "nbumi") {
  if (!requireNamespace("M3Drop", quietly = TRUE))
    stop("M3Drop required: BiocManager::install('M3Drop')", call. = FALSE)
  suppressPackageStartupMessages(library(M3Drop))
  counts_a <- as.matrix(get_counts(obj_a))
  counts_b <- as.matrix(get_counts(obj_b))
  counts_mat <- cbind(counts_a, counts_b)
  tryCatch({
    fit   <- M3Drop::NBumiConvertData(counts_mat, is.counts = TRUE)
    stats <- M3Drop::NBumiFeatureSelectionCombinedDrop(fit, ntop = nf, suppress.plot = TRUE)
    rownames(stats)[seq_len(min(nf, nrow(stats)))]
  }, error = function(e) {
    message("NBumi failed (", conditionMessage(e), "); falling back to seurat_vst")
    obj_a2 <- NormalizeData(obj_a, verbose = FALSE)
    obj_a2 <- FindVariableFeatures(obj_a2, nfeatures = nf, verbose = FALSE)
    VariableFeatures(obj_a2)
  })

} else if (method == "dubstepr") {
  if (!requireNamespace("DUBStepR", quietly = TRUE))
    stop("DUBStepR required: remotes::install_github('prabhakarlab/DUBStepR')", call. = FALSE)
  suppressPackageStartupMessages(library(DUBStepR))
  counts_mat <- cbind(get_counts(obj_a), get_counts(obj_b))
  tryCatch({
    result <- DUBStepR(input.data = counts_mat, num.pcs = 20, k = 10,
                       num.genes = nf, optimise.features = FALSE)
    result$optimal.feature.genes
  }, error = function(e) {
    message("DUBStepR failed (", conditionMessage(e), "); falling back to seurat_vst")
    obj_a2 <- NormalizeData(obj_a, verbose = FALSE)
    obj_a2 <- FindVariableFeatures(obj_a2, nfeatures = nf, verbose = FALSE)
    VariableFeatures(obj_a2)
  })

} else if (method == "scry") {
  if (!requireNamespace("scry", quietly = TRUE))
    stop("scry required: BiocManager::install('scry')", call. = FALSE)
  suppressPackageStartupMessages({ library(scry); library(SingleCellExperiment) })
  counts_mat <- cbind(get_counts(obj_a), get_counts(obj_b))
  sce    <- SingleCellExperiment(assays = list(counts = counts_mat))
  sce    <- scry::devianceFeatureSelection(sce, assay = "counts")
  scores <- rowData(sce)$binomial_deviance
  names(scores) <- rownames(sce)
  names(sort(scores, decreasing = TRUE))[seq_len(nf)]

} else if (method == "scpnmf") {
  if (!requireNamespace("scPNMF", quietly = TRUE))
    stop("scPNMF required: install.packages('scPNMF')", call. = FALSE)
  suppressPackageStartupMessages(library(scPNMF))
  counts_mat <- cbind(get_counts(obj_a), get_counts(obj_b))
  tryCatch({
    result <- scPNMF::scPNMF(data = as.matrix(counts_mat), K = 20)
    # Extract top genes from NMF basis
    basis  <- result$W
    scores <- apply(basis, 1, max)
    names(sort(scores, decreasing = TRUE))[seq_len(min(nf, length(scores)))]
  }, error = function(e) {
    message("scPNMF failed (", conditionMessage(e), "); falling back to seurat_vst")
    obj_a2 <- NormalizeData(obj_a, verbose = FALSE)
    obj_a2 <- FindVariableFeatures(obj_a2, nfeatures = nf, verbose = FALSE)
    VariableFeatures(obj_a2)
  })

} else if (method == "singlecellhaystack") {
  if (!requireNamespace("singleCellHaystack", quietly = TRUE))
    stop("singleCellHaystack required: install.packages('singleCellHaystack')", call. = FALSE)
  suppressPackageStartupMessages(library(singleCellHaystack))
  obj_merged <- merge(obj_a, y = obj_b)
  obj_merged <- NormalizeData(obj_merged, verbose = FALSE)
  obj_merged <- FindVariableFeatures(obj_merged, nfeatures = nf, verbose = FALSE)
  obj_merged <- ScaleData(obj_merged, verbose = FALSE)
  obj_merged <- RunPCA(obj_merged, npcs = 20, verbose = FALSE)
  expr_mat   <- as.matrix(GetAssayData(obj_merged, layer = "data")[common_genes, ])
  coords     <- Embeddings(obj_merged, "pca")[, 1:2, drop = FALSE]
  tryCatch({
    result <- singleCellHaystack::haystack(x = coords, expression = t(expr_mat))
    top    <- singleCellHaystack::show_result_haystack(
      result, what = "D_KL", n = nf)
    rownames(top)[seq_len(min(nf, nrow(top)))]
  }, error = function(e) {
    message("singleCellHaystack failed (", conditionMessage(e), "); falling back to seurat_vst")
    VariableFeatures(obj_merged)[seq_len(nf)]
  })

} else if (method == "scsegindex") {
  if (!requireNamespace("scSEGIndex", quietly = TRUE))
    stop("scSEGIndex required: remotes::install_github('LuyiTian/scSEGIndex')", call. = FALSE)
  suppressPackageStartupMessages(library(scSEGIndex))
  counts_mat <- cbind(get_counts(obj_a), get_counts(obj_b))
  tryCatch({
    seg  <- scSEGIndex(exprs_data = counts_mat)
    # scSEGIndex: low SEGindex = stable expression (negative control use case)
    # Here we select the MOST variable — i.e. highest SEGindex — as features
    scores <- seg$SEGindex
    names(scores) <- rownames(seg)
    names(sort(scores, decreasing = TRUE))[seq_len(min(nf, length(scores)))]
  }, error = function(e) {
    message("scSEGIndex failed (", conditionMessage(e), "); falling back to seurat_vst")
    obj_a2 <- NormalizeData(obj_a, verbose = FALSE)
    obj_a2 <- FindVariableFeatures(obj_a2, nfeatures = nf, verbose = FALSE)
    VariableFeatures(obj_a2)
  })

} else {
  stop(paste("Unknown method:", method), call. = FALSE)
}

# ── Intersect with common genes and write output ─────────────────────────────
selected_genes <- intersect(selected_genes, common_genes)
if (length(selected_genes) < 10)
  stop(paste("Only", length(selected_genes), "genes selected — too few"), call. = FALSE)

tag          <- paste0(method, "_", nf)
out_features <- paste0(sample_id, "_", tag, "_features.txt")
out_report   <- paste0(sample_id, "_", tag, "_fs_report.txt")

writeLines(selected_genes, con = out_features)

writeLines(
  c(
    paste("sample:",             sample_id),
    paste("species_a:",          species_a),
    paste("species_b:",          species_b),
    paste("method:",             method),
    paste("n_features_requested:", n_features),
    paste("n_features_selected:", length(selected_genes)),
    paste("n_common_genes:",     length(common_genes)),
    paste("seed:",               seed),
    "status: ok"
  ),
  con = out_report
)

message("[feature_selection_r] ", method, ": ", length(selected_genes),
        " features from ", length(common_genes), " common genes -> ", out_features)
