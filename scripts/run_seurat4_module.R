suppressPackageStartupMessages({
  library(Seurat)
})

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

input_a       <- get_arg("input_a")
input_b       <- get_arg("input_b")
sample_id     <- get_arg("sample_id")
species_a     <- get_arg("species_a")
species_b     <- get_arg("species_b")
normalization <- get_arg("normalization", required = FALSE)
if (is.null(normalization)) normalization <- "norm_data"
features_file <- get_arg("features_file", required = FALSE)

# Helper: return integration feature list from file or SelectIntegrationFeatures.
get_integration_features <- function(obj_list, features_file, n = 2000) {
  if (!is.null(features_file)) {
    genes   <- readLines(features_file)
    common  <- Reduce(intersect, lapply(obj_list, rownames))
    genes   <- intersect(genes, common)
    if (length(genes) < 10)
      stop("features_file has fewer than 10 valid genes after intersection", call. = FALSE)
    genes
  } else {
    nf <- min(n, min(sapply(obj_list, nrow)))
    SelectIntegrationFeatures(object.list = obj_list, nfeatures = nf)
  }
}

obj_a <- readRDS(input_a)
obj_b <- readRDS(input_b)

if (!inherits(obj_a, "Seurat") || !inherits(obj_b, "Seurat")) {
  stop("Both inputs must be Seurat objects in .rds format", call. = FALSE)
}

obj_a$batch <- paste0(sample_id, "_", species_a)
obj_b$batch <- paste0(sample_id, "_", species_b)

if (!("celltype" %in% colnames(obj_a@meta.data))) {
  obj_a$celltype <- "unknown"
}
if (!("celltype" %in% colnames(obj_b@meta.data))) {
  obj_b$celltype <- "unknown"
}

obj_list <- list(obj_a, obj_b)

if (normalization == "sctransform") {
  obj_list <- lapply(obj_list, function(x) SCTransform(x, verbose = FALSE))
  features <- get_integration_features(obj_list, features_file, n = 3000)
  obj_list <- PrepSCTIntegration(object.list = obj_list, anchor.features = features)
  anchors  <- FindIntegrationAnchors(object.list = obj_list,
                                     normalization.method = "SCT",
                                     anchor.features = features)
  combined <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
  DefaultAssay(combined) <- "integrated"
  combined <- ScaleData(combined, verbose = FALSE)
  combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
} else if (normalization == "scran") {
  obj_list <- lapply(obj_list, function(x) scran_lognorm(x))
  features <- get_integration_features(obj_list, features_file)
  anchors  <- FindIntegrationAnchors(object.list = obj_list, anchor.features = features)
  combined <- IntegrateData(anchorset = anchors)
  DefaultAssay(combined) <- "integrated"
  combined <- ScaleData(combined, verbose = FALSE)
  combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
} else if (normalization == "pre_normalized") {
  # Input normalized externally; skip normalization, go straight to integration.
  features <- get_integration_features(obj_list, features_file)
  anchors  <- FindIntegrationAnchors(object.list = obj_list, anchor.features = features)
  combined <- IntegrateData(anchorset = anchors)
  DefaultAssay(combined) <- "integrated"
  combined <- ScaleData(combined, verbose = FALSE)
  combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
} else {
  obj_list <- lapply(obj_list, NormalizeData)
  features <- get_integration_features(obj_list, features_file)
  anchors  <- FindIntegrationAnchors(object.list = obj_list, anchor.features = features)
  combined <- IntegrateData(anchorset = anchors)
  DefaultAssay(combined) <- "integrated"
  combined <- ScaleData(combined, verbose = FALSE)
  combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
}

seurat_res <- as.data.frame(Embeddings(combined, reduction = "pca"))
seurat_res$cell <- rownames(seurat_res)
seurat_res <- seurat_res[, c("cell", setdiff(colnames(seurat_res), "cell"))]

pca_out <- paste0(sample_id, "_seurat4_embedding.tsv")
rds_out <- paste0(sample_id, "_seurat4_integration.rds")
report_out <- paste0(sample_id, "_seurat4_report.txt")

write.table(seurat_res, file = pca_out, sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(combined, file = rds_out)

writeLines(
  c(
    paste("sample:", sample_id),
    paste("species_a:", species_a),
    paste("species_b:", species_b),
    paste("normalization:", normalization),
    paste("features_file:", if (!is.null(features_file)) features_file else "none"),
    paste("n_integration_features:", length(features)),
    paste("cells:", ncol(combined)),
    paste("genes:", nrow(combined)),
    "status: ok"
  ),
  con = report_out
)
