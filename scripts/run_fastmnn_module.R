suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratWrappers)
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

apply_features <- function(obj, features_file, n = 2000) {
  if (!is.null(features_file)) {
    genes <- readLines(features_file)
    genes <- intersect(genes, rownames(obj))
    if (length(genes) < 10)
      stop("features_file has fewer than 10 valid genes after intersection", call. = FALSE)
    VariableFeatures(obj) <- genes
  } else {
    obj <- FindVariableFeatures(obj, selection.method = "vst",
                                nfeatures = min(n, nrow(obj)), verbose = FALSE)
  }
  obj
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

common_genes <- intersect(rownames(obj_a), rownames(obj_b))
if (length(common_genes) < 50) {
  stop("Insufficient shared genes between inputs for fastMNN", call. = FALSE)
}
obj_a <- obj_a[common_genes, ]
obj_b <- obj_b[common_genes, ]

merged <- merge(obj_a, y = obj_b, add.cell.ids = c("a", "b"))

if (normalization == "sctransform") {
  merged <- SCTransform(merged, verbose = FALSE)
  DefaultAssay(merged) <- "SCT"
  merged <- apply_features(merged, features_file)
} else if (normalization == "scran") {
  merged <- scran_lognorm(merged)
  merged <- apply_features(merged, features_file)
  merged <- ScaleData(merged, verbose = FALSE)
} else if (normalization == "pre_normalized") {
  merged <- JoinLayers(merged, assay = "RNA")
  merged <- apply_features(merged, features_file)
  merged <- ScaleData(merged, verbose = FALSE)
} else {
  merged <- NormalizeData(merged)
  merged <- apply_features(merged, features_file)
  merged <- ScaleData(merged, verbose = FALSE)
}

set.seed(42)
merged <- RunFastMNN(SplitObject(merged, split.by = "batch"))

emb <- as.data.frame(Embeddings(merged, reduction = "mnn"))
if (ncol(emb) < 2) {
  stop("fastMNN returned fewer than 2 embedding dimensions", call. = FALSE)
}
colnames(emb) <- paste0("fastmnn_", seq_len(ncol(emb)))
emb$cell <- rownames(emb)
emb <- emb[, c("cell", setdiff(colnames(emb), "cell"))]

if ("batch" %in% colnames(merged@meta.data)) {
  emb$batch <- merged@meta.data[emb$cell, "batch"]
}
if ("celltype" %in% colnames(merged@meta.data)) {
  emb$celltype <- merged@meta.data[emb$cell, "celltype"]
}

embedding_out <- paste0(sample_id, "_fastmnn_embedding.tsv")
rds_out <- paste0(sample_id, "_fastmnn_integration.rds")
report_out <- paste0(sample_id, "_fastmnn_report.txt")

write.table(emb, file = embedding_out, sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(merged, file = rds_out)

writeLines(
  c(
    paste("sample:", sample_id),
    paste("species_a:", species_a),
    paste("species_b:", species_b),
    paste("normalization:", normalization),
    paste("features_file:", if (!is.null(features_file)) features_file else "none"),
    paste("n_genes_used:", length(VariableFeatures(merged))),
    paste("cells:", nrow(emb)),
    paste("dims:", ncol(Embeddings(merged, reduction = "mnn"))),
    "status: ok"
  ),
  con = report_out
)
