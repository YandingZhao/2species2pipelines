suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
})

scran_lognorm <- function(obj) {
  suppressPackageStartupMessages({
    library(SingleCellExperiment); library(scran); library(scuttle)
  })
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

obj_a <- readRDS(input_a)
obj_b <- readRDS(input_b)

if (!inherits(obj_a, "Seurat") || !inherits(obj_b, "Seurat")) {
  stop("Both inputs must be Seurat objects in .rds format", call. = FALSE)
}

# Force deterministic per-dataset batch annotations for Harmony.
obj_a$batch <- paste0(sample_id, "_", species_a)
obj_b$batch <- paste0(sample_id, "_", species_b)

merged <- merge(obj_a, y = obj_b, add.cell.ids = c("a", "b"))

if (!("celltype" %in% colnames(merged@meta.data))) {
  merged$celltype <- "unknown"
}

if (normalization == "sctransform") {
  merged <- SCTransform(merged, verbose = FALSE)
  merged <- RunPCA(merged, assay = "SCT", npcs = 30, verbose = FALSE)
  merged <- RunHarmony(merged, group.by.vars = "batch")
} else if (normalization == "scran") {
  merged <- scran_lognorm(merged)
  nfeatures <- min(2000, nrow(merged))
  merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = nfeatures)
  merged <- ScaleData(merged)
  merged <- RunPCA(merged, npcs = 30, verbose = FALSE)
  merged <- RunHarmony(merged, group.by.vars = "batch", theta = 2)
} else {
  merged <- NormalizeData(merged)
  nfeatures <- min(2000, nrow(merged))
  merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = nfeatures)
  merged <- ScaleData(merged)
  merged <- RunPCA(merged, npcs = 30, verbose = FALSE)
  merged <- RunHarmony(merged, group.by.vars = "batch", theta = 2)
}

harmony_res <- as.data.frame(Embeddings(merged, reduction = "harmony"))
harmony_res$cell <- rownames(harmony_res)
harmony_res <- harmony_res[, c("cell", setdiff(colnames(harmony_res), "cell"))]

pca_out <- paste0(sample_id, "_harmony_embedding.tsv")
rds_out <- paste0(sample_id, "_harmony_integration.rds")
report_out <- paste0(sample_id, "_harmony_report.txt")

write.table(harmony_res, file = pca_out, sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(merged, file = rds_out)

writeLines(
  c(
    paste("sample:", sample_id),
    paste("species_a:", species_a),
    paste("species_b:", species_b),
    paste("cells:", ncol(merged)),
    paste("genes:", nrow(merged)),
    "status: ok"
  ),
  con = report_out
)
