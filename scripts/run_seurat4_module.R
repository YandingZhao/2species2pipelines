suppressPackageStartupMessages({
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name) {
  key <- paste0("--", name)
  idx <- which(args == key)
  if (length(idx) == 0 || idx == length(args)) {
    stop(paste("Missing required argument:", key), call. = FALSE)
  }
  args[idx + 1]
}

input_a <- get_arg("input_a")
input_b <- get_arg("input_b")
sample_id <- get_arg("sample_id")
species_a <- get_arg("species_a")
species_b <- get_arg("species_b")

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

obj_list <- lapply(obj_list, function(x) {
  x <- NormalizeData(x)
  nfeature <- min(2000, nrow(x))
  FindVariableFeatures(x, selection.method = "vst", nfeatures = nfeature)
})

features <- SelectIntegrationFeatures(object.list = obj_list)
anchors <- FindIntegrationAnchors(object.list = obj_list, anchor.features = features)
combined <- IntegrateData(anchorset = anchors)

DefaultAssay(combined) <- "integrated"
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)

seurat_res <- as.data.frame(Embeddings(combined, reduction = "pca"))
seurat_res$cell <- rownames(seurat_res)
seurat_res <- seurat_res[, c("cell", setdiff(colnames(seurat_res), "cell"))]

pca_out <- paste0(sample_id, "_seurat4_pca.tsv")
rds_out <- paste0(sample_id, "_seurat4_integration.rds")
report_out <- paste0(sample_id, "_seurat4_report.txt")

write.table(seurat_res, file = pca_out, sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(combined, file = rds_out)

writeLines(
  c(
    paste("sample:", sample_id),
    paste("species_a:", species_a),
    paste("species_b:", species_b),
    paste("cells:", ncol(combined)),
    paste("genes:", nrow(combined)),
    "status: ok"
  ),
  con = report_out
)
