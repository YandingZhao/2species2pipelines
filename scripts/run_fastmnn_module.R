suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratWrappers)
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

common_genes <- intersect(rownames(obj_a), rownames(obj_b))
if (length(common_genes) < 50) {
  stop("Insufficient shared genes between inputs for fastMNN", call. = FALSE)
}
obj_a <- obj_a[common_genes, ]
obj_b <- obj_b[common_genes, ]

merged <- merge(obj_a, y = obj_b, add.cell.ids = c("a", "b"))
merged <- NormalizeData(merged)
nfeature <- min(2000, nrow(merged))
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = nfeature)
merged <- ScaleData(merged, verbose = FALSE)

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
    paste("cells:", nrow(emb)),
    paste("dims:", ncol(Embeddings(merged, reduction = "mnn"))),
    "status: ok"
  ),
  con = report_out
)
