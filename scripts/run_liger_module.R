suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratWrappers)
  library(rliger)
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

input_a       <- get_arg("input_a")
input_b       <- get_arg("input_b")
sample_id     <- get_arg("sample_id")
species_a     <- get_arg("species_a")
species_b     <- get_arg("species_b")
features_file <- get_arg("features_file", required = FALSE)

obj_a <- readRDS(input_a)
obj_b <- readRDS(input_b)

if (!inherits(obj_a, "Seurat") || !inherits(obj_b, "Seurat")) {
  stop("Both inputs must be Seurat objects in .rds format", call. = FALSE)
}

obj_a$batch <- paste0(sample_id, "_", species_a)
obj_b$batch <- paste0(sample_id, "_", species_b)

merged <- merge(obj_a, y = obj_b, add.cell.ids = c("a", "b"))

if (!("celltype" %in% colnames(merged@meta.data))) {
  merged$celltype <- "unknown"
}

# LIGER normalizes raw counts internally — always use counts layer.
# Join split Seurat v5 layers before extracting the count matrix.
merged <- JoinLayers(merged, assay = "RNA")

# Feature selection: external list takes priority over internal HVG.
if (!is.null(features_file)) {
  genes <- readLines(features_file)
  genes <- intersect(genes, rownames(merged))
  if (length(genes) < 10)
    stop("features_file has fewer than 10 valid genes after intersection", call. = FALSE)
  VariableFeatures(merged) <- genes
} else {
  merged <- FindVariableFeatures(merged, selection.method = "vst",
                                 nfeatures = min(2000L, nrow(merged)), verbose = FALSE)
}

# LIGER pipeline via SeuratWrappers (wraps rliger iNMF factorization).
# k=20 rank, lambda=5 regularization — standard defaults from Welch et al. 2019.
merged <- RunOptimizeALS(merged, k = 20, lambda = 5, split.by = "batch")
merged <- RunQuantileNorm(merged, split.by = "batch")

liger_emb <- as.data.frame(Embeddings(merged, reduction = "iNMF"))
liger_emb$cell <- rownames(liger_emb)
liger_emb <- liger_emb[, c("cell", setdiff(colnames(liger_emb), "cell"))]

emb_out    <- paste0(sample_id, "_liger_embedding.tsv")
rds_out    <- paste0(sample_id, "_liger_integration.rds")
report_out <- paste0(sample_id, "_liger_report.txt")

write.table(liger_emb, file = emb_out, sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(merged, file = rds_out)

writeLines(
  c(
    paste("sample:", sample_id),
    paste("species_a:", species_a),
    paste("species_b:", species_b),
    paste("features_file:", if (!is.null(features_file)) features_file else "none"),
    paste("rliger_version:", as.character(packageVersion("rliger"))),
    paste("n_genes_used:", length(VariableFeatures(merged))),
    paste("cells:", ncol(merged)),
    paste("genes:", nrow(merged)),
    "status: ok"
  ),
  con = report_out
)
