suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(orthogene)
})

input_a <- "Dog_pt15_Immune_Lymphoid_diet.rds"
input_b <- "Human_X00004_Immune_Lymphoid_diet.rds"
sample_id <- "task4_demo"
species_a <- "dog"
species_b <- "human"

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name) {
  key <- paste0("--", name)
  idx <- which(args == key)
  if (length(idx) == 0 || idx == length(args)) {
    stop(paste("Missing required argument:", key), call. = FALSE)
  }
  args[idx + 1]
}

extract_counts <- function(obj) {
  DefaultAssay(obj) <- "RNA"
  counts <- tryCatch(
    LayerData(obj, assay = "RNA", layer = "counts"),
    error = function(e) NULL
  )

  if (is.null(counts)) {
    counts <- GetAssayData(obj, assay = "RNA", slot = "counts")
  }

  if (!inherits(counts, "dgCMatrix")) {
    counts <- as(counts, "dgCMatrix")
  }
  counts
}

normalize_species_label <- function(name) {
  lower <- tolower(name)
  if (lower == "canine") {
    return("dog")
  }
  lower
}

input_a <- get_arg("input_a")
input_b <- get_arg("input_b")
sample_id <- get_arg("sample_id")
species_a <- normalize_species_label(get_arg("species_a"))
species_b <- normalize_species_label(get_arg("species_b"))

obj_a <- readRDS(input_a)
obj_b <- readRDS(input_b)

if (!inherits(obj_a, "Seurat") || !inherits(obj_b, "Seurat")) {
  stop("Both inputs must be Seurat objects in .rds format", call. = FALSE)
}

counts_a <- extract_counts(obj_a)
counts_b <- extract_counts(obj_b)

# Convert species_a genes into species_b ortholog symbols.
counts_a_converted <- orthogene::convert_orthologs(
  gene_df = counts_a,
  input_species = species_a,
  output_species = species_b,
  non121_strategy = "drop_both_species",
  method = "gprofiler"
)

if (!inherits(counts_a_converted, "dgCMatrix")) {
  counts_a_converted <- as(counts_a_converted, "dgCMatrix")
}

rownames(counts_a_converted) <- make.unique(rownames(counts_a_converted))
rownames(counts_b) <- make.unique(rownames(counts_b))

common_genes <- intersect(rownames(counts_a_converted), rownames(counts_b))
if (length(common_genes) < 50) {
  stop("Insufficient shared genes after ortholog conversion", call. = FALSE)
}

counts_a_final <- counts_a_converted[common_genes, , drop = FALSE]
counts_b_final <- counts_b[common_genes, , drop = FALSE]

obj_a_out <- CreateSeuratObject(
  counts = counts_a_final,
  meta.data = obj_a@meta.data[colnames(counts_a_final), , drop = FALSE],
  assay = "RNA",
  min.cells = 0,
  min.features = 0
)
obj_b_out <- CreateSeuratObject(
  counts = counts_b_final,
  meta.data = obj_b@meta.data[colnames(counts_b_final), , drop = FALSE],
  assay = "RNA",
  min.cells = 0,
  min.features = 0
)

output_a <- paste0(sample_id, "_a_ortholog.rds")
output_b <- paste0(sample_id, "_b_ortholog.rds")
report_out <- paste0(sample_id, "_ortholog_report.txt")

saveRDS(obj_a_out, file = output_a)
saveRDS(obj_b_out, file = output_b)

writeLines(
  c(
    paste("sample:", sample_id),
    paste("species_a_input:", species_a),
    paste("species_b_target:", species_b),
    paste("genes_a_before:", nrow(counts_a)),
    paste("genes_a_after_conversion:", nrow(counts_a_converted)),
    paste("genes_b_before:", nrow(counts_b)),
    paste("shared_genes_after_conversion:", length(common_genes)),
    "status: ok"
  ),
  con = report_out
)