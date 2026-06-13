#!/usr/bin/env Rscript
# run_ortholog_convert.R
#
# Unified ortholog conversion script supporting multiple gene-matching strategies.
#
# Strategies
# ----------
# ensembl_1to1   : Ensembl/gProfiler strict 1-to-1 orthologs (drop multi-mappers)
# ensembl_1tomany: Ensembl/gProfiler all orthologs; aggregate multi-mappers by summing
# symbol         : Intersect gene symbols directly (no database required)
# blast_bbh_60   : BLAST bidirectional best hits >=60% AA identity (pre-computed map)
# blast_bbh_80   : BLAST bidirectional best hits >=80% AA identity
# blast_bbh_90   : BLAST bidirectional best hits >=90% AA identity
# orthofinder    : OrthoFinder one-to-one orthologs (pre-computed map)
#
# Usage
# -----
# Rscript run_ortholog_convert.R \
#   --strategy   ensembl_1to1 \
#   --input_a    Dog.rds \
#   --input_b    Human.rds \
#   --species_a  dog \
#   --species_b  human \
#   --output_a   dog_converted.rds \
#   --output_b   human_converted.rds \
#   [--map_file  ortholog_map.tsv]   # required for blast_* and orthofinder
#
# Map file format (TSV, no header):
#   col1 = species_a gene symbol, col2 = species_b gene symbol

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
})

# ── Argument parsing ─────────────────────────────────────────────────────────

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, required = TRUE, default = NULL) {
  key <- paste0("--", name)
  idx <- which(args == key)
  if (length(idx) == 0) {
    if (required) stop(paste("Missing required argument:", key), call. = FALSE)
    return(default)
  }
  if (idx == length(args)) stop(paste("No value for argument:", key), call. = FALSE)
  args[idx + 1]
}

STRATEGY    <- get_arg("strategy")
INPUT_A     <- get_arg("input_a")
INPUT_B     <- get_arg("input_b")
SPECIES_A   <- tolower(get_arg("species_a"))
SPECIES_B   <- tolower(get_arg("species_b"))
OUTPUT_A    <- get_arg("output_a")
OUTPUT_B    <- get_arg("output_b")
MAP_FILE    <- get_arg("map_file", required = FALSE)

VALID_STRATEGIES <- c("ensembl_1to1", "ensembl_1tomany", "symbol",
                       "blast_bbh_60", "blast_bbh_80", "blast_bbh_90",
                       "orthofinder")
if (!STRATEGY %in% VALID_STRATEGIES) {
  stop(paste("Unknown strategy:", STRATEGY,
             "\nChoose from:", paste(VALID_STRATEGIES, collapse = ", ")),
       call. = FALSE)
}

# BLAST and OrthoFinder strategies require a pre-computed map file
if (grepl("^blast_|^orthofinder$", STRATEGY) && is.null(MAP_FILE)) {
  stop(paste("--map_file required for strategy:", STRATEGY), call. = FALSE)
}

# ── Helpers ──────────────────────────────────────────────────────────────────

extract_counts <- function(obj) {
  DefaultAssay(obj) <- "RNA"
  counts <- tryCatch(
    LayerData(obj, assay = "RNA", layer = "counts"),
    error = function(e) GetAssayData(obj, assay = "RNA", slot = "counts")
  )
  if (!inherits(counts, "dgCMatrix")) counts <- as(counts, "dgCMatrix")
  counts
}

make_seurat <- function(counts, meta) {
  CreateSeuratObject(
    counts   = counts,
    meta.data = meta[colnames(counts), , drop = FALSE],
    assay    = "RNA",
    min.cells    = 0,
    min.features = 0
  )
}

aggregate_duplicate_rows <- function(mat) {
  # Sum counts of genes that map to the same output symbol (1-to-many case)
  dup_genes <- unique(rownames(mat)[duplicated(rownames(mat))])
  if (length(dup_genes) == 0) return(mat)

  unique_part <- mat[!rownames(mat) %in% dup_genes, , drop = FALSE]
  agg_list <- lapply(dup_genes, function(g) {
    sub <- mat[rownames(mat) == g, , drop = FALSE]
    Matrix::colSums(sub)
  })
  agg_mat <- do.call(rbind, agg_list)
  rownames(agg_mat) <- dup_genes
  rbind(unique_part, agg_mat)
}

apply_external_map <- function(counts_a, counts_b, map_file) {
  # Read pre-computed gene mapping: col1 = species_a gene, col2 = species_b gene
  map <- read.table(map_file, sep = "\t", header = FALSE,
                    stringsAsFactors = FALSE, col.names = c("gene_a", "gene_b"))
  map <- map[map$gene_a %in% rownames(counts_a) &
             map$gene_b %in% rownames(counts_b), ]
  if (nrow(map) < 50) {
    stop(paste("Too few mappable genes in map file:", nrow(map)), call. = FALSE)
  }

  # Build converted species_a matrix (rename rows to species_b symbols)
  counts_a_sub <- counts_a[map$gene_a, , drop = FALSE]
  rownames(counts_a_sub) <- map$gene_b
  counts_a_sub <- aggregate_duplicate_rows(counts_a_sub)

  common <- intersect(rownames(counts_a_sub), rownames(counts_b))
  list(a = counts_a_sub[common, , drop = FALSE],
       b = counts_b[common, , drop = FALSE])
}

# ── Load data ────────────────────────────────────────────────────────────────

cat("Loading RDS files...\n")
obj_a <- readRDS(INPUT_A)
obj_b <- readRDS(INPUT_B)

if (!inherits(obj_a, "Seurat") || !inherits(obj_b, "Seurat")) {
  stop("Both inputs must be Seurat objects (.rds)", call. = FALSE)
}

counts_a <- extract_counts(obj_a)
counts_b <- extract_counts(obj_b)

cat(sprintf("Species A (%s): %d genes x %d cells\n",
            SPECIES_A, nrow(counts_a), ncol(counts_a)))
cat(sprintf("Species B (%s): %d genes x %d cells\n",
            SPECIES_B, nrow(counts_b), ncol(counts_b)))

# ── Strategy dispatch ────────────────────────────────────────────────────────

if (STRATEGY == "symbol") {
  # ── Symbol overlap ──────────────────────────────────────────────────────
  cat("Strategy: symbol — intersecting gene symbols directly\n")
  rownames(counts_a) <- make.unique(rownames(counts_a))
  rownames(counts_b) <- make.unique(rownames(counts_b))
  common <- intersect(rownames(counts_a), rownames(counts_b))
  if (length(common) < 50) stop("Too few shared gene symbols", call. = FALSE)
  counts_a_final <- counts_a[common, , drop = FALSE]
  counts_b_final <- counts_b[common, , drop = FALSE]

} else if (STRATEGY %in% c("blast_bbh_60", "blast_bbh_80",
                             "blast_bbh_90", "orthofinder")) {
  # ── External map (BLAST / OrthoFinder) ─────────────────────────────────
  cat(sprintf("Strategy: %s — applying pre-computed gene map from %s\n",
              STRATEGY, MAP_FILE))
  result <- apply_external_map(counts_a, counts_b, MAP_FILE)
  counts_a_final <- result$a
  counts_b_final <- result$b

} else {
  # ── Ensembl via orthogene (ensembl_1to1 or ensembl_1tomany) ────────────
  if (!requireNamespace("orthogene", quietly = TRUE)) {
    stop("Package 'orthogene' required. Install with: BiocManager::install('orthogene')",
         call. = FALSE)
  }

  non121 <- if (STRATEGY == "ensembl_1to1") "drop_both_species" else "keep_both_species"
  cat(sprintf("Strategy: %s — orthogene::convert_orthologs (non121='%s')\n",
              STRATEGY, non121))

  counts_a_converted <- orthogene::convert_orthologs(
    gene_df          = counts_a,
    input_species    = SPECIES_A,
    output_species   = SPECIES_B,
    non121_strategy  = non121,
    method           = "gprofiler"
  )
  if (!inherits(counts_a_converted, "dgCMatrix")) {
    counts_a_converted <- as(counts_a_converted, "dgCMatrix")
  }

  # Aggregate duplicate output genes (only relevant for ensembl_1tomany)
  counts_a_converted <- aggregate_duplicate_rows(counts_a_converted)
  rownames(counts_b) <- make.unique(rownames(counts_b))

  common <- intersect(rownames(counts_a_converted), rownames(counts_b))
  if (length(common) < 50) stop("Too few shared genes after ortholog conversion", call. = FALSE)

  counts_a_final <- counts_a_converted[common, , drop = FALSE]
  counts_b_final <- counts_b[common, , drop = FALSE]
}

# ── Build output Seurat objects ──────────────────────────────────────────────

cat(sprintf("Shared genes after conversion: %d\n", nrow(counts_a_final)))

obj_a_out <- make_seurat(counts_a_final,
                          obj_a@meta.data[colnames(counts_a_final), , drop = FALSE])
obj_b_out <- make_seurat(counts_b_final,
                          obj_b@meta.data[colnames(counts_b_final), , drop = FALSE])

saveRDS(obj_a_out, file = OUTPUT_A)
saveRDS(obj_b_out, file = OUTPUT_B)

# ── Report ───────────────────────────────────────────────────────────────────

report_file <- sub("\\.rds$", "_ortholog_report.txt", OUTPUT_A)
writeLines(c(
  paste("strategy:", STRATEGY),
  paste("species_a:", SPECIES_A),
  paste("species_b:", SPECIES_B),
  paste("genes_a_input:", nrow(counts_a)),
  paste("genes_b_input:", nrow(counts_b)),
  paste("shared_genes:", nrow(counts_a_final)),
  paste("cells_a:", ncol(counts_a_final)),
  paste("cells_b:", ncol(counts_b_final)),
  if (!is.null(MAP_FILE)) paste("map_file:", MAP_FILE) else NULL,
  "status: ok"
), con = report_file)

cat(sprintf("Done. Outputs: %s, %s\n", OUTPUT_A, OUTPUT_B))
