suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
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

write_h5ad_with_python <- function(mtx_file, obs_file, var_file, output_h5ad) {
  py_file <- tempfile(fileext = ".py")
  py_code <- c(
    "import sys",
    "import anndata as ad",
    "import pandas as pd",
    "from scipy.io import mmread",
    "",
    "mtx_file, obs_file, var_file, out_file = sys.argv[1:5]",
    "X = mmread(mtx_file).tocsr().transpose()",
    "obs = pd.read_csv(obs_file, sep='\\t', index_col=0)",
    "var = pd.read_csv(var_file, sep='\\t', index_col=0)",
    "obs.index = pd.Index(obs.index.astype(str), dtype='object')",
    "var.index = pd.Index(var.index.astype(str), dtype='object')",
    "for df in (obs, var):",
    "    for col in df.columns:",
    "        if pd.api.types.is_string_dtype(df[col].dtype):",
    "            df[col] = df[col].astype('object')",
    "adata = ad.AnnData(X=X, obs=obs, var=var)",
    "adata.write_h5ad(out_file)"
  )
  writeLines(py_code, con = py_file)

  py_exec <- Sys.getenv("RETICULATE_PYTHON")
  if (!nzchar(py_exec)) {
    py_exec <- "python3"
  }

  status <- system2(py_exec, c(py_file, mtx_file, obs_file, var_file, output_h5ad))
  unlink(py_file)

  if (!identical(status, 0L)) {
    stop("Python conversion to .h5ad failed", call. = FALSE)
  }
}

convert_one <- function(input_rds, output_h5ad) {
  obj <- readRDS(input_rds)
  if (!inherits(obj, "Seurat")) {
    stop(paste("Input is not a Seurat object:", input_rds), call. = FALSE)
  }

  counts <- extract_counts(obj)

  obs <- obj@meta.data
  if (!("batch" %in% colnames(obs))) {
    obs$batch <- "unknown"
  }
  if (!("celltype" %in% colnames(obs))) {
    obs$celltype <- "unknown"
  }
  obs <- obs[colnames(counts), , drop = FALSE]

  var <- data.frame(gene = rownames(counts), row.names = rownames(counts), stringsAsFactors = FALSE)

  mtx_file <- tempfile(fileext = ".mtx")
  obs_file <- tempfile(fileext = ".tsv")
  var_file <- tempfile(fileext = ".tsv")

  writeMM(counts, file = mtx_file)
  write.table(obs, file = obs_file, sep = "\t", quote = FALSE, col.names = NA)
  write.table(var, file = var_file, sep = "\t", quote = FALSE, col.names = NA)

  write_h5ad_with_python(mtx_file, obs_file, var_file, output_h5ad)

  unlink(c(mtx_file, obs_file, var_file))
}

input_a <- get_arg("input_a")
input_b <- get_arg("input_b")
sample_id <- get_arg("sample_id")

output_a <- paste0(sample_id, "_a.h5ad")
output_b <- paste0(sample_id, "_b.h5ad")

convert_one(input_a, output_a)
convert_one(input_b, output_b)
