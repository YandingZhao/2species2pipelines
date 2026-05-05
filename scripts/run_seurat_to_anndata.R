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

  # Detect available layers; post-integration objects may have counts.1, counts.2, etc.
  all_layers <- tryCatch(Layers(obj[["RNA"]]), error = function(e) character(0))
  count_layers <- grep("^counts", all_layers, value = TRUE)

  if (length(count_layers) > 1) {
    # Multiple split layers (e.g. after merge of two species) — combine by cbind
    matrices <- lapply(count_layers, function(lyr) {
      m <- tryCatch(LayerData(obj, assay = "RNA", layer = lyr), error = function(e) NULL)
      if (!is.null(m) && !inherits(m, "dgCMatrix")) m <- as(m, "dgCMatrix")
      m
    })
    matrices <- Filter(Negate(is.null), matrices)
    counts <- do.call(cbind, matrices)
    # Restore cell ordering to match the Seurat object
    counts <- counts[, colnames(obj), drop = FALSE]
  } else if (length(count_layers) == 1) {
    counts <- tryCatch(
      LayerData(obj, assay = "RNA", layer = count_layers),
      error = function(e) NULL
    )
    if (is.null(counts)) {
      counts <- GetAssayData(obj, assay = "RNA", slot = "counts")
    }
  } else {
    counts <- GetAssayData(obj, assay = "RNA", slot = "counts")
  }

  if (!inherits(counts, "dgCMatrix")) {
    counts <- as(counts, "dgCMatrix")
  }
  counts
}


extract_reductions <- function(obj, cell_names) {
  reduction_names <- Reductions(obj)
  reductions <- list()

  for (reduction_name in reduction_names) {
    embeddings <- tryCatch(
      Embeddings(obj[[reduction_name]]),
      error = function(e) NULL
    )

    if (is.null(embeddings) || nrow(embeddings) == 0) {
      next
    }

    embeddings <- embeddings[cell_names, , drop = FALSE]
    reductions[[reduction_name]] <- embeddings
  }

  reductions
}


write_h5ad_with_python <- function(mtx_file, obs_file, var_file, reduction_files, output_h5ad) {
  py_file <- tempfile(fileext = ".py")
  py_code <- c(
    "import sys",
    "import numpy as np",
    "import anndata as ad",
    "import pandas as pd",
    "from scipy.io import mmread",
    "",
    "mtx_file, obs_file, var_file, reductions_file, out_file = sys.argv[1:6]",
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
    "try:",
    "    reductions = pd.read_csv(reductions_file, sep='\\t')",
    "except pd.errors.EmptyDataError:",
    "    reductions = pd.DataFrame()",
    "for reduction_name, reduction_path in reductions.items():",
    "    reduction = pd.read_csv(reduction_path.iloc[0], sep='\\t', index_col=0)",
    "    reduction = reduction.loc[obs.index]",
    "    adata.obsm[f'X_{reduction_name}'] = np.asarray(reduction)",
    "adata.write_h5ad(out_file)"
  )
  writeLines(py_code, con = py_file)

  py_exec <- Sys.getenv("RETICULATE_PYTHON")
  if (!nzchar(py_exec)) {
    py_exec <- "python3"
  }

  status <- system2(py_exec, c(py_file, mtx_file, obs_file, var_file, reduction_files, output_h5ad))
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
  reductions <- extract_reductions(obj, colnames(counts))

  mtx_file <- tempfile(fileext = ".mtx")
  obs_file <- tempfile(fileext = ".tsv")
  var_file <- tempfile(fileext = ".tsv")
  reduction_manifest <- tempfile(fileext = ".tsv")
  reduction_files <- character(0)

  writeMM(counts, file = mtx_file)
  write.table(obs, file = obs_file, sep = "\t", quote = FALSE, col.names = NA)
  write.table(var, file = var_file, sep = "\t", quote = FALSE, col.names = NA)

  if (length(reductions) > 0) {
    reduction_files <- vapply(names(reductions), function(reduction_name) {
      reduction_file <- tempfile(pattern = paste0(reduction_name, "_"), fileext = ".tsv")
      write.table(reductions[[reduction_name]], file = reduction_file, sep = "\t", quote = FALSE, col.names = NA)
      reduction_file
    }, character(1))
  }

  write.table(
    as.data.frame(as.list(reduction_files), check.names = FALSE),
    file = reduction_manifest,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )

  write_h5ad_with_python(mtx_file, obs_file, var_file, reduction_manifest, output_h5ad)

  unlink(c(mtx_file, obs_file, var_file, reduction_manifest, unname(reduction_files)))
}

input_a <- get_arg("input_a")
output_a <- paste0(gsub("\\.rds$", "", basename(input_a)), ".h5ad")
convert_one(input_a, output_a)

tryCatch({
  input_b <- get_arg("input_b")
  output_b <- paste0(gsub("\\.rds$", "", basename(input_b)), ".h5ad")
  convert_one(input_b, output_b)
}, error = function(e) {
  message("No second input provided, skipping conversion for input_b.")
})
