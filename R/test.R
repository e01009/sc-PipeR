#!/usr/bin/env Rscript

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))) y else x
}

abort <- function(message) {
  stop(message, call. = FALSE)
}

required_packages <- c(
  "yaml",
  "jsonlite",
  "Seurat",
  "SeuratObject",
  "fs",
  "glue",
  "readr",
  "dplyr",
  "tibble",
  "sessioninfo"
)

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  abort(paste0(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    ". Install/restore these in renv before running the pipeline."
  ))
}

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else file.path("config", "params.yaml")

if (!file.exists(config_path)) {
  abort(paste0("Config file not found: ", normalizePath(config_path, mustWork = FALSE)))
}

params <- yaml::read_yaml(config_path)

get_param <- function(path, default = NULL, required = FALSE) {
  value <- params
  for (key in strsplit(path, "\\.", fixed = FALSE)[[1]]) {
    if (is.null(value) || is.null(value[[key]])) {
      if (required) abort(paste0("Missing required config value: ", path))
      return(default)
    }
    value <- value[[key]]
  }
  value %||% default
}

project_name <- get_param("input.project_name", get_param("project.name", "scPipeR_sample"))
input_type <- tolower(get_param("input.type", "10x"))
input_path <- get_param("input.path", required = TRUE)
organism <- tolower(get_param("project.organism", "human"))
seed <- as.integer(get_param("project.seed", 1234))

min_cells <- as.integer(get_param("input.min_cells", 3))
min_features <- as.integer(get_param("qc.min_features", 200))
max_features <- as.integer(get_param("qc.max_features", 6000))
max_percent_mt <- as.numeric(get_param("qc.max_percent_mt", 10))

normalization_method <- get_param("normalization.method", "LogNormalize")
scale_factor <- as.numeric(get_param("normalization.scale_factor", 10000))
variable_features_method <- get_param("normalization.variable_features", "vst")
n_variable_features <- as.integer(get_param("normalization.n_variable_features", 2000))

n_pcs <- as.integer(get_param("reduction.n_pcs", 30))
umap_dims <- as.integer(get_param("reduction.umap_dims", n_pcs))
clustering_resolution <- as.numeric(get_param("clustering.resolution", 0.5))

marker_logfc_threshold <- as.numeric(get_param("markers.logfc_threshold", 0.25))
marker_min_pct <- as.numeric(get_param("markers.min_pct", 0.1))
marker_only_pos <- isTRUE(get_param("markers.only_pos", TRUE))

output_base_dir <- get_param("output.base_dir", "Output")
run_name <- get_param(
  "output.run_name",
  paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_", project_name)
)
save_seurat_object <- isTRUE(get_param("output.save_seurat_object", TRUE))

resolve_10x_dir <- function(path) {
  expected_files <- c("matrix.mtx", "barcodes.tsv", "genes.tsv")
  expected_files_gz <- paste0(expected_files, ".gz")

  has_10x_files <- function(candidate) {
    all(file.exists(file.path(candidate, expected_files))) ||
      all(file.exists(file.path(candidate, expected_files_gz))) ||
      (
        file.exists(file.path(candidate, "matrix.mtx")) &&
          file.exists(file.path(candidate, "barcodes.tsv")) &&
          file.exists(file.path(candidate, "features.tsv"))
      ) ||
      (
        file.exists(file.path(candidate, "matrix.mtx.gz")) &&
          file.exists(file.path(candidate, "barcodes.tsv.gz")) &&
          file.exists(file.path(candidate, "features.tsv.gz"))
      )
  }

  if (has_10x_files(path)) {
    return(path)
  }

  nested_matrices <- fs::dir_ls(path, recurse = TRUE, regexp = "matrix\\.mtx(\\.gz)?$")
  candidate_dirs <- unique(dirname(nested_matrices))
  matching_dirs <- candidate_dirs[vapply(candidate_dirs, has_10x_files, logical(1))]

  if (length(matching_dirs) == 1) {
    return(matching_dirs[[1]])
  }

  if (length(matching_dirs) > 1) {
    abort(paste0(
      "Multiple 10x matrix directories found under input.path. Set input.path to one of: ",
      paste(normalizePath(matching_dirs, mustWork = FALSE), collapse = "; ")
    ))
  }

  abort(paste0(
    "No 10x matrix directory found at or below input.path: ",
    normalizePath(path, mustWork = FALSE),
    ". Expected matrix.mtx, barcodes.tsv, and genes.tsv/features.tsv files."
  ))
}

if (!input_type %in% c("10x", "rds")) {
  abort("input.type must be either '10x' or 'rds'.")
}

if (!file.exists(input_path)) {
  abort(paste0("Input path not found: ", normalizePath(input_path, mustWork = FALSE)))
}

if (min_features > max_features) {
  abort("qc.min_features cannot be greater than qc.max_features.")
}

if (n_pcs < 2) {
  abort("reduction.n_pcs must be at least 2.")
}

set.seed(seed)

run_dir <- file.path(output_base_dir, run_name)
dirs <- c(
  run_dir,
  file.path(run_dir, "objects"),
  file.path(run_dir, "tables"),
  file.path(run_dir, "logs"),
  file.path(run_dir, "config")
)
invisible(lapply(dirs, fs::dir_create))

log_file <- file.path(run_dir, "logs", "run.log")
log_message <- function(...) {
  line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", glue::glue(...))
  message(line)
  cat(line, "\n", file = log_file, append = TRUE)
}

writeLines(
  c(
    "sc-PipeR ver0 run",
    paste0("Started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0("R version: ", getRversion()),
    paste0("Config: ", normalizePath(config_path, mustWork = FALSE)),
    paste0("Input: ", normalizePath(input_path, mustWork = FALSE)),
    paste0("Output: ", normalizePath(run_dir, mustWork = FALSE))
  ),
  con = file.path(run_dir, "logs", "summary.txt")
)

file.copy(config_path, file.path(run_dir, "config", basename(config_path)), overwrite = TRUE)

qc_summary <- function(object, stage) {
  tibble::tibble(
    stage = stage,
    cells = ncol(object),
    genes = nrow(object),
    median_features = stats::median(object$nFeature_RNA),
    median_counts = stats::median(object$nCount_RNA),
    median_percent_mt = stats::median(object$percent.mt),
    mean_percent_mt = mean(object$percent.mt)
  )
}

log_message("Loading input as {input_type}.")
if (input_type == "10x") {
  input_path <- resolve_10x_dir(input_path)
  log_message("Resolved 10x data directory: {normalizePath(input_path, mustWork = FALSE)}")
  counts <- Seurat::Read10X(data.dir = input_path)
  if (is.list(counts)) {
    if ("Gene Expression" %in% names(counts)) {
      counts <- counts[["Gene Expression"]]
    } else {
      counts <- counts[[1]]
    }
  }
  seu <- Seurat::CreateSeuratObject(
    counts = counts,
    project = project_name,
    min.cells = min_cells,
    min.features = 0
  )
} else {
  seu <- readRDS(input_path)
  if (!inherits(seu, "Seurat")) {
    abort("The .rds input must contain a Seurat object.")
  }
}

mt_pattern <- switch(
  organism,
  human = "^MT-",
  mouse = "^mt-",
  abort("Unsupported project.organism. Use 'human' or 'mouse' for ver0.")
)

log_message("Calculating QC metrics.")
seu[["percent.mt"]] <- Seurat::PercentageFeatureSet(seu, pattern = mt_pattern)

qc_before <- qc_summary(seu, "before_filtering")

log_message(
  "Filtering cells with nFeature_RNA >= {min_features}, nFeature_RNA <= {max_features}, percent.mt <= {max_percent_mt}."
)
cell_metadata <- seu@meta.data
keep_cells <- rownames(cell_metadata)[
  cell_metadata$nFeature_RNA >= min_features &
    cell_metadata$nFeature_RNA <= max_features &
    cell_metadata$percent.mt <= max_percent_mt
]
seu <- subset(seu, cells = keep_cells)

if (ncol(seu) == 0) {
  abort("No cells remain after filtering. Loosen QC thresholds and rerun.")
}

qc_after <- qc_summary(seu, "after_filtering")
qc_table <- dplyr::bind_rows(qc_before, qc_after)

log_message("Normalizing data.")
seu <- Seurat::NormalizeData(
  object = seu,
  normalization.method = normalization_method,
  scale.factor = scale_factor,
  verbose = FALSE
)

log_message("Finding variable features.")
seu <- Seurat::FindVariableFeatures(
  object = seu,
  selection.method = variable_features_method,
  nfeatures = n_variable_features,
  verbose = FALSE
)

log_message("Scaling data.")
seu <- Seurat::ScaleData(seu, verbose = FALSE)

available_pcs <- min(n_pcs, ncol(seu) - 1, length(Seurat::VariableFeatures(seu)))
if (available_pcs < 2) {
  abort("Not enough cells or variable features to run PCA. Check input data and QC thresholds.")
}
if (available_pcs < n_pcs) {
  log_message("Reducing requested PCs from {n_pcs} to {available_pcs} based on available data.")
}
pc_dims <- seq_len(available_pcs)
umap_dims <- min(umap_dims, available_pcs)

log_message("Running PCA.")
seu <- Seurat::RunPCA(seu, npcs = available_pcs, verbose = FALSE)

log_message("Running neighbor graph, clustering, and UMAP.")
seu <- Seurat::FindNeighbors(seu, dims = pc_dims, verbose = FALSE)
seu <- Seurat::FindClusters(seu, resolution = clustering_resolution, verbose = FALSE)
seu <- Seurat::RunUMAP(seu, dims = seq_len(umap_dims), verbose = FALSE)

log_message("Finding cluster markers.")
markers <- tryCatch(
  Seurat::FindAllMarkers(
    object = seu,
    only.pos = marker_only_pos,
    logfc.threshold = marker_logfc_threshold,
    min.pct = marker_min_pct
  ),
  error = function(err) {
    log_message("Marker detection failed: {conditionMessage(err)}")
    tibble::tibble()
  }
)

metadata <- tibble::rownames_to_column(seu@meta.data, var = "cell_barcode")
cluster_counts <- metadata |>
  dplyr::count(seurat_clusters, name = "n_cells") |>
  dplyr::arrange(seurat_clusters)

umap_embeddings <- Seurat::Embeddings(seu, reduction = "umap") |>
  as.data.frame() |>
  tibble::rownames_to_column(var = "cell_barcode") |>
  dplyr::left_join(
    metadata |> dplyr::select(cell_barcode, seurat_clusters),
    by = "cell_barcode"
  )

run_metadata <- list(
  project_name = project_name,
  input_type = input_type,
  input_path = normalizePath(input_path, mustWork = FALSE),
  organism = organism,
  seed = seed,
  parameters = list(
    min_cells = min_cells,
    min_features = min_features,
    max_features = max_features,
    max_percent_mt = max_percent_mt,
    normalization_method = normalization_method,
    scale_factor = scale_factor,
    variable_features_method = variable_features_method,
    n_variable_features = n_variable_features,
    requested_pcs = n_pcs,
    used_pcs = available_pcs,
    used_umap_dims = umap_dims,
    clustering_resolution = clustering_resolution,
    marker_logfc_threshold = marker_logfc_threshold,
    marker_min_pct = marker_min_pct,
    marker_only_pos = marker_only_pos
  )
)

log_message("Writing outputs.")
readr::write_csv(qc_table, file.path(run_dir, "tables", "qc_summary.csv"))
readr::write_csv(metadata, file.path(run_dir, "tables", "cell_metadata.csv"))
readr::write_csv(cluster_counts, file.path(run_dir, "tables", "cluster_counts.csv"))
readr::write_csv(markers, file.path(run_dir, "tables", "markers.csv"))
readr::write_csv(umap_embeddings, file.path(run_dir, "tables", "umap_embeddings.csv"))
jsonlite::write_json(
  run_metadata,
  file.path(run_dir, "logs", "run_metadata.json"),
  pretty = TRUE,
  auto_unbox = TRUE,
  null = "null"
)

if (save_seurat_object) {
  saveRDS(seu, file.path(run_dir, "objects", "seurat_processed.rds"))
}

writeLines(
  capture.output(sessioninfo::session_info()),
  con = file.path(run_dir, "logs", "sessioninfo.txt")
)

summary_lines <- c(
  "sc-PipeR ver0 run complete",
  paste0("Completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("Output directory: ", normalizePath(run_dir, mustWork = FALSE)),
  paste0("Input directory: ", normalizePath(input_path, mustWork = FALSE)),
  paste0("Cells removed by filtering: ", qc_before$cells - qc_after$cells),
  paste0("Cells before filtering: ", qc_before$cells),
  paste0("Cells after filtering: ", qc_after$cells),
  paste0("Genes/features: ", qc_after$genes),
  paste0("PCs used: ", available_pcs),
  paste0("UMAP dims used: ", umap_dims),
  paste0("Clusters: ", dplyr::n_distinct(metadata$seurat_clusters)),
  paste0("Marker rows: ", nrow(markers))
)
writeLines(summary_lines, con = file.path(run_dir, "logs", "summary.txt"))

log_message("Run complete. Output directory: {normalizePath(run_dir, mustWork = FALSE)}")
