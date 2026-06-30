#!/usr/bin/env Rscript

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))) y else x
}

abort <- function(message) {
  stop(message, call. = FALSE)
}

if (!requireNamespace("yaml", quietly = TRUE)) {
  abort("Missing required package: yaml. Install/restore it in renv before running the pipeline.")
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

validate_config <- function(config) {
  errors <- character()

  add_error <- function(message) {
    errors <<- c(errors, message)
  }

  get_config_value <- function(path) {
    value <- config
    for (key in strsplit(path, "\\.", fixed = FALSE)[[1]]) {
      if (!is.list(value) || is.null(value[[key]])) {
        return(NULL)
      }
      value <- value[[key]]
    }
    value
  }

  is_scalar_value <- function(value) {
    !is.null(value) && length(value) == 1 && !is.na(value)
  }

  validate_section <- function(path) {
    value <- get_config_value(path)
    if (is.null(value) || !is.list(value)) {
      add_error(paste0(path, " must be a mapping."))
    }
  }

  validate_string <- function(path, required = TRUE, allowed = NULL, allow_null = FALSE) {
    value <- get_config_value(path)

    if (is.null(value)) {
      if (required && !allow_null) add_error(paste0(path, " is required."))
      return(invisible(NULL))
    }

    if (!is_scalar_value(value) || !is.character(value) || !nzchar(trimws(value))) {
      if (allow_null && is.null(value)) return(invisible(NULL))
      add_error(paste0(path, " must be a non-empty string."))
      return(invisible(NULL))
    }

    if (!is.null(allowed) && !tolower(trimws(value)) %in% allowed) {
      add_error(paste0(path, " must be one of: ", paste(allowed, collapse = ", "), "."))
    }
  }

  validate_number <- function(path, required = TRUE, integer = FALSE, min = NULL, max = NULL, positive = FALSE) {
    value <- get_config_value(path)

    if (is.null(value)) {
      if (required) add_error(paste0(path, " is required."))
      return(invisible(NULL))
    }

    if (!is_scalar_value(value) || !is.numeric(value)) {
      add_error(paste0(path, " must be a numeric value."))
      return(invisible(NULL))
    }

    if (integer && value != as.integer(value)) {
      add_error(paste0(path, " must be an integer."))
    }
    if (positive && value <= 0) {
      add_error(paste0(path, " must be > 0."))
    }
    if (!is.null(min) && value < min) {
      add_error(paste0(path, " must be >= ", min, "."))
    }
    if (!is.null(max) && value > max) {
      add_error(paste0(path, " must be <= ", max, "."))
    }
  }

  validate_bool <- function(path, required = TRUE) {
    value <- get_config_value(path)

    if (is.null(value)) {
      if (required) add_error(paste0(path, " is required."))
      return(invisible(NULL))
    }

    if (!is_scalar_value(value) || !is.logical(value)) {
      add_error(paste0(path, " must be true or false."))
    }
  }

  if (!is.list(config)) {
    abort("Config validation failed:\n- Config file must contain a YAML mapping.")
  }

  for (section in c("project", "input", "qc", "normalization", "reduction", "clustering", "markers", "output")) {
    validate_section(section)
  }

  validate_string("project.name", required = FALSE)
  validate_string("project.organism", allowed = c("human", "mouse"))
  validate_number("project.seed", integer = TRUE)

  validate_string("input.type", allowed = c("10x", "rds"))
  validate_string("input.path")
  validate_string("input.project_name", required = FALSE)
  validate_number("input.min_cells", integer = TRUE, min = 0)

  validate_number("qc.min_features", integer = TRUE, min = 0)
  validate_number("qc.max_features", integer = TRUE, min = 0)
  validate_number("qc.max_percent_mt", min = 0, max = 100)

  validate_string("normalization.method", allowed = c("lognormalize", "clr", "rc"))
  validate_number("normalization.scale_factor", positive = TRUE)
  validate_string("normalization.variable_features", allowed = c("vst", "mean.var.plot", "dispersion"))
  validate_number("normalization.n_variable_features", integer = TRUE, min = 1)

  validate_number("reduction.n_pcs", integer = TRUE, min = 2)
  validate_number("reduction.umap_dims", integer = TRUE, min = 2)

  validate_number("clustering.resolution", min = 0)

  validate_number("markers.logfc_threshold", min = 0)
  validate_number("markers.min_pct", min = 0, max = 1)
  validate_bool("markers.only_pos")

  validate_string("output.base_dir")
  validate_bool("output.save_seurat_object")

  run_name_value <- get_config_value("output.run_name")
  if (!is.null(run_name_value) && (!is_scalar_value(run_name_value) || !is.character(run_name_value))) {
    add_error("output.run_name must be null or a string.")
  }

  min_features_value <- get_config_value("qc.min_features")
  max_features_value <- get_config_value("qc.max_features")
  if (
    is.numeric(min_features_value) &&
      is.numeric(max_features_value) &&
      length(min_features_value) == 1 &&
      length(max_features_value) == 1 &&
      !is.na(min_features_value) &&
      !is.na(max_features_value) &&
      min_features_value > max_features_value
  ) {
    add_error("qc.min_features cannot be greater than qc.max_features.")
  }

  input_type_value <- get_config_value("input.type")
  input_path_value <- get_config_value("input.path")
  if (is_scalar_value(input_path_value) && is.character(input_path_value)) {
    if (!file.exists(input_path_value)) {
      add_error(paste0("input.path does not exist: ", normalizePath(input_path_value, mustWork = FALSE)))
    }
    if (
      is_scalar_value(input_type_value) &&
        is.character(input_type_value) &&
        tolower(trimws(input_type_value)) == "rds"
    ) {
      if (dir.exists(input_path_value)) {
        add_error("For input.type 'rds', input.path must point to a .rds file, not a directory.")
      }
      if (!grepl("\\.rds$", input_path_value, ignore.case = TRUE)) {
        add_error("For input.type 'rds', input.path must point to a file ending in .rds.")
      }
    }
  }

  if (length(errors) > 0) {
    abort(paste0("Config validation failed:\n- ", paste(errors, collapse = "\n- ")))
  }

  invisible(TRUE)
}

validate_config(params)

required_packages <- c(
  "jsonlite",
  "Seurat",
  "SeuratObject",
  "Matrix",
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

project_name <- get_param("input.project_name", get_param("project.name", "scPipeR_sample"))
input_type <- tolower(trimws(as.character(get_param("input.type", "10x"))))
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
feature_count_column <- NULL
count_column <- NULL

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
    median_features = stats::median(object[[feature_count_column]][, 1]),
    median_counts = stats::median(object[[count_column]][, 1]),
    median_percent_mt = stats::median(object$percent.mt),
    mean_percent_mt = mean(object$percent.mt)
  )
}

get_counts_matrix <- function(object, assay) {
  tryCatch(
    SeuratObject::GetAssayData(object = object, assay = assay, layer = "counts"),
    error = function(layer_error) {
      tryCatch(
        SeuratObject::GetAssayData(object = object, assay = assay, slot = "counts"),
        error = function(slot_error) {
          abort(paste0(
            "Unable to read counts from assay '",
            assay,
            "'. The .rds input must contain counts data for QC filtering."
          ))
        }
      )
    }
  )
}

prepare_qc_columns <- function(object) {
  assay <- SeuratObject::DefaultAssay(object)
  feature_col <- paste0("nFeature_", assay)
  count_col <- paste0("nCount_", assay)

  if (!all(c(feature_col, count_col) %in% colnames(object@meta.data))) {
    counts <- get_counts_matrix(object, assay)
    object[[feature_col]] <- Matrix::colSums(counts > 0)
    object[[count_col]] <- Matrix::colSums(counts)
  }

  list(
    object = object,
    feature_col = feature_col,
    count_col = count_col,
    assay = assay
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
  seu <- tryCatch(
    readRDS(input_path),
    error = function(err) {
      abort(paste0("Unable to read .rds input: ", conditionMessage(err)))
    }
  )
  if (!inherits(seu, "Seurat")) {
    abort("The .rds input must contain a Seurat object.")
  }
  if (ncol(seu) == 0 || nrow(seu) == 0) {
    abort("The .rds Seurat object must contain at least one cell and one feature.")
  }
}

qc_prep <- prepare_qc_columns(seu)
seu <- qc_prep$object
feature_count_column <- qc_prep$feature_col
count_column <- qc_prep$count_col
active_assay <- qc_prep$assay

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
  "Filtering cells with {feature_count_column} >= {min_features}, {feature_count_column} <= {max_features}, percent.mt <= {max_percent_mt}."
)
cell_metadata <- seu@meta.data
keep_cells <- rownames(cell_metadata)[
  cell_metadata[[feature_count_column]] >= min_features &
    cell_metadata[[feature_count_column]] <= max_features &
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
  active_assay = active_assay,
  organism = organism,
  seed = seed,
  parameters = list(
    min_cells = min_cells,
    feature_count_column = feature_count_column,
    count_column = count_column,
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
  paste0("Input path: ", normalizePath(input_path, mustWork = FALSE)),
  paste0("Input type: ", input_type),
  paste0("Active assay: ", active_assay),
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
