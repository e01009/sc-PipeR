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

  for (section in c("project", "input", "qc", "normalization", "modules", "reduction", "clustering", "markers", "plots", "report", "output")) {
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

  validate_string("normalization.method", allowed = c("lognormalize", "clr", "rc", "sctransform", "scran"))
  validate_number("normalization.scale_factor", positive = TRUE)
  validate_string("normalization.variable_features", allowed = c("vst", "mean.var.plot", "dispersion"))
  validate_number("normalization.n_variable_features", integer = TRUE, min = 1)

  validate_bool("modules.run_scater_qc")
  validate_bool("modules.run_doublet_detection")
  validate_bool("modules.run_annotation")
  validate_string(
    "modules.annotation_reference",
    allowed = c("hpca", "blueprint_encode", "dice", "monaco", "mouse_rnaseq", "immgen")
  )

  validate_number("reduction.n_pcs", integer = TRUE, min = 2)
  validate_number("reduction.umap_dims", integer = TRUE, min = 2)

  validate_number("clustering.resolution", min = 0)

  validate_number("markers.logfc_threshold", min = 0)
  validate_number("markers.min_pct", min = 0, max = 1)
  validate_bool("markers.only_pos")

  validate_bool("plots.enabled")
  validate_bool("plots.export_pdf")
  validate_number("plots.width", positive = TRUE)
  validate_number("plots.height", positive = TRUE)
  validate_number("plots.dpi", integer = TRUE, min = 72)
  validate_number("plots.top_marker_count", integer = TRUE, min = 1)

  validate_bool("report.enabled")

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
  "SingleCellExperiment",
  "SummarizedExperiment",
  "Matrix",
  "fs",
  "glue",
  "readr",
  "dplyr",
  "tibble",
  "sessioninfo"
)
if (isTRUE(get_param("plots.enabled", TRUE))) {
  required_packages <- c(required_packages, "ggplot2")
}
if (isTRUE(get_param("report.enabled", TRUE))) {
  required_packages <- c(required_packages, "rmarkdown")
}
normalization_method_for_packages <- tolower(trimws(as.character(get_param("normalization.method", "LogNormalize"))))
if (
  normalization_method_for_packages == "scran" ||
    isTRUE(get_param("modules.run_scater_qc", TRUE)) ||
    isTRUE(get_param("modules.run_doublet_detection", TRUE)) ||
    isTRUE(get_param("modules.run_annotation", TRUE))
) {
  required_packages <- c(required_packages, "SingleCellExperiment", "SummarizedExperiment")
}
if (normalization_method_for_packages == "scran" || isTRUE(get_param("modules.run_scater_qc", TRUE))) {
  required_packages <- c(required_packages, "scater")
}
if (normalization_method_for_packages == "scran") {
  required_packages <- c(required_packages, "scran")
}
if (normalization_method_for_packages == "sctransform") {
  required_packages <- c(required_packages, "sctransform")
}
if (isTRUE(get_param("modules.run_doublet_detection", TRUE))) {
  required_packages <- c(required_packages, "scDblFinder")
}
if (isTRUE(get_param("modules.run_annotation", TRUE))) {
  required_packages <- c(required_packages, "SingleR", "celldex")
}
required_packages <- unique(required_packages)

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
normalization_method_lower <- tolower(trimws(as.character(normalization_method)))
scale_factor <- as.numeric(get_param("normalization.scale_factor", 10000))
variable_features_method <- get_param("normalization.variable_features", "vst")
n_variable_features <- as.integer(get_param("normalization.n_variable_features", 2000))

run_scater_qc <- isTRUE(get_param("modules.run_scater_qc", TRUE))
run_doublet_detection <- isTRUE(get_param("modules.run_doublet_detection", TRUE))
run_annotation <- isTRUE(get_param("modules.run_annotation", TRUE))
annotation_reference <- tolower(trimws(as.character(get_param("modules.annotation_reference", "hpca"))))

n_pcs <- as.integer(get_param("reduction.n_pcs", 30))
umap_dims <- as.integer(get_param("reduction.umap_dims", n_pcs))
clustering_resolution <- as.numeric(get_param("clustering.resolution", 0.5))

marker_logfc_threshold <- as.numeric(get_param("markers.logfc_threshold", 0.25))
marker_min_pct <- as.numeric(get_param("markers.min_pct", 0.1))
marker_only_pos <- isTRUE(get_param("markers.only_pos", TRUE))

plots_enabled <- isTRUE(get_param("plots.enabled", TRUE))
plots_export_pdf <- isTRUE(get_param("plots.export_pdf", FALSE))
plot_width <- as.numeric(get_param("plots.width", 8))
plot_height <- as.numeric(get_param("plots.height", 6))
plot_dpi <- as.integer(get_param("plots.dpi", 300))
top_marker_count <- as.integer(get_param("plots.top_marker_count", 10))

report_enabled <- isTRUE(get_param("report.enabled", TRUE))

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
  file.path(run_dir, "plots"),
  file.path(run_dir, "plots", "qc"),
  file.path(run_dir, "plots", "reduction"),
  file.path(run_dir, "plots", "markers"),
  file.path(run_dir, "report"),
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
    "sc-PipeR version 0.2 run",
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

seurat_to_sce <- function(object, assay = SeuratObject::DefaultAssay(object)) {
  tryCatch(
    Seurat::as.SingleCellExperiment(object, assay = assay),
    error = function(err) {
      abort(paste0("Unable to convert Seurat object to SingleCellExperiment: ", conditionMessage(err)))
    }
  )
}

coldata_to_data_frame <- function(sce_object) {
  as.data.frame(SummarizedExperiment::colData(sce_object))
}

add_sce_columns_to_seurat <- function(object, sce_object, columns, prefix) {
  coldata <- coldata_to_data_frame(sce_object)
  matching_columns <- intersect(columns, colnames(coldata))

  for (column in matching_columns) {
    object[[paste0(prefix, column)]] <- coldata[colnames(object), column]
  }

  object
}

set_assay_data <- function(object, assay, data_layer) {
  tryCatch(
    SeuratObject::SetAssayData(object, assay = assay, layer = "data", new.data = data_layer),
    error = function(layer_error) {
      SeuratObject::SetAssayData(object, assay = assay, slot = "data", new.data = data_layer)
    }
  )
}

run_scater_qc_module <- function(object) {
  log_message("Running scater QC metrics.")
  sce_object <- seurat_to_sce(object, assay = active_assay)
  mitochondrial_features <- grepl(mt_pattern, rownames(sce_object))
  sce_object <- scater::addPerCellQC(sce_object, subsets = list(mito = mitochondrial_features))
  scater_columns <- c("sum", "detected", "subsets_mito_sum", "subsets_mito_detected", "subsets_mito_percent")
  object <- add_sce_columns_to_seurat(object, sce_object, scater_columns, "scater_")
  scater_qc_table <- coldata_to_data_frame(sce_object)[colnames(object), scater_columns, drop = FALSE]
  scater_qc_table <- tibble::rownames_to_column(scater_qc_table, var = "cell_barcode")

  list(object = object, sce = sce_object, table = scater_qc_table)
}

run_doublet_detection_module <- function(object, sce_object = NULL) {
  log_message("Running scDblFinder doublet detection.")
  if (is.null(sce_object)) {
    sce_object <- seurat_to_sce(object, assay = active_assay)
  }

  sce_object <- scDblFinder::scDblFinder(sce_object)
  doublet_columns <- grep("^scDblFinder\\.", colnames(coldata_to_data_frame(sce_object)), value = TRUE)
  object <- add_sce_columns_to_seurat(object, sce_object, doublet_columns, "")
  doublet_table <- coldata_to_data_frame(sce_object)[colnames(object), doublet_columns, drop = FALSE]
  doublet_table <- tibble::rownames_to_column(doublet_table, var = "cell_barcode")

  list(object = object, sce = sce_object, table = doublet_table)
}

run_scran_normalization <- function(object) {
  log_message("Normalizing data with scran.")
  sce_object <- seurat_to_sce(object, assay = active_assay)
  clusters <- scran::quickCluster(sce_object)
  sce_object <- scran::computeSumFactors(sce_object, clusters = clusters)
  sce_object <- scater::logNormCounts(sce_object)
  object <- set_assay_data(object, active_assay, SummarizedExperiment::assay(sce_object, "logcounts"))

  size_factor_table <- tibble::tibble(
    cell_barcode = colnames(object),
    scran_size_factor = SingleCellExperiment::sizeFactors(sce_object)
  )
  object[["scran_size_factor"]] <- size_factor_table$scran_size_factor

  list(object = object, sce = sce_object, table = size_factor_table)
}

load_celldex_reference <- function(reference_name) {
  switch(
    reference_name,
    hpca = celldex::HumanPrimaryCellAtlasData(),
    blueprint_encode = celldex::BlueprintEncodeData(),
    dice = celldex::DatabaseImmuneCellExpressionData(),
    monaco = celldex::MonacoImmuneData(),
    mouse_rnaseq = celldex::MouseRNAseqData(),
    immgen = celldex::ImmGenData(),
    abort(paste0("Unsupported annotation reference: ", reference_name))
  )
}

reference_label_column <- function(reference) {
  reference_coldata <- coldata_to_data_frame(reference)
  if ("label.main" %in% colnames(reference_coldata)) {
    return("label.main")
  }
  if ("label.fine" %in% colnames(reference_coldata)) {
    return("label.fine")
  }
  abort("The selected celldex reference does not contain label.main or label.fine.")
}

run_annotation_module <- function(object) {
  log_message("Running SingleR annotation.")
  annotation_sce <- seurat_to_sce(object, assay = SeuratObject::DefaultAssay(object))
  reference <- load_celldex_reference(annotation_reference)
  label_column <- reference_label_column(reference)
  reference_labels <- coldata_to_data_frame(reference)[[label_column]]
  predictions <- SingleR::SingleR(test = annotation_sce, ref = reference, labels = reference_labels)
  prediction_table <- as.data.frame(predictions)
  prediction_table <- tibble::rownames_to_column(prediction_table, var = "cell_barcode")

  object[["SingleR_label"]] <- prediction_table$labels[match(colnames(object), prediction_table$cell_barcode)]
  if ("pruned.labels" %in% colnames(prediction_table)) {
    object[["SingleR_pruned_label"]] <- prediction_table$pruned.labels[match(colnames(object), prediction_table$cell_barcode)]
  }

  list(object = object, table = prediction_table, reference = annotation_reference, label_column = label_column)
}

save_plot <- function(plot, path_without_ext, width = plot_width, height = plot_height) {
  png_path <- paste0(path_without_ext, ".png")
  ggplot2::ggsave(
    filename = png_path,
    plot = plot,
    width = width,
    height = height,
    dpi = plot_dpi,
    units = "in"
  )

  output_paths <- png_path
  if (plots_export_pdf) {
    pdf_path <- paste0(path_without_ext, ".pdf")
    ggplot2::ggsave(
      filename = pdf_path,
      plot = plot,
      width = width,
      height = height,
      units = "in"
    )
    output_paths <- c(output_paths, pdf_path)
  }

  output_paths
}

write_plot_outputs <- function(object, marker_table) {
  plot_paths <- character()

  if (!plots_enabled) {
    return(plot_paths)
  }

  log_message("Writing plots.")

  qc_plot <- Seurat::VlnPlot(
    object,
    features = c(feature_count_column, count_column, "percent.mt"),
    ncol = 3,
    pt.size = 0.1
  )
  plot_paths <- c(plot_paths, save_plot(qc_plot, file.path(run_dir, "plots", "qc", "qc_violin"), width = 12, height = 5))

  scatter_feature_count <- Seurat::FeatureScatter(
    object,
    feature1 = count_column,
    feature2 = feature_count_column
  )
  plot_paths <- c(plot_paths, save_plot(scatter_feature_count, file.path(run_dir, "plots", "qc", "feature_count_scatter")))

  scatter_mt_count <- Seurat::FeatureScatter(
    object,
    feature1 = count_column,
    feature2 = "percent.mt"
  )
  plot_paths <- c(plot_paths, save_plot(scatter_mt_count, file.path(run_dir, "plots", "qc", "mitochondrial_count_scatter")))

  elbow_plot <- Seurat::ElbowPlot(object, ndims = available_pcs)
  plot_paths <- c(plot_paths, save_plot(elbow_plot, file.path(run_dir, "plots", "reduction", "pca_elbow")))

  umap_cluster_plot <- Seurat::DimPlot(object, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
  plot_paths <- c(plot_paths, save_plot(umap_cluster_plot, file.path(run_dir, "plots", "reduction", "umap_clusters")))

  umap_qc_plot <- Seurat::FeaturePlot(
    object,
    features = c(feature_count_column, count_column, "percent.mt"),
    reduction = "umap",
    ncol = 3
  )
  plot_paths <- c(plot_paths, save_plot(umap_qc_plot, file.path(run_dir, "plots", "reduction", "umap_qc_metrics"), width = 12, height = 5))

  if (nrow(marker_table) > 0 && "cluster" %in% colnames(marker_table) && "gene" %in% colnames(marker_table)) {
    marker_score_column <- if ("avg_log2FC" %in% colnames(marker_table)) "avg_log2FC" else "avg_logFC"
    if (marker_score_column %in% colnames(marker_table)) {
      ordered_markers <- marker_table[order(marker_table$cluster, -marker_table[[marker_score_column]]), , drop = FALSE]
      top_markers <- dplyr::bind_rows(lapply(split(ordered_markers, ordered_markers$cluster), utils::head, top_marker_count))
      top_marker_genes <- unique(top_markers$gene)

      if (length(top_marker_genes) > 0) {
        marker_plot <- Seurat::DotPlot(object, features = top_marker_genes) +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
        plot_paths <- c(plot_paths, save_plot(marker_plot, file.path(run_dir, "plots", "markers", "top_marker_dotplot"), width = 12, height = 7))
      }
    }
  }

  plot_paths
}

relative_path <- function(path) {
  normalized_path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  normalized_run_dir <- normalizePath(run_dir, winslash = "/", mustWork = FALSE)
  sub(paste0("^", gsub("([\\^$.|?*+(){}\\[\\]\\\\])", "\\\\\\1", normalized_run_dir), "/?"), "", normalized_path)
}

write_report <- function(plot_paths) {
  if (!report_enabled) {
    return(NULL)
  }

  log_message("Writing HTML report.")

  report_rmd <- file.path(run_dir, "report", "sc_piper_report.Rmd")
  report_html <- file.path(run_dir, "report", "sc_piper_report.html")
  png_plot_paths <- plot_paths[grepl("\\.png$", plot_paths, ignore.case = TRUE)]
  plot_sections <- if (length(png_plot_paths) > 0) {
    paste0(
      "### ",
      tools::file_path_sans_ext(basename(png_plot_paths)),
      "\n\n![](../",
      vapply(png_plot_paths, relative_path, character(1)),
      ")\n"
    )
  } else {
    "No plots were exported for this run.\n"
  }
  module_sections <- character()
  if (nrow(scater_qc_table) > 0) {
    module_sections <- c(
      module_sections,
      "## scater QC Metrics",
      "",
      "```{r}",
      "readr::read_csv('../tables/scater_qc_metrics.csv', show_col_types = FALSE)",
      "```",
      ""
    )
  }
  if (nrow(doublet_table) > 0) {
    module_sections <- c(
      module_sections,
      "## scDblFinder Doublet Calls",
      "",
      "```{r}",
      "readr::read_csv('../tables/doublet_calls.csv', show_col_types = FALSE)",
      "```",
      ""
    )
  }
  if (nrow(scran_size_factor_table) > 0) {
    module_sections <- c(
      module_sections,
      "## scran Size Factors",
      "",
      "```{r}",
      "readr::read_csv('../tables/scran_size_factors.csv', show_col_types = FALSE)",
      "```",
      ""
    )
  }
  if (nrow(annotation_table) > 0) {
    module_sections <- c(
      module_sections,
      "## SingleR Annotation",
      "",
      "```{r}",
      "readr::read_csv('../tables/singleR_annotation.csv', show_col_types = FALSE)",
      "```",
      ""
    )
  }

  report_lines <- c(
    "---",
    "title: \"sc-PipeR version 0.2 report\"",
    "output:",
    "  html_document:",
    "    toc: true",
    "    toc_depth: 2",
    "---",
    "",
    "```{r setup, include=FALSE}",
    "knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = FALSE)",
    "```",
    "",
    "## Run Summary",
    "",
    paste0("- Project: ", project_name),
    paste0("- Input type: ", input_type),
    paste0("- Input path: ", normalizePath(input_path, mustWork = FALSE)),
    paste0("- Active assay: ", active_assay),
    paste0("- Cells before filtering: ", qc_before$cells),
    paste0("- Cells after filtering: ", qc_after$cells),
    paste0("- Genes/features: ", qc_after$genes),
    paste0("- PCs used: ", available_pcs),
    paste0("- UMAP dimensions used: ", umap_dims),
    paste0("- Clusters: ", dplyr::n_distinct(metadata$seurat_clusters)),
    paste0("- Marker rows: ", nrow(markers)),
    paste0("- scater QC metrics: ", if (nrow(scater_qc_table) > 0) "enabled" else "not run"),
    paste0("- Doublet calls: ", if (nrow(doublet_table) > 0) nrow(doublet_table) else "not run"),
    paste0("- SingleR annotations: ", if (nrow(annotation_table) > 0) nrow(annotation_table) else "not run"),
    "",
    "## QC Summary",
    "",
    "```{r}",
    "readr::read_csv('../tables/qc_summary.csv', show_col_types = FALSE)",
    "```",
    "",
    "## Cluster Counts",
    "",
    "```{r}",
    "readr::read_csv('../tables/cluster_counts.csv', show_col_types = FALSE)",
    "```",
    "",
    module_sections,
    "## Plots",
    "",
    plot_sections
  )

  writeLines(report_lines, con = report_rmd)
  rmarkdown::render(
    input = report_rmd,
    output_file = basename(report_html),
    output_dir = dirname(report_html),
    quiet = TRUE,
    envir = new.env(parent = globalenv())
  )

  report_html
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
    abort("Unsupported project.organism. Use 'human' or 'mouse'.")
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

scater_qc_table <- tibble::tibble()
doublet_table <- tibble::tibble()
scran_size_factor_table <- tibble::tibble()
annotation_table <- tibble::tibble()

sce <- seurat_to_sce(seu, assay = active_assay)

if (run_scater_qc) {
  scater_result <- run_scater_qc_module(seu)
  seu <- scater_result$object
  sce <- scater_result$sce
  scater_qc_table <- scater_result$table
}

if (run_doublet_detection) {
  doublet_result <- run_doublet_detection_module(seu, sce)
  seu <- doublet_result$object
  sce <- doublet_result$sce
  doublet_table <- doublet_result$table
}

if (normalization_method_lower %in% c("lognormalize", "clr", "rc")) {
  log_message("Normalizing data with Seurat {normalization_method}.")
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
} else if (normalization_method_lower == "sctransform") {
  log_message("Normalizing data with SCTransform.")
  seu <- Seurat::SCTransform(
    object = seu,
    variable.features.n = n_variable_features,
    verbose = FALSE
  )
  active_assay <- SeuratObject::DefaultAssay(seu)
} else if (normalization_method_lower == "scran") {
  scran_result <- run_scran_normalization(seu)
  seu <- scran_result$object
  sce <- scran_result$sce
  scran_size_factor_table <- scran_result$table

  log_message("Finding variable features.")
  seu <- Seurat::FindVariableFeatures(
    object = seu,
    selection.method = variable_features_method,
    nfeatures = n_variable_features,
    verbose = FALSE
  )

  log_message("Scaling data.")
  seu <- Seurat::ScaleData(seu, verbose = FALSE)
} else {
  abort(paste0("Unsupported normalization.method: ", normalization_method))
}

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

if (run_annotation) {
  annotation_result <- run_annotation_module(seu)
  seu <- annotation_result$object
  annotation_table <- annotation_result$table
}

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
    marker_only_pos = marker_only_pos,
    run_scater_qc = run_scater_qc,
    run_doublet_detection = run_doublet_detection,
    run_annotation = run_annotation,
    annotation_reference = annotation_reference,
    plots_enabled = plots_enabled,
    plots_export_pdf = plots_export_pdf,
    plot_width = plot_width,
    plot_height = plot_height,
    plot_dpi = plot_dpi,
    top_marker_count = top_marker_count,
    report_enabled = report_enabled
  )
)

sce <- seurat_to_sce(seu, assay = SeuratObject::DefaultAssay(seu))

log_message("Writing outputs.")
readr::write_csv(qc_table, file.path(run_dir, "tables", "qc_summary.csv"))
readr::write_csv(metadata, file.path(run_dir, "tables", "cell_metadata.csv"))
readr::write_csv(cluster_counts, file.path(run_dir, "tables", "cluster_counts.csv"))
readr::write_csv(markers, file.path(run_dir, "tables", "markers.csv"))
readr::write_csv(umap_embeddings, file.path(run_dir, "tables", "umap_embeddings.csv"))
if (nrow(scater_qc_table) > 0) {
  readr::write_csv(scater_qc_table, file.path(run_dir, "tables", "scater_qc_metrics.csv"))
}
if (nrow(doublet_table) > 0) {
  readr::write_csv(doublet_table, file.path(run_dir, "tables", "doublet_calls.csv"))
}
if (nrow(scran_size_factor_table) > 0) {
  readr::write_csv(scran_size_factor_table, file.path(run_dir, "tables", "scran_size_factors.csv"))
}
if (nrow(annotation_table) > 0) {
  readr::write_csv(annotation_table, file.path(run_dir, "tables", "singleR_annotation.csv"))
}
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
saveRDS(sce, file.path(run_dir, "objects", "single_cell_experiment.rds"))

plot_paths <- write_plot_outputs(seu, markers)
if (length(plot_paths) > 0) {
  readr::write_csv(
    tibble::tibble(path = vapply(plot_paths, relative_path, character(1))),
    file.path(run_dir, "plots", "plot_manifest.csv")
  )
}
report_path <- write_report(plot_paths)

writeLines(
  capture.output(sessioninfo::session_info()),
  con = file.path(run_dir, "logs", "sessioninfo.txt")
)

summary_lines <- c(
  "sc-PipeR version 0.2 run complete",
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
  paste0("Marker rows: ", nrow(markers)),
  paste0("scater QC metrics: ", if (nrow(scater_qc_table) > 0) "written" else "not run"),
  paste0(
    "Doublets detected: ",
    if (nrow(doublet_table) > 0 && "scDblFinder.class" %in% colnames(doublet_table)) {
      sum(doublet_table[["scDblFinder.class"]] == "doublet", na.rm = TRUE)
    } else {
      "not run"
    }
  ),
  paste0("SingleR annotations: ", if (nrow(annotation_table) > 0) nrow(annotation_table) else "not run"),
  paste0("Plot files: ", length(plot_paths)),
  paste0("HTML report: ", report_path %||% "not generated")
)
writeLines(summary_lines, con = file.path(run_dir, "logs", "summary.txt"))

log_message("Run complete. Output directory: {normalizePath(run_dir, mustWork = FALSE)}")
