#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(jsonlite))

args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag) {
  positions <- which(args == flag)
  if (length(positions) != 1L || positions + 1L > length(args)) {
    stop(sprintf("Expected exactly one %s VALUE argument", flag))
  }
  args[[positions + 1L]]
}

method <- value_after("--method")
counts_path <- value_after("--counts")
metadata_path <- value_after("--metadata")
factors_path <- value_after("--output-factors")
receipt_path <- value_after("--output-receipt")
supported <- c("tmm_master_dhs_v1", "tmm_background_10kb_v1")
if (!(method %in% supported)) {
  stop(sprintf("Unsupported TMM method: %s", method))
}

metadata <- read.delim(
  metadata_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
expected_metadata <- c(
  "library_id", "assay", "context", "total_units", "signal_sha256",
  "count_source", "count_source_sha256"
)
if (!identical(colnames(metadata), expected_metadata)) {
  stop("Unexpected TMM metadata columns")
}
if (
  nrow(metadata) < 4L || anyDuplicated(metadata$library_id) ||
  any(!(metadata$assay %in% c("atac", "h3k27ac"))) ||
  any(!is.finite(metadata$total_units)) || any(metadata$total_units <= 0)
) {
  stop("Invalid TMM library metadata")
}

counts <- read.delim(
  gzfile(counts_path),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
base_columns <- c("feature_id", "chrom", "start", "end", "width_bp")
if (
  nrow(counts) < 1L ||
  !identical(colnames(counts)[seq_along(base_columns)], base_columns) ||
  !identical(colnames(counts)[-(seq_along(base_columns))], metadata$library_id)
) {
  stop("TMM count matrix columns do not match metadata")
}
count_matrix <- as.matrix(counts[, metadata$library_id, drop = FALSE])
storage.mode(count_matrix) <- "double"
if (
  any(!is.finite(count_matrix)) || any(count_matrix < 0) ||
  any(count_matrix != floor(count_matrix))
) {
  stop("TMM count matrix must contain non-negative integers")
}

factor_rows <- list()
assay_receipts <- list()
for (assay_name in c("atac", "h3k27ac")) {
  positions <- which(metadata$assay == assay_name)
  if (length(positions) < 2L) {
    stop(sprintf("TMM requires at least two %s libraries", assay_name))
  }
  assay_counts <- count_matrix[, positions, drop = FALSE]
  retained_features <- rowSums(assay_counts) > 0
  feature_count <- sum(retained_features)
  if (feature_count < 1L) {
    stop(sprintf("TMM has no nonzero %s features", assay_name))
  }
  assay_metadata <- metadata[positions, , drop = FALSE]
  factors <- edgeR::calcNormFactors(
    assay_counts[retained_features, , drop = FALSE],
    lib.size = assay_metadata$total_units,
    method = "TMM",
    logratioTrim = 0.30,
    sumTrim = 0.05,
    doWeighting = TRUE,
    Acutoff = -1e10
  )
  effective_sizes <- assay_metadata$total_units * factors
  factor_rows[[assay_name]] <- data.frame(
    library_id = assay_metadata$library_id,
    assay = assay_metadata$assay,
    context = assay_metadata$context,
    total_units = assay_metadata$total_units,
    feature_count = rep(feature_count, length(positions)),
    tmm_normalization_factor = factors,
    effective_library_size = effective_sizes,
    normalization_method = rep(method, length(positions)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  assay_receipts[[assay_name]] <- list(
    library_count = length(positions),
    nonzero_feature_count = feature_count,
    normalization_factor_product = unname(prod(factors))
  )
}
factor_table <- do.call(rbind, factor_rows)

write_atomic_table <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  old_digits <- getOption("digits")
  on.exit(options(digits = old_digits), add = TRUE)
  options(digits = 17)
  write.table(
    value,
    file = temporary,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )
  if (!file.rename(temporary, path)) {
    stop(sprintf("Could not atomically replace %s", path))
  }
}

write_atomic_json <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  jsonlite::write_json(
    value,
    path = temporary,
    auto_unbox = TRUE,
    pretty = TRUE,
    digits = 17,
    null = "null"
  )
  cat("\n", file = temporary, append = TRUE)
  if (!file.rename(temporary, path)) {
    stop(sprintf("Could not atomically replace %s", path))
  }
}

write_atomic_table(factor_table, factors_path)
write_atomic_json(
  list(
    status = "ok",
    schema_version = 1L,
    normalization_method = method,
    implementation = "edgeR::calcNormFactors",
    edgeR_version = as.character(packageVersion("edgeR")),
    R_version = R.version.string,
    parameters = list(
      method = "TMM",
      logratioTrim = 0.30,
      sumTrim = 0.05,
      doWeighting = TRUE,
      Acutoff = -1e10,
      library_size = "total filtered assay units",
      assay_stratification = TRUE
    ),
    assays = assay_receipts
  ),
  receipt_path
)
