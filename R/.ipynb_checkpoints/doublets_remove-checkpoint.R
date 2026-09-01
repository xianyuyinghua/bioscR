#!/usr/bin/env Rscript

# Typical usage:
# seu <- doublet_filter(seu)
# seu <- doublet_filter(seu, assay = "RNA", sample_col = "sample",
#                       remove_doublets = FALSE)

doublets_remove <- function(seu,
                            assay = "RNA",
                            sample_col = "sample",
                            remove_doublets = TRUE,
                            seed = 123,
                            expected_doublet_rate = NULL,
                            verbose = TRUE,
                            ...
                           ){

  if (!inherits(seu, "Seurat")) {
    stop("'seu' must be a Seurat object.", call. = FALSE)
  }
  if (!is.character(assay) || length(assay) != 1L || !assay %in% names(seu@assays)) {
    stop("Assay '", assay, "' is absent from the Seurat object.", call. = FALSE)
  }
  if (!is.logical(remove_doublets) || length(remove_doublets) != 1L ||
      is.na(remove_doublets)) {
    stop("'remove_doublets' must be TRUE or FALSE.", call. = FALSE)
  }

  required <- c("SeuratObject", "SingleCellExperiment", "SummarizedExperiment",
                "scDblFinder")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    stop(
      "Missing R package(s): ", paste(missing, collapse = ", "), ".\n",
      "Install them with:\n",
      "  if (!requireNamespace('BiocManager', quietly = TRUE)) ",
      "install.packages('BiocManager')\n",
      "  BiocManager::install(c('SingleCellExperiment', 'scDblFinder'))",
      call. = FALSE
    )
  }

  # Raw UMI counts: Seurat v5, with a v4 fallback.
  counts <- tryCatch(
    SeuratObject::LayerData(seu, assay = assay, layer = "counts"),
    error = function(e) NULL
  )
  if (is.null(counts)) {
    counts <- tryCatch(
      SeuratObject::GetAssayData(seu, assay = assay, slot = "counts"),
      error = function(e) NULL
    )
  }
  if (is.null(counts) || nrow(counts) == 0L || ncol(counts) == 0L) {
    stop("Assay '", assay, "' has no non-empty raw counts layer.", call. = FALSE)
  }

  # sample_col identifies independent libraries; use NULL for a single library.
  samples <- NULL
  if (!is.null(sample_col)) {
    if (!is.character(sample_col) || length(sample_col) != 1L ||
        !sample_col %in% colnames(seu[[]])) {
      stop("Metadata column '", sample_col, "' is absent from the Seurat object.",
           call. = FALSE)
    }
    sample_values <- as.character(seu[[sample_col, drop = TRUE]])
    if (anyNA(sample_values) || any(sample_values == "")) {
      stop("Metadata column '", sample_col, "' contains missing/empty values.",
           call. = FALSE)
    }
    samples <- factor(sample_values)
  }

  if (verbose) {
    message(
      "Running scDblFinder: ", ncol(seu), " cells; ",
      if (is.null(samples)) 1L else nlevels(samples), " sample(s); assay = '",
      assay, "'; sample_col = ",
      if (is.null(sample_col)) "NULL" else paste0("'", sample_col, "'"), "."
    )
  }

  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
  args <- list(sce, samples = samples, verbose = verbose)
  if (!is.null(expected_doublet_rate)) {
    args$dbr <- expected_doublet_rate
  }
  args <- c(args, list(...))

  set.seed(seed)
  sce <- do.call(scDblFinder::scDblFinder, args)

  # Write predictions back into Seurat metadata in matching cell order.
  result <- as.data.frame(SummarizedExperiment::colData(sce))
  result <- result[colnames(seu), , drop = FALSE]
  result_columns <- grep("^scDblFinder\\.", colnames(result), value = TRUE)
  for (column in result_columns) {
    seu[[column]] <- result[[column]]
  }
  seu$is_doublet <- seu$scDblFinder.class == "doublet"
  seu$is_singlet <- seu$scDblFinder.class == "singlet"

  if (verbose) {
    message("Doublet prediction:")
    print(table(seu$scDblFinder.class, useNA = "ifany"))
  }

  # FALSE only annotates; TRUE returns a Seurat object containing singlets only.
  if (!remove_doublets) {
    return(seu)
  }

  keep <- !is.na(seu$is_singlet) & seu$is_singlet
  filtered <- subset(seu, cells = colnames(seu)[keep])
  if (verbose) {
    message("Retained ", ncol(filtered), " singlets; removed ", sum(!keep),
            " doublet/uncalled cells.")
  }
        
  return(filtered)
}
