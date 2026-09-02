#' Calculate observed-to-expected cell-distribution preference
#'
#' Calculates the observed-to-expected ratio (Ro/e) for every cell-type and
#' group combination in a Seurat object or cell-level metadata table. The
#' calculation follows the pooled-cell STARTRAC-dist approach and returns both
#' numerical results and an annotated heatmap.
#'
#' @param object A Seurat object or a `data.frame` containing cell-level
#'   metadata.
#' @param group_col A character scalar naming the metadata column that defines
#'   the biological groups, conditions, tissues, or other categories whose
#'   cell distributions are compared.
#' @param celltype_col A character scalar naming the metadata column containing
#'   cell-type or cell-cluster annotations.
#' @param p_adjust_method A multiple-testing correction method passed to
#'   [stats::p.adjust()]. The default is `"BH"`, matching the current STARTRAC
#'   implementation. For example, use `"bonferroni"` for Bonferroni correction.
#' @param na_action How cells with missing or blank values in `group_col` or
#'   `celltype_col` are handled. With `"omit"` (the default), these cells are
#'   excluded. With `"fail"`, the function stops and reports their number.
#' @param min_expected A non-negative numeric scalar used only as a warning
#'   threshold for expected cell counts. Combinations below this value are not
#'   removed. Defaults to `5`; use `0` to disable the warning in practice.
#' @param roe_limits A numeric vector of length two giving the heatmap colour
#'   limits. Outlying values are squished for colour display only; calculated
#'   values and labels remain unchanged. Defaults to `c(0, 2)`.
#' @param nudge_y A numeric scalar controlling the vertical position of the
#'   Ro/e text in each heatmap tile. Significance symbols are placed at
#'   `nudge_y - 0.4`. Defaults to `0.15`.
#'
#' @details
#' The function constructs a cell-type by group contingency table. For cell
#' type *i* and group *j*, the expected count under independence is
#' \deqn{E_{ij} = \frac{N_{i+}N_{+j}}{N},}
#' and the distribution-preference score is
#' \deqn{R_{o/e,ij} = \frac{O_{ij}}{E_{ij}}.}
#' A value greater than 1 indicates relative enrichment, a value below 1
#' indicates relative depletion, and a value equal to 1 agrees with the random
#' expectation determined by the marginal cell counts.
#'
#' For each cell-type/group combination, a two-sided Fisher exact test compares
#' that combination against all other cells. P values across all combinations
#' are adjusted together using `p_adjust_method`. Heatmap symbols are based on
#' adjusted P values: `*` for P < 0.05, `**` for P < 0.01, and `***` for
#' P < 0.001.
#'
#' This pooled cell-level analysis does not account for donor, sample, or other
#' biological-replicate structure. It is therefore best treated as a
#' descriptive measure of distribution preference. Sparse combinations should
#' also be interpreted cautiously.
#'
#' @return An object of class `"roe_result"`, containing:
#' \describe{
#'   \item{observed}{Observed cell-count matrix.}
#'   \item{expected}{Expected-count matrix under independence.}
#'   \item{roe}{Observed-to-expected ratio matrix.}
#'   \item{p_value}{Two-sided Fisher exact-test P-value matrix.}
#'   \item{p_adjusted}{Multiple-testing-adjusted P-value matrix.}
#'   \item{odds_ratio}{Fisher exact-test odds-ratio matrix.}
#'   \item{long}{Long-format results including counts, Ro/e, P values,
#'   preferences, and significance symbols.}
#'   \item{plot}{A `ggplot` heatmap annotated with Ro/e values and significance
#'   symbols.}
#'   \item{diagnostics}{Cell counts, adjustment method, global chi-squared
#'   P value, and low-expected-count diagnostics.}
#' }
#'
#' @examples
#' \dontrun{
#' seu <- readRDS("03_Seurat_Celltype.rds")
#' roe_result <- roe_calculate(
#'   object = seu,
#'   group_col = "group",
#'   celltype_col = "Celltype"
#' )
#'
#' roe_result$roe
#' roe_result$plot
#'
#' ggplot2::ggsave(
#'   "RoE_heatmap.pdf",
#'   plot = roe_result$plot,
#'   width = 7,
#'   height = 10
#' )
#' }
#'
#' @seealso [stats::chisq.test()], [stats::fisher.test()], [stats::p.adjust()]
#' @export
roe_calculate <- function(object,
                          group_col,
                          celltype_col,
                          p_adjust_method = "BH",
                          na_action = c("omit", "fail"), # 控制缺失注释的处理方式
                          min_expected = 5, # 低期望细胞数的警告阈值
                          roe_limits = c(0, 2),
                          nudge_y = 0.15
                         ){
  na_action <- match.arg(na_action)

  if (!is.character(group_col) || length(group_col) != 1L || is.na(group_col)) {
    stop("`group_col` must be one non-missing character string.", call. = FALSE)
  }
  if (!is.character(celltype_col) || length(celltype_col) != 1L || is.na(celltype_col)) {
    stop("`celltype_col` must be one non-missing character string.", call. = FALSE)
  }
  if (!is.numeric(min_expected) || length(min_expected) != 1L ||
      is.na(min_expected) || min_expected < 0) {
    stop("`min_expected` must be one non-negative number.", call. = FALSE)
  }
  if (!is.numeric(roe_limits) || length(roe_limits) != 2L ||
      anyNA(roe_limits) || roe_limits[1] >= roe_limits[2]) {
    stop("`roe_limits` must contain two increasing numeric values.",
         call. = FALSE)
  }
  valid_adjust <- stats::p.adjust.methods
  if (!is.character(p_adjust_method) || length(p_adjust_method) != 1L ||
      !p_adjust_method %in% valid_adjust) {
    stop("`p_adjust_method` must be one of: ",
         paste(valid_adjust, collapse = ", "), ".", call. = FALSE)
  }

  meta <- if (inherits(object, "Seurat")) {
    object[[]]
  } else if (is.data.frame(object)) {
    object
  } else {
    stop("`object` must be a Seurat object or a data.frame.", call. = FALSE)
  }

  missing_cols <- setdiff(c(group_col, celltype_col), colnames(meta))
  if (length(missing_cols)) {
    stop("Metadata column(s) not found: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  dat <- data.frame(
    celltype = as.character(meta[[celltype_col]]),
    group = as.character(meta[[group_col]]),
    stringsAsFactors = FALSE
  )
  bad <- is.na(dat$celltype) | is.na(dat$group) |
    !nzchar(trimws(dat$celltype)) | !nzchar(trimws(dat$group))
  n_omitted <- sum(bad)
  if (n_omitted && na_action == "fail") {
    stop(n_omitted, " cells have missing or blank group/cell-type values.",
         call. = FALSE)
  }
  dat <- dat[!bad, , drop = FALSE]
  if (!nrow(dat)) stop("No annotated cells remain for analysis.", call. = FALSE)
  if (length(unique(dat$celltype)) < 2L || length(unique(dat$group)) < 2L) {
    stop("Ro/e requires at least two cell types and two groups.", call. = FALSE)
  }

  observed <- table(dat$celltype, dat$group, dnn = c("celltype", "group"))
  observed <- unclass(observed)
  expected <- outer(rowSums(observed), colSums(observed)) / sum(observed)
  dimnames(expected) <- dimnames(observed)
  roe <- observed / expected

  p_value <- odds_ratio <- matrix(
    NA_real_, nrow(observed), ncol(observed), dimnames = dimnames(observed)
  )
  total <- sum(observed)
  row_total <- rowSums(observed)
  col_total <- colSums(observed)
  for (i in seq_len(nrow(observed))) {
    for (j in seq_len(ncol(observed))) {
      a <- observed[i, j]
      tab_2x2 <- matrix(c(
        a,
        row_total[i] - a,
        col_total[j] - a,
        total - row_total[i] - col_total[j] + a
      ), nrow = 2L, byrow = TRUE)
      ft <- stats::fisher.test(tab_2x2, alternative = "two.sided")
      p_value[i, j] <- ft$p.value
      odds_ratio[i, j] <- unname(ft$estimate)
    }
  }
  p_adjusted <- matrix(
    stats::p.adjust(as.vector(p_value), method = p_adjust_method),
    nrow = nrow(p_value), ncol = ncol(p_value), dimnames = dimnames(p_value)
  )

  long <- expand.grid(
    celltype = rownames(observed), group = colnames(observed),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  long$observed <- as.vector(observed)
  long$expected <- as.vector(expected)
  long$roe <- as.vector(roe)
  long$odds_ratio <- as.vector(odds_ratio)
  long$p_value <- as.vector(p_value)
  long$p_adjusted <- as.vector(p_adjusted)
  long$preference <- ifelse(long$roe > 1, "enriched",
                            ifelse(long$roe < 1, "depleted", "as_expected"))
  long$significance <- ifelse(
    long$p_adjusted < 0.001, "***",
    ifelse(long$p_adjusted < 0.01, "**",
           ifelse(long$p_adjusted < 0.05, "*", ""))
  )
  long$heatmap_label <- paste0(sprintf("%.2f", long$roe),"\n", long$significance)
  long$group <- factor(long$group,levels = levels(meta[[group_col]]))
    
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("scales", quietly = TRUE)) {
    stop("Packages `ggplot2` and `scales` are required for `result$plot`.",
         call. = FALSE)
  }
    
  # Heatmap
  heatmap_plot <- ggplot2::ggplot(long, ggplot2::aes(x = group, y = celltype, fill = roe)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
    #ggplot2::geom_text(ggplot2::aes(label = heatmap_label),size = 4,lineheight = 0.8) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", roe)),size = 4,colour = "black",nudge_y = nudge_y) +
    ggplot2::geom_text(ggplot2::aes(label = significance),size = 5,colour = "black",nudge_y = nudge_y - 0.4)+
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                  midpoint = 1, limits = roe_limits, oob = scales::squish,
                                  name = "Ro/e") +
    ggplot2::scale_y_discrete(limits = rev(rownames(observed))) +
    ggplot2::labs(#title = "Cell distribution preference (Ro/e)",
                  #subtitle = paste0("Fisher's exact test; ", p_adjust_method,"-adjusted P values"),
                  x = "",y = "") +
    ggplot2::theme_classic(base_size = 18) +
    ggplot2::theme(axis.line = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   axis.text = ggplot2::element_text(face = "bold",colour = "black"),
                   plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5,face = "italic",size = 15),
                   plot.margin = ggplot2::margin(b = 10,l = 10,t = 20,r = 10)
                  )

  low_expected <- expected < min_expected
  if (any(low_expected)) {
    warning(sum(low_expected), " of ", length(expected),
            " expected counts are below ", min_expected,
            "; interpret sparse combinations cautiously.", call. = FALSE)
  }

  structure(list(
    observed = observed,
    expected = expected,
    roe = roe,
    p_value = p_value,
    p_adjusted = p_adjusted,
    odds_ratio = odds_ratio,
    long = long,
    plot = heatmap_plot,
    diagnostics = list(
      n_cells_used = nrow(dat),
      n_cells_omitted = n_omitted,
      p_adjust_method = p_adjust_method,
      global_chisq_p_value = suppressWarnings(stats::chisq.test(observed)$p.value),
      n_expected_below_threshold = sum(low_expected),
      min_expected = min(expected)
    )
  ), class = "roe_result")
}

#' Print a Ro/e analysis result
#'
#' Prints the number of cells, contingency-table dimensions, global
#' chi-squared P value, and rounded Ro/e matrix.
#'
#' @param x An object of class `"roe_result"` returned by [roe_calculate()].
#' @param ... Additional arguments, currently unused.
#'
#' @return `x`, invisibly.
#' @method print roe_result
#' @export
print.roe_result <- function(x, ...) {
  cat("Ro/e cell-distribution analysis\n")
  cat("  cells used:", x$diagnostics$n_cells_used, "\n")
  cat("  cell types x groups:", nrow(x$roe), "x", ncol(x$roe), "\n")
  cat("  global chi-square P:",
      format.pval(x$diagnostics$global_chisq_p_value), "\n")
  print(round(x$roe, 3))
  invisible(x)
}
