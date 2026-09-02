% Manually maintained documentation for R/roe_calculate.R.
\name{roe_calculate}
\alias{roe_calculate}
\title{Calculate observed-to-expected cell-distribution preference}
\usage{
roe_calculate(
  object,
  group_col,
  celltype_col,
  p_adjust_method = "BH",
  na_action = c("omit", "fail"),
  min_expected = 5,
  roe_limits = c(0, 2),
  nudge_y = 0.15
)
}
\arguments{
\item{object}{A Seurat object or a \code{data.frame} containing cell-level
metadata.}

\item{group_col}{A character scalar naming the metadata column that defines
the biological groups, conditions, tissues, or other categories whose cell
distributions are compared.}

\item{celltype_col}{A character scalar naming the metadata column containing
cell-type or cell-cluster annotations.}

\item{p_adjust_method}{A multiple-testing correction method passed to
\code{stats::p.adjust()}. The default is \code{"BH"}, matching the current
STARTRAC implementation. For example, use \code{"bonferroni"} for Bonferroni
correction.}

\item{na_action}{How cells with missing or blank values in \code{group_col} or
\code{celltype_col} are handled. With \code{"omit"} (the default), these cells
are excluded. With \code{"fail"}, the function stops and reports their number.}

\item{min_expected}{A non-negative numeric scalar used only as a warning
threshold for expected cell counts. The default is \code{5}. Combinations below
this threshold are not removed and their Ro/e values and tests are still
returned. Set to \code{0} to disable this warning in practice.}

\item{roe_limits}{A numeric vector of length two giving the lower and upper
limits of the heatmap colour scale. Values outside these limits are squished
for colour display only; calculated Ro/e values and text labels are unchanged.
Defaults to \code{c(0, 2)}.}

\item{nudge_y}{A numeric scalar controlling the vertical position of the Ro/e
text inside each heatmap tile. The significance symbols are placed at
\code{nudge_y - 0.4}. Defaults to \code{0.15}.}
}
\value{
An object of class \code{"roe_result"}, implemented as a list with components:
\describe{
  \item{\code{observed}}{Matrix of observed cell counts for every cell-type and
  group combination.}
  \item{\code{expected}}{Matrix of expected counts under independence between
  cell type and group.}
  \item{\code{roe}}{Matrix of observed-to-expected ratios.}
  \item{\code{p_value}}{Matrix of two-sided Fisher exact-test P values.}
  \item{\code{p_adjusted}}{Matrix of multiple-testing-adjusted P values.}
  \item{\code{odds_ratio}}{Matrix of Fisher exact-test odds-ratio estimates.}
  \item{\code{long}}{A long-format data frame containing counts, Ro/e values,
  P values, distribution preference, significance symbols, and heatmap labels.}
  \item{\code{plot}}{A \code{ggplot} heatmap. Tile colours represent Ro/e;
  every tile contains its Ro/e value and, when significant, a symbol.}
  \item{\code{diagnostics}}{A list containing the number of cells used and
  omitted, adjustment method, global chi-squared P value, and diagnostics for
  low expected counts.}
}
}
\description{
Calculates the observed-to-expected ratio (Ro/e) for every cell-type and group
combination in a Seurat object or cell-level metadata table. The calculation
follows the pooled-cell STARTRAC-dist approach and returns both numerical
results and an annotated heatmap.
}
\details{
The function first constructs a cell-type by group contingency table. For cell
type \eqn{i} and group \eqn{j}, the expected count under independence is

\deqn{E_{ij} = \frac{N_{i+} N_{+j}}{N},}

and the distribution-preference score is

\deqn{R_{o/e,ij} = \frac{O_{ij}}{E_{ij}}.}

A value greater than 1 indicates relative enrichment, a value below 1 indicates
relative depletion, and a value equal to 1 indicates agreement with the random
expectation determined by the marginal cell counts.

For each cell-type and group combination, a two-by-two table contrasting that
combination against all other cells is evaluated with a two-sided Fisher exact
test. P values across all combinations are adjusted together using
\code{p_adjust_method}. Heatmap significance symbols use adjusted P values:
\code{*} for P < 0.05, \code{**} for P < 0.01, and \code{***} for P < 0.001.

This is a pooled cell-level analysis. It does not account for donor, sample, or
other biological replicate structure. Consequently, the result is best treated
as a descriptive measure of distribution preference rather than donor-level
statistical inference. Sparse cell-type/group combinations can also produce
unstable Ro/e values and should be interpreted cautiously.
}
\examples{
\dontrun{
seu <- readRDS("03_Seurat_Celltype.rds")

roe_result <- roe_calculate(
  object = seu,
  group_col = "group",
  celltype_col = "Celltype"
)

# Ro/e matrix and complete long-format results
roe_result$roe
head(roe_result$long)

# Display and save the annotated heatmap
roe_result$plot
ggplot2::ggsave(
  "RoE_heatmap.pdf",
  plot = roe_result$plot,
  width = 7,
  height = 10
)

# Stop instead of silently omitting cells with missing annotations
strict_result <- roe_calculate(
  object = seu,
  group_col = "group",
  celltype_col = "Celltype",
  na_action = "fail",
  p_adjust_method = "bonferroni"
)
}
}
\seealso{
\code{stats::chisq.test()}, \code{stats::fisher.test()},
\code{stats::p.adjust()}
}
