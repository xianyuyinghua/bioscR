% Documentation for R/or_heatmap.R; keep in sync with its opening help comments.
\name{or_heatmap}
\alias{or_heatmap}
\title{
Plot cell-type enrichment odds ratios across groups
}
\description{
Calculates Fisher exact-test odds ratios for cell-type and group combinations and returns an annotated heatmap with its underlying odds ratios and significance symbols.
}
\usage{
or_heatmap(
  seurat = seurat,
  group = "group",
  celltype = "Celltype",
  row_col = "celltype_group",
  colors = c("#0073C2FF", "white", "#CD534CFF")
)
}
\arguments{
\item{seurat}{A Seurat object containing the group and cell-type metadata columns. Supply this argument explicitly.}
\item{group}{Name of the metadata column defining groups.}
\item{celltype}{Name of the metadata column containing cell-type annotations.}
\item{row_col}{With "celltype_group", cell types form rows and groups form columns. Any other value transposes the displayed matrix and returned data.}
\item{colors}{Three colours for the lower endpoint, log2(OR + 1) = 1, and upper endpoint of the heatmap scale.}
}
\details{
The function builds a cell-type by group count table and replaces zero counts with one before testing. For each combination, a two-sided Fisher exact test compares that combination with the remaining cells. P values are adjusted together using the Benjamini-Hochberg method. Symbols are shown only when OR > 1.5 or OR < 0.5: * for adjusted P < 0.05, ** for adjusted P < 0.01, and *** for adjusted P < 0.001.

Heatmap colours encode log2(OR + 1), with endpoints at its minimum and 99th percentile and a middle colour at 1. Rows and columns are not clustered. Factor levels determine their order; non-factor cell types use first-occurrence order, while non-factor groups are ordered by group_new_levels(). Use at least two observed cell types and groups, without unused factor levels.

The implementation uses data.table, plyr, tidyr, circlize, ComplexHeatmap, grid, gridExtra, and cowplot. Unqualified grid and heatmap drawing functions must be available; group_new_levels() is required when group is not a factor. This is a pooled cell-level comparison and does not model biological replicates.
}
\value{
A list with plot, a cowplot/ggplot heatmap, and data, a list containing or_df (untransformed odds ratios) and p_df (a matrix of significance strings, not numeric P values). The orientation of both data components follows row_col.
}
\examples{
\dontrun{
library(ComplexHeatmap)
library(grid)
seu <- readRDS("seurat.rds")
seu$group <- droplevels(factor(seu$group))
seu$Celltype <- droplevels(factor(seu$Celltype))
result <- or_heatmap(seurat = seu)
result$plot
result$data$or_df
result$data$p_df
}
}
