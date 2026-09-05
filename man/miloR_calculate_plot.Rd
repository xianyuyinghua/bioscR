% Documentation for R/miloR_calculate_plot.R; keep in sync with its opening help comments.
\name{miloR_calculate_plot}
\alias{miloR_calculate_plot}
\title{
Calculate and plot pairwise Milo neighbourhood differential abundance
}
\description{
Performs neighbourhood differential abundance analysis for every pair of group factor levels in a Seurat object and saves neighbourhood graphs, cell-type beeswarm plots with boxplots, and annotated result tables.
}
\usage{
miloR_calculate_plot(
  seurat,
  cell_type_col = "Celltype",
  group_col = "group",
  reduction = "umap.harmony",
  celltype_colors = basicR::get_colors(number = 50.1),
  continuous_colors = c("#070091", "lightgrey", "#910000"),
  formats = c("pdf", "png"),
  output_dir = "03_Celltype_Marker"
)
}
\arguments{
\item{seurat}{A Seurat object with expression data suitable for conversion by Seurat::as.SingleCellExperiment(), metadata columns for groups and cell types, and a factor column named sample identifying biological samples. Each sample must belong to one group; supply biological replication for the differential abundance model.}

\item{cell_type_col}{Name of a factor metadata column containing cell-type annotations. Defaults to "Celltype". The current plotting code also accesses the annotated result column Celltype explicitly, so use "Celltype" for the supported workflow.}

\item{group_col}{Name of a factor metadata column defining comparison groups, with at least two observed levels. Remove unused levels before calling. The name must be usable directly in an R formula because the model is constructed as ~ group_col without quoting the column name.}

\item{reduction}{Name of the embedding used for neighbourhood graph visualisation. Its uppercase form is passed as the layout name in the converted SingleCellExperiment object. Defaults to "umap.harmony". This argument does not select the dimensional reduction used to construct the nearest-neighbour graph.}

\item{celltype_colors}{Vector of colours for the cell-type boxplots, with enough colours for the annotated cell types. Defaults to basicR::get_colors(number = 50.1).}

\item{continuous_colors}{Three colours for negative, central, and positive log fold changes, respectively.}

\item{formats}{Character vector of figure formats passed to basicR::save_figure(). Defaults to PDF and PNG.}

\item{output_dir}{Writable directory for figures and CSV tables. Defaults to "03_Celltype_Marker". Create the directory before calling; the function contains no explicit directory-creation step.}
}
\details{
Comparisons follow the order returned by combn() on the levels of the group factor and run sequentially. Each pair is subsetted, unused group, sample, and cell-type levels are dropped, and the subset is converted to a SingleCellExperiment and then a Milo object.

The nearest-neighbour graph uses k = 30 and d = 30. Neighbourhood sampling uses prop = 0.2, k = 30, d = 30, and refined = TRUE. Neighbourhood distances also use d = 30. These calls use the default reduced-dimension selection of miloR; provide input suitable for its 30-dimensional workflow. The reduction argument is used only for graph visualisation. No random seed is set internally; set one before calling if needed.

Cells are counted by the fixed metadata column sample. The design contains one row per sample and uses only the grouping variable, with no additional covariates or explicit contrast supplied to testNhoods(). Under standard treatment contrasts, the second retained group level is compared with the first. Neighbourhoods are annotated using cell_type_col, and the plotted Celltype order is reversed relative to the input factor levels.

Both plotNhoodGraphDA() and plotDAbeeswarm() receive alpha = 1. The plots therefore do not apply a conventional 0.05 significance threshold. The graph fill scale is limited to [-5, 5], with outlying colours squished to the endpoints; the beeswarm plot includes reference lines at log fold changes -1 and 1.

For each group pair, the function saves a neighbourhood graph with prefix 07_Celltype_Cell_OR_milo_by_group_ and a beeswarm/boxplot figure and CSV table with prefix 07_Celltype_Cell_Beeswarm_and_Boxplots_by_group_. Each prefix is followed by the two group levels joined by an underscore. Despite the OR text in the graph filename, the displayed values are log fold changes, not odds ratios. Graph figures are 6.5 by 5.5 inches and beeswarm figures are 6.5 by 6 inches, all at 600 dpi. CSV files contain the annotated testNhoods() results. There is no result cache or explicit save of the Milo object.

The function attaches dplyr, Seurat, miloR, scran, scater, SeuratWrappers, SingleCellExperiment, patchwork, ggbeeswarm, scales, forcats, data.table, stringr, and basicR. It also calls ggnewscale::new_scale_color() and uses unqualified ggplot2 plotting functions; these must be available in the evaluation environment.
}
\value{
A list with one element per group pair, containing the return values of the final basicR::save_figure() call for the beeswarm/boxplot figure. The primary outputs are the saved figures and CSV tables; the function does not directly return the Milo objects or differential abundance tables.
}
\examples{
\dontrun{
library(Seurat)
library(ggplot2)
seu <- readRDS("seurat.rds")
seu$group <- droplevels(factor(seu$group))
seu$sample <- droplevels(factor(seu$sample))
seu$Celltype <- droplevels(factor(seu$Celltype))
dir.create("milo_output", recursive = TRUE, showWarnings = FALSE)

# Supply an object suitable for the 30-dimensional graph workflow
# and containing the requested embedding.
set.seed(123)
miloR_calculate_plot(
  seurat = seu,
  cell_type_col = "Celltype",
  group_col = "group",
  reduction = "umap.harmony",
  output_dir = "milo_output"
)
}
}
