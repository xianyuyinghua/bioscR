% Documentation for R/augur_calculate_plot.R; keep in sync with its opening help comments.
\name{augur_calculate_plot}
\alias{augur_calculate_plot}
\title{
Calculate and plot pairwise Augur cell-type prioritisation
}
\description{
Runs Augur for every pair of group factor levels in a Seurat object, caches the results, and saves AUC tables and cell-type lollipop plots.
}
\usage{
augur_calculate_plot(seurat,
                                 cell_type_col = "Celltype",
                                 group_col = "group",
                                 colors = basicR::get_colors(number = 50.1),
                                 formats = c("pdf", "png"),
                                 output_dir = "03_Celltype_Marker"
                                )
}
\arguments{
\item{seurat}{A Seurat object with expression data and metadata columns specified by cell_type_col and group_col. A factor metadata column named sample is also required.}
\item{cell_type_col}{Name of the metadata column containing cell-type annotations.}
\item{group_col}{Name of a factor metadata column defining the groups. Supply at least two observed levels and remove unused levels before calling.}
\item{colors}{Vector of plotting colours, with enough colours for the cell types. Defaults to basicR::get_colors(number = 50.1).}
\item{formats}{Character vector of output formats passed to basicR::save_figure(), by default PDF and PNG.}
\item{output_dir}{Existing writable directory for RDS caches, CSV tables, and figures. The function does not create this directory before saving RDS files.}
}
\details{
Each pair is analysed with Augur::calculate_auc(), using min_cells = 10, subsample_size = 10, and n_threads = 1. Comparisons run through parallel::mclapply() with one core per pair; multiple pairs require a platform supporting forked parallelism. No random seed is set internally.

Files use the prefix 06_Augur_AUC_Celltype_by_ followed by the two group levels joined with an underscore. Existing RDS files are reused without checking whether the input data have changed. CSV tables and figures are regenerated from these caches. AUC tables are sorted in ascending order; plots include a reference line at AUC 0.5 and are saved at 6 by 6 inches and 600 dpi.

The implementation uses Augur, Seurat, dplyr, ggplot2, and basicR. Unqualified helpers, including sym(), the pipe operator, and plotting functions, must be available in the evaluation environment.
}
\value{
A list with one element per group pair, containing the return values of basicR::save_figure(). The Augur results themselves are saved as RDS files; sorted AUC tables are saved as CSV files.
}
\examples{
\dontrun{
library(Seurat)
library(dplyr)
library(ggplot2)
seu <- readRDS("seurat.rds")
seu$group <- droplevels(factor(seu$group))
seu$sample <- factor(seu$sample)
dir.create("augur_output", recursive = TRUE, showWarnings = FALSE)
augur_calculate_plot(seurat = seu, output_dir = "augur_output")
}
}
