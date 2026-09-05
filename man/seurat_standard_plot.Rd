% Documentation for R/seurat_standard_plot.R; keep in sync with its opening help comments.
\name{seurat_standard_plot}
\alias{seurat_standard_plot}
\title{
Save standard Seurat embedding and cell-composition plots
}
\description{
Creates embedding plots for combinations of grouping and splitting metadata columns, together with cell-count and cell-proportion plots for each split.
}
\usage{
seurat_standard_plot(seurat = seurat,
                                 reduction = "umap.harmony",
                                 group_by = c("Cluster","Celltype"),
                                 split_by = c("sample","group"),
                                 colors = basicR::get_colors(100.1),
                                 formats = c("pdf","png"),
                                 output_dir = file.path("Output","03_Celltype_Marker")
                                )
}
\arguments{
\item{seurat}{A Seurat object with the requested reduction and metadata columns. Supply this argument explicitly. Grouping columns should be factors with the desired cell-type order.}
\item{reduction}{Name of an existing dimensional reduction, by default "umap.harmony".}
\item{group_by}{Character vector of metadata column prefixes used for colouring cells. Prefixes are combined into a regular expression anchored at the start of the column name.}
\item{split_by}{Character vector of metadata column prefixes used for splitting plots, matched in the same way as group_by. An unsplit embedding is also generated.}
\item{colors}{Vector of colours passed to the embedding and composition plots. Defaults to basicR::get_colors(100.1).}
\item{formats}{Character vector of output formats passed to basicR::save_figure(), by default PDF and PNG.}
\item{output_dir}{Directory passed to basicR::save_figure() for saving plots. Defaults to Output/03_Celltype_Marker. Create a writable directory before use.}
}
\details{
For every matched grouping column, the function creates split and unsplit Seurat::DimPlot() figures, each with either a legend or cell-group labels. Point size is calculated from the total cell count and clipped to [0.01, 8]. Facet layout and figure dimensions are selected automatically.

For each splitting column, cell-count and within-group cell-proportion figures include alluvial flows. Composition file names are shared between the legend and label iterations, so these files are saved again during the loop. Embedding file names include a grouping index, grouping column, reduction (dots replaced by underscores), Legend or Label, and splitting column (NULL for unsplit plots). Composition file names contain Cell_Number or Cell_Proportion. Figures are saved at 600 dpi.

The function attaches dplyr, Seurat, ggplot2, and basicR, and additionally uses stringr, tibble, tidyr, ggalluvial, and scales. It uses group_new_levels() for non-factor splitting columns.
}
\value{
A list with one element per grouping, splitting, and legend combination. Elements are the return values of the final composition save for split combinations, and NULL for unsplit combinations. The primary output is the set of saved figure files.
}
\examples{
\dontrun{
seu <- readRDS("seurat.rds")
seu$Celltype <- factor(seu$Celltype)
dir.create("standard_plots", recursive = TRUE, showWarnings = FALSE)
seurat_standard_plot(
  seurat = seu, reduction = "umap",
  group_by = "Celltype", split_by = c("sample", "group"),
  output_dir = "standard_plots"
)
}
}
