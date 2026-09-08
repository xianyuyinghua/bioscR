% Documentation for R/inferCNV_analysis.R; keep in sync with its opening help comments.
\name{inferCNV_analysis}
\alias{inferCNV_analysis}
\title{
Run inferCNV analysis and plot cell-type CNV profiles
}
\description{
Runs inferCNV on raw counts from a Seurat object, caches the analysis, and saves a chromosome heatmap and cell-type CNV score boxplot.
}
\usage{
inferCNV_analysis(seurat,
                              celltype_col = "Celltype",
                              sample_col = "sample",
                              species = "human",
                              ref_name = NULL,
                              assay = "RNA",
                              heatmap_downsample = NULL,
                              heatmap_cnv_probs = c(0.025,0.975),
                              celltype_col_colors = basicR::get_colors(number = 50.1),
                              chr_colors = c("#E64B35", "#4DBBD5"),
                              ncor = 10,
                              output_dir = ".",
                              formats = c("pdf","png")
                             )
}
\arguments{
\item{seurat}{A Seurat object with gene symbols as feature names, a counts layer in the selected assay, and the metadata columns specified below.}
\item{celltype_col}{Name of a factor metadata column containing cell-type annotations. Remove unused factor levels before calling. Defaults to "Celltype".}
\item{sample_col}{Name of the metadata column identifying samples. Used together with celltype_col for heatmap downsampling. Defaults to "sample".}
\item{species}{Species passed to AnnoProbe::annoGene(), for example "human" or "mouse". Defaults to "human". Chromosome ordering and annotation are fixed to chr1 through chr22 in the current implementation.}
\item{ref_name}{A single reference cell-type name present in celltype_col. Required despite the NULL default. Use a name without hyphens, because metadata labels are converted from hyphens to underscores internally, while ref_name is not converted.}
\item{assay}{Assay containing the counts layer. Defaults to "RNA".}
\item{heatmap_downsample}{Optional positive integer giving the maximum number of cells per sample and cell-type combination in the heatmap. NULL uses all cells. Analysis and CNV scoring always use all retained cells.}
\item{heatmap_cnv_probs}{Two quantile probabilities used to set the heatmap colour range symmetrically around 1. Defaults to c(0.025, 0.975).}
\item{celltype_col_colors}{Vector of plotting colours, with enough colours for all cell-type factor levels. Defaults to basicR::get_colors(number = 50.1).}
\item{chr_colors}{Colours repeated across the 22 chromosome annotation levels. Defaults to c("#E64B35", "#4DBBD5").}
\item{ncor}{Number of threads passed to infercnv::run(). Defaults to 10.}
\item{output_dir}{Existing writable directory for input tables, RDS caches, CSV tables, and figures. Defaults to the working directory. The function does not create this directory. The inferCNV run directory is separately fixed to infercnv_output in the working directory.}
\item{formats}{Character vector of figure formats passed to basicR::save_figure(). Defaults to c("pdf", "png").}
}
\details{
Gene positions are obtained with AnnoProbe::annoGene(). Mitochondrial and sex chromosomes (chrM, chrX, and chrY) are excluded; duplicated gene symbols are removed and genes are ordered by chromosome and start position. Hyphens in cell names, sample labels, and cell-type labels are replaced by underscores.

The inferCNV run uses cutoff = 0.1, cluster_by_groups = FALSE, hclust_method = "ward.D2", analysis_mode = "subclusters", denoise = TRUE, HMM = FALSE, and leiden_resolution = "auto". Intermediate and final inferCNV plotting are disabled. The function sets options(scipen = 100) when running inferCNV and does not restore the previous value.

Input files are expFile.txt, geneFile.txt, groupFiles.txt, and cell_meta.rds. The initial and analysed inferCNV objects are saved as infercnv_obj.rds and infercnv_obj_analysis.rds. Existing cell_meta.rds and infercnv_obj_analysis.rds files are reused without checking the input or parameters. For a fresh analysis, use a new output directory and a working directory without a previous infercnv_output run. A partial cache containing cell_meta.rds but no analysed object cannot be resumed reliably because input file path variables are only assigned during input preparation.

The heatmap matrix is saved to 01_inferCNV_Heatmap.csv with cells in rows and genes in columns. Heatmap figures use the same filename stem. Rows are split by cell type and columns by chromosome, without clustering. Downsampling sets the random seed to 123 and affects only this heatmap and its CSV table.

CNV score thresholds are derived from the reference expression matrix: the lower and upper bounds are the mean of gene-wise means minus or plus twice the mean of gene-wise standard deviations. Their difference defines one copy interval. Each gene contributes 0 within the central interval, 1 in either adjacent interval, and 2 beyond them. Scores are summed per cell. Metadata and scores are saved to 02_Celltype_TotalCNV.csv; the boxplot uses the same filename stem and includes a dashed line at the reference median score.

The function loads dplyr, ggplot2, Seurat, infercnv, AnnoProbe, ComplexHeatmap, circlize, and grid, and also uses tibble and basicR. These packages must be installed. If ref_name is NULL or absent from the original cell-type metadata, the current implementation calls q("no"), which terminates the R session.
}
\value{
The return value of the final basicR::save_figure() call. Analysis objects, plotting data, and figures are saved to disk; the function does not return a modified Seurat object.
}
\examples{
\dontrun{
seu <- readRDS("seurat.rds")
seu$Celltype <- droplevels(factor(seu$Celltype))
# Replace T_cells with a reference level present in seu$Celltype.
dir.create("cnv_results", recursive = TRUE, showWarnings = FALSE)
inferCNV_analysis(
  seurat = seu,
  ref_name = "T_cells",
  heatmap_downsample = 200,
  output_dir = "cnv_results"
)
}
}
