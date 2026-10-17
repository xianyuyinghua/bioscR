% Documentation for R/Trajectory_Monocle3.R; keep in sync with its opening help comments.
\name{Trajectory_Monocle3}
\alias{Trajectory_Monocle3}
\title{
Infer a Monocle3 trajectory using an existing Seurat embedding
}
\description{
Builds and caches a Monocle3 cell_data_set from RNA counts and an existing
Seurat embedding, orders cells from selected roots, tests graph-associated
genes, and exports pseudotime tables and trajectory figures.
}
\usage{
Trajectory_Monocle3(
  seurat,
  reduction = "umap.harmony",
  celltype_col = "Celltype",
  celltype_start = NULL,
  point_start_use = TRUE,
  celltype_colors = basicR::get_colors(number = 50.1),
  output_dir = "./04_Pseudotime",
  output_rds = "./Rdata/04_Monocle3.cds",
  formats = c("pdf", "png")
)
}
\arguments{
\item{seurat}{A Seurat object, required with no default. A fresh analysis uses
SeuratObject::LayerData() with assay = "RNA" and layer = "counts",
cell metadata with matching cell names, and the requested embedding.
Supply a compatible counts layer and enough cells and features for the
fixed 30-dimensional preprocessing. All counts cells must occur in the
metadata and embedding. Ignored when output_rds exists.}
\item{reduction}{Single name of a Seurat reduction. Defaults to "umap.harmony".
It must contain at least two columns; only the first two are used,
reordered to the CDS cell order and stored under the name "UMAP".
No new UMAP is computed. Only used for a fresh analysis.}
\item{celltype_col}{Single name of a character or factor cell-type metadata
column. Defaults to "Celltype". Used for candidate-root annotation and
cell-type plotting, including plotting a cached object. With
point_start_use = FALSE, actual root-cell selection always uses the
literal Celltype column, irrespective of this argument. Use "Celltype"
for that workflow, or ensure both columns have consistent labels.}
\item{celltype_start}{Single starting cell-type label. Defaults to NULL,
which triggers readline() during a fresh analysis. Provide an observed,
non-missing label. It identifies candidate leaf nodes for display; with
point_start_use = FALSE, all cells with this label in the literal
Celltype column are passed as root_cells. No explicit label validation
is performed. Ignored when output_rds exists.}
\item{point_start_use}{Single non-missing logical value. Defaults to TRUE.
For a fresh analysis, TRUE always prompts for comma-separated principal
node IDs, for example "Y_1,Y_2", and passes them as root_pr_nodes.
There is no argument for supplying these IDs programmatically. Input
is split on commas without trimming whitespace or checking candidates.
FALSE selects root cells by celltype_start without a node prompt.
For non-interactive use, set FALSE and supply celltype_start explicitly.
Ignored when output_rds exists.}
\item{celltype_colors}{Vector of valid colours with enough entries for the
cell-type categories in the final labelled trajectory plot. Defaults to
basicR::get_colors(number = 50.1). Named colours can match category
labels. This palette is not applied to the candidate-root figure or
the pseudotime-coloured figures.}
\item{output_dir}{Existing writable directory for CSV tables and figures.
Defaults to "./04_Pseudotime". The function does not create it. Relative
paths use the working directory. On successful completion, files with
names containing "02_Monocle3_RootNodes" are deleted from this directory,
including pre-existing matching files.}
\item{output_rds}{Single path to the CDS cache, independent of output_dir.
Defaults to "./Rdata/04_Monocle3.cds". Despite the extension, this is an
RDS file written with saveRDS(). Its parent directory must exist and be
writable for a fresh analysis. An existing file must contain a readable,
compatible ordered cell_data_set with UMAP coordinates, graph information,
and the requested cell-type metadata; it is reused without validation
against the supplied Seurat object or trajectory parameters.}
\item{formats}{Character vector of formats accepted by basicR::save_figure(),
such as "pdf", "png", "tiff", or "svg". Defaults to c("pdf", "png").
Used for all figures, including the temporary root-selection figure.
The corresponding graphics devices and fonts must be available.}
}
\details{
For a fresh analysis, counts come from the fixed RNA assay and counts
layer. Metadata are reordered to counts columns, and gene_short_name is
set to the counts row names. The wrapper performs no explicit gene or
cell filtering before constructing the CDS. It calls preprocess_cds()
with num_dim = 30 and norm_method = "log", copies the first two embedding
coordinates into UMAP, calls cluster_cells() with reduction_method =
"UMAP" and resolution = 0.1, then learn_graph() with use_partition = TRUE,
close_loop = FALSE, and learn_graph_control = list(minimal_branch_len = 30,
prune_graph = TRUE). There is no explicit random seed in this wrapper.

Root candidates are degree-one principal-graph nodes with at least one
assigned cell of celltype_start, using the closest-vertex mapping and the
column named by celltype_col. Counts and within-node proportions are
computed by node and cell type. No minimum proportion, majority rule,
or highest-proportion selection is applied. The temporary figure
02_Monocle3_RootNodes highlights these candidates. Manually entered nodes
are not restricted to the displayed candidates by this wrapper. With
point_start_use = FALSE, root cells are selected from the literal Celltype
column across the whole CDS, not only from candidate leaf nodes.

Pseudotime is obtained from order_cells() using the selected roots; the
wrapper does not reverse it, select an endpoint, or verify root coverage
of every partition. Missing or non-finite pseudotime values are not
explicitly filtered before CSV export. There are no sample/group
comparisons, treatment contrasts, or fold-change directions. Clustering
and graph learning use Monocle3 partitions; cell-type labels guide root
selection and colouring rather than separate analyses.

The ordered CDS is saved before gene testing. If output_rds already
exists, construction, preprocessing, graph learning, root prompts, and
root selection are skipped. The cached object replaces the supplied
input without checking parameter consistency. Gene testing, CSV export,
final plotting, and root-figure cleanup still run on every call. Use an
unused cache path to change inputs or roots. A failure after caching can
leave a usable cache alongside incomplete outputs.

All CDS genes are split in row order into batches of at most 5000 genes,
not 1000 as stated in an internal comment. Each batch is passed to
graph_test() with neighbor_graph = "principal_graph" and cores = 1.
The batch results are row-bound; q_value values are taken from each
separate call without a new correction across all batches. The exported
subset requires status == "OK", q_value < 0.05, and morans_I > 0.1,
and is sorted by decreasing morans_I. Missing filter values are dropped
by dplyr::filter(). These are graph-associated genes, not a pairwise
differential-expression contrast or a fitted pseudotime regression.

The file 02_Monocle3_Deg_PseudoTime_sig.csv contains a gene column derived
from the combined result row names and all graph_test() result columns,
including status, q_value, and morans_I, for retained genes. The full
unfiltered result is neither saved separately nor added to the cached CDS.
The file 02_Monocle3_pseudotime.csv contains one column named
monocle3_pseudotime with cell IDs as row names. Both CSV files include
row names through write.csv(). The cache retains counts, cell and gene
metadata, preprocessing results, UMAP coordinates, clusters/partitions,
and principal-graph and ordering information in the cell_data_set.
Retrieve it with readRDS(output_rds).

Four final figure stems are written under output_dir:
03_Monocle3_pseudotime_noNode_noTrajectory (pseudotime colours, no graph
or node labels); 03_Monocle3_pseudotime_noNode_Trajectory (pseudotime
colours and graph, no node labels); 03_Monocle3_pseudotime_Node_Trajectory
(pseudotime colours, graph, and root labels only); and
03_Monocle3_pseudotime_Label_Trajectory (cell-type colours and graph,
with principal-node IDs for branch, leaf, and root nodes, not cell-type
text labels or every principal node). All figures, including the root
candidate figure, pass width = 7, height = 6, and res = 600 to
basicR::save_figure(). Existing output tables and figures may be overwritten.

Plotting uses a locally modified plot_cells() helper, not a direct call
to monocle3::plot_cells(). The numeric colour scale uses viridis option
"C"; the helper adds overlapping point layers. It reads Monocle3 internal
slots, so compatibility with the installed Monocle3 version is required.
On successful completion, list.files() with the unanchored pattern
"02_Monocle3_RootNodes" identifies files to remove. This also runs with
a cached CDS and can remove unrelated files containing that substring.
Earlier errors prevent this final cleanup.

The function attaches dplyr, basicR, monocle3, Seurat, ggplot2, and
SingleCellExperiment even when reading a cache. Executed paths also use
SeuratObject, SummarizedExperiment, assertthat, methods, igraph, tibble,
viridis, ggrepel, stats, and utils. Packages must be installed with their
runtime dependencies. The local plotting helper also references Matrix,
reshape2, and ggrastr in gene-expression or rasterisation branches, but
those branches are not used by this wrapper. DESCRIPTION currently omits
monocle3, assertthat, methods, igraph, tibble, viridis, and ggrepel, as well
as the unused-branch references Matrix, reshape2, and ggrastr. No external
script is sourced or launched. Input is the in-memory Seurat object or
an existing RDS cache. A fresh run with default root-selection settings
requires an interactive console and valid root input; invalid or empty
input is not handled by this wrapper.
}
\value{
An invisible logical vector returned by file.remove() if matching root
figure files are found: one element per file, TRUE for successful removal
and FALSE for failure. If no matching files exist, the final conditional
returns invisible NULL. The function does not return the CDS, gene-test
table, pseudotime vector, or plots. Its analysis results are the RDS cache,
two CSV tables, and figures described above.
}
\examples{
\dontrun{
library(Seurat)
library(bioscR)
seu <- readRDS("seurat.rds")
# Requires RNA counts, at least two embedding columns, and Celltype metadata.
# Replace Progenitor with an observed starting cell-type label.
stopifnot("Celltype" \%in\% colnames(seu@meta.data))
stopifnot("Progenitor" \%in\% seu$Celltype)
stopifnot(ncol(Embeddings(seu, "umap.harmony")) >= 2L)
dir.create("monocle3_results", recursive = TRUE, showWarnings = FALSE)
dir.create("Rdata", recursive = TRUE, showWarnings = FALSE)
cache <- "Rdata/Monocle3.cds"
# Use an unused cache path for a fresh trajectory with these roots.
set.seed(123)
Trajectory_Monocle3(
  seurat = seu,
  celltype_col = "Celltype",
  celltype_start = "Progenitor",
  point_start_use = FALSE,
  output_dir = "monocle3_results",
  output_rds = cache
)
cds <- readRDS(cache)
# For manual node selection in an interactive console, use
# point_start_use = TRUE and a new cache path; inspect the root figure
# when prompted and enter valid node IDs separated by commas, without spaces.
}
}
