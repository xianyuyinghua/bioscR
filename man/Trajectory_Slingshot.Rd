% Documentation for R/Trajectory_Slingshot.R; keep in sync with its opening help comments.
\name{Trajectory_Slingshot}
\alias{Trajectory_Slingshot}
\title{
Infer Slingshot lineages and plot dynamic gene modules
}
\description{
Fits Slingshot trajectories on a two-dimensional Seurat embedding, tests
gene dynamics with tradeSeq, clusters branch expression profiles, and
exports trajectory plots and a GO BP/KEGG-annotated module heatmap.
}
\usage{
Trajectory_Slingshot(
  seurat,
  reduction = "umap.harmony",
  celltype_col = "Celltype",
  celltype_start = NULL,
  celltype_end = NULL,
  celltype_colors = basicR::get_colors(number = 50.1),
  species = "human",
  npoints = 100,
  num_clusters = 6,
  min_pct = 0.1,
  deg_pname = "p.adj",
  deg_pvalue = 0.05,
  topn_genes = 1000,
  topn_show = 15,
  enrich_pname = "p.value",
  enrich_pvalue = 0.05,
  enrich_topn = 10,
  expr_source = "fitted",
  ht_colors = c("#19475E", "#2E8B8B", "#F2E33B"),
  ncors = 4,
  output_dir = "./04_Pseudotime",
  formats = c("pdf", "png")
)
}
\arguments{
\item{seurat}{A Seurat object, required with no default, convertible using
as.SingleCellExperiment(assay = "RNA"). Supply accessible RNA counts and
normalised data layers with matching genes and cells, the requested
embedding, and cell-type metadata. Both layers are read even with
expr_source = "fitted". Gene identifiers must be compatible with the
installed ClusterGVis enrichment workflow and its annotation database;
this wrapper performs no identifier conversion itself. Retained cells
need positive library sizes for the logarithmic offset.}
\item{reduction}{Single Seurat reduction name. Defaults to "umap.harmony".
The embedding must have exactly two columns, all finite coordinates,
and the same set of cell names as the converted object. Rows are
reordered to the converted cell order, and coordinates are stored under
toupper(reduction). No dimensional reduction is computed.}
\item{celltype_col}{Single name of a metadata column containing non-missing
cluster/cell-type labels. Defaults to "Celltype". Converted to a factor
with unused levels dropped, and used as Slingshot clusterLabels and for
plotting. Labels define trajectory clusters, not treatment comparisons.}
\item{celltype_start}{Starting cluster label passed as start.clus to Slingshot.
Defaults to NULL, leaving start selection to Slingshot. Any supplied
label must occur among the retained cell-type factor levels. There is
no interactive prompt.}
\item{celltype_end}{Character vector of terminal cluster labels passed as
end.clus to Slingshot. Defaults to NULL, leaving endpoint constraints
unspecified. All supplied labels must be observed factor levels.}
\item{celltype_colors}{Vector of valid cell-type colours, defaulting to
basicR::get_colors(number = 50.1). An unnamed vector must contain at
least as many colours as observed cell-type factor levels; the first
colours are assigned in factor-level order. A named vector must cover
all observed cell-type labels and is reordered by those labels. Extra
colours are ignored. Too few unnamed colours or missing required labels
in a named vector cause a stopifnot() error before trajectory fitting.
Use unique names matching cell-type labels to make the mapping explicit.}
\item{species}{Single string, "human" or "mouse"; defaults to "human".
Selects organism = "hsa" with OrgDb = org.Hs.eg.db for human, or
organism = "mmu" with OrgDb = org.Mm.eg.db for mouse, in both GO BP and
KEGG enrichment calls. Supply gene identifiers appropriate for the
selected species. Other values are not validated and leave the local
organism and orgdb variables undefined.}
\item{npoints}{Integer number of heatmap columns per lineage, at least 2;
defaults to 100. Used as nPoints for fitted curves or as the number of
equal-width observed pseudotime bins. Lineages need non-empty, finite
pseudotime ranges; observed-mode interpolation requires at least two
non-missing bins when filling gaps. No range validation is implemented.}
\item{num_clusters}{Positive integer number of k-means gene modules. Defaults
to 6. Must not exceed the number of displayed genes with usable,
distinct expression profiles or the available module palette. The
explicit check only compares the total significant gene count with
num_clusters, before truncation by topn_genes.}
\item{min_pct}{Detection fraction in [0, 1], defaulting to 0.1. A gene must
have counts >= 1 in at least ceiling(min_pct * retained_cell_count)
cells. The default is 10 percent. The count threshold is fixed; no
min_count argument exists. The range is not explicitly validated.}
\item{deg_pname}{Single choice "p.adj" (default) or "p.value". Selects the
global association-test padj or pvalue column for gene filtering and
ranking. Choices are validated only after model fitting and export
of the full test table; lineage-specific P values do not select genes.}
\item{deg_pvalue}{Numeric significance threshold in (0, 1], defaulting to
0.05. Only converged genes with finite selected P values strictly
below this threshold are retained. The range is not validated.}
\item{topn_genes}{Positive integer maximum number of significant genes used
for the heatmap, modules, and enrichment. Defaults to 1000. Genes are
ranked by the selected global P value, then decreasing waldStat.
Supply at least num_clusters; fewer genes are used if fewer pass the
filter. Also controls heatmap height, even if fewer genes are displayed.}
\item{topn_show}{Non-negative integer maximum number of gene labels per
module. Defaults to 15. Selects the first rows per cluster in
st_data$wide.res as returned by clusterData(); there is no subsequent
P-value sort within each module. No range validation is implemented.}
\item{enrich_pname}{Single choice "p.value" (default) or "p.adj". Selects
pvalue or p.adjust for post-enrichment filtering and ranking. If the
selected column is absent, the code silently falls back to pvalue.
Validation occurs after enrichment; "p.adj" also emits a warning
about the earlier enrichCluster() cutoff.}
\item{enrich_pvalue}{Numeric threshold in (0, 1], defaulting to 0.05. Passed
as pvalueCutoff to both enrichCluster() calls, then used as a strict
upper bound on the finite values selected by enrich_pname. Earlier
filtering is controlled by the installed enrichment implementation.
The wrapper does not validate the range or recalculate enrichment P values.}
\item{enrich_topn}{Positive integer maximum number of GO BP terms and KEGG
pathways shown per module, separately for each database. Defaults to
10. Display terms are ordered by the selected enrichment P value.
No range validation is implemented.}
\item{expr_source}{Single string specifying branch expression values.
Defaults to "fitted", which uses tradeSeq::predictSmooth(). Use
"observed" for unweighted bin means from the RNA data layer. Every
value other than "fitted" enters the observed branch; spelling is
not validated. Both modes still fit tradeSeq and select genes using
its association tests.}
\item{ht_colors}{Three valid colours corresponding to the fixed heatmap
Z-score breaks -2, 0, and 2. Defaults to
c("#19475E", "#2E8B8B", "#F2E33B"). Passed in htColList to
ClusterGVis::visCluster(); no length validation is implemented.}
\item{ncors}{Positive integer number of workers for
BiocParallel::MulticoreParam(), defaulting to 4. fitGAM() always receives
parallel = TRUE and RNGseed = 123. Use an environment supporting this
backend, normally Unix-like systems; no alternative backend is exposed.}
\item{output_dir}{Writable output directory, defaulting to "./04_Pseudotime".
Created recursively if necessary. Relative paths use the working
directory. Existing files may be overwritten; no saved result is reused
as a cache, and files from previous runs are not removed.}
\item{formats}{Character vector of formats accepted by basicR::save_figure(),
such as "pdf", "png", "tiff", and "svg". Defaults to c("pdf", "png").
Used for all figures. Suitable graphics devices, fonts, and support for
the heatmap object returned by visCluster() are required.}
}
\details{
The RNA assay is converted to SingleCellExperiment, embedding rows are
matched to cells, and cell-type labels are converted to a factor. Slingshot
receives start.clus and end.clus as supplied, omega = FALSE, times = NULL,
shrink = TRUE, extend = "y", reweight = TRUE, reassign = TRUE,
thresh = 0.001, maxit = 50, and stretch = 1. All inferred lineages are
retained. Lineage indices follow the returned ordering, not predefined
biological endpoints. Consult the saved lineage paths for direction.
No sample or treatment comparisons, fold changes, or endpoint contrasts
are performed. Smooth curves are plotted in their stored curve order;
the separate node plot connects cell-type coordinate means, not medians.

Per-lineage pseudotime is exported with na = TRUE, preserving missing
values outside each lineage, together with curve weights. For gene
modelling, only cells with a positive sum of lineage weights are retained.
No explicit assertion checks agreement of pseudotime NA positions and
zero weights, despite an internal comment. Detection filtering uses these
retained cells. Complete pseudotime from na = FALSE and curve weights are
passed to fitGAM(). Counts for all genes are retained in the input matrix,
while genes = which(keep_genes) limits the genes fitted. The offset is
log(colSums(counts_heat) * edgeR::calcNormFactors(counts_heat)), computed
before gene filtering. fitGAM() uses nknots = 6.

associationTest() runs with global = TRUE, lineages = TRUE, and
inverse = "eigen". Global and available lineage-specific P values are
recomputed with the chi-square upper tail where the Wald statistic and
degrees of freedom are finite and degrees of freedom are positive.
BH adjustment is applied only to finite global P values of converged
genes. The full test table is sorted by padj and decreasing waldStat,
regardless of deg_pname. Significant genes are separately selected and
ranked by deg_pname; only the first topn_genes enter subsequent analysis.
There is no separate export of the complete significant-gene subset.

In fitted mode, predictSmooth(tidy = FALSE) supplies the branch matrix.
Columns are assumed to occur in lineage blocks of npoints columns. The
pseudotime grid for each branch is reconstructed from the t/j and l/j
columns (named t1, l1, and so on) in colData(gam_sce)$tradeSeq$dm, using
cells with lineage assignment equal to 1. This depends on tradeSeq's
internal design-matrix structure. In observed mode, cells with non-missing
pseudotime in each lineage are put into equal-width bins. The matrix
contains unweighted means of normalised RNA data, with empty bins filled
by within-lineage linear interpolation and constant endpoint extension
via approx(rule = 2). No values are borrowed from other lineages. The
observed-mode exported grid is an evenly spaced sequence including the
range endpoints, not the bin midpoints. Degenerate ranges or too few
occupied bins are not handled explicitly.

clusterData() receives the unscaled branch matrix with scaleData = TRUE,
clusterMethod = "kmeans", clusterNum = num_clusters, and
kmeansParamsList = list(nstart = 50, iter.max = 100). The returned object
is labelled type = "monocle" and geneType = "branched" for presentation;
no Monocle analysis is performed. Heatmap columns are split by lineage,
with an additional actual-pseudotime colour annotation on a shared range.
Gene labels follow the returned cluster table order rather than an
explicit within-module significance ranking.

GO BP and KEGG enrichment are performed separately for each module
formed from the selected topn_genes genes (at most 1000 by default).
Both calls use the species-matched OrgDb and organism, topn = NULL,
add.gene = TRUE, and seed = 123. The detected-gene vector named
background_genes is only used in a message; it is not passed as an
enrichment universe and does not expand the selected module gene lists.
Background selection and identifier handling therefore follow enrichCluster()
defaults. The wrapper passes enrich_pvalue as pvalueCutoff before its
own selected-column filter. It warns that "p.adj" is a later filter,
but does not inspect or override upstream filtering. If p.adjust is
missing, raw pvalue is used even in "p.adj" mode. KEGG errors are caught
for up to three attempts, with a five-second sleep after each failed
attempt, including the last. A NULL result emits a warning and allows
execution to continue; GO BP errors are not caught.

Output tables and RDS files under output_dir are as follows. The files
02_Slingshot_pseudotime.csv and 02_Slingshot_weights.csv have cells in
rows and lineages in columns, with row names included; pseudotime columns
are renamed slingPseudotime_1, slingPseudotime_2, and so on.
02_Slingshot_lineages.csv has lineage indices and cluster paths joined
by arrows. 02_Slingshot_sds.rds stores the object returned by slingshot()
on the SingleCellExperiment input, including embedding and trajectory
information; it is not the separate SlingshotDataSet conversion.
04_Slingshot_tradeSeq_models.rds stores the fitted model object.
05_Slingshot_dynamic_gene_test.csv contains all association-test columns,
gene identifiers, convergence flags, and global padj values.

05_Slingshot_branch_expression_matrix.csv contains displayed genes in
rows and lineage/grid-point columns, before clusterData() scaling, with
gene row names included. 05_Slingshot_branch_pseudotime_grid.csv maps
column names to pseudotime values. 06_Slingshot_gene_cluster.csv is the
returned st_data$wide.res table, including gene/module assignments and
expression-profile columns. The optional 07_Slingshot_enrich_BP_all.csv
and 07_Slingshot_enrich_KEGG_all.csv retain enrichment result columns
after filtering; "all" does not mean unfiltered tests. Corresponding
files ending in _display.csv contain group, Description, pvalue, and
ratio for up to enrich_topn terms per module. Here pvalue contains the
selected statistic (possibly adjusted P), and ratio is GeneRatio
converted from a numerator/denominator string to a percentage, or NA if
that column is absent. These tables, the lineage-path table, the test
table, the grid table, and the cluster table omit CSV row names. NULL
enrichment results are not written, so stale files may remain from
previous runs.

Figure stems are 03_Slingshot_celltype (all cells and smooth lineages,
7 by 6 inches), 03_Slingshot_pseudotime_nolegend_LineageN (one lineage
at a time, excluding its missing pseudotimes, 6 by 6 inches),
03_Slingshot_pseudotime_lineage (combined panels, 4 times the panel-column
count by 5 times the panel-row count), and 03_Slingshot_node (cell-type
means and cluster paths, 8 by 6 inches), each at res = 600. The combined
layout uses bioBulkR::facet_nrow_ncol(). 08_Slingshot_Branch_Heatmap
combines expression, module trends, and enrichment annotations at
res = 300. Its height is 5 + 0.01 * topn_genes. If the longest displayed
term exceeds 130 characters, width is 10 + 2 * number_of_lineages;
otherwise it is 4 + 3 * number_of_lineages + 0.1 * longest_term_length.
The returned visCluster() object is passed directly to save_figure(),
not wrapped in a function as suggested by an internal comment.

If either enrichment display table is empty, max(nchar(...)) can warn;
if both are empty, the computed heatmap width is non-finite and export
can fail. No fallback handles this case. Empty retained cell/gene sets,
constant gene profiles, too few displayed genes for k-means, and invalid
parameter ranges likewise have no comprehensive checks. Lineage colours
use the ggsci JCO palette through basicR and are limited by that palette;
module colours use the first num_clusters entries of the fixed basicR
palette. The wrapper resets the random seed to 123 before Slingshot,
fitGAM(), and clustering and passes seed 123 to enrichment. It changes
repr.plot.width, repr.plot.height, and repr.plot.res options without
restoring them. Intermediate outputs may remain if later stages fail.

The function attaches dplyr, basicR, Seurat, ggplot2, slingshot, Matrix,
tradeSeq, ClusterGVis, clusterProfiler, org.Hs.eg.db, org.Mm.eg.db,
BiocParallel, and SingleCellExperiment. Both annotation packages are
required even for human analysis. It also calls edgeR,
SummarizedExperiment, viridis, ggnewscale, bioBulkR, patchwork,
ComplexHeatmap, circlize, grid, stats, utils, and
bioscR::calculate_cell_point_size(); basicR uses ggsci for lineage
colours. Compatible package versions must support the accessed layers,
tradeSeq design matrix, and ClusterGVis plotting interface. KEGG access
may require network connectivity. No external scripts or input files
are read by this wrapper, and no saved results are loaded on rerun.
DESCRIPTION currently omits slingshot, Matrix, tradeSeq, ClusterGVis,
clusterProfiler, org.Hs.eg.db, org.Mm.eg.db, BiocParallel, edgeR,
viridis, bioBulkR, ComplexHeatmap, circlize, and grid; ggsci is also
undeclared there as an indirect palette dependency.
}
\value{
The return value of the final save_figure() call, normally resolved from
basicR; its type and visibility depend on that implementation. No
structured result list, modified Seurat object, fitted model, or heatmap
is explicitly returned. Retrieve trajectory and fitted-model objects
from the two RDS files and tabular results from the CSV files described
above. The module object and combined plot are not saved as R objects.
}
\examples{
\dontrun{
library(Seurat)
library(bioscR)
# Requires the listed analysis packages, a compatible parallel environment,
# and working enrichment services. This example uses human data; for mouse
# data, set species = "mouse" to use the corresponding annotation database.
seu <- readRDS("human_seurat.rds")
# Supply RNA counts and normalised data, plus a finite two-column embedding.
stopifnot("Celltype" \%in\% colnames(seu@meta.data))
stopifnot(!anyNA(seu$Celltype))
emb <- Embeddings(seu, "umap.harmony")
stopifnot(ncol(emb) == 2L, all(is.finite(emb)))
# Replace Progenitor with a starting label present in your data.
stopifnot("Progenitor" \%in\% seu$Celltype)
Trajectory_Slingshot(
  seurat = seu,
  celltype_start = "Progenitor",
  species = "human",
  expr_source = "fitted",
  output_dir = "slingshot_results"
)
sds <- readRDS("slingshot_results/02_Slingshot_sds.rds")
models <- readRDS("slingshot_results/04_Slingshot_tradeSeq_models.rds")
# Outputs depend on sufficient significant genes and non-empty enrichment.
}
}
