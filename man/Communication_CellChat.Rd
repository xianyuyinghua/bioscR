% Documentation for R/Communication_CellChat.R; keep in sync with its opening help comments.
\name{Communication_CellChat}
\alias{Communication_CellChat}
\title{
Infer and compare CellChat communication networks from Seurat data
}
\description{
Runs group-wise CellChat inference, optionally merges groups for comparison,
and exports cached objects, network tables, and communication figures.
}
\usage{
Communication_CellChat(seurat,
  celltype = "Celltype",
  group_col = "group",
  group_levels = NULL,
  analysis_mode = "all",
  smooth_use = TRUE,
  key_cell = NULL,
  source_cell = NULL,
  target_cell = NULL,
  hierarchy_top = NULL,
  plot_top_num = c(0.1, 0.25, 0.5, 1),
  plot_count_weight = c("count", "weight"),
  plot_source_target = c("source", "target"),
  plot_allcirc_allheat_allbub_singcirc_singbub = c(TRUE, TRUE, TRUE, TRUE, TRUE),
  plot_pathway_hier_circ_heat_chord_bar_violin_singRoleheat = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
  plot_allsignalingRole_scatt_outheat_inheat = c(TRUE, TRUE, TRUE),
  plot_allcommunication_outgo_income = c(TRUE, TRUE),
  plot_Compute_Aggregate_Analysis_Pvalue = 0.05,
  plot_Net_netP_Pvalue = 0.05,
  plot_prob_Quantile = 0,
  plot_geneExp_minpct = 0.1,
  plot_geneExp_Pvalue = 0.05,
  plot_geneExp_logFC = 0.25,
  ncore = 10,
  formats = c("pdf","png"),
  output_rds_dir = "./Rdata",
  output_plot_dir = "./05_Regulation"
)
}
\arguments{
\item{seurat}{A Seurat object, required with no default. Must contain
normalised expression accessible as seurat@assays$RNA$data and the metadata
columns specified below. The wrapper does not normalise counts or join
layers. An RNA Assay5 without a data slot is not supported by this access.
Use human or mouse gene symbols, including ACTB or Actb respectively, and
matching expression column names and metadata row names.}
\item{celltype}{Single cell-type metadata column name. Defaults to "Celltype".
Supply a factor with observed, non-missing labels and no unused levels;
some plots inspect factor levels. Missing annotations are removed before
createCellChat(). One hierarchy-plot sizing expression explicitly reads
cellchat@meta$Celltype even when a different column is selected.}
\item{group_col}{Single grouping metadata column name, defaulting to "group".
NULL searches column names for the case-insensitive substring "group" and
requires exactly one match. An absent specified column or an ambiguous
search terminates the R session via quit(). Supply non-missing group labels.}
\item{group_levels}{Character vector specifying all observed groups in order,
or NULL (default). With NULL, the first observed label equal to "control"
or "normal", ignoring case, is moved first and other labels keep their
observed order; otherwise labels are sorted. Supply each observed group
exactly once, with no extras. This is not a group-subsetting argument;
subset seurat beforehand. Use exactly two groups for "compare" or "all".
The second group is compared with the first. In "single" mode, all listed
groups are analysed, but filtering with multiple groups uses only the
first two groups to determine shared cell types.}
\item{analysis_mode}{Single, case-sensitive string: "all" (default), "single",
or "compare". "single" infers separate objects per group; "compare" merges
group-wise objects; "all" performs both inference workflows independently.
Although membership is checked, vectors are not supported by later if()
conditions. "all" plots only the first two single-group objects and assumes
exactly two groups when indexing its comparison output directory.}
\item{smooth_use}{Single non-missing logical, default TRUE. Calls smoothData()
with the species-matched PPI network before computeCommunProb(). The latter
receives population.size = FALSE but no explicit raw.use argument; whether
smoothed data are used therefore depends on the installed CellChat default.}
\item{key_cell}{Non-empty character vector of cell-type labels to highlight
in individual circle and bubble plots. Defaults to NULL, but NULL or an
empty vector terminates the R session, even if these plots are disabled.
Labels must occur in the input and remain available after shared-cell-type
filtering. This parameter does not restrict the inference to these cells.}
\item{source_cell}{Optional character vector of sender cell-type labels,
default NULL. Must occur in the input and retained networks. Extra
source-to-target bubble plots, and comparison rank plots, require both
source_cell and target_cell to be non-empty. Each individual label and the
whole supplied vector are combined with the corresponding target choices.}
\item{target_cell}{Optional character vector of receiver cell-type labels,
default NULL. Same requirements as source_cell. Supplying only one of these
two arguments does not enable the extra source-to-target plots.}
\item{hierarchy_top}{Non-empty character vector of cell-type labels passed
as vertex.receiver for single-group hierarchy plots. Defaults to NULL, but
is required even when hierarchy plotting is disabled; NULL or absent labels
terminate the R session. Labels must survive filtering and be suitable for
the installed CellChat hierarchy layout; there is no layout-specific check.}
\item{plot_top_num}{Numeric vector of edge fractions in (0, 1], defaulting to
c(0.1, 0.25, 0.5, 1). Passed as top to overall circle and differential-circle
plots. For single-group CSV matrices, missing entries become zero and
entries below quantile(matrix, 1 - top) become zero, retaining ties.
Does not control individual-cell circle plots. Range is not validated.}
\item{plot_count_weight}{Character vector drawn from "count" and "weight",
defaulting to c("count", "weight"). Selects interaction counts or summed
weights for circle plots, their matrices, and single-group heatmaps.
Comparison bar plots and heatmaps always include both measures regardless
of this argument. Choices are not validated.}
\item{plot_source_target}{Character vector drawn from "source" and "target",
defaulting to c("source", "target"). Controls sender/receiver individual
circle plots and key-cell comparison bubbles. NULL disables these branches.
Single-group key-cell bubbles always include both directions; explicit
source-to-target bubbles are not controlled by this argument.}
\item{plot_allcirc_allheat_allbub_singcirc_singbub}{Logical vector of length
five with no missing values, default c(TRUE, TRUE, TRUE, TRUE, TRUE).
Positions enable overall circles, overall heatmaps, overall bubbles,
individual-cell circles, and individual-cell bubbles, respectively, in
both single and comparison plotting. The fifth position also gates extra
source-to-target bubbles. Comparison interaction bars and rank plots are
unconditional. Length and missing values are not checked.}
\item{plot_pathway_hier_circ_heat_chord_bar_violin_singRoleheat}{Logical
vector of length seven with no missing values, default
c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE). Positions enable single-group
pathway hierarchy, circle, heatmap, chord, ligand-receptor contribution
bar, gene-expression violin, and signaling-role network heatmap plots.
Also controls creation of pathway directories, including in comparison
mode, but produces no comparison pathway plots. The fourth branch refers
to undefined vertex.receiver; disable it unless that binding is available.
Length and missing values are not checked.}
\item{plot_allsignalingRole_scatt_outheat_inheat}{Logical vector of length
three with no missing values, default c(TRUE, TRUE, TRUE). Enables
single-group signaling-role scatter, outgoing heatmap, and incoming
heatmap, in that order. Ignored by comparison plotting. Length and missing
values are not checked.}
\item{plot_allcommunication_outgo_income}{Logical vector of length two with
no missing values, default c(TRUE, TRUE). Enables outgoing and incoming
single-group communication-pattern plots, respectively. selectK() is
plotted, but identifyCommunicationPatterns() always uses k = 3; the data
must support this rank. Ignored by comparison plotting. Not validated.}
\item{plot_Compute_Aggregate_Analysis_Pvalue}{Numeric significance threshold
in (0, 1], default 0.05. Passed as thresh to computeCommunProbPathway(),
aggregateNet(), and comparison-inference netAnalysis_computeCentrality().
Centrality recomputed during single-group plotting uses CellChat defaults.
This does not set the inference gene-detection thresholds. Not validated.}
\item{plot_Net_netP_Pvalue}{Numeric significance threshold in (0, 1], default
0.05. Passed to netMappingDEG() and subsetCommunication() as thresh, and
to comparison rankNet() as both thresh and cutoff.pvalue. Most plotting
calls otherwise retain their own defaults. Not validated.}
\item{plot_prob_Quantile}{Single quantile probability in [0, 1], default 0.
Bubble candidates must have prob strictly greater than this quantile of
the currently selected communication table. Thus 0 excludes the minimum
and ties at the minimum; 1 leaves no candidates. Recomputed for each cell
selection, with na.rm = TRUE. Not validated.}
\item{plot_geneExp_minpct}{Numeric detection fraction in [0, 1], default 0.1.
The code multiplies this by 100 before passing ligand.pct.1 and
receptor.pct.1 to subsetCommunication(); comparison mode also passes it
for both pct.2 arguments. This assumes percentage-scale CellChat columns;
compatibility with the installed implementation must be checked. Does not
change identifyOverExpressedGenes() inference thresholds. Not validated.}
\item{plot_geneExp_Pvalue}{Numeric gene-expression significance threshold
in (0, 1], default 0.05. Passed for both ligand.pvalues and receptor.pvalues
in comparison bubble-candidate filtering. Single-group filtering passes
NULL for these arguments, so this parameter has no effect there. The
wrapper performs no additional multiple-testing correction. Not validated.}
\item{plot_geneExp_logFC}{Non-negative numeric threshold, default 0.25.
Comparison bubble candidates require both abs(ligand.logFC) and
abs(receptor.logFC) to be at least this value. Opposite signs are allowed.
The log scale is inherited from CellChat. Single-group filtering does not
use this threshold. Not validated.}
\item{ncore}{Positive integer worker count, default 10, passed to
future::plan("multisession", workers = ncore) for fresh inference.
Group-wise future_lapply() uses future.seed = TRUE; no fixed seed is set
by this wrapper. Ignored when cached inference is reused. Not validated.}
\item{formats}{Character vector of formats accepted by basicR::save_figure(),
default c("pdf", "png"). Applies to figures, not CSV, RDS, or log files.
Requires suitable graphics devices and support for the returned plot,
heatmap, grob, and plotting-function objects. Not validated here.}
\item{output_rds_dir}{Writable cache directory, default "./Rdata". Created
recursively. Stores or reads 05_Cellchat_Single_obj.rds and
05_Cellchat_Compare_obj.rds as needed. Existing caches are reused without
checking data, groups, smoothing, thresholds, or package versions. Use a
separate directory for each analysis configuration; see cache limitations.}
\item{output_plot_dir}{Writable figure/table/log directory, default
"./05_Regulation". Created recursively. Relative paths use the working
directory. Existing output files can be overwritten; obsolete files are
not removed. Group and cell labels are used in filenames without sanitising.}
}
\details{
For fresh inference with more than one group, the wrapper first retains
only cell types shared by the first two groups, applying this restriction
to every group, including in "single" mode. It removes cells with missing
cell-type annotations when preparing expression and metadata. It then
creates one CellChat object per group, chooses CellChatDB.human when ACTB
is detected among Seurat feature names, or CellChatDB.mouse when Actb is
detected, with human taking precedence. No other species or identifier
mapping is implemented; absence of both markers leaves the database
undefined. The matching PPI.human or PPI.mouse is used when smoothing.

Each group runs subsetData(), identifyOverExpressedGenes(), and
identifyOverExpressedInteractions(), optional smoothData(), then
computeCommunProb(population.size = FALSE), filterCommunication(min.cells
= 10), computeCommunProbPathway(), and aggregateNet(). The complete chosen
CellChat database is assigned; no signaling category is selected. There
is no sample-level replicate model or covariate adjustment in this wrapper.

Comparison inference additionally computes centrality before mergeCellChat().
Merged differential expression uses group.dataset = "datasets",
pos.dataset = group_levels[2], features.name = "<second>_vs_<first>",
only.pos = FALSE, thresh.pc = 0, thresh.fc = 0, and thresh.p = 1.
Plot comparisons use datasets 1 and 2, and exported differences are the
second matrix minus the first. With extra groups, the positive dataset
can be tested against other datasets while plots still compare only the
first two; use exactly two groups. Increased/decreased bubble labels refer
to max.dataset = 2/1 relative to the second group, not the sign of the gene
log-fold changes. No additional test of a communication difference is
performed by the wrapper's bubble filtering.

Before bubble selection, netMappingDEG(variable.all = TRUE) uses "features"
for single objects or the comparison feature name for merged objects.
Missing mapped P values become 1; missing log-fold changes and detection
percentages become 0. subsetCommunication() receives the thresholds
described above, followed by the absolute-fold-change filter in comparison
mode and the strict probability-quantile filter in both modes. Unique
interaction_name values form pairLR.use. Consequently, an interaction
passing in one selected row can be displayed for other plotted cell pairs
or datasets; these filters do not mask individual bubble rows or alter
the networks used for circle, heatmap, or pathway plots. Candidate mapping
and filtering run even when all bubble switches are FALSE.

The function reads no external script or raw input file. Supply seurat
in memory; cached RDS files are the only file inputs. It attaches Seurat,
CellChat, dplyr, repr, ggsci, ggplot2, patchwork, viridis, NMF, ggalluvial,
grid, and basicR with library(), and attaches future.apply for fresh
inference. It also calls future and, in some comparison circles, cowplot.
Compatible CellChat APIs must provide smoothData(), netMappingDEG(), the
ligand/receptor filtering arguments, and the expected object slots.
CellChat, repr, ggsci, viridis, NMF, ggalluvial, grid, future.apply, future,
and cowplot are not declared in the current DESCRIPTION. Package attachment
and availability of unqualified functions, including plan(), are required.
No automatic installation or external download is implemented.

Fresh inference sets future.globals.maxSize to 100 * 1024^3 bytes and does
not restore it. On success it calls plan("sequential") rather than restoring
the previous plan; errors can leave workers active. Packages remain attached
and some graphics parameters change. Several validation failures call quit()
and terminate the whole R session rather than raising a catchable R error.
Most plotting messages and warnings are suppressed, and many bubble and
selected comparison-rank errors are silently caught, so a completed call
does not guarantee that all requested files were written.

Cache detection uses unanchored regular-expression filename matches but
reading uses exact filenames. In "all" mode both existence checks mistakenly
look for 05_Cellchat_Single_obj.rds. A single cache alone therefore skips
inference and then fails when reading the missing comparison cache. Neither
cache contents nor input compatibility are validated. Use a fresh directory
or a matching pair of caches for "all"; plotting changes are not saved back
into cached objects.

Known plotting limitations include the undefined vertex.receiver in the
single-group chord branch, hard-coded Celltype in hierarchy sizing, and
use of @data on ggplot results in many bubble and comparison-rank branches.
The latter requires a compatible plot representation and can fail for
objects without that slot; some such errors are uncaught. The default
arguments therefore do not guarantee a successful run. The default NULL
key_cell and hierarchy_top fail initial checks, even with all plots disabled.
Numeric ranges, switch-vector lengths, and post-filter label availability
are not comprehensively validated.

Single-group outputs are under 01_Single/<group>; comparison outputs are
under 02_Compare. Single-group pathway figures and hierarchy probability
matrices are under 01_Signa_Pathway within each group directory. Files
include 01_CellChat_all_<measure>_full.csv and top-filtered matrices,
04_CellChat_*_BublePlot.csv, 10_CellChat_SignalingRole_Scatter.csv,
11_CellChat_SignalingRole_Heatmap_<direction>.csv, comparison
02_CellChat_Circle_*_diff.csv, 03_CellChat_Bar_weight.csv (pathway information
flow), and 05_CellChat_Buble_L-R_pairs_*.csv. CSVs include row names.
Figures share descriptive stems and the requested format extensions; they
include interaction bars, circles, heatmaps, bubbles, pathway plots, and
single-group pattern diagnostics, pattern heatmaps, river plots, and dot
plots as enabled. Not every figure has a CSV counterpart. The resolved
options are written to CellChat_Mode_<mode>_Smooth_<logical>_Arguments.log
in output_plot_dir before inference or cache reuse.
}
\value{
No stable analysis-result list is returned. "single" finishes with
invisible NULL from its plotting loop. "compare" and "all" invisibly return
the incidental result of the final comparison-plot branch, usually NULL
or a list of plotting side-effect results when explicit source/target
bubbles run. Do not use this value as a CellChat object.

The analysis results are saved in output_rds_dir. The Single RDS contains
a named list of CellChat objects, one per group in group_levels order,
including expression/metadata, the selected database, inferred
ligand-receptor networks in net, and pathway networks in netP. The Compare
RDS contains one merged CellChat object with group-specific networks,
joint identities/metadata, and differential-expression features named
"<second>_vs_<first>". These are the pre-plotting objects; centrality or
communication patterns added locally during single-group plotting are not
persisted. Read the appropriate RDS explicitly to retrieve results.
}
\examples{
\dontrun{
# Requires real normalised RNA data, compatible CellChat, and all dependencies.
# Use a fresh output directory and labels present in both groups.
seu <- readRDS("normalised_seurat.rds")
seu <- subset(seu, subset = group \%in\% c("Control", "Treatment"))
seu$Celltype <- droplevels(factor(seu$Celltype))
Communication_CellChat(
  seurat = seu,
  group_levels = c("Control", "Treatment"),
  analysis_mode = "single",
  key_cell = "T cells",
  hierarchy_top = "T cells",
  # Disable the chord branch, which refers to undefined vertex.receiver.
  plot_pathway_hier_circ_heat_chord_bar_violin_singRoleheat =
    c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE),
  ncore = 2,
  output_rds_dir = "cellchat_single_cache",
  output_plot_dir = "cellchat_single_plots"
)
single_results <- readRDS("cellchat_single_cache/05_Cellchat_Single_obj.rds")
control_network <- single_results[["Control"]]

# A two-group comparison, subject to the plotting compatibility notes above.
Communication_CellChat(
  seurat = seu,
  group_levels = c("Control", "Treatment"),
  analysis_mode = "compare",
  key_cell = "T cells",
  hierarchy_top = "T cells",
  ncore = 2,
  output_rds_dir = "cellchat_compare_cache",
  output_plot_dir = "cellchat_compare_plots"
)
comparison <- readRDS("cellchat_compare_cache/05_Cellchat_Compare_obj.rds")
}
}
