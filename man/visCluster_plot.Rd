% Documentation for R/visCluster_plot.R; keep in sync with its opening help comments.
\name{visCluster_plot}
\alias{visCluster_plot}
\title{
Plot cell-type expression profiles with GO and KEGG enrichment
}
\description{
Identifies differential genes in a Seurat object, selects top genes per cell type, and combines expression profiles with GO biological process and KEGG enrichment annotations using ClusterGVis.
}
\usage{
visCluster_plot(seurat,
                            celltype = "Celltype",
                            group_by = "group",
                            deg_pname = "p.value",
                            deg_pvalue = 0.05,
                            logfc_threshold = 0.5,
                            min_pct = 0.1,
                            test_use = "wilcox",
                            only_pos = FALSE,
                            species = 9606,
                            genes_topn_enrich = 50,
                            genes_topn_show = 20,
                            pathway_topn = 10,
                            enrich_pname = "p.value",
                            enrich_pvalue = 0.05
)
}
\arguments{
\item{seurat}{A Seurat object with expression data suitable for Seurat differential expression analysis and gene identifiers compatible with the selected organism annotation database. Uses the active assay and the metadata columns specified below.}
\item{celltype}{Name of a factor metadata column containing cell-type annotations. Remove unused levels and missing values before calling. Defaults to "Celltype".}
\item{group_by}{Name of a factor metadata column defining comparisons. If different from celltype, exactly two observed levels are required, with both groups represented in each cell type; the second level is compared against the first. If identical to celltype, calls Seurat::FindAllMarkers(). Defaults to "group".}
\item{deg_pname}{Differential expression P-value column to filter and rank: "p.value" uses p_val and "p.adj" uses p_val_adj. Defaults to "p.value".}
\item{deg_pvalue}{Strict upper bound for the selected differential expression P value. Defaults to 0.05.}
\item{logfc_threshold}{Passed to logfc.threshold in Seurat::FindMarkers() or Seurat::FindAllMarkers(). Defaults to 0.5.}
\item{min_pct}{Passed to min.pct in the Seurat marker test. Defaults to 0.1.}
\item{test_use}{Differential expression test passed to test.use in the Seurat marker test. Defaults to "wilcox". Results must contain avg_log2FC for downstream ranking.}
\item{only_pos}{Logical passed to only.pos in the Seurat marker test. Defaults to FALSE, retaining both directions. In two-group comparisons, positive fold changes indicate higher expression in the second group level.}
\item{species}{Numeric taxonomy identifier: 9606 for human or 10090 for mouse. No other species are supported. Defaults to 9606.}
\item{genes_topn_enrich}{Positive integer giving the maximum number of differential genes per cluster used for both expression profiles and enrichment. Defaults to 50. NULL is not supported by the current implementation because the selected-gene table would remain undefined.}
\item{genes_topn_show}{Positive integer giving the maximum number of labelled genes per cluster, selected from the genes_topn_enrich subset. Defaults to 20. Increasing this value cannot add genes outside that subset.}
\item{pathway_topn}{Maximum number of enrichment terms per group passed as topn to enrichCluster() for plot annotations. Defaults to 10. Does not limit the returned full enrichment tables.}
\item{enrich_pname}{Enrichment filtering option, "p.value" (default) or "p.adj". Plot annotations always use raw P values. With "p.adj", a warning is issued and the returned full enrichment tables are additionally filtered by p.adjust < enrich_pvalue.}
\item{enrich_pvalue}{Enrichment cutoff passed as pvalueCutoff to enrichCluster(). Also used as the strict adjusted P-value cutoff for returned full tables when enrich_pname = "p.adj". Defaults to 0.05.}
}
\details{
Both metadata columns must be factors whose levels exactly match their observed values. Their level order determines cell-type ordering and the direction of two-group comparisons. When group_by equals celltype, marker testing uses FindAllMarkers() before the function resets active identities; set Idents(seurat) to the intended cell-type factor before calling for compatibility with the installed Seurat version.

Genes passing the differential expression cutoff are ranked within each cluster by increasing selected P value, breaking ties by decreasing absolute avg_log2FC. Up to genes_topn_enrich genes per cluster are passed to prepareDataFromscRNA() with showAverage = TRUE, after setting active identities to celltype. Thus expression profiles use the selected subset, not every significant gene. Labels are selected from the same subset, and duplicate labels are removed. An error is raised if no genes pass the differential expression thresholds.

GO biological process (type = "BP") and KEGG enrichment are run separately with seed = 5201314. Plot annotations use topn = pathway_topn. Additional enrichment calls with topn = NULL produce the full tables returned after the enrichment cutoff. With enrich_pname = "p.adj", these full tables receive an additional adjusted P-value filter, whereas plot annotations retain raw P-value filtering; the displayed terms may therefore differ from the returned tables.

The combined plot uses plot.type = "both", left-side gene labels and expression profiles, GO and KEGG annotation blocks, and enrichment bars. Colours are obtained from basicR::get_colors(). The function does not explicitly save figures or result tables; plot contains the value returned by visCluster(), which may also draw on the current graphics device.

The implementation attaches dplyr, Seurat, ClusterGVis, scRNAtoolVis, org.Hs.eg.db, org.Mm.eg.db, and ggplot2, and also uses basicR. Both organism annotation packages must be installed regardless of species. Every call sources ~/scripts/Function/01_Bio_Data_and_Plot_Function.R; this external script must exist and its dependencies must be available. It is not supplied by bioscR.
}
\value{
A list with plot (the return value of visCluster()) and result, a nested list containing DEGs_all (all genes passing the differential expression cutoff), DEGs_enrich (the selected genes used for expression profiles and enrichment), TopGenes (the subset used for labels), enrichBP (the full GO biological process enrichment table), and enrichKEGG (the full KEGG enrichment table). Enrichment tables may be NULL when no result is available. Full enrichment tables are not restricted by pathway_topn. The function does not return a modified Seurat object.
}
\examples{
\dontrun{
library(Seurat)
# Requires ~/scripts/Function/01_Bio_Data_and_Plot_Function.R
# and the packages listed in Details.
seu <- readRDS("seurat.rds")
seu$Celltype <- droplevels(factor(seu$Celltype))
# Replace control and treatment with the two observed group names.
seu$group <- factor(seu$group, levels = c("control", "treatment"))
res <- visCluster_plot(
  seurat = seu, celltype = "Celltype", group_by = "group",
  deg_pname = "p.adj", species = 9606,
  genes_topn_enrich = 50, genes_topn_show = 10
)
res$plot
res$result$DEGs_all
res$result$enrichBP
res$result$enrichKEGG

# Find cell-type markers instead of comparing two groups within cell types.
Idents(seu) <- seu$Celltype
markers <- visCluster_plot(
  seurat = seu, celltype = "Celltype", group_by = "Celltype",
  only_pos = TRUE
)
}
}
