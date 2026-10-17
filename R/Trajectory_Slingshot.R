#' Infer Slingshot lineages and plot dynamic gene modules
#'
#' Fits Slingshot trajectories on a two-dimensional Seurat embedding, tests
#' gene dynamics with tradeSeq, clusters branch expression profiles, and
#' exports trajectory plots and a GO BP/KEGG-annotated module heatmap.
#'
#' @param seurat A Seurat object, required with no default, convertible using
#' as.SingleCellExperiment(assay = "RNA"). Supply accessible RNA counts and
#' normalised data layers with matching genes and cells, the requested
#' embedding, and cell-type metadata. Both layers are read even with
#' expr_source = "fitted". Gene identifiers must be compatible with the
#' installed ClusterGVis enrichment workflow and its annotation database;
#' this wrapper performs no identifier conversion itself. Retained cells
#' need positive library sizes for the logarithmic offset.
#' @param reduction Single Seurat reduction name. Defaults to "umap.harmony".
#' The embedding must have exactly two columns, all finite coordinates,
#' and the same set of cell names as the converted object. Rows are
#' reordered to the converted cell order, and coordinates are stored under
#' toupper(reduction). No dimensional reduction is computed.
#' @param celltype_col Single name of a metadata column containing non-missing
#' cluster/cell-type labels. Defaults to "Celltype". Converted to a factor
#' with unused levels dropped, and used as Slingshot clusterLabels and for
#' plotting. Labels define trajectory clusters, not treatment comparisons.
#' @param celltype_start Starting cluster label passed as start.clus to Slingshot.
#' Defaults to NULL, leaving start selection to Slingshot. Any supplied
#' label must occur among the retained cell-type factor levels. There is
#' no interactive prompt.
#' @param celltype_end Character vector of terminal cluster labels passed as
#' end.clus to Slingshot. Defaults to NULL, leaving endpoint constraints
#' unspecified. All supplied labels must be observed factor levels.
#' @param celltype_colors Vector of valid cell-type colours, defaulting to
#' basicR::get_colors(number = 50.1). An unnamed vector must contain at
#' least as many colours as observed cell-type factor levels; the first
#' colours are assigned in factor-level order. A named vector must cover
#' all observed cell-type labels and is reordered by those labels. Extra
#' colours are ignored. Too few unnamed colours or missing required labels
#' in a named vector cause a stopifnot() error before trajectory fitting.
#' Use unique names matching cell-type labels to make the mapping explicit.
#' @param species Single string, "human" or "mouse"; defaults to "human".
#' Selects organism = "hsa" with OrgDb = org.Hs.eg.db for human, or
#' organism = "mmu" with OrgDb = org.Mm.eg.db for mouse, in both GO BP and
#' KEGG enrichment calls. Supply gene identifiers appropriate for the
#' selected species. Other values are not validated and leave the local
#' organism and orgdb variables undefined.
#' @param npoints Integer number of heatmap columns per lineage, at least 2;
#' defaults to 100. Used as nPoints for fitted curves or as the number of
#' equal-width observed pseudotime bins. Lineages need non-empty, finite
#' pseudotime ranges; observed-mode interpolation requires at least two
#' non-missing bins when filling gaps. No range validation is implemented.
#' @param num_clusters Positive integer number of k-means gene modules. Defaults
#' to 6. Must not exceed the number of displayed genes with usable,
#' distinct expression profiles or the available module palette. The
#' explicit check only compares the total significant gene count with
#' num_clusters, before truncation by topn_genes.
#' @param min_pct Detection fraction in [0, 1], defaulting to 0.1. A gene must
#' have counts >= 1 in at least ceiling(min_pct * retained_cell_count)
#' cells. The default is 10 percent. The count threshold is fixed; no
#' min_count argument exists. The range is not explicitly validated.
#' @param deg_pname Single choice "p.adj" (default) or "p.value". Selects the
#' global association-test padj or pvalue column for gene filtering and
#' ranking. Choices are validated only after model fitting and export
#' of the full test table; lineage-specific P values do not select genes.
#' @param deg_pvalue Numeric significance threshold in (0, 1], defaulting to
#' 0.05. Only converged genes with finite selected P values strictly
#' below this threshold are retained. The range is not validated.
#' @param topn_genes Positive integer maximum number of significant genes used
#' for the heatmap, modules, and enrichment. Defaults to 1000. Genes are
#' ranked by the selected global P value, then decreasing waldStat.
#' Supply at least num_clusters; fewer genes are used if fewer pass the
#' filter. Also controls heatmap height, even if fewer genes are displayed.
#' @param topn_show Non-negative integer maximum number of gene labels per
#' module. Defaults to 15. Selects the first rows per cluster in
#' st_data$wide.res as returned by clusterData(); there is no subsequent
#' P-value sort within each module. No range validation is implemented.
#' @param enrich_pname Single choice "p.value" (default) or "p.adj". Selects
#' pvalue or p.adjust for post-enrichment filtering and ranking. If the
#' selected column is absent, the code silently falls back to pvalue.
#' Validation occurs after enrichment; "p.adj" also emits a warning
#' about the earlier enrichCluster() cutoff.
#' @param enrich_pvalue Numeric threshold in (0, 1], defaulting to 0.05. Passed
#' as pvalueCutoff to both enrichCluster() calls, then used as a strict
#' upper bound on the finite values selected by enrich_pname. Earlier
#' filtering is controlled by the installed enrichment implementation.
#' The wrapper does not validate the range or recalculate enrichment P values.
#' @param enrich_topn Positive integer maximum number of GO BP terms and KEGG
#' pathways shown per module, separately for each database. Defaults to
#' 10. Display terms are ordered by the selected enrichment P value.
#' No range validation is implemented.
#' @param expr_source Single string specifying branch expression values.
#' Defaults to "fitted", which uses tradeSeq::predictSmooth(). Use
#' "observed" for unweighted bin means from the RNA data layer. Every
#' value other than "fitted" enters the observed branch; spelling is
#' not validated. Both modes still fit tradeSeq and select genes using
#' its association tests.
#' @param ht_colors Three valid colours corresponding to the fixed heatmap
#' Z-score breaks -2, 0, and 2. Defaults to
#' c("#19475E", "#2E8B8B", "#F2E33B"). Passed in htColList to
#' ClusterGVis::visCluster(); no length validation is implemented.
#' @param ncors Positive integer number of workers for
#' BiocParallel::MulticoreParam(), defaulting to 4. fitGAM() always receives
#' parallel = TRUE and RNGseed = 123. Use an environment supporting this
#' backend, normally Unix-like systems; no alternative backend is exposed.
#' @param output_dir Writable output directory, defaulting to "./04_Pseudotime".
#' Created recursively if necessary. Relative paths use the working
#' directory. Existing files may be overwritten; no saved result is reused
#' as a cache, and files from previous runs are not removed.
#' @param formats Character vector of formats accepted by basicR::save_figure(),
#' such as "pdf", "png", "tiff", and "svg". Defaults to c("pdf", "png").
#' Used for all figures. Suitable graphics devices, fonts, and support for
#' the heatmap object returned by visCluster() are required.
#'
#' @details
#' The RNA assay is converted to SingleCellExperiment, embedding rows are
#' matched to cells, and cell-type labels are converted to a factor. Slingshot
#' receives start.clus and end.clus as supplied, omega = FALSE, times = NULL,
#' shrink = TRUE, extend = "y", reweight = TRUE, reassign = TRUE,
#' thresh = 0.001, maxit = 50, and stretch = 1. All inferred lineages are
#' retained. Lineage indices follow the returned ordering, not predefined
#' biological endpoints. Consult the saved lineage paths for direction.
#' No sample or treatment comparisons, fold changes, or endpoint contrasts
#' are performed. Smooth curves are plotted in their stored curve order;
#' the separate node plot connects cell-type coordinate means, not medians.
#'
#' Per-lineage pseudotime is exported with na = TRUE, preserving missing
#' values outside each lineage, together with curve weights. For gene
#' modelling, only cells with a positive sum of lineage weights are retained.
#' No explicit assertion checks agreement of pseudotime NA positions and
#' zero weights, despite an internal comment. Detection filtering uses these
#' retained cells. Complete pseudotime from na = FALSE and curve weights are
#' passed to fitGAM(). Counts for all genes are retained in the input matrix,
#' while genes = which(keep_genes) limits the genes fitted. The offset is
#' log(colSums(counts_heat) * edgeR::calcNormFactors(counts_heat)), computed
#' before gene filtering. fitGAM() uses nknots = 6.
#'
#' associationTest() runs with global = TRUE, lineages = TRUE, and
#' inverse = "eigen". Global and available lineage-specific P values are
#' recomputed with the chi-square upper tail where the Wald statistic and
#' degrees of freedom are finite and degrees of freedom are positive.
#' BH adjustment is applied only to finite global P values of converged
#' genes. The full test table is sorted by padj and decreasing waldStat,
#' regardless of deg_pname. Significant genes are separately selected and
#' ranked by deg_pname; only the first topn_genes enter subsequent analysis.
#' There is no separate export of the complete significant-gene subset.
#'
#' In fitted mode, predictSmooth(tidy = FALSE) supplies the branch matrix.
#' Columns are assumed to occur in lineage blocks of npoints columns. The
#' pseudotime grid for each branch is reconstructed from the t/j and l/j
#' columns (named t1, l1, and so on) in colData(gam_sce)$tradeSeq$dm, using
#' cells with lineage assignment equal to 1. This depends on tradeSeq's
#' internal design-matrix structure. In observed mode, cells with non-missing
#' pseudotime in each lineage are put into equal-width bins. The matrix
#' contains unweighted means of normalised RNA data, with empty bins filled
#' by within-lineage linear interpolation and constant endpoint extension
#' via approx(rule = 2). No values are borrowed from other lineages. The
#' observed-mode exported grid is an evenly spaced sequence including the
#' range endpoints, not the bin midpoints. Degenerate ranges or too few
#' occupied bins are not handled explicitly.
#'
#' clusterData() receives the unscaled branch matrix with scaleData = TRUE,
#' clusterMethod = "kmeans", clusterNum = num_clusters, and
#' kmeansParamsList = list(nstart = 50, iter.max = 100). The returned object
#' is labelled type = "monocle" and geneType = "branched" for presentation;
#' no Monocle analysis is performed. Heatmap columns are split by lineage,
#' with an additional actual-pseudotime colour annotation on a shared range.
#' Gene labels follow the returned cluster table order rather than an
#' explicit within-module significance ranking.
#'
#' GO BP and KEGG enrichment are performed separately for each module
#' formed from the selected topn_genes genes (at most 1000 by default).
#' Both calls use the species-matched OrgDb and organism, topn = NULL,
#' add.gene = TRUE, and seed = 123. The detected-gene vector named
#' background_genes is only used in a message; it is not passed as an
#' enrichment universe and does not expand the selected module gene lists.
#' Background selection and identifier handling therefore follow enrichCluster()
#' defaults. The wrapper passes enrich_pvalue as pvalueCutoff before its
#' own selected-column filter. It warns that "p.adj" is a later filter,
#' but does not inspect or override upstream filtering. If p.adjust is
#' missing, raw pvalue is used even in "p.adj" mode. KEGG errors are caught
#' for up to three attempts, with a five-second sleep after each failed
#' attempt, including the last. A NULL result emits a warning and allows
#' execution to continue; GO BP errors are not caught.
#'
#' Output tables and RDS files under output_dir are as follows. The files
#' 02_Slingshot_pseudotime.csv and 02_Slingshot_weights.csv have cells in
#' rows and lineages in columns, with row names included; pseudotime columns
#' are renamed slingPseudotime_1, slingPseudotime_2, and so on.
#' 02_Slingshot_lineages.csv has lineage indices and cluster paths joined
#' by arrows. 02_Slingshot_sds.rds stores the object returned by slingshot()
#' on the SingleCellExperiment input, including embedding and trajectory
#' information; it is not the separate SlingshotDataSet conversion.
#' 04_Slingshot_tradeSeq_models.rds stores the fitted model object.
#' 05_Slingshot_dynamic_gene_test.csv contains all association-test columns,
#' gene identifiers, convergence flags, and global padj values.
#'
#' 05_Slingshot_branch_expression_matrix.csv contains displayed genes in
#' rows and lineage/grid-point columns, before clusterData() scaling, with
#' gene row names included. 05_Slingshot_branch_pseudotime_grid.csv maps
#' column names to pseudotime values. 06_Slingshot_gene_cluster.csv is the
#' returned st_data$wide.res table, including gene/module assignments and
#' expression-profile columns. The optional 07_Slingshot_enrich_BP_all.csv
#' and 07_Slingshot_enrich_KEGG_all.csv retain enrichment result columns
#' after filtering; "all" does not mean unfiltered tests. Corresponding
#' files ending in _display.csv contain group, Description, pvalue, and
#' ratio for up to enrich_topn terms per module. Here pvalue contains the
#' selected statistic (possibly adjusted P), and ratio is GeneRatio
#' converted from a numerator/denominator string to a percentage, or NA if
#' that column is absent. These tables, the lineage-path table, the test
#' table, the grid table, and the cluster table omit CSV row names. NULL
#' enrichment results are not written, so stale files may remain from
#' previous runs.
#'
#' Figure stems are 03_Slingshot_celltype (all cells and smooth lineages,
#' 7 by 6 inches), 03_Slingshot_pseudotime_nolegend_LineageN (one lineage
#' at a time, excluding its missing pseudotimes, 6 by 6 inches),
#' 03_Slingshot_pseudotime_lineage (combined panels, 4 times the panel-column
#' count by 5 times the panel-row count), and 03_Slingshot_node (cell-type
#' means and cluster paths, 8 by 6 inches), each at res = 600. The combined
#' layout uses bioBulkR::facet_nrow_ncol(). 08_Slingshot_Branch_Heatmap
#' combines expression, module trends, and enrichment annotations at
#' res = 300. Its height is 5 + 0.01 * topn_genes. If the longest displayed
#' term exceeds 130 characters, width is 10 + 2 * number_of_lineages;
#' otherwise it is 4 + 3 * number_of_lineages + 0.1 * longest_term_length.
#' The returned visCluster() object is passed directly to save_figure(),
#' not wrapped in a function as suggested by an internal comment.
#'
#' If either enrichment display table is empty, max(nchar(...)) can warn;
#' if both are empty, the computed heatmap width is non-finite and export
#' can fail. No fallback handles this case. Empty retained cell/gene sets,
#' constant gene profiles, too few displayed genes for k-means, and invalid
#' parameter ranges likewise have no comprehensive checks. Lineage colours
#' use the ggsci JCO palette through basicR and are limited by that palette;
#' module colours use the first num_clusters entries of the fixed basicR
#' palette. The wrapper resets the random seed to 123 before Slingshot,
#' fitGAM(), and clustering and passes seed 123 to enrichment. It changes
#' repr.plot.width, repr.plot.height, and repr.plot.res options without
#' restoring them. Intermediate outputs may remain if later stages fail.
#'
#' The function attaches dplyr, basicR, Seurat, ggplot2, slingshot, Matrix,
#' tradeSeq, ClusterGVis, clusterProfiler, org.Hs.eg.db, org.Mm.eg.db,
#' BiocParallel, and SingleCellExperiment. Both annotation packages are
#' required even for human analysis. It also calls edgeR,
#' SummarizedExperiment, viridis, ggnewscale, bioBulkR, patchwork,
#' ComplexHeatmap, circlize, grid, stats, utils, and
#' bioscR::calculate_cell_point_size(); basicR uses ggsci for lineage
#' colours. Compatible package versions must support the accessed layers,
#' tradeSeq design matrix, and ClusterGVis plotting interface. KEGG access
#' may require network connectivity. No external scripts or input files
#' are read by this wrapper, and no saved results are loaded on rerun.
#' DESCRIPTION currently omits slingshot, Matrix, tradeSeq, ClusterGVis,
#' clusterProfiler, org.Hs.eg.db, org.Mm.eg.db, BiocParallel, edgeR,
#' viridis, bioBulkR, ComplexHeatmap, circlize, and grid; ggsci is also
#' undeclared there as an indirect palette dependency.
#'
#' @return The return value of the final save_figure() call, normally resolved from
#' basicR; its type and visibility depend on that implementation. No
#' structured result list, modified Seurat object, fitted model, or heatmap
#' is explicitly returned. Retrieve trajectory and fitted-model objects
#' from the two RDS files and tabular results from the CSV files described
#' above. The module object and combined plot are not saved as R objects.
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(bioscR)
#' # Requires the listed analysis packages, a compatible parallel environment,
#' # and working enrichment services. This example uses human data; for mouse
#' # data, set species = "mouse" to use the corresponding annotation database.
#' seu <- readRDS("human_seurat.rds")
#' # Supply RNA counts and normalised data, plus a finite two-column embedding.
#' stopifnot("Celltype" %in% colnames(seu@meta.data))
#' stopifnot(!anyNA(seu$Celltype))
#' emb <- Embeddings(seu, "umap.harmony")
#' stopifnot(ncol(emb) == 2L, all(is.finite(emb)))
#' # Replace Progenitor with a starting label present in your data.
#' stopifnot("Progenitor" %in% seu$Celltype)
#' Trajectory_Slingshot(
#'   seurat = seu,
#'   celltype_start = "Progenitor",
#'   species = "human",
#'   expr_source = "fitted",
#'   output_dir = "slingshot_results"
#' )
#' sds <- readRDS("slingshot_results/02_Slingshot_sds.rds")
#' models <- readRDS("slingshot_results/04_Slingshot_tradeSeq_models.rds")
#' # Outputs depend on sufficient significant genes and non-empty enrichment.
#' }
#' @export
Trajectory_Slingshot <- function(seurat,
                                 reduction = "umap.harmony",
                                 celltype_col = "Celltype",
                                 celltype_start = NULL,
                                 celltype_end = NULL,
                                 celltype_colors = basicR::get_colors(number = 50.1),
                                 species = "human",     # mouse / human
                                 npoints = 100,      # 每条分支的拟时序区间数（= 热图列数/分支；多分支版式按每分支 100 列）
                                 num_clusters = 6,   # 基因模块数
                                 min_pct = 0.1,      # 表达过滤：至少在 10% 细胞中检出
                                 deg_pname = "p.adj",    # 动态基因 p.adj/p.value
                                 deg_pvalue = 0.05,      # 动态基因 q 值阈值
                                 topn_genes = 1000,      # 参与热图的动态基因数（按 q 值排序取前 N）
                                 topn_show = 15,         # 每个模块标注的基因数
                                 enrich_pname = "p.value",  # 富集p p.adj/p.value
                                 enrich_pvalue = 0.05,      # 富集原始 p 值阈值
                                 enrich_topn = 10,         # 每个模块展示的通路数
                                 expr_source = "fitted",    # 热图表达值来源：fitted = tradeSeq 平滑拟合曲线（更平滑）；observed = 分箱观测均值
                                 ht_colors = c("#19475E", "#2E8B8B", "#F2E33B"),  # Z-score 配色：低-中-高（青→黄）
                                 ncors = 4,
                                 output_dir = "./04_Pseudotime",
                                 formats = c("pdf","png")
                                ){

    
    suppressPackageStartupMessages({
        library(dplyr)
        library(basicR)
        library(Seurat)
        library(ggplot2)
        library(slingshot)
        library(Matrix)
        library(tradeSeq)
        library(ClusterGVis)
        library(clusterProfiler)
        library(org.Hs.eg.db)
        library(org.Mm.eg.db)
        library(BiocParallel)      
        library(SingleCellExperiment)
    })
    
    # ANSI 颜色代码
    red <- "\033[31m"
    green <- "\033[32m"
    yellow <- "\033[33m"
    blue <- "\033[34m"
    magenta <- "\033[35m"
    cyan <- "\033[36m"
    reset <- "\033[0m"

    #  Create dir  
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    ###################################################################################################################
    #################################################   Main   ########################################################
    ###################################################################################################################

    #-------------------------------------------------   Slingshot   -------------------------------------------------
    # as SingleCellExperiment
    sce <- as.SingleCellExperiment(seurat,assay = "RNA")
    
    # 显式传递降维坐标，核对细胞顺序，避免依赖转换后的降维名称
    embedding <- Embeddings(seurat, reduction = reduction)
    stopifnot(setequal(rownames(embedding), colnames(sce)))
    embedding <- embedding[colnames(sce), , drop = FALSE]
    stopifnot(ncol(embedding) == 2L, all(is.finite(embedding)))
    reducedDim(sce, toupper(reduction)) <- embedding
    sce[[celltype_col]] <- droplevels(factor(sce[[celltype_col]]))
    stopifnot(!anyNA(sce[[celltype_col]]),all(c(celltype_start, celltype_end) %in% levels(sce[[celltype_col]])))

    celltypes <- levels(sce[[celltype_col]])
    
    if (is.null(names(celltype_colors))) {
        # 未命名配色：按照细胞类型的 factor levels 顺序匹配
        stopifnot(length(celltype_colors) >= length(celltypes))
        celltype_colors <- setNames(
            celltype_colors[seq_along(celltypes)],
            celltypes
        )
    } else {
        # 已命名配色：按照细胞类型名称匹配
        stopifnot(all(celltypes %in% names(celltype_colors)))
        celltype_colors <- celltype_colors[celltypes]
    }
    
    # Slingshot
    set.seed(123)
    sds <- slingshot(sce, 
                    reducedDim = toupper(reduction), 
                    clusterLabels = celltype_col, 
                    start.clus = celltype_start,
                    end.clus = celltype_end,
                    omega = FALSE, # 不启用长距离连边限制
                    times = NULL,
                    shrink = TRUE,
                    extend = "y",
                    reweight = TRUE,
                    reassign = TRUE,
                    thresh = 0.001,
                    maxit = 50,
                    stretch = 1
                    )
    # 查看真实路径；分支编号不代表特定生物学终点
    lineage_list <- slingLineages(sds)

    # 提取 UMAP 和 lineage 数据
    dim_data <- as.data.frame(reducedDims(sds)[[toupper(reduction)]])
    colnames(dim_data) <- c("dim_1", "dim_2")
    dim_data[[celltype_col]] <- sce[[celltype_col]]
    
    lineages <- SlingshotDataSet(sds)

    # 按曲线顺序连接投影点，保留全部 lineage
    curves <- slingCurves(sds)
    lineage_df <- lapply(seq_along(curves),function(i){
      curve <- curves[[i]]$s[curves[[i]]$ord, , drop = FALSE]
      curve_df <- as.data.frame(curve)
      colnames(curve_df) <- c("dim_1", "dim_2")
      curve_df$lineage <- as.factor(i)
      return(curve_df)
    })
    lineage_df <- bind_rows(lineage_df)
    
    # 保留原变量名；使用均值中心，与默认 use.median = FALSE 一致
    # 节点图为簇中心连接示意，不是平滑曲线
    dim_median <- dim_data %>% group_by(!!sym(celltype_col)) %>%
          dplyr::summarize(dim_1_median = mean(dim_1),dim_2_median = mean(dim_2),.groups = "drop")

    # 提取 Lineage 的连接点并创建路径数据框
    lineage_paths <- lapply(seq_along(lineage_list),function(i){
      points <- lineage_list[[i]]
      coords <- dim_median %>%
        dplyr::filter(!!sym(celltype_col) %in% points) %>%
        dplyr::arrange(match(!!sym(celltype_col), points))
      coords$lineage <- as.factor(i)
      return(coords)  
    })
    lineage_paths <- bind_rows(lineage_paths)


    # 保留每条 lineage 的拟时序及权重，不跨分支填补 NA
    slingPseudotime_df <- as.data.frame(slingPseudotime(sds, na = TRUE))
    stopifnot(identical(rownames(slingPseudotime_df), rownames(dim_data)))
    slingPseudotime_names <- paste0("slingPseudotime_", seq_len(ncol(slingPseudotime_df)))
    colnames(slingPseudotime_df) <- slingPseudotime_names
    dim_data <- cbind(dim_data, slingPseudotime_df)
    lineage_weights <- slingCurveWeights(sds)
    lineage_colors <- setNames(basicR::get_colors(number = as.numeric(paste0(length(curves),".1")),package = "ggsci",name = "jco"), 
                               as.character(seq_along(curves))
                              )
    write.csv(slingPseudotime_df, file.path(output_dir, "02_Slingshot_pseudotime.csv"))
    write.csv(lineage_weights, file.path(output_dir, "02_Slingshot_weights.csv"))
    saveRDS(sds, file.path(output_dir, "02_Slingshot_sds.rds"))
    df_Slingshot_lineages <- data.frame(lineage = seq_along(lineage_list),path = vapply(lineage_list, paste, character(1), collapse = " -> "))
    write.csv(df_Slingshot_lineages,file.path(output_dir, "02_Slingshot_lineages.csv"), row.names = FALSE)

    
    #-------------------------------------------------   Heatmap   ---------------------------------------------------
    #####  1. 准备 counts、拟时序与分支权重
    # 沿用上方 Slingshot 的 sds；细胞顺序与 counts 列顺序逐一核对
    counts_heat <- Seurat::GetAssayData(seurat, assay = "RNA", layer = "counts")[, colnames(sds), drop = FALSE]
    norm_heat <- Seurat::GetAssayData(seurat, assay = "RNA", layer = "data")[, colnames(sds), drop = FALSE]
    pt_na <- slingPseudotime(sds)              # NA 版：细胞不在该分支上为 NA
    pt_heat <- slingPseudotime(sds, na = FALSE)  # 补全版：tradeSeq 拟合用
    weights_heat <- slingCurveWeights(sds)
    stopifnot(identical(colnames(counts_heat), rownames(pt_heat)))
    
    # 只用有分支权重的细胞；并核对 NA 结构与权重一致（NA 即该分支权重为 0）
    keep_cells <- rowSums(weights_heat) > 0
    counts_heat <- counts_heat[, keep_cells, drop = FALSE]
    norm_heat <- norm_heat[, keep_cells, drop = FALSE]
    pt_heat <- pt_heat[keep_cells, , drop = FALSE]
    pt_na <- pt_na[keep_cells, , drop = FALSE]
    weights_heat <- weights_heat[keep_cells, , drop = FALSE]
    
    # 基因表达过滤：至少在 min_pct 比例的细胞中达到 min_count 条 counts
    keep_genes <- Matrix::rowSums(counts_heat >= 1) >= ceiling(min_pct * ncol(counts_heat))
    background_genes <- rownames(counts_heat)[keep_genes]
    
    lineage_labels <- names(slingLineages(sds))
    n_lineage <- length(lineage_labels)
    stopifnot(n_lineage == ncol(pt_heat))
    message("受检基因 ", length(background_genes), " 个；细胞 ", ncol(counts_heat), " 个")


    #####  2. 动态基因：NB-GAM + associationTest
    # 对应 monocle2 的 differentialGeneTest(~sm.ns(Pseudotime))：
    # 以 slingshot 的拟时序为自变量、分支权重为细胞归属，拟合负二项 GAM 后做全局关联检验
    offset_heat <- log(Matrix::colSums(counts_heat) * edgeR::calcNormFactors(counts_heat))
    set.seed(123)
    gam_sce <- tradeSeq::fitGAM(counts = counts_heat,
                                pseudotime = pt_heat,
                                cellWeights = weights_heat,
                                genes = which(keep_genes),
                                offset = offset_heat,
                                nknots = 6,
                                parallel = TRUE,
                                BPPARAM = BiocParallel::MulticoreParam(workers = ncors,RNGseed = 123),
                                verbose = TRUE
                               )
    saveRDS(gam_sce, file.path(output_dir, "04_Slingshot_tradeSeq_models.rds"))
    
    # 该数据用默认 Cholesky 求逆会产生大量 NA，统一改用特征分解求逆
    assoc_res <- tradeSeq::associationTest(gam_sce, global = TRUE, lineages = TRUE, inverse = "eigen")
    for(suffix in c("", paste0("_", seq_len(n_lineage)))) {
        sc <- paste0("waldStat", suffix); dc <- paste0("df", suffix); pc <- paste0("pvalue", suffix)
        if (is.null(assoc_res[[sc]])) next
        ok <- is.finite(assoc_res[[sc]]) & is.finite(assoc_res[[dc]]) & assoc_res[[dc]] > 0
        # 直接算卡方上尾，避免 1 - pchisq 在极小 P 值处丢精度
        assoc_res[[pc]][ok] <- pchisq(assoc_res[[sc]][ok], df = assoc_res[[dc]][ok], lower.tail = FALSE)
    }
    assoc_res$gene <- rownames(assoc_res)
    assoc_res$converged <- rowData(gam_sce)$tradeSeq$converged
    assoc_res$padj <- NA_real_
    valid <- assoc_res$converged %in% TRUE & is.finite(assoc_res$pvalue)
    assoc_res$padj[valid] <- p.adjust(assoc_res$pvalue[valid], method = "BH")
    assoc_res <- assoc_res[order(assoc_res$padj, -assoc_res$waldStat, na.last = TRUE), ]
    write.csv(assoc_res, file.path(output_dir, "05_Slingshot_dynamic_gene_test.csv"), row.names = FALSE)
    
    # ---- 动态基因筛选：deg_pname 决定用哪一种 P 值 ----
    #   "p.adj"   -> BH 校正 P（assoc_res$padj）
    #   "p.value" -> 原始 P   （assoc_res$pvalue）
    if (!deg_pname %in% c("p.value", "p.adj")) {
        stop("deg_pname 只能是 'p.value' 或 'p.adj'")
    }
    deg_col <- if (deg_pname == "p.value") "pvalue" else "padj"
    
    # 只保留收敛、且所选 P 值可用的基因；按所选 P 值升序（并列时 waldStat 降序）排序
    deg_ok <- assoc_res$converged %in% TRUE & is.finite(assoc_res[[deg_col]])
    assoc_sig <- assoc_res[deg_ok & assoc_res[[deg_col]] < deg_pvalue, , drop = FALSE]
    assoc_sig <- assoc_sig[order(assoc_sig[[deg_col]], -assoc_sig$waldStat), , drop = FALSE]
    diff_genes <- assoc_sig$gene
    if (length(diff_genes) < num_clusters) stop("动态基因不足，请检查拟合与检验结果。")
    
    # 展示基因：按所选 P 值取前 topn_genes 个
    # （branch_mat 的行序 = 这里的顺序，所以每模块标注的也是该模块此指标最小的 topn_show 个）
    show_genes <- head(diff_genes, topn_genes)
    message("动态基因（", deg_pname, " < ", deg_pvalue, "）：",length(diff_genes), " 个；用于热图 ", length(show_genes), " 个")

    ##### 3. 分支动态表达矩阵与聚类
    # 表达矩阵有两种来源，用 expr_source 切换：
    #   fitted   —— tradeSeq 对每条分支的平滑拟合曲线（predictSmooth），平滑、接近 monocle2 分支热图的观感
    #   observed —— 每条分支按拟时序等宽分箱后的箱内平均表达，不含模型假设，但噪声更明显
    # pseudotime_grid 记录每个热图列对应的「实际拟时序」，后面画热图顶部的拟时序色条要用
    pseudotime_grid <- vector("list", n_lineage)
    
    if (expr_source == "fitted") {
        branch_mat <- tradeSeq::predictSmooth(gam_sce, gene = show_genes,
                                              nPoints = npoints, tidy = FALSE)
        stopifnot(ncol(branch_mat) == n_lineage * npoints)
        colnames(branch_mat) <- paste0(rep(lineage_labels, each = npoints), "_",
                                       sprintf("%03d", rep(seq_len(npoints), n_lineage)))
        # 与 tradeSeq 内部一致：分支网格取该分支被指派细胞(l_j==1)的拟时序范围
        dm_pt <- SummarizedExperiment::colData(gam_sce)$tradeSeq$dm
        for (j in seq_len(n_lineage)) {
            tt <- dm_pt[[paste0("t", j)]][dm_pt[[paste0("l", j)]] == 1]
            pseudotime_grid[[j]] <- seq(min(tt), max(tt), length.out = npoints)
        }
    } else {
        branch_mat <- matrix(NA_real_, nrow = length(show_genes), ncol = n_lineage * npoints,
                             dimnames = list(show_genes,
                                 paste0(rep(lineage_labels, each = npoints), "_",
                                        sprintf("%03d", rep(seq_len(npoints), n_lineage)))))
        for (j in seq_len(n_lineage)) {
            idx <- (j - 1) * npoints + seq_len(npoints)
            cells_j <- which(!is.na(pt_na[, j]))          # 落在该分支上的细胞
            tt <- pt_na[cells_j, j]
            pseudotime_grid[[j]] <- seq(min(tt), max(tt), length.out = npoints)
            brk <- seq(min(tt), max(tt), length.out = npoints + 1)
            bin <- findInterval(tt, brk, rightmost.closed = TRUE)
            bin[bin < 1] <- 1; bin[bin > npoints] <- npoints
            for (b in seq_len(npoints)) {
                cs <- cells_j[bin == b]
                if (length(cs) > 0) branch_mat[, idx[b]] <- Matrix::rowMeans(norm_heat[show_genes, cs, drop = FALSE])
            }
            # 没有细胞的区间在该分支内线性插值，不跨分支借值
            for (i in seq_len(nrow(branch_mat))) {
                v <- branch_mat[i, idx]
                if (anyNA(v)) {
                    ok <- which(!is.na(v))
                    branch_mat[i, idx] <- approx(ok, v[ok], xout = seq_len(npoints), rule = 2)$y
                }
            }
        }
    }
    write.csv(branch_mat, file.path(output_dir, "05_Slingshot_branch_expression_matrix.csv"))
    
    # 顶部拟时序色条用的向量：长度 = n_lineage * nPoints，与热图列一一对应
    pseudotime_per_col <- unlist(pseudotime_grid, use.names = FALSE)
    stopifnot(length(pseudotime_per_col) == ncol(branch_mat))
    write.csv(data.frame(column = colnames(branch_mat), pseudotime = pseudotime_per_col),
              file.path(output_dir, "05_Slingshot_branch_pseudotime_grid.csv"), row.names = FALSE)
    message("各分支拟时序范围：", paste(vapply(pseudotime_grid, function(x)
            sprintf("%.2f-%.2f", min(x), max(x)), character(1)), collapse = " / "))
    
    # 聚类：scaleData = TRUE 表示先按基因做 Z-score 再 kmeans
    set.seed(123)
    st_data <- ClusterGVis::clusterData(obj = branch_mat,
                                        scaleData = TRUE,
                                        clusterMethod = "kmeans",
                                        clusterNum = num_clusters,
                                        kmeansParamsList = list(nstart = 50, iter.max = 100)
                                       )
    # clusterData 不写 pseudotime/geneType，这里补上，visCluster 才会按「多分支」版式绘制
    st_data$type <- "monocle"
    st_data$geneType <- "branched"
    st_data$pseudotime <- factor(rep(lineage_labels, each = npoints), levels = lineage_labels)
    write.csv(st_data$wide.res, file.path(output_dir, "06_Slingshot_gene_cluster.csv"), row.names = FALSE)
    message("模块大小：", paste(table(st_data$wide.res$cluster), collapse = " / "))


    #####  4. GO BP 与 KEGG 富集
    if(species == "human"){
        organism = "hsa"
        orgdb = org.Hs.eg.db
    }else if(species == "mouse"){
        organism = "mmu"
        orgdb = org.Mm.eg.db
    }
    
    enrich_bp <- ClusterGVis::enrichCluster(object = st_data, OrgDb = orgdb,
                                            type = "BP", organism = organism,
                                            pvalueCutoff = enrich_pvalue,
                                            topn = NULL, add.gene = TRUE, seed = 123)
    
    enrich_kegg <- NULL
    for (attempt in seq_len(3)) {
        enrich_kegg <- tryCatch(
            ClusterGVis::enrichCluster(object = st_data, OrgDb = orgdb,
                                       type = "KEGG", organism = organism,
                                       pvalueCutoff = enrich_pvalue,
                                       topn = NULL, add.gene = TRUE, seed = 123),
            error = function(e) { message("KEGG 第 ", attempt, " 次失败：", conditionMessage(e)); NULL })
        if (!is.null(enrich_kegg)) break
        Sys.sleep(5)
    }
    if (is.null(enrich_kegg)) warning("KEGG 富集未成功（通常为网络问题），将只保留 GO BP 注释。", call. = FALSE)
    
    # ---- 富集显著性用哪一种 P 值：enrich_pname ----
    #   "p.value" -> 原始 P（$pvalue）；"p.adj" -> BH 校正 P（$p.adjust）
    # 注意：enrichCluster() 内部的 pvalueCutoff 始终按「原始 p 值」筛（与 visCluster_plot.R 一致），
    #       p.adj 模式是在已通过原始 P 的条目里再做一次校正 P 过滤 + 排序。
    if (!enrich_pname %in% c("p.value", "p.adj")) {
        stop("enrich_pname 只能是 'p.value' 或 'p.adj'")
    }
    if (enrich_pname == "p.adj") {
        warning("enrichCluster() 内部按原始 p 值筛选；p.adj 仅用于后续过滤与排序，注意口径。", call. = FALSE)
    }
    enrich_col <- if (enrich_pname == "p.adj") "p.adjust" else "pvalue"
    
    # 全量结果表按所选 P 值过滤后再保存
    filter_enrich_all <- function(x) {
        if (is.null(x) || nrow(x) == 0) return(x)
        pcol <- if (enrich_col %in% names(x)) enrich_col else "pvalue"
        x[is.finite(x[[pcol]]) & x[[pcol]] < enrich_pvalue, , drop = FALSE]
    }
    enrich_bp_all   <- filter_enrich_all(enrich_bp)
    enrich_kegg_all <- filter_enrich_all(enrich_kegg)
    if (!is.null(enrich_bp_all))   write.csv(enrich_bp_all,   file.path(output_dir, "07_Slingshot_enrich_BP_all.csv"),   row.names = FALSE)
    if (!is.null(enrich_kegg_all)) write.csv(enrich_kegg_all, file.path(output_dir, "07_Slingshot_enrich_KEGG_all.csv"), row.names = FALSE)
    
    # 每个模块按所选 P 值取前 enrich_topn 条用于绘图，并换算 GeneRatio 用于柱状图
    top_enrich_for_plot <- function(x, topn) {
        if (is.null(x) || nrow(x) == 0) return(NULL)
        pcol <- if (enrich_col %in% names(x)) enrich_col else "pvalue"
        x <- x[is.finite(x[[pcol]]) & x[[pcol]] < enrich_pvalue, , drop = FALSE]
        x <- x %>% dplyr::group_by(group) %>% dplyr::arrange(.data[[pcol]], .by_group = TRUE) %>%
            dplyr::slice_head(n = topn) %>% dplyr::ungroup()
        x$ratio <- if ("GeneRatio" %in% names(x))
            vapply(strsplit(as.character(x$GeneRatio), "/"),
                   function(z) as.numeric(z[1]) / as.numeric(z[2]) * 100, numeric(1)) else NA_real_
        # 交给 visCluster 画柱子/着色：p.adj 模式下第三列即为校正 P
        data.frame(group = x$group, Description = x$Description,
                   pvalue = x[[pcol]], ratio = x$ratio, check.names = FALSE)
    }
    bp_plot   <- top_enrich_for_plot(enrich_bp,   enrich_topn)
    kegg_plot <- top_enrich_for_plot(enrich_kegg, enrich_topn)
    if (!is.null(bp_plot))   write.csv(bp_plot,   file.path(output_dir, "07_Slingshot_enrich_BP_display.csv"),   row.names = FALSE)
    if (!is.null(kegg_plot)) write.csv(kegg_plot, file.path(output_dir, "07_Slingshot_enrich_KEGG_display.csv"), row.names = FALSE)
    
    # 通路条目与所属模块同色：模块颜色按该模块的条目数重复
    module_cols <- basicR::get_colors(number = 50.1)[seq_len(num_clusters)]
    path_cols <- function(x) {
        if (is.null(x) || nrow(x) == 0) return(NULL)
        rep(module_cols, as.numeric(table(factor(x$group, levels = paste0("C", seq_len(num_clusters))))))
    }
    bp_col   <- path_cols(bp_plot)
    kegg_col <- path_cols(kegg_plot)
    message("GO BP 条目：", if (is.null(bp_plot)) 0 else nrow(bp_plot),"；KEGG 条目：", if (is.null(kegg_plot)) 0 else nrow(kegg_plot))
                                      
    
    ###################################################################################################################
    ##################################################   Plot   #######################################################
    ###################################################################################################################

    #-------------------------------------------------   Linenage   ---------------------------------------------------
    # 所有 Linenage
    ggplot() +
      geom_point(data = dim_data, aes(x = dim_1, y = dim_2, color = !!sym(celltype_col)), 
                 size = bioscR::calculate_cell_point_size(cell_count = nrow(dim_data))) +
      scale_color_manual(values = celltype_colors) +
      guides(color = guide_legend(override.aes = list(size = 4),title = NULL))+
      ggnewscale::new_scale_color()+
      geom_path(data = lineage_df, aes(x = dim_1, y = dim_2, group = lineage, color = lineage), linewidth = 1, show.legend = FALSE,
                arrow = grid::arrow(length = grid::unit(0.3, "cm"),type = "closed",ends = "last",angle = 25)) +
      scale_color_manual(values = lineage_colors) +
      labs(x = paste0(toupper(reduction),1),y = paste0(toupper(reduction),2))+
      theme_minimal(base_size = 16) +
      theme(legend.position = "right",
            panel.grid = element_blank(),
            axis.text = element_blank(),
            axis.title = element_blank()
           ) -> p1
    basicR::save_figure(obj = p1,
                        filename = file.path(output_dir, "03_Slingshot_celltype"),
                        width = 7, height = 6, res = 600,
                        formats = formats
                       )

    # 单个 Linenage
    plot_list <- lapply(slingPseudotime_names,function(spt){
        spt_num <- strsplit(spt,"_")[[1]][2] %>% as.numeric()
        line_df <- lineage_df[lineage_df$lineage == spt_num, , drop = FALSE]
        lineage_color <- lineage_colors[as.character(spt_num)]
        df <- dim_data[!is.na(dim_data[[spt]]),]
        ggplot() +
          geom_point(data = df, aes(x = dim_1, y = dim_2, color = !!sym(spt)), size = bioscR::calculate_cell_point_size(cell_count = nrow(dim_data))) +
          viridis::scale_color_viridis()+
          ggnewscale::new_scale_color()+
          geom_path(data = line_df, aes(x = dim_1, y = dim_2, group = lineage, color = lineage), linewidth = 1, show.legend = FALSE,
                    arrow = grid::arrow(length = grid::unit(0.3, "cm"),type = "closed",ends = "last",angle = 25)) +
          scale_color_manual(values = lineage_color) +
          labs(x = "",y = "", title = paste0("Lineage ", spt_num))+
          theme_minimal(base_size = 16) +
          theme(legend.position = "none",
                panel.grid = element_blank(),
                axis.text = element_blank(),
                plot.title = element_text(hjust = 0.5)
               )
    })
    names(plot_list) <- slingPseudotime_names
    
    # 单个图
    lapply(names(plot_list),function(pname){
        p <- plot_list[[pname]]
        linage_number <- pname %>% strsplit("_") %>% sapply(.,"[",2) 
        basicR::save_figure(obj = p,
                            filename = file.path(output_dir, paste0("03_Slingshot_pseudotime_nolegend_Lineage",linage_number)),
                            width = 6,height = 6,res = 600, formats = formats
                           )
    })
    
    spt_split <- grep("_",slingPseudotime_names,value = T)
    nrow_ncol <- bioBulkR::facet_nrow_ncol(facet_number = length(spt_split))
    combined_plot <- patchwork::wrap_plots(plot_list[spt_split], ncol = nrow_ncol[2])
    
    basicR::save_figure(obj = combined_plot,
                        filename = file.path(output_dir, "03_Slingshot_pseudotime_lineage"),
                        width = 4 * nrow_ncol[2],
                        height = 5 * nrow_ncol[1],
                        res = 600, formats = formats
                       )

    # 节点图保留全部亚群和路径，不删除其他分析遗留的细胞标签
    ggplot() +
      geom_point(data = dim_data, aes(x = dim_1, y = dim_2, color = !!sym(celltype_col)), 
                 size = bioscR::calculate_cell_point_size(cell_count = nrow(dim_data))) +
      geom_point(data = dim_median, aes(x = dim_1_median, y = dim_2_median), color = 'red', size = 3,show.legend = FALSE) +
      geom_path(data = lineage_paths, aes(x = dim_1_median, y = dim_2_median, group = lineage), color = 'red', linewidth = 1,show.legend = FALSE) +
      scale_color_manual(values = celltype_colors) +
      guides(color = guide_legend(override.aes = list(size = 4),title = NULL))+
      theme_minimal(base_size = 16) +
      labs(x = "", y = "") + 
      theme(legend.position = "right",
            panel.grid = element_blank(),
            axis.text = element_blank()
            ) -> p2
    basicR::save_figure(obj = p2,
                        filename = file.path(output_dir, "03_Slingshot_node"),
                        width = 8, height = 6, res = 600,
                        formats = formats
                       ) 

            
    #-------------------------------------------------   Heatmap   ---------------------------------------------------
    ##### 合并热图：表达 + 模块趋势 + GO BP + KEGG
            
    # 每个模块标注前 topn_show 个基因（模块内按聚类顺序）
    markGenes <- st_data$wide.res %>% dplyr::group_by(cluster) %>%
        dplyr::slice_head(n = topn_show) %>% dplyr::pull(gene)
    
    # 分支配色与上方 Slingshot 图保持一致
    branch_cols <- unname(lineage_colors)
    
    # ---- 顶部双行注释：上行 = 分支（Lineage），下行 = 各分支的实际拟时序 ----
    # pseudotime_per_col 来自上一步：每列对应的真实拟时序值，三条分支各自的范围不同
    # （本数据 Lineage1/2/3 的拟时序上限分别约为 16.8 / 18.4 / 8.2，所以三段的颜色深浅不同）
    column_lineage <- factor(rep(lineage_labels, each = npoints), levels = lineage_labels)
    pt_range <- range(pseudotime_per_col)
    heat_anno <- ComplexHeatmap::HeatmapAnnotation(
        Lineage = column_lineage,
        Pseudotime = pseudotime_per_col,
        col = list(Lineage = stats::setNames(branch_cols, lineage_labels),
                   Pseudotime = circlize::colorRamp2(seq(pt_range[1], pt_range[2], length.out = 3),c("#F7FBFF", "#6BAED6", "#08306B"))
                  ),
        annotation_name_gp = grid::gpar(fontsize = 8)
    )
    
    # 图幅（内联显示用）；导出尺寸见下方 save_figure

    nchr_max_bp <- max(nchar(bp_plot$Description))
    nchr_max_kegg <- max(nchar(kegg_plot$Description))
    nchr_max_bp_kegg <- max(nchr_max_bp,nchr_max_kegg)
    if(nchr_max_bp_kegg > 130){ wordwrap = TRUE}else{ wordwrap = FALSE } 
    if(wordwrap){        
        heatmap_width <-  10 + 2 * length(lineage_labels) 
    }else{
        heatmap_width <-  4 + 3 * length(lineage_labels) + nchr_max_bp_kegg * 0.1
    }
    heatmap_height <- 5 + 0.01 * topn_genes
            
    options(repr.plot.width = heatmap_width, repr.plot.height = heatmap_height, repr.plot.res = 100)
    
    ClusterGVis::visCluster(object = st_data,
                            plotType = "both",      # 热图 + 模块趋势线
                            column_names_rot = 45,
                            htColList = list(col_range = c(-2, 0, 2),col_color = ht_colors),  # Z-score 配色：低-中-高（青→黄），与之前的图一致
                            show_row_dend = FALSE, 
                            markGenes = markGenes,
                            markGenesSide = "left",
                            pseudotimeCol = branch_cols,    # 分支配色
                            heatmapAnnotation = heat_anno,  # 顶部注释：分支 + 相对拟时序（替代默认单行注释）
                            columnSplit = column_lineage,   # 与上面注释对应的分列（传 heatmapAnnotation 时必须给）
                            annoTermData = bp_plot,         # GO BP 注释
                            goCol = bp_col,
                            byGo = "anno_block",
                            panelArg = c(2, 0.25, 3, "grey90", NA),          # 左侧趋势面板
                            ggplotPanelArg = c(10, 0.5, 6, "grey90", NA),    # 右侧富集文字块（收窄给热图让位）
                            width = grid::unit(3 * length(lineage_labels), "cm"),
                            annoKeggData = kegg_plot,       # KEGG 注释
                            keggCol = kegg_col,
                            byKegg = "anno_block",
                            clusterOrder = seq_len(num_clusters),
                            raster_quality = 5,  # 热图主体栅格分辨率（否则 IRkernel 下会糊）
                            wordWrap = wordwrap,  # 条目文字换行，避免文字块过宽
                            addNewLine = TRUE,
                            addBar = TRUE
                           ) -> p
    
    # 必须用「函数形式」交给 save_figure：visCluster 返回的是 HeatmapList
    save_figure(obj = p,
                filename = file.path(output_dir, "08_Slingshot_Branch_Heatmap"),
                width = heatmap_width, height = heatmap_height, res = 300, formats = formats
               )        
            
}