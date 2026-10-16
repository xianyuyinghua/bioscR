#' Plot cell-type expression profiles with GO and KEGG enrichment
#'
#' Identifies differential genes in a Seurat object, selects top genes per cell
#' type, and combines expression profiles with GO biological process and KEGG
#' enrichment annotations using ClusterGVis.
#'
#' @param seurat A Seurat object with expression data suitable for Seurat
#'   differential expression analysis and gene identifiers compatible with the
#'   selected organism annotation database. Uses the active assay and the
#'   metadata columns specified below.
#' @param celltype Name of a factor metadata column containing cell-type
#'   annotations. Remove unused levels and missing values before calling.
#'   Defaults to "Celltype".
#' @param group_by Name of a factor metadata column defining comparisons. If
#'   different from celltype, exactly two observed levels are required, with
#'   both groups represented in each cell type; the second level is compared
#'   against the first. If identical to celltype, calls
#'   Seurat::FindAllMarkers(). Defaults to "group".
#' @param deg_pname Differential expression P-value column to filter and rank:
#'   "p.value" uses p_val and "p.adj" uses p_val_adj. Defaults to "p.value".
#' @param deg_pvalue Strict upper bound for the selected differential
#'   expression P value. Defaults to 0.05.
#' @param logfc_threshold Passed to logfc.threshold in Seurat::FindMarkers() or
#'   Seurat::FindAllMarkers(). Defaults to 0.5.
#' @param min_pct Passed to min.pct in the Seurat marker test. Defaults to 0.1.
#' @param test_use Differential expression test passed to test.use in the
#'   Seurat marker test. Defaults to "wilcox". Results must contain avg_log2FC
#'   for downstream ranking.
#' @param only_pos Logical passed to only.pos in the Seurat marker test.
#'   Defaults to FALSE, retaining both directions. In two-group comparisons,
#'   positive fold changes indicate higher expression in the second group
#'   level.
#' @param species Numeric taxonomy identifier: 9606 for human or 10090 for
#'   mouse. No other species are supported. Defaults to 9606.
#' @param genes_topn_enrich Positive integer giving the maximum number of
#'   differential genes per cluster used for both expression profiles and
#'   enrichment. Defaults to 50. NULL is not supported by the current
#'   implementation because the selected-gene table would remain undefined.
#' @param genes_topn_show Positive integer giving the maximum number of
#'   labelled genes per cluster, selected from the genes_topn_enrich subset.
#'   Defaults to 20. Increasing this value cannot add genes outside that
#'   subset.
#' @param pathway_topn Maximum number of enrichment terms per group passed as
#'   topn to enrichCluster() for plot annotations. Defaults to 10. Does not
#'   limit the returned full enrichment tables.
#' @param enrich_pname Enrichment filtering option, "p.value" (default) or
#'   "p.adj". Plot annotations always use raw P values. With "p.adj", a warning
#'   is issued and the returned full enrichment tables are additionally
#'   filtered by p.adjust < enrich_pvalue.
#' @param enrich_pvalue Enrichment cutoff passed as pvalueCutoff to
#'   enrichCluster(). Also used as the strict adjusted P-value cutoff for
#'   returned full tables when enrich_pname = "p.adj". Defaults to 0.05.
#'
#' @details
#' Both metadata columns must be factors whose levels exactly match their
#' observed values. Their level order determines cell-type ordering and the
#' direction of two-group comparisons. When group_by equals celltype, marker
#' testing uses FindAllMarkers() before the function resets active identities;
#' set Idents(seurat) to the intended cell-type factor before calling for
#' compatibility with the installed Seurat version.
#'
#' Genes passing the differential expression cutoff are ranked within each
#' cluster by increasing selected P value, breaking ties by decreasing absolute
#' avg_log2FC. Up to genes_topn_enrich genes per cluster are passed to
#' prepareDataFromscRNA() with showAverage = TRUE, after setting active
#' identities to celltype. Thus expression profiles use the selected subset,
#' not every significant gene. Labels are selected from the same subset, and
#' duplicate labels are removed. An error is raised if no genes pass the
#' differential expression thresholds.
#'
#' GO biological process (type = "BP") and KEGG enrichment are run separately
#' with seed = 5201314. Plot annotations use topn = pathway_topn. Additional
#' enrichment calls with topn = NULL produce the full tables returned after the
#' enrichment cutoff. With enrich_pname = "p.adj", these full tables receive an
#' additional adjusted P-value filter, whereas plot annotations retain raw
#' P-value filtering; the displayed terms may therefore differ from the
#' returned tables.
#'
#' The combined plot uses plot.type = "both", left-side gene labels and
#' expression profiles, GO and KEGG annotation blocks, and enrichment bars.
#' Colours are obtained from basicR::get_colors(). The function does not
#' explicitly save figures or result tables; plot contains the value returned
#' by visCluster(), which may also draw on the current graphics device.
#'
#' The implementation attaches dplyr, Seurat, ClusterGVis, scRNAtoolVis,
#' org.Hs.eg.db, org.Mm.eg.db, and ggplot2, and also uses basicR. Both organism
#' annotation packages must be installed regardless of species. Every call
#' sources ~/scripts/Function/01_Bio_Data_and_Plot_Function.R; this external
#' script must exist and its dependencies must be available. It is not supplied
#' by bioscR.
#'
#' @return A list with plot (the return value of visCluster()) and result, a
#'   nested list containing DEGs_all (all genes passing the differential
#'   expression cutoff), DEGs_enrich (the selected genes used for expression
#'   profiles and enrichment), TopGenes (the subset used for labels), enrichBP
#'   (the full GO biological process enrichment table), and enrichKEGG (the
#'   full KEGG enrichment table). Enrichment tables may be NULL when no result
#'   is available. Full enrichment tables are not restricted by pathway_topn.
#'   The function does not return a modified Seurat object.
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' # Requires ~/scripts/Function/01_Bio_Data_and_Plot_Function.R
#' # and the packages listed in Details.
#' seu <- readRDS("seurat.rds")
#' seu$Celltype <- droplevels(factor(seu$Celltype))
#' # Replace control and treatment with the two observed group names.
#' seu$group <- factor(seu$group, levels = c("control", "treatment"))
#' res <- visCluster_plot(
#'   seurat = seu, celltype = "Celltype", group_by = "group",
#'   deg_pname = "p.adj", species = 9606,
#'   genes_topn_enrich = 50, genes_topn_show = 10
#' )
#' res$plot
#' res$result$DEGs_all
#' res$result$enrichBP
#' res$result$enrichKEGG
#'
#' # Find cell-type markers instead of comparing two groups within cell types.
#' Idents(seu) <- seu$Celltype
#' markers <- visCluster_plot(
#'   seurat = seu, celltype = "Celltype", group_by = "Celltype",
#'   only_pos = TRUE
#' )
#' }
#' @export
visCluster_plot <- function(seurat,
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
                           ){

    suppressPackageStartupMessages({
        library(dplyr)
        library(Seurat)
        library(ClusterGVis)
        library(scRNAtoolVis)
        library(org.Hs.eg.db)
        library(org.Mm.eg.db)
        library(ggplot2)
    })

    red <- "\033[31m"
    reset <- "\033[0m"

    source("~/scripts/Function/01_Bio_Data_and_Plot_Function.R")

    ########################################
    # 参数检查
    ########################################
    if(!deg_pname %in% c("p.value", "p.adj")){
        stop("deg_pname must be 'p.value' or 'p.adj'")
    }

    if(!enrich_pname %in% c("p.value", "p.adj")){
        stop("enrich_pname must be 'p.value' or 'p.adj'")
    }

    if(deg_pname == "p.value"){ degpname <- "p_val"}

    if(deg_pname == "p.adj"){ degpname <- "p_val_adj"}

    ########################################
    # 检查 celltype factor
    ########################################

    celltype_levels <- levels(seurat@meta.data[[celltype]])
    celltype_unique <- unique(as.character(seurat@meta.data[[celltype]]))

    if(!setequal(celltype_levels, celltype_unique)){
        cat("\n",red,"Celltype factor levels do not match unique values. Please check factor setting.",reset,"\n")
        flush.console()
        stop(".........")
    }

    ########################################
    # 检查 group factor
    ########################################
    group_levels <- levels(seurat@meta.data[[group_by]])
    group_unique <- unique(as.character(seurat@meta.data[[group_by]]))

    if(!setequal(group_levels, group_unique)){
        cat("\n",red,"group factor levels do not match unique values. Please check factor setting.",reset,"\n")
        flush.console()
        stop(".........")
    }

    ########################################
    # FindMarkers / FindAllMarkers
    ########################################
    if(group_by == celltype){
        seurat_markers_all <- Seurat::FindAllMarkers(seurat,
                                                     group.by = group_by,
                                                     logfc.threshold = logfc_threshold,
                                                     min.pct = min_pct,
                                                     test.use = test_use,
                                                     only.pos = only_pos
                                                    )

        seurat_markers_all <- seurat_markers_all %>% dplyr::filter(.data[[degpname]] < deg_pvalue)
        
    }else{
        if(length(group_levels) != 2){
            cat("\n",red,"groups number is not 2 !",reset,"\n")
            flush.console()
            stop(".........")
        }


        seurat_markers_list <-
            lapply(
                celltype_levels,
                function(ctype){
                    subseurat <- seurat[,seurat@meta.data[[celltype]] == ctype]
                    subseurat@meta.data[[celltype]] <- droplevels(subseurat@meta.data[[celltype]])

                    if("sample" %in% colnames(subseurat@meta.data)){
                        if(is.factor(subseurat@meta.data$sample)){
                            subseurat@meta.data$sample <- droplevels(subseurat@meta.data$sample)
                        }
                    }

                    if(is.factor(subseurat@meta.data[[group_by]])){
                        subseurat@meta.data[[group_by]] <- droplevels(subseurat@meta.data[[group_by]])
                    }

                    result <- Seurat::FindMarkers(subseurat,
                                                  group.by = group_by,
                                                  ident.1 = group_levels[2],
                                                  ident.2 = group_levels[1],
                                                  logfc.threshold = logfc_threshold,
                                                  min.pct = min_pct,
                                                  test.use = test_use,
                                                  only.pos = only_pos
                                                 )

                    result$cluster <- ctype
                    result$gene <- rownames(result)

                    result <- result %>%  dplyr::filter(.data[[degpname]] < deg_pvalue)
                    return(result)
                }
            )


        seurat_markers_all <- dplyr::bind_rows( seurat_markers_list)
    }

    ########################################
    # 检查 DEG
    ########################################

    if(nrow(seurat_markers_all) == 0){
        stop("No DEGs passed the specified thresholds.")
    }
    if(!is.null(genes_topn_enrich)){
        seurat_markers_enrich <- seurat_markers_all %>% dplyr::group_by(cluster) %>%
            dplyr::arrange(.data[[degpname]],dplyr::desc(abs(avg_log2FC)),.by_group = TRUE) %>%
            dplyr::slice_head(n = genes_topn_enrich) %>%
            dplyr::ungroup()
    }
    ########################################
    # Top genes
    #
    # 仅用于热图标签
    # 不用于 enrichment
    ########################################

    seurat_markers_top <- seurat_markers_enrich %>% dplyr::group_by(cluster) %>%
        dplyr::arrange(.data[[degpname]],dplyr::desc(abs(avg_log2FC)),.by_group = TRUE) %>%
        dplyr::slice_head(n = genes_topn_show) %>%
        dplyr::ungroup()

    markGenes <- unique( seurat_markers_top$gene)


    ########################################
    # species
    ########################################

    if(species == 9606){
        organism_name <- "hsa"
        OrgDb_Data <- org.Hs.eg.db
    }else if(species == 10090){
        organism_name <- "mmu"
        OrgDb_Data <- org.Mm.eg.db
    }else{
        stop("species only supports 9606 (human) or 10090 (mouse).")
    }


    ########################################
    # prepareData
    #
    # 这里使用：
    # 所有满足 DEG 阈值的基因
    ########################################

    Idents(seurat) <- seurat@meta.data[[celltype]]
    st_data <- prepareDataFromscRNA(object = seurat,
                                    diffData = seurat_markers_enrich,
                                    showAverage = TRUE
                                   )


    ########################################
    # enrichment
    #
    # 重要：
    #
    # enrichCluster 自己负责
    # pvalueCutoff + topn
    #
    # 不再手工 group_by / slice
    # 防止破坏 C1/C2/... mapping
    ########################################


    if(enrich_pname == "p.value"){

        enrich <- enrichCluster(object = st_data,
                                OrgDb = OrgDb_Data,
                                type = "BP",
                                organism = organism_name,
                                pvalueCutoff = enrich_pvalue,
                                topn = pathway_topn,
                                seed = 5201314
                                )


        enrich.kegg <- enrichCluster(object = st_data,
                                    OrgDb = OrgDb_Data,
                                    type = "KEGG",
                                    organism = organism_name,
                                    pvalueCutoff = enrich_pvalue,
                                    topn = pathway_topn,
                                    seed = 5201314
                                    )

    }else{
        warning(paste0(
                "ClusterGVis::enrichCluster() internally filters ",
                "using raw pvalue. ",
                "enrich_pname = 'p.adj' is not directly supported ",
                "without rewriting enrichCluster(). ",
                "Current enrichment therefore uses p.value."
        ))


        enrich <- enrichCluster(object = st_data,
                                OrgDb = OrgDb_Data,
                                type = "BP",
                                organism = organism_name,
                                pvalueCutoff = enrich_pvalue,
                                topn = pathway_topn,
                                seed = 5201314
                                )


        enrich.kegg <-  enrichCluster(object = st_data,
                                        OrgDb = OrgDb_Data,
                                        type = "KEGG",
                                        organism = organism_name,
                                        pvalueCutoff = enrich_pvalue,
                                        topn = pathway_topn,
                                        seed = 5201314
                                        )
    }


    # ==================== 完整富集结果 ====================
    
    # GO BP：全部显著通路
    enrich_all <- enrichCluster(object = st_data,
                                OrgDb = OrgDb_Data,
                                type = "BP",
                                organism = organism_name,
                                pvalueCutoff = enrich_pvalue,
                                topn = NULL,
                                seed = 5201314
                                )
    
    # KEGG：全部显著通路
    enrich.kegg_all <- enrichCluster(object = st_data,
                                    OrgDb = OrgDb_Data,
                                    type = "KEGG",
                                    organism = organism_name,
                                    pvalueCutoff = enrich_pvalue,
                                    topn = NULL,
                                    seed = 5201314
                                    )
    
    # 如果要求按照调整后 P 值筛选
    if (enrich_pname == "p.adj") {
        if (!is.null(enrich_all)) {
            enrich_all <- enrich_all %>% dplyr::filter(p.adjust < enrich_pvalue)
        }
    
        if (!is.null(enrich.kegg_all)) {
            enrich.kegg_all <- enrich.kegg_all %>% dplyr::filter(p.adjust < enrich_pvalue)
        }
    }
    
    ########################################
    # colors
    ########################################

    colors_all <- c(basicR::get_colors(number = 10.1,package = "ggsci",name = "jco"),basicR::get_colors(number = 100.1))


    if(!is.null(enrich) &&
       nrow(enrich) > 0){
        colors_go <- colors_all[seq_along(unique(enrich$group))]
    }else{
        colors_go <- NULL
    }


    if(!is.null(enrich.kegg) &&
       nrow(enrich.kegg) > 0){
        colors_kegg <- colors_all[seq_along(unique(enrich.kegg$group))]
    }else{
        colors_kegg <- NULL
    }


    ########################################
    # visCluster
    ########################################

    p <- visCluster(object = st_data,
                    plot.type = "both",
                    column_names_rot = 45,
                    show_row_dend = FALSE,
                    markGenes = markGenes, # Top marker labels
                    markGenes.side = "left",
                    annoTerm.data = enrich, # Go
                    go.col = if(!is.null(enrich) && nrow(enrich) > 0){rep(colors_go,times =as.numeric(table(factor(enrich$group,levels =unique(enrich$group)))))}else{NULL},
                    by.go = "anno_block",
                    line.side = "left",
                    cluster.order = seq_len(length(levels(st_data$long.res$cell_type))), # cluster order
                    ggplot.panel.arg =c(10,0.5,12,"grey90",NA),
                    annoKegg.data = enrich.kegg, # KEGG
                    kegg.col = if(!is.null(enrich.kegg) && nrow(enrich.kegg) > 0){ rep(colors_kegg,times =as.numeric(table(factor(enrich.kegg$group,levels =unique(enrich.kegg$group)))))}else{NULL},
                    by.kegg = "anno_block",
                    add.bar = TRUE, 
                    wordWrap = FALSE,
                    addNewLine = FALSE
                   )


    ########################################
    # return
    ########################################

    return(
        list(plot = p, 
             result = list(DEGs_all = seurat_markers_all,DEGs_enrich = seurat_markers_enrich,TopGenes = seurat_markers_top,enrichBP = enrich_all,enrichKEGG = enrich.kegg_all)
            )
    )
}