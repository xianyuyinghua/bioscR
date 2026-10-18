#' Infer and compare CellChat communication networks from Seurat data
#'
#' Runs group-wise CellChat inference, optionally merges groups for comparison,
#' and exports cached objects, network tables, and communication figures.
#'
#' @param seurat A Seurat object, required with no default. Must contain
#' normalised expression accessible as seurat@assays$RNA$data and the metadata
#' columns specified below. The wrapper does not normalise counts or join
#' layers. An RNA Assay5 without a data slot is not supported by this access.
#' Use human or mouse gene symbols, including ACTB or Actb respectively, and
#' matching expression column names and metadata row names.
#' @param celltype Single cell-type metadata column name. Defaults to "Celltype".
#' Supply a factor with observed, non-missing labels and no unused levels;
#' some plots inspect factor levels. Missing annotations are removed before
#' createCellChat(). One hierarchy-plot sizing expression explicitly reads
#' cellchat@meta$Celltype even when a different column is selected.
#' @param group_col Single grouping metadata column name, defaulting to "group".
#' NULL searches column names for the case-insensitive substring "group" and
#' requires exactly one match. An absent specified column or an ambiguous
#' search terminates the R session via quit(). Supply non-missing group labels.
#' @param group_levels Character vector specifying all observed groups in order,
#' or NULL (default). With NULL, the first observed label equal to "control"
#' or "normal", ignoring case, is moved first and other labels keep their
#' observed order; otherwise labels are sorted. Supply each observed group
#' exactly once, with no extras. This is not a group-subsetting argument;
#' subset seurat beforehand. Use exactly two groups for "compare" or "all".
#' The second group is compared with the first. In "single" mode, all listed
#' groups are analysed, but filtering with multiple groups uses only the
#' first two groups to determine shared cell types.
#' @param analysis_mode Single, case-sensitive string: "all" (default), "single",
#' or "compare". "single" infers separate objects per group; "compare" merges
#' group-wise objects; "all" performs both inference workflows independently.
#' Although membership is checked, vectors are not supported by later if()
#' conditions. "all" plots only the first two single-group objects and assumes
#' exactly two groups when indexing its comparison output directory.
#' @param smooth_use Single non-missing logical, default TRUE. Calls smoothData()
#' with the species-matched PPI network before computeCommunProb(). The latter
#' receives population.size = FALSE but no explicit raw.use argument; whether
#' smoothed data are used therefore depends on the installed CellChat default.
#' @param key_cell Non-empty character vector of cell-type labels to highlight
#' in individual circle and bubble plots. Defaults to NULL, but NULL or an
#' empty vector terminates the R session, even if these plots are disabled.
#' Labels must occur in the input and remain available after shared-cell-type
#' filtering. This parameter does not restrict the inference to these cells.
#' @param source_cell Optional character vector of sender cell-type labels,
#' default NULL. Must occur in the input and retained networks. Extra
#' source-to-target bubble plots, and comparison rank plots, require both
#' source_cell and target_cell to be non-empty. Each individual label and the
#' whole supplied vector are combined with the corresponding target choices.
#' @param target_cell Optional character vector of receiver cell-type labels,
#' default NULL. Same requirements as source_cell. Supplying only one of these
#' two arguments does not enable the extra source-to-target plots.
#' @param hierarchy_top Non-empty character vector of cell-type labels passed
#' as vertex.receiver for single-group hierarchy plots. Defaults to NULL, but
#' is required even when hierarchy plotting is disabled; NULL or absent labels
#' terminate the R session. Labels must survive filtering and be suitable for
#' the installed CellChat hierarchy layout; there is no layout-specific check.
#' @param plot_top_num Numeric vector of edge fractions in (0, 1], defaulting to
#' c(0.1, 0.25, 0.5, 1). Passed as top to overall circle and differential-circle
#' plots. For single-group CSV matrices, missing entries become zero and
#' entries below quantile(matrix, 1 - top) become zero, retaining ties.
#' Does not control individual-cell circle plots. Range is not validated.
#' @param plot_count_weight Character vector drawn from "count" and "weight",
#' defaulting to c("count", "weight"). Selects interaction counts or summed
#' weights for circle plots, their matrices, and single-group heatmaps.
#' Comparison bar plots and heatmaps always include both measures regardless
#' of this argument. Choices are not validated.
#' @param plot_source_target Character vector drawn from "source" and "target",
#' defaulting to c("source", "target"). Controls sender/receiver individual
#' circle plots and key-cell comparison bubbles. NULL disables these branches.
#' Single-group key-cell bubbles always include both directions; explicit
#' source-to-target bubbles are not controlled by this argument.
#' @param plot_allcirc_allheat_allbub_singcirc_singbub Logical vector of length
#' five with no missing values, default c(TRUE, TRUE, TRUE, TRUE, TRUE).
#' Positions enable overall circles, overall heatmaps, overall bubbles,
#' individual-cell circles, and individual-cell bubbles, respectively, in
#' both single and comparison plotting. The fifth position also gates extra
#' source-to-target bubbles. Comparison interaction bars and rank plots are
#' unconditional. Length and missing values are not checked.
#' @param plot_pathway_hier_circ_heat_chord_bar_violin_singRoleheat Logical
#' vector of length seven with no missing values, default
#' c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE). Positions enable single-group
#' pathway hierarchy, circle, heatmap, chord, ligand-receptor contribution
#' bar, gene-expression violin, and signaling-role network heatmap plots.
#' Also controls creation of pathway directories, including in comparison
#' mode, but produces no comparison pathway plots. The fourth branch refers
#' to undefined vertex.receiver; disable it unless that binding is available.
#' Length and missing values are not checked.
#' @param plot_allsignalingRole_scatt_outheat_inheat Logical vector of length
#' three with no missing values, default c(TRUE, TRUE, TRUE). Enables
#' single-group signaling-role scatter, outgoing heatmap, and incoming
#' heatmap, in that order. Ignored by comparison plotting. Length and missing
#' values are not checked.
#' @param plot_allcommunication_outgo_income Logical vector of length two with
#' no missing values, default c(TRUE, TRUE). Enables outgoing and incoming
#' single-group communication-pattern plots, respectively. selectK() is
#' plotted, but identifyCommunicationPatterns() always uses k = 3; the data
#' must support this rank. Ignored by comparison plotting. Not validated.
#' @param plot_Compute_Aggregate_Analysis_Pvalue Numeric significance threshold
#' in (0, 1], default 0.05. Passed as thresh to computeCommunProbPathway(),
#' aggregateNet(), and comparison-inference netAnalysis_computeCentrality().
#' Centrality recomputed during single-group plotting uses CellChat defaults.
#' This does not set the inference gene-detection thresholds. Not validated.
#' @param plot_Net_netP_Pvalue Numeric significance threshold in (0, 1], default
#' 0.05. Passed to netMappingDEG() and subsetCommunication() as thresh, and
#' to comparison rankNet() as both thresh and cutoff.pvalue. Most plotting
#' calls otherwise retain their own defaults. Not validated.
#' @param plot_prob_Quantile Single quantile probability in [0, 1], default 0.
#' Bubble candidates must have prob strictly greater than this quantile of
#' the currently selected communication table. Thus 0 excludes the minimum
#' and ties at the minimum; 1 leaves no candidates. Recomputed for each cell
#' selection, with na.rm = TRUE. Not validated.
#' @param plot_geneExp_minpct Numeric detection fraction in [0, 1], default 0.1.
#' The code multiplies this by 100 before passing ligand.pct.1 and
#' receptor.pct.1 to subsetCommunication(); comparison mode also passes it
#' for both pct.2 arguments. This assumes percentage-scale CellChat columns;
#' compatibility with the installed implementation must be checked. Does not
#' change identifyOverExpressedGenes() inference thresholds. Not validated.
#' @param plot_geneExp_Pvalue Numeric gene-expression significance threshold
#' in (0, 1], default 0.05. Passed for both ligand.pvalues and receptor.pvalues
#' in comparison bubble-candidate filtering. Single-group filtering passes
#' NULL for these arguments, so this parameter has no effect there. The
#' wrapper performs no additional multiple-testing correction. Not validated.
#' @param plot_geneExp_logFC Non-negative numeric threshold, default 0.25.
#' Comparison bubble candidates require both abs(ligand.logFC) and
#' abs(receptor.logFC) to be at least this value. Opposite signs are allowed.
#' The log scale is inherited from CellChat. Single-group filtering does not
#' use this threshold. Not validated.
#' @param ncore Positive integer worker count, default 10, passed to
#' future::plan("multisession", workers = ncore) for fresh inference.
#' Group-wise future_lapply() uses future.seed = TRUE; no fixed seed is set
#' by this wrapper. Ignored when cached inference is reused. Not validated.
#' @param formats Character vector of formats accepted by basicR::save_figure(),
#' default c("pdf", "png"). Applies to figures, not CSV, RDS, or log files.
#' Requires suitable graphics devices and support for the returned plot,
#' heatmap, grob, and plotting-function objects. Not validated here.
#' @param output_rds_dir Writable cache directory, default "./Rdata". Created
#' recursively. Stores or reads 05_Cellchat_Single_obj.rds and
#' 05_Cellchat_Compare_obj.rds as needed. Existing caches are reused without
#' checking data, groups, smoothing, thresholds, or package versions. Use a
#' separate directory for each analysis configuration; see cache limitations.
#' @param output_plot_dir Writable figure/table/log directory, default
#' "./05_Regulation". Created recursively. Relative paths use the working
#' directory. Existing output files can be overwritten; obsolete files are
#' not removed. Group and cell labels are used in filenames without sanitising.
#'
#' @details
#' For fresh inference with more than one group, the wrapper first retains
#' only cell types shared by the first two groups, applying this restriction
#' to every group, including in "single" mode. It removes cells with missing
#' cell-type annotations when preparing expression and metadata. It then
#' creates one CellChat object per group, chooses CellChatDB.human when ACTB
#' is detected among Seurat feature names, or CellChatDB.mouse when Actb is
#' detected, with human taking precedence. No other species or identifier
#' mapping is implemented; absence of both markers leaves the database
#' undefined. The matching PPI.human or PPI.mouse is used when smoothing.
#'
#' Each group runs subsetData(), identifyOverExpressedGenes(), and
#' identifyOverExpressedInteractions(), optional smoothData(), then
#' computeCommunProb(population.size = FALSE), filterCommunication(min.cells
#' = 10), computeCommunProbPathway(), and aggregateNet(). The complete chosen
#' CellChat database is assigned; no signaling category is selected. There
#' is no sample-level replicate model or covariate adjustment in this wrapper.
#'
#' Comparison inference additionally computes centrality before mergeCellChat().
#' Merged differential expression uses group.dataset = "datasets",
#' pos.dataset = group_levels[2], features.name = "<second>_vs_<first>",
#' only.pos = FALSE, thresh.pc = 0, thresh.fc = 0, and thresh.p = 1.
#' Plot comparisons use datasets 1 and 2, and exported differences are the
#' second matrix minus the first. With extra groups, the positive dataset
#' can be tested against other datasets while plots still compare only the
#' first two; use exactly two groups. Increased/decreased bubble labels refer
#' to max.dataset = 2/1 relative to the second group, not the sign of the gene
#' log-fold changes. No additional test of a communication difference is
#' performed by the wrapper's bubble filtering.
#'
#' Before bubble selection, netMappingDEG(variable.all = TRUE) uses "features"
#' for single objects or the comparison feature name for merged objects.
#' Missing mapped P values become 1; missing log-fold changes and detection
#' percentages become 0. subsetCommunication() receives the thresholds
#' described above, followed by the absolute-fold-change filter in comparison
#' mode and the strict probability-quantile filter in both modes. Unique
#' interaction_name values form pairLR.use. Consequently, an interaction
#' passing in one selected row can be displayed for other plotted cell pairs
#' or datasets; these filters do not mask individual bubble rows or alter
#' the networks used for circle, heatmap, or pathway plots. Candidate mapping
#' and filtering run even when all bubble switches are FALSE.
#'
#' The function reads no external script or raw input file. Supply seurat
#' in memory; cached RDS files are the only file inputs. It attaches Seurat,
#' CellChat, dplyr, repr, ggsci, ggplot2, patchwork, viridis, NMF, ggalluvial,
#' grid, and basicR with library(), and attaches future.apply for fresh
#' inference. It also calls future and, in some comparison circles, cowplot.
#' Compatible CellChat APIs must provide smoothData(), netMappingDEG(), the
#' ligand/receptor filtering arguments, and the expected object slots.
#' CellChat, repr, ggsci, viridis, NMF, ggalluvial, grid, future.apply, future,
#' and cowplot are not declared in the current DESCRIPTION. Package attachment
#' and availability of unqualified functions, including plan(), are required.
#' No automatic installation or external download is implemented.
#'
#' Fresh inference sets future.globals.maxSize to 100 * 1024^3 bytes and does
#' not restore it. On success it calls plan("sequential") rather than restoring
#' the previous plan; errors can leave workers active. Packages remain attached
#' and some graphics parameters change. Several validation failures call quit()
#' and terminate the whole R session rather than raising a catchable R error.
#' Most plotting messages and warnings are suppressed, and many bubble and
#' selected comparison-rank errors are silently caught, so a completed call
#' does not guarantee that all requested files were written.
#'
#' Cache detection uses unanchored regular-expression filename matches but
#' reading uses exact filenames. In "all" mode both existence checks mistakenly
#' look for 05_Cellchat_Single_obj.rds. A single cache alone therefore skips
#' inference and then fails when reading the missing comparison cache. Neither
#' cache contents nor input compatibility are validated. Use a fresh directory
#' or a matching pair of caches for "all"; plotting changes are not saved back
#' into cached objects.
#'
#' Known plotting limitations include the undefined vertex.receiver in the
#' single-group chord branch, hard-coded Celltype in hierarchy sizing, and
#' use of @data on ggplot results in many bubble and comparison-rank branches.
#' The latter requires a compatible plot representation and can fail for
#' objects without that slot; some such errors are uncaught. The default
#' arguments therefore do not guarantee a successful run. The default NULL
#' key_cell and hierarchy_top fail initial checks, even with all plots disabled.
#' Numeric ranges, switch-vector lengths, and post-filter label availability
#' are not comprehensively validated.
#'
#' Single-group outputs are under 01_Single/<group>; comparison outputs are
#' under 02_Compare. Single-group pathway figures and hierarchy probability
#' matrices are under 01_Signa_Pathway within each group directory. Files
#' include 01_CellChat_all_<measure>_full.csv and top-filtered matrices,
#' 04_CellChat_*_BublePlot.csv, 10_CellChat_SignalingRole_Scatter.csv,
#' 11_CellChat_SignalingRole_Heatmap_<direction>.csv, comparison
#' 02_CellChat_Circle_*_diff.csv, 03_CellChat_Bar_weight.csv (pathway information
#' flow), and 05_CellChat_Buble_L-R_pairs_*.csv. CSVs include row names.
#' Figures share descriptive stems and the requested format extensions; they
#' include interaction bars, circles, heatmaps, bubbles, pathway plots, and
#' single-group pattern diagnostics, pattern heatmaps, river plots, and dot
#' plots as enabled. Not every figure has a CSV counterpart. The resolved
#' options are written to CellChat_Mode_<mode>_Smooth_<logical>_Arguments.log
#' in output_plot_dir before inference or cache reuse.
#'
#' @return No stable analysis-result list is returned. "single" finishes with
#' invisible NULL from its plotting loop. "compare" and "all" invisibly return
#' the incidental result of the final comparison-plot branch, usually NULL
#' or a list of plotting side-effect results when explicit source/target
#' bubbles run. Do not use this value as a CellChat object.
#'
#' The analysis results are saved in output_rds_dir. The Single RDS contains
#' a named list of CellChat objects, one per group in group_levels order,
#' including expression/metadata, the selected database, inferred
#' ligand-receptor networks in net, and pathway networks in netP. The Compare
#' RDS contains one merged CellChat object with group-specific networks,
#' joint identities/metadata, and differential-expression features named
#' "<second>_vs_<first>". These are the pre-plotting objects; centrality or
#' communication patterns added locally during single-group plotting are not
#' persisted. Read the appropriate RDS explicitly to retrieve results.
#'
#' @examples
#' \dontrun{
#' # Requires real normalised RNA data, compatible CellChat, and all dependencies.
#' # Use a fresh output directory and labels present in both groups.
#' seu <- readRDS("normalised_seurat.rds")
#' seu <- subset(seu, subset = group %in% c("Control", "Treatment"))
#' seu$Celltype <- droplevels(factor(seu$Celltype))
#' Communication_CellChat(
#'   seurat = seu,
#'   group_levels = c("Control", "Treatment"),
#'   analysis_mode = "single",
#'   key_cell = "T cells",
#'   hierarchy_top = "T cells",
#'   # Disable the chord branch, which refers to undefined vertex.receiver.
#'   plot_pathway_hier_circ_heat_chord_bar_violin_singRoleheat =
#'     c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE),
#'   ncore = 2,
#'   output_rds_dir = "cellchat_single_cache",
#'   output_plot_dir = "cellchat_single_plots"
#' )
#' single_results <- readRDS("cellchat_single_cache/05_Cellchat_Single_obj.rds")
#' control_network <- single_results[["Control"]]
#'
#' # A two-group comparison, subject to the plotting compatibility notes above.
#' Communication_CellChat(
#'   seurat = seu,
#'   group_levels = c("Control", "Treatment"),
#'   analysis_mode = "compare",
#'   key_cell = "T cells",
#'   hierarchy_top = "T cells",
#'   ncore = 2,
#'   output_rds_dir = "cellchat_compare_cache",
#'   output_plot_dir = "cellchat_compare_plots"
#' )
#' comparison <- readRDS("cellchat_compare_cache/05_Cellchat_Compare_obj.rds")
#' }
#' @export
Communication_CellChat <- function(seurat,
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
                                  ){

# CellChat Analysis


# ANSI 颜色代码
red <- "\033[31m"
green <- "\033[32m"
yellow <- "\033[33m"
blue <- "\033[34m"
magenta <- "\033[35m"
cyan <- "\033[36m"
reset <- "\033[0m"

# colors
colors <- c('#4DBBD5FF','#E64B35FF','#00A087FF','#3C5488FF','#925E9FFF','#91D1C2FF',
            '#8491B4FF','#7E6148FF','#0072B5FF','#E18727FF','#B09C85FF','#20854EFF',
            '#6F99ADFF','#FFDC91FF','#00468BFF','#FDAF91FF','#B24745FF','#6699FFFF',
            '#99991EFF','#FFCCCCFF','#358000FF','#99CCFFFF','#FFCC00FF','#8C564BFF',
            '#BCBD22FF','#996600FF','#5CB85CFF','#F39B7FFF','#CE3D32FF','#749B58FF',
            '#466983FF','#F0E685FF','#D595A7FF','#924822FF','#7A65A5FF','#C75127FF',
            '#FFA319FF','#8A9045FF','#8F3931FF','#00AF66FF','#748AA6FF','#D0DFE6FF',
            '#C71000FF','#008EA0FF','#8A4198FF','#D5E4A2FF','#5A9599FF','#FF6348FF',
            '#B7E4F9FF','#FF95A8FF','#526E2DFF','#FB6467FF','#E89242FF','#69C8ECFF',
            '#917C5DFF','#FED439FF','#709AE1FF','#D2AF81FF','#FD7446FF')



suppressPackageStartupMessages({
    library(Seurat)
    library(CellChat)
    library(dplyr)
    library(repr)
    library(ggsci)
    library(ggplot2)
    library(patchwork)
    library(viridis)
    library(NMF)
    library(ggalluvial)
    library(grid)
    library(basicR)
})

############################################################################  function  #########################################################################
# Group levels set
group_new_levels <- function(vector,group_levels = NULL){
    # 设置分组水平
    
    # colors
    red <- "\033[31m"
    green <- "\033[32m"
    yellow <- "\033[33m"
    blue <- "\033[34m"
    magenta <- "\033[35m"
    cyan <- "\033[36m"
    reset <- "\033[0m"
    
    # 获取列名
    vector <- unique(vector) %>% as.character()
    if(length(group_levels) == 0){
        if(tolower(vector) %in% c("control", "normal") %>% any()){

            # 找出小写化后等于 "control" 或 "normal" 的列名
            target_col <- vector[tolower(vector) %in% c("control", "normal")]
            
            # 如果存在符合条件的列
            if (length(target_col) > 0) {
              # 将目标列放在第一位，其余列保持原顺序但去掉目标列
              new_order <- c(target_col[1], setdiff(vector, target_col[1]))
              return(new_order)  
            }
        }else{
            return(sort(vector))  
        }
    }else{
        if(all(vector %in% group_levels)){
            return(group_levels)
        }else{
            cat("\n",yellow,"group_levels Not Match the vector!",reset,"\n")
        }
    }
}




# 把函数参数收进 opt：后续正文与内部函数沿用脚本里 opt$xxx 的写法，
# 这样内部函数只需传一个 opt，不必逐个列出十几个参数
opt <- list(celltype = celltype,
            group_col = group_col,
            group_levels = group_levels,
            analysis_mode = analysis_mode,
            smooth_use = smooth_use,
            key_cell = key_cell,
            source_cell = source_cell,
            target_cell = target_cell,
            hierarchy_top = hierarchy_top,
            plot_top_num = plot_top_num,
            plot_count_weight = plot_count_weight,
            plot_source_target = plot_source_target,
            plot_allcirc_allheat_allbub_singcirc_singbub = plot_allcirc_allheat_allbub_singcirc_singbub,
            plot_pathway_hier_circ_heat_chord_bar_violin_singRoleheat = plot_pathway_hier_circ_heat_chord_bar_violin_singRoleheat,
            plot_allsignalingRole_scatt_outheat_inheat = plot_allsignalingRole_scatt_outheat_inheat,
            plot_allcommunication_outgo_income = plot_allcommunication_outgo_income,
            plot_Compute_Aggregate_Analysis_Pvalue = plot_Compute_Aggregate_Analysis_Pvalue,
            plot_Net_netP_Pvalue = plot_Net_netP_Pvalue,
            plot_prob_Quantile = plot_prob_Quantile,
            plot_geneExp_minpct = plot_geneExp_minpct,
            plot_geneExp_Pvalue = plot_geneExp_Pvalue,
            plot_geneExp_logFC = plot_geneExp_logFC,
            ncore = ncore,
            formats = formats,
            output_rds_dir = output_rds_dir,
            output_plot_dir = output_plot_dir
           )

# group
seurat_group_name <- colnames(seurat@meta.data)[grepl("group",tolower(colnames(seurat@meta.data)))]
if(is.null(opt$group_col)){
    if(length(seurat_group_name) == 1 ){
        opt$group_col <- seurat_group_name
    }else if(length(seurat_group_name) == 0){
        cat("\n",red,"The Seurat object does not contain the group column.",reset,"\n")
        flush.console()
        quit(save = "no", status = 0)
    }else{
        cat("\n",yellow,"The Seurat object group column: ",reset,seurat_group_name,"\n")
        flush.console()
        quit(save = "no", status = 0)        
    }
}else{
    if(!(opt$group_col %in% colnames(seurat@meta.data))){
        cat("\n",yellow,"The Seurat object group column: ",reset,seurat_group_name,"\n")
        flush.console()
        quit(save = "no", status = 0) 
    }
}

# group_levels
opt$group_levels <- group_new_levels(vector = unique(seurat@meta.data[[opt$group_col]]),group_levels = opt$group_levels)
# analysis_mode
if(!all(opt$analysis_mode %in% c("all","single","compare")) ){
    cat("\n",red,opt$analysis_mode," not in 'all,single,compare'!",reset,"\n")
    quit(save = "no", status = 0)
}
# key_cell
if( length(opt$key_cell) == 0 ){
    cat("\n",yellow,"Please set the key_cell.",reset,"\n")
    flush.console()
    quit(save = "no", status = 0) 
}else{
    if( !all(opt$key_cell %in% seurat@meta.data[[opt$celltype]] ) ){
        key_cells_setdiff <- setdiff(opt$key_cell,unique(seurat@meta.data[[opt$celltype]])) 
        cat("\n",red,key_cells_setdiff,"Not in Celltype!",reset,"\n")
        flush.console()
        quit(save = "no", status = 0)
    }
}

# hierarchy_top
if( length(opt$hierarchy_top) == 0 ){
    cat("\n",yellow,"Please set the hierarchy_top.",reset,"\n")
    flush.console()
    quit(save = "no", status = 0) 
}else{
    if( !all(opt$hierarchy_top %in% seurat@meta.data[[opt$celltype]] ) ){
        hierarchy_top_setdiff <- setdiff(opt$hierarchy_top,unique(seurat@meta.data[[opt$celltype]])) 
        cat("\n",red,hierarchy_top_setdiff,"Not in Celltype!",reset,"\n")
        flush.console()
        quit(save = "no", status = 0)
    }
}
# Source Cell
if( length(opt$source_cell) != 0 ){
    if( !all(opt$source_cell %in% seurat@meta.data[[opt$celltype]] ) ){
        source_cell_setdiff <- setdiff(opt$source_cell,unique(seurat@meta.data[[opt$celltype]])) 
        cat("\n",red,source_cell_setdiff,"Not in Celltype!",reset,"\n")
        flush.console()
        quit(save = "no", status = 0)
    }
}
# Target Cell
if( length(opt$target_cell) != 0 ){
    if( !all(opt$target_cell %in% seurat@meta.data[[opt$celltype]] ) ){
        target_cell_setdiff <- setdiff(opt$target_cell,unique(seurat@meta.data[[opt$celltype]])) 
        cat("\n",red,target_cell_setdiff,"Not in Celltype!",reset,"\n")
        flush.console()
        quit(save = "no", status = 0)
    }
}


# print opt 
cat("\n\n")
print(str(opt))
cat("\n\n")


# dir
if (!dir.exists(opt$output_rds_dir)){dir.create(opt$output_rds_dir, recursive = T)}
if (!dir.exists(opt$output_plot_dir)){dir.create(opt$output_plot_dir,recursive = T)}

# 保存参数设置
sink(paste0(opt$output_plot_dir,"/CellChat_Mode_",opt$analysis_mode,"_Smooth_",opt$smooth_use,"_Arguments.log"))
str(opt,list.len = Inf,max.level = Inf,vec.len = 1000)
#dput(opt)
sink()

#####################################################################################################################
##################################################  Main Function  ##################################################
#####################################################################################################################
# Compare模式运行函数
run_compare_analysis <- function(seurat_all, group_levels, opt) {
    library(future.apply) 
    
    Idents(seurat_all) <- opt$group_col
    seurat_all_tmp <- SplitObject(seurat_all)
    seurat_all_tmp <- seurat_all_tmp[match(group_levels, names(seurat_all_tmp))]
    
    compare_list <- future_lapply(names(seurat_all_tmp), function(i){
        seurat <- seurat_all_tmp[[i]]
        data <- seurat@assays$RNA$data
        metadata <- seurat@meta.data[!is.na(seurat@meta.data[[opt$celltype]]),]
        data <- data[,rownames(metadata)]
        
        cellchat <- createCellChat(data, meta = metadata, group.by = opt$celltype)
        cellchat <- addMeta(cellchat, meta = metadata)
        
        if (any(grepl("\\bACTB\\b", rownames(seurat)))) {
            CellChatDB <- CellChatDB.human
        }else if(any(grepl("\\bActb\\b", rownames(seurat)))){
            CellChatDB <- CellChatDB.mouse
        }
        cellchat@DB <- CellChatDB
        
        cellchat <- subsetData(cellchat)
        #future::plan("multisession", workers = ncor)
        cellchat <- identifyOverExpressedGenes(cellchat)
        cellchat <- identifyOverExpressedInteractions(cellchat)
        
        if (opt$smooth_use) {
            if (any(grepl("\\bACTB\\b", rownames(seurat)))) {
                cellchat <- smoothData(cellchat, adj = PPI.human)
            }else if(any(grepl("\\bActb\\b", rownames(seurat)))){
                cellchat <- smoothData(cellchat, adj = PPI.mouse)
            }
        }
        
        cellchat <- computeCommunProb(cellchat, population.size = FALSE)
        cellchat <- filterCommunication(cellchat, min.cells = 10)
        cellchat <- computeCommunProbPathway(cellchat,thresh = opt$plot_Compute_Aggregate_Analysis_Pvalue)
        cellchat <- aggregateNet(cellchat,thresh = opt$plot_Compute_Aggregate_Analysis_Pvalue)
        cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP",thresh = opt$plot_Compute_Aggregate_Analysis_Pvalue)
        
        return(cellchat)
    }, future.seed = TRUE)
    
    names(compare_list) <- names(seurat_all_tmp)
    
    cellchat <- mergeCellChat(compare_list, add.names = names(compare_list))
    cellchat <- identifyOverExpressedGenes(cellchat, 
                                           group.dataset = "datasets",
                                           pos.dataset = group_levels[2],
                                           features.name = paste0(group_levels[2],"_vs_",group_levels[1]),
                                           only.pos = FALSE,
                                           thresh.pc = 0,
                                           thresh.fc = 0,
                                           thresh.p = 1
                                          )
    return(cellchat)
}

# Single模式运行函数
run_single_analysis <- function(seurat_all, group_levels, opt) {
    library(future.apply)   
    single_list <- future_lapply(group_levels, function(group_select){
        seurat <- subset(seurat_all, !!sym(opt$group_col) %in% group_select)
        seurat@active.ident <- seurat@meta.data[[opt$celltype]]
        
        data <- seurat@assays$RNA$data
        metadata <- seurat@meta.data[!is.na(seurat@meta.data[[opt$celltype]]),]
        data <- data[,rownames(metadata)]
        
        cellchat <- createCellChat(data, meta = metadata, group.by = opt$celltype)
        cellchat <- addMeta(cellchat, meta = metadata)
        
        if (any(grepl("\\bACTB\\b", rownames(seurat)))) {
            CellChatDB <- CellChatDB.human
        } else if (any(grepl("\\bActb\\b", rownames(seurat)))) {
            CellChatDB <- CellChatDB.mouse
        }
        cellchat@DB <- CellChatDB
        
        cellchat <- subsetData(cellchat)
        #future::plan("multisession", workers = ncor)
        cellchat <- identifyOverExpressedGenes(cellchat)
        cellchat <- identifyOverExpressedInteractions(cellchat)
        
        if (opt$smooth_use) {
            if (any(grepl("\\bACTB\\b", rownames(seurat)))) {
                cellchat <- smoothData(cellchat, adj = PPI.human)
            }else if(any(grepl("\\bActb\\b", rownames(seurat)))) {
                cellchat <- smoothData(cellchat, adj = PPI.mouse)
            }
        }
        
        cellchat <- computeCommunProb(cellchat, population.size = FALSE)
        cellchat <- filterCommunication(cellchat, min.cells = 10)
        cellchat <- computeCommunProbPathway(cellchat,thresh = opt$plot_Compute_Aggregate_Analysis_Pvalue)
        cellchat <- aggregateNet(cellchat,thresh = opt$plot_Compute_Aggregate_Analysis_Pvalue)
        
        return(cellchat)
    }, future.seed = TRUE)
    
    names(single_list) <- group_levels
    return(single_list)
}

# Single模式出图函数
single_cellchat_plot <- function( cellchat, opt, sub_output_dir ){
    # ANSI 颜色代码
    red <- "\033[31m"
    green <- "\033[32m"
    yellow <- "\033[33m"
    blue <- "\033[34m"
    magenta <- "\033[35m"
    cyan <- "\033[36m"
    reset <- "\033[0m"
    
    # 根据信号通路显著性和基因表达过滤信号通路
    net.pct <- netMappingDEG(cellchat, features.name = "features",variable.all = TRUE,thresh = opt$plot_Net_netP_Pvalue)
    # 将其中的NA替换为值
    pvalue_cols <- grep("\\.pvalues$", colnames(net.pct), value = TRUE)
    logFC_cols  <- grep("\\.logFC$",   colnames(net.pct), value = TRUE)
    pct_cols    <- grep("\\.pct\\.",   colnames(net.pct), value = TRUE)
    # 替换NA
    net.pct[pvalue_cols] <- lapply(net.pct[pvalue_cols],\(x) replace(x, is.na(x), 1))
    net.pct[logFC_cols] <- lapply(net.pct[logFC_cols], \(x) replace(x, is.na(x), 0))
    net.pct[pct_cols] <- lapply(net.pct[pct_cols],\(x) replace(x, is.na(x), 0))
    
    # net
    df_net <- subsetCommunication(cellchat,
                                  net = net.pct,
                                  sources.use = NULL,
                                  slot.name = "net",
                                  thresh = opt$plot_Net_netP_Pvalue,   
                                  ligand.pvalues = NULL,
                                  ligand.logFC = NULL,
                                  ligand.pct.1 = opt$plot_geneExp_minpct*100,
                                  ligand.pct.2 = NULL,
                                  receptor.pvalues = NULL,
                                  receptor.logFC = NULL,
                                  receptor.pct.1 = opt$plot_geneExp_minpct*100,
                                  receptor.pct.2 = NULL
                                 )
    df_net <- df_net %>% dplyr::filter(prob > quantile(df_net$prob,probs = opt$plot_prob_Quantile,na.rm = TRUE))
    pairLR_use <- data.frame(interaction_name = unique(df_net$interaction_name))
    
    ### plot
    if(opt$plot_allcirc_allheat_allbub_singcirc_singbub[1]){
        cat("\n\t",yellow," Plot Single:  allcir",reset,"\n")
        flush.console() 
        # 整体圆环图
        groupSize <- as.numeric(table(cellchat@idents)) %>% log()
        lapply(opt$plot_count_weight,function(type){
            
            # 原始完整矩阵
            mat <- cellchat@net[[type]]
            write.csv(mat,file.path(opt$output_plot_dir,sub_output_dir, paste0("01_CellChat_all_", type, "_full.csv")))
            
            lapply(opt$plot_top_num,function(i){
                if (type == "count"){title_name = "Number of interactions"}else{title_name = "Interaction weights/strength"}

                # 保存数据，与 CellChat netVisual_circle 的 top 过滤逻辑一致
                mat_top <- mat
                mat_top[is.na(mat_top)] <- 0
                threshold <- stats::quantile(mat_top,probs = 1 - i,na.rm = TRUE)
                mat_top[mat_top < threshold] <- 0
                # 保存矩阵格式
                write.csv(mat_top,file.path(opt$output_plot_dir,sub_output_dir,paste0("01_CellChat_all_", type, "_", i, ".csv")))
                
                p1 <- function(){
                    par(bg = "white",mar = c(1, 1, 1, 1),cex = 1.5)
                    netVisual_circle(cellchat@net[[type]],
                                     vertex.weight = groupSize,
                                     weight.scale = T, 
                                     label.edge= F,
                                     title.name = paste0(title_name,"  (",i,")"),
                                     margin = c(0.3, 0.1, 0.4, 0.1),
                                     arrow.size = 0.8,
                                     vertex.label.cex = 0.8,
                                     top = i)
                }
                save_figure(obj = p1,
                            filename =  file.path(opt$output_plot_dir,sub_output_dir, paste0("01_CellChat_all_",type,"_",i)),
                            width = 6.5,
                            height = 6.5,
                            res = 600,
                            formats = opt$formats
                           )
            })
            
        })
    } 
    
    if( opt$plot_allcirc_allheat_allbub_singcirc_singbub[2] ){
        cat("\n\t",yellow," Plot Single:  allheat",reset,"\n")
        flush.console() 
        # 整体热图
        lapply(opt$plot_count_weight,function(type){
            p1 <- netVisual_heatmap(cellchat,
                                    measure = type,
                                    font.size = 12,
                                    font.size.title = 13
                                   )
            #par(mar = c(2, 1, 1, 1),cex = 1.5)
            save_figure(obj = p1,
                        filename =  file.path(opt$output_plot_dir,sub_output_dir,paste0("03_CellChat_all_Heatmap_",type)),
                        width = 8.5,
                        height = 8,
                        res = 600,
                        formats = opt$formats
                       )
        })
    }



    if( opt$plot_allcirc_allheat_allbub_singcirc_singbub[3] ){
        cat("\n\t",yellow," Plot Single:  allbub",reset,"\n")
        flush.console() 
        # 所有细胞类型气泡图
        if(nrow(pairLR_use) != 0){
            tryCatch({
                p1 <- netVisual_bubble(cellchat, remove.isolate = FALSE,pairLR.use = pairLR_use,angle.x = 45)+
                    theme(plot.margin = margin(t = 5,r = 5,b = 5,l = 50,unit = "pt"))                
                width = 4 + length(unique(as.character(p1@data$source.target))) * 0.2
                height = 5 + length(unique(as.character(p1@data$interaction_name_2))) * 0.1

                # 保存对应的气泡图数据
                write.csv(p1@data,file = file.path(opt$output_plot_dir,sub_output_dir,"04_CellChat_all_BublePlot.csv"))
                
                save_figure(obj = p1,
                            filename = file.path(opt$output_plot_dir, sub_output_dir, paste0("04_CellChat_all_BublePlot")),
                            width = min(50,width),
                            height = min(50,height),
                            res = 300,
                            formats = opt$formats
                           )     
            }, error = function(e) {NULL})
        
        }else{
            cat("\n\t",red,"In Single mode, no signaling pathways meet the current threshold criteria.",reset,"\n")
        }

    }
    
    if( opt$plot_allcirc_allheat_allbub_singcirc_singbub[4] ){
        cat("\n\t",yellow," Plot Single:  singcirc",reset,"\n")
        flush.console() 
        # 单个 循环图
        if( !is.null(opt$plot_source_target) ){
            if( "source" %in% opt$plot_source_target ){
                    # 信号源（配体）
                    groupSize <- as.numeric(table(cellchat@idents)) %>% log()
                    for(type in opt$plot_count_weight){
                        if (type == "count"){title_name = "Number of interactions"}else{title_name = "Interaction weights/strength"}
                        if(is.null(opt$key_cell)){
                            # 全部运行一遍
                            lapply(rownames(cellchat@net[[type]]),function(cellselect){
                                p1 <- function(){
                                    par(bg = "white",mar = c(1, 1, 1, 1),cex = 1.4)
                                    netVisual_circle(cellchat@net[[type]],
                                                     vertex.weight = groupSize,
                                                     weight.scale = T,
                                                     sources.use = cellselect,
                                                     edge.weight.max = max(cellchat@net[[type]]),
                                                     arrow.size = 1,
                                                     vertex.label.cex = 0.8,
                                                     margin = c(0.3, 0.1, 0.4, 0.1),
                                                     title.name = paste0(cellselect,"  (",title_name,")")
                                                    )
                                }
                                save_figure(obj = p1,
                                            filename = file.path(opt$output_plot_dir,sub_output_dir, paste0("02_CellChat_Single_Source_",type,"_",cellselect)),
                                            width = 6.5,
                                            height = 6.5,
                                            res = 600,
                                            formats = opt$formats
                                           )
                            })
                        }else{
                            if(all(opt$key_cell %in% rownames(cellchat@net[[type]]) )){
                                lapply(opt$key_cell,function(cellselect){
                                    # 只运行指定的
                                    p1 <- function(){
                                        par(bg = "white",mar = c(1, 1, 1, 1),cex = 1.4)    
                                        netVisual_circle(cellchat@net[[type]],
                                                         vertex.weight = groupSize,
                                                         weight.scale = T,
                                                         sources.use = cellselect,
                                                         edge.weight.max = max(cellchat@net[[type]]),
                                                         arrow.size = 1,
                                                         vertex.label.cex = 0.8,
                                                         margin = c(0.3, 0.1, 0.4, 0.1),
                                                         title.name = paste0(cellselect,"  (",title_name,")")
                                                        )
                                    }
                                    save_figure(obj = p1,
                                                filename = file.path(opt$output_plot_dir,sub_output_dir, paste0("02_CellChat_Single_Source_",type,"_",cellselect)),
                                                width = 6.5,
                                                height = 6.5,
                                                res = 600,
                                                formats = opt$formats
                                               )                            
                                })                             
                            }else{
                                cat("\n",blue,setdiff(opt$key_cell,rownames(cellchat@net[[type]])),"Not in Celltype",reset,"\n")
                                quit(save = "no")
                            }

                        }
                    }
            }
            if(  "target" %in% opt$plot_source_target  ){
                # 接收者（受体）
                groupSize <- as.numeric(table(cellchat@idents)) %>% log()
                for(type in opt$plot_count_weight){
                    if (type == "count"){title_name = "Number of interactions"}else{title_name = "Interaction weights/strength"}
                        if( is.null(opt$key_cell) ){
                            # 全部运行一遍
                            lapply(rownames(cellchat@net[[type]]),function(cell_name){
                                p1 <- function(){
                                    par(bg = "white",mar = c(1, 1, 1, 1),cex = 1.4)   
                                    netVisual_circle(cellchat@net[[type]],
                                                     vertex.weight = groupSize,
                                                     weight.scale = T, 
                                                     label.edge= F,
                                                     title.name = paste0(cell_name,"  (",title_name,")"),
                                                     margin = c(0.3, 0.1, 0.4, 0.1),
                                                     arrow.size = 1,
                                                     vertex.label.cex = 0.8,
                                                     targets.use = cell_name,
                                                     top = 1
                                                    )
                                }
                                save_figure(obj = p1,
                                            filename = file.path(opt$output_plot_dir,sub_output_dir, paste0("02_CellChat_Single_Target_",type,"_",cell_name)),
                                            width = 6.5,
                                            height = 6.5,
                                            res = 600,
                                            formats = opt$formats
                                           )  
                            })
                        }else{
                            # 只出指定细胞的
                            if(all(opt$key_cell %in% rownames(cellchat@net[[type]]) )){
                                lapply(opt$key_cell,function(cell_name){
                                    p1 <- function(){
                                        par(bg = "white",mar = c(1, 1, 1, 1),cex = 1.4)   
                                        netVisual_circle(cellchat@net[[type]],
                                                         vertex.weight = groupSize,
                                                         weight.scale = T, 
                                                         label.edge= F,
                                                         title.name = paste0(cell_name,"  (",title_name,")"),
                                                         margin = c(0.3, 0.1, 0.4, 0.1),
                                                         arrow.size = 1,
                                                         vertex.label.cex = 0.8,
                                                         targets.use = cell_name,
                                                         top = 1)
                                    }
                                    save_figure(obj = p1,
                                                filename = file.path(opt$output_plot_dir,sub_output_dir, paste0("02_CellChat_Single_Target_",type,"_",cell_name)),
                                                width = 6.5,
                                                height = 6.5,
                                                res = 600,
                                                formats = opt$formats
                                               ) 
                                })
                            }else{
                                cat("\n",blue,setdiff(opt$key_cell,rownames(cellchat@net[[type]])),"Not in Celltype",reset,"\n")
                                quit(save = "no")                      
                            }
                    
                        }
                }
            }
        }
    }

    if( opt$plot_allcirc_allheat_allbub_singcirc_singbub[5] ){
        cat("\n\t",yellow," Plot Single:  singbub",reset,"\n")
        flush.console() 
        # 单个气泡图
        # 指定细胞类型为源的气泡图
        
        if( is.na(opt$key_cell) %>% all() ){
            # 特定细胞
            sources_uses <- 1:ceiling(length(levels(cellchat@idents)) * 0.2)
        }else{
            sources_uses <- match(opt$key_cell,levels(cellchat@idents))
        }


        # 拆分成数字向量
        sources_use_list <- as.list(sources_uses)
        sources_use_list[[length(sources_use_list) + 1]] <- sources_uses
        
        lapply(sources_use_list,function(sources_use){
            # 
            df_net <- subsetCommunication(cellchat,
                                          net = net.pct,
                                          sources.use = sources_use,
                                          slot.name = "net",
                                          thresh = opt$plot_Net_netP_Pvalue,   
                                          ligand.pvalues = NULL,
                                          ligand.logFC = NULL,
                                          ligand.pct.1 = opt$plot_geneExp_minpct*100,
                                          ligand.pct.2 = NULL,
                                          receptor.pvalues = NULL,
                                          receptor.logFC = NULL,
                                          receptor.pct.1 = opt$plot_geneExp_minpct*100,
                                          receptor.pct.2 = NULL
                                         )
            df_net <- df_net %>% dplyr::filter(prob > quantile(df_net$prob,probs = opt$plot_prob_Quantile,na.rm = TRUE))
            pairLR_use <- data.frame(interaction_name = unique(df_net$interaction_name))
            
            # Plot
            if(nrow(pairLR_use) != 0){
                tryCatch({
                    p1 <- netVisual_bubble(cellchat, sources.use = sources_use, remove.isolate = FALSE,pairLR.use = pairLR_use,angle.x = 45)+
                        theme(plot.margin = margin(t = 5,r = 5,b = 5,l = 50,unit = "pt"))

                    # 保存气泡图数据
                    write.csv(p1@data,file = file.path(opt$output_plot_dir,sub_output_dir,paste0("04_CellChat_",paste(levels(cellchat@idents)[sources_use],collapse = "_"),"_Source_BublePlot.csv")))
                    
                    width = 4 + length(unique(as.character(p1@data$source.target))) * 0.2
                    height = 5 + length(unique(as.character(p1@data$interaction_name_2))) * 0.1
                    save_figure(obj = p1,
                                filename = file.path(opt$output_plot_dir,sub_output_dir,paste0("04_CellChat_",paste(levels(cellchat@idents)[sources_use],collapse = "_"),"_Source_BublePlot")),
                                width = min(50,width),
                                height = min(50,height),
                                res = 300,
                                formats = opt$formats
                               )    
                }, error = function(e) {NULL})

            }else{
                cat("\n\t",red,"In Single mode, no signaling pathways meet the current threshold criteria.",reset,"\n")
                flush.console()
            }


            
            df_net <- subsetCommunication(cellchat,
                                          net = net.pct,
                                          targets.use = sources_use,
                                          slot.name = "net",
                                          thresh = opt$plot_Net_netP_Pvalue,   
                                          ligand.pvalues = NULL,
                                          ligand.logFC = NULL,
                                          ligand.pct.1 = opt$plot_geneExp_minpct*100,
                                          ligand.pct.2 = NULL,
                                          receptor.pvalues = NULL,
                                          receptor.logFC = NULL,
                                          receptor.pct.1 = opt$plot_geneExp_minpct*100,
                                          receptor.pct.2 = NULL
                                         )
            df_net <- df_net %>% dplyr::filter(prob > quantile(df_net$prob,probs = opt$plot_prob_Quantile,na.rm = TRUE))
            pairLR_use <- data.frame(interaction_name = unique(df_net$interaction_name))
            
            if(nrow(pairLR_use) != 0){
                tryCatch({
                    p1 <- netVisual_bubble(cellchat, targets.use = sources_use, remove.isolate = FALSE,pairLR.use = pairLR_use,angle.x = 45)+
                        theme(plot.margin = margin(t = 5,r = 5,b = 5,l = 50,unit = "pt")) 
                    
                    # 保存气泡图数据
                    write.csv(p1@data,file = file.path(opt$output_plot_dir,sub_output_dir,paste0("04_CellChat_",paste(levels(cellchat@idents)[sources_use],collapse = "_"),"_Target_BublePlot.csv")))
                              
                    width = 4 + length(unique(as.character(p1@data$source.target))) * 0.2
                    height = 5 + length(unique(as.character(p1@data$interaction_name_2))) * 0.1
                    save_figure(obj = p1,
                                filename = file.path(opt$output_plot_dir,sub_output_dir,paste0("04_CellChat_",paste(levels(cellchat@idents)[sources_use],collapse = "_"),"_Target_BublePlot")),
                                width = min(50,width),
                                height = min(50,height),
                                res = 300,
                                formats = opt$formats
                               )     
                }, error = function(e) {NULL})
            }else{
                cat("\n\t",red,"In Single mode, no signaling pathways meet the current threshold criteria.",reset,"\n")
                flush.console()
            }
        })

        if(length(opt$source_cell) != 0 & length(opt$target_cell) != 0){
            cat("\n\t",yellow," Plot Single:  Source & Target Bub",reset,"\n")
            flush.console() 
            # 设置了源和受细胞
            source_cells <- match(opt$source_cell,levels(cellchat@idents))
            target_cells <- match(opt$target_cell,levels(cellchat@idents))
            # 拆分成数字向量
            source_cells_list <- as.list(source_cells)
            source_cells_list[[length(source_cells_list) + 1]] <- source_cells
        
            target_cells_list <- as.list(target_cells)
            target_cells_list[[length(target_cells_list) + 1]] <- target_cells
            
            cell_pairs <- expand.grid(source = source_cells_list,target = target_cells_list,stringsAsFactors = FALSE)
            
            lapply(1:nrow(cell_pairs),function(i){
                source_cell_use <- cell_pairs$source[[i]]
                target_cell_use <- cell_pairs$target[[i]]

                df_net <- subsetCommunication(cellchat,
                                              net = net.pct,
                                              sources.use = source_cell_use,
                                              targets.use = target_cell_use,
                                              slot.name = "net",
                                              thresh = opt$plot_Net_netP_Pvalue,   
                                              ligand.pvalues = NULL,
                                              ligand.logFC = NULL,
                                              ligand.pct.1 = opt$plot_geneExp_minpct*100,
                                              ligand.pct.2 = NULL,
                                              receptor.pvalues = NULL,
                                              receptor.logFC = NULL,
                                              receptor.pct.1 = opt$plot_geneExp_minpct*100,
                                              receptor.pct.2 = NULL
                                             )
                df_net <- df_net %>% dplyr::filter(prob > quantile(df_net$prob,probs = opt$plot_prob_Quantile,na.rm = TRUE))
                pairLR_use <- data.frame(interaction_name = unique(df_net$interaction_name))
                if(nrow(pairLR_use) != 0){
                    tryCatch({
                        p1 <- netVisual_bubble(cellchat, sources.use = source_cell_use,targets.use = target_cell_use, remove.isolate = FALSE,pairLR.use = pairLR_use,angle.x = 45)+
                            theme(plot.margin = margin(t = 5,r = 5,b = 5,l = 50,unit = "pt"))
                        
                        # 保存气泡图数据
                        write.csv(p1@data,file = file.path(opt$output_plot_dir,sub_output_dir,
                                                         paste0("04_CellChat_",paste(levels(cellchat@idents)[source_cell_use],collapse = "_"),"_to_",paste(levels(cellchat@idents)[target_cell_use],collapse = "_"),"_BublePlot.csv")
                                                        ))
                                  
                        width = 4 + length(unique(as.character(p1@data$source.target))) * 0.2
                        height = 5 + length(unique(as.character(p1@data$interaction_name_2))) * 0.1
                        save_figure(obj = p1,
                                    filename = file.path(opt$output_plot_dir,sub_output_dir,
                                                         paste0("04_CellChat_",paste(levels(cellchat@idents)[source_cell_use],collapse = "_"),"_to_",paste(levels(cellchat@idents)[target_cell_use],collapse = "_"),"_BublePlot")
                                                        ),
                                    width = min(50,width),
                                    height = min(50,height),
                                    res = 300,
                                    formats = opt$formats
                                   )     
                    }, error = function(e) {NULL})
                                   
                }else{
                    cat("\n\t",red,"In Single mode, no signaling pathways meet the current threshold criteria.",reset,"\n")
                    flush.console()                    
                }
            })
        } 


    }

    
    if(opt$plot_pathway_hier_circ_heat_chord_bar_violin_singRoleheat[1]){
        cat("\n\t",yellow," Plot Single Pathway:  Hierarchy",reset,"\n")
        flush.console() 
        # 目标通路层次图
        if( is.na(opt$hierarchy_top) %>% all() ){
            # 特定细胞
            sources_use <- 1:ceiling(length(levels(cellchat@idents)) * 0.2)
        }else{
            sources_use <- match(opt$hierarchy_top,levels(cellchat@idents))
        }
        # 层次图
        lapply(cellchat@netP$pathways,function(pathway){
            old_par <- par(no.readonly = TRUE)
            on.exit(par(old_par), add = TRUE)

            # 提取当前通路的通讯概率矩阵
            mat <- cellchat@netP$prob[, , pathway]
            # 保存矩阵
            write.csv(mat,file = file.path(opt$output_plot_dir,sub_output_dir,"01_Signa_Pathway",paste0("01_CellChat_Hierarchy_Pathway_", pathway, ".csv")))            

            width_add = levels(cellchat@meta$Celltype) %>% nchar() %>% max()
            width = 10 + width_add*0.15
            
            p <- function(){
                par(bg = "white",mar = c(0.2, 2, 5, 2),cex = 1.9)   
                netVisual_aggregate(cellchat, 
                                    signaling = pathway,  
                                    vertex.receiver = sources_use,
                                    layout = "hierarchy",
                                    title.space = 4
                                    )
            }
            save_figure(obj = p,
                        filename = file.path(opt$output_plot_dir, sub_output_dir,"01_Signa_Pathway", paste0("01_CellChat_Hierarchy_Pathway_",pathway)),
                        width = width,
                        height = 6,
                        res = 600,
                        formats = opt$formats
                       )
        })
    }
                                   
    if(opt$plot_pathway_hier_circ_heat_chord_bar_violin_singRoleheat[2]){
        cat("\n\t",yellow," Plot Single Pathway:  Circle",reset,"\n")
        flush.console() 
        # 目标通路圆环图
        if( is.na(opt$key_cell) %>% all() ){
            # 特定细胞
            sources_use <- 1:ceiling(length(levels(cellchat@idents)) * 0.2)
        }else{
            sources_use <- match(opt$key_cell,levels(cellchat@idents))
        }
        # 圆环图
        lapply(cellchat@netP$pathways,function(pathway){
            p <- function(){
                par(bg = "white",mar = c(1, 1, 1, 1),cex = 1.5)   
                netVisual_aggregate(cellchat, 
                                    signaling = pathway,  
                                    vertex.receiver = sources_use,
                                    layout = "circle",
                                    vertex.label.cex = 0.8,
                                    margin = c(0.3, 0.1, 0.4, 0.1),
                                    title.space = 4
                                   )
            }
            save_figure(obj = p,
                        filename = file.path(opt$output_plot_dir, sub_output_dir,"01_Signa_Pathway", paste0("01_CellChat_Circle_Pathway_",pathway)),
                        width = 6.5,
                        height = 6.5,
                        res = 600,
                        formats = opt$formats
                       )
        })
        
    } 


    if(opt$plot_pathway_hier_circ_heat_chord_bar_violin_singRoleheat[3]){
        cat("\n\t",yellow," Plot Single Pathway:  Heatmap",reset,"\n")
        flush.console()
        # 目标通路热图
        lapply(cellchat@netP$pathways,function(pathway){
            p1 <- netVisual_heatmap(cellchat, 
                                    signaling = pathway, 
                                    color.heatmap = "Reds",
                                    font.size = 14,
                                    font.size.title = 16
                                   )
            save_figure(obj = p1,
                        filename = file.path(opt$output_plot_dir, sub_output_dir,"01_Signa_Pathway",paste0("02_CellChat_Heatmap_Pathway_",pathway)),
                        width = 6.5,
                        height = 6.5,
                        res = 600,
                        formats = opt$formats
                       )
        })
    }

    if(opt$plot_pathway_hier_circ_heat_chord_bar_violin_singRoleheat[4]){
        cat("\n\t",yellow," Plot Single Pathway:  Chord",reset,"\n")
        flush.console()
        # 目标通路和和弦图
        groupSize <- as.numeric(table(cellchat@idents)) %>% log()
        old_par <- par(no.readonly = TRUE)
        on.exit(par(old_par), add = TRUE)
        lapply(cellchat@netP$pathways,function(pathway){
            p1 <- function(){
                par(bg = "white",mar = c(1, 1, 1, 1),cex = 1.5)   
                netVisual_aggregate(cellchat, 
                                          signaling = pathway, 
                                          layout = "chord",
                                          vertex.receiver = vertex.receiver, 
                                          vertex.weight  = groupSize
                                         )
            }
            save_figure(obj = p1,
                        filename = file.path(opt$output_plot_dir,sub_output_dir,"01_Signa_Pathway",paste0("03_CellChat_Chord_Pathway_",pathway)),
                        width = 6.5,
                        height = 6.5,
                        res = 600,
                        formats = opt$formats
                       )      
        })  
    }

    if( opt$plot_pathway_hier_circ_heat_chord_bar_violin_singRoleheat[5] ){
        cat("\n\t",yellow," Plot Single Pathway:  Bar",reset,"\n")
        flush.console()
        # 每个信号通路下的配受体对贡献
        par(mar = c(1, 1, 1, 1))
        lapply(cellchat@netP$pathways,function(pathway){
            p1 <- netAnalysis_contribution(cellchat, signaling = pathway,font.size = 11,font.size.title = 14)
            save_figure(obj = p1,
                        filename = file.path(opt$output_plot_dir,sub_output_dir,"01_Signa_Pathway",paste0("05_CellChat_contribution_Pathway_",pathway)),
                        width = 6.5,
                        height = 6.5,
                        res = 600,
                        formats = opt$formats
                       ) 
        })
    }
                                   
    if( opt$plot_pathway_hier_circ_heat_chord_bar_violin_singRoleheat[6] ){
        cat("\n\t",yellow," Plot Single Pathway:  Violin",reset,"\n")
        flush.console()
        # 信号通路小提琴图
        lapply(cellchat@netP$pathways,function(pathway){
            pathway_interaction <- cellchat@DB$interaction[cellchat@DB$interaction$pathway_name == pathway, c("ligand", "receptor")]
            genes <- unique(c(pathway_interaction$ligand,pathway_interaction$receptor))
            
            p1 <- plotGeneExpression(cellchat, signaling = pathway, enriched.only = TRUE, type = "violin")
            save_figure(obj = p1,
                        filename = file.path(opt$output_plot_dir,sub_output_dir,"01_Signa_Pathway",paste0("09_CellChat_Violin_Sig_Expression_Pathway_",pathway)),
                        width = length(unique(cellchat@meta[[opt$celltype]])) * 0.5,
                        height = floor(sqrt(length(genes))) * 0.5 + 4,
                        res = 600,
                        formats = opt$formats
                       )     
        })
    }
                                   
    # 信号的主要发送和接收
    if( opt$plot_pathway_hier_circ_heat_chord_bar_violin_singRoleheat[7] ){
        cat("\n\t",yellow," Plot Single Pathway:  SignalingRole",reset,"\n")
        flush.console()
        # 整个的信号模式
        cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular comm
        lapply(cellchat@netP$pathways,function(pathway){
            p1 <- function(){
                netAnalysis_signalingRole_network(cellchat,signaling = pathway, width = length(unique(cellchat@meta[[opt$celltype]])), height = 5, font.size = 16,font.size.title = 18)
            }
            save_figure(obj = p1,
                        filename = file.path(opt$output_plot_dir,sub_output_dir,"01_Signa_Pathway",paste0("10_CellChat_SignalingRole_Pathway_",pathway)),
                        width = length(unique(cellchat@meta[[opt$celltype]]))*0.6,
                        height = 5,
                        res = 600,
                        formats = opt$formats
                       )   
        })
    }  

    if( opt$plot_allsignalingRole_scatt_outheat_inheat[1] ){
        cat("\n\t",yellow," Plot Single:  SignalingRole Scatter",reset,"\n")
        flush.console()
        
        cellchat <- netAnalysis_computeCentrality(cellchat)
	# 全局Scatter
        gg1 <- netAnalysis_signalingRole_scatter(cellchat,font.size = 18,label.size = 5)
        # 保存全局散点图数据
        write.csv(gg1$data,file = file.path(opt$output_plot_dir,sub_output_dir,"10_CellChat_SignalingRole_Scatter.csv"))
        save_figure(obj = gg1,
                    filename = file.path(opt$output_plot_dir,sub_output_dir,paste0("10_CellChat_SignalingRole_Scatter")),
                    width = 6,
                    height = 5.5,
                    res = 600,
                    formats = opt$formats
                   )  
        
    }
                                   
    if( opt$plot_allsignalingRole_scatt_outheat_inheat[2] ){
        cat("\n\t",yellow," Plot Single:  SignalingRole Heatmap Outgoing",reset,"\n")
        flush.console()
        cellchat <- netAnalysis_computeCentrality(cellchat)
        # 全局Scatter Out
        ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
        ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing",width = 4 + ncol(ht1@matrix)*0.2,height = 4 + nrow(ht1@matrix)*0.25)
        # 保存热图对应矩阵
        write.csv(ht1@matrix,file = file.path(opt$output_plot_dir,sub_output_dir,"11_CellChat_SignalingRole_Heatmap_outgoing.csv"))
        save_figure(obj = ht1,
                    filename = file.path(opt$output_plot_dir,sub_output_dir,paste0("11_CellChat_SignalingRole_Heatmap_outgoing")),
                    width = 4 + ncol(ht1@matrix)*0.15,
                    height = 4 + nrow(ht1@matrix)*0.1,
                    res = 600,
                    formats = opt$formats
                   ) 
        
    }
                                   
    if( opt$plot_allsignalingRole_scatt_outheat_inheat[3] ){
        cat("\n\t",yellow," Plot Single:  SignalingRole Heatmap Incoming",reset,"\n")
        flush.console()
        cellchat <- netAnalysis_computeCentrality(cellchat)
        # 全局Scatter In
        ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
        ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",width = 4 + ncol(ht2@matrix)*0.2,height = 4 + nrow(ht2@matrix)*0.25)
        # 保存热图对应矩阵
        write.csv(ht2@matrix,file = file.path(opt$output_plot_dir,sub_output_dir,"11_CellChat_SignalingRole_Heatmap_incoming.csv"))
        save_figure(obj = ht2,
                    filename = file.path(opt$output_plot_dir,sub_output_dir,paste0("11_CellChat_SignalingRole_Heatmap_incoming")),
                    width = 4 + ncol(ht2@matrix)*0.15,
                    height = 4 + nrow(ht2@matrix)*0.1,
                    res = 600,
                    formats = opt$formats
                   ) 
    }

    if( opt$plot_allcommunication_outgo_income[1] ){
        cat("\n\t",yellow," Plot Single: Outgoing Global Communication Patterns",reset,"\n")
        flush.console()
        # 发出的全局通讯模式    
        pattern = "outgoing"
        # 推断模式数量  
        plot_k <- selectK(cellchat, pattern = pattern)
        save_figure(obj = plot_k,
                    filename = file.path(opt$output_plot_dir,sub_output_dir,paste0("12_CellChat_global_communication_patterns_selectK_",pattern)),
                    width = 8,
                    height = 4,
                    res = 600,
                    formats = opt$formats
                   )
        
        grid::grid.newpage()
        cellchat <- identifyCommunicationPatterns(cellchat,pattern = pattern, k = 3, height = length(cellchat@netP$pathways) * 0.28,width = length(unique(cellchat@meta[[opt$celltype]])))
        
        # 模式图
        plot_grob <- grid::grid.grab()
        save_figure(obj = plot_grob,
                    filename = file.path(opt$output_plot_dir,sub_output_dir,paste0("12_CellChat_global_communication_patterns_Show_",pattern)),
                    width = 2*(4 + length(unique(cellchat@meta[[opt$celltype]]))*0.3),
                    height = length(cellchat@netP$pathways) * 0.25,
                    res = 600,
                    formats = opt$formats
                   )
        
        # 河流图
        plot_river <- netAnalysis_river(cellchat, pattern = pattern)
        save_figure(obj = plot_river,
                    filename = file.path(opt$output_plot_dir,sub_output_dir,paste0("12_CellChat_global_communication_patterns_River_",pattern)),
                    width = 6,
                    height = 4 + length(cellchat@netP$pathways)*0.15,
                    res = 600,
                    formats = opt$formats
                   )
        
        # dot plot
        pd <- netAnalysis_dot(cellchat, pattern = pattern)
        save_figure(obj = pd,
                    filename = file.path(opt$output_plot_dir,sub_output_dir,paste0("12_CellChat_global_communication_patterns_DotPlot_",pattern)),
                    width = 4 + length(cellchat@netP$pathways)*0.15,
                    height = 4 + length(unique(cellchat@meta[[opt$celltype]]))*0.1,
                    res = 600,
                    formats = opt$formats
                   )
    }
                                   
   if( opt$plot_allcommunication_outgo_income[2] ){
        cat("\n\t",yellow," Plot Single: Incoming Global Communication Patterns",reset,"\n")
        flush.console()
       # 接受的全局通讯模式  

        pattern = "incoming"
        # 推断模式数量   
        plot_k <- selectK(cellchat, pattern = pattern)
        save_figure(obj = plot_k,
                    filename = file.path(opt$output_plot_dir,sub_output_dir,paste0("13_CellChat_global_communication_patterns_selectK_",pattern)),
                    width = 8,
                    height = 4,
                    res = 600,
                    formats = opt$formats
                   )        
        
        grid::grid.newpage()
        cellchat <- identifyCommunicationPatterns(cellchat,pattern = pattern, k = 3,height = length(cellchat@netP$pathways) * 0.28, width = length(unique(cellchat@meta[[opt$celltype]])))
        plot_grob <- grid::grid.grab()
        
        # 模式图
        save_figure(obj = plot_grob,
                    filename = file.path(opt$output_plot_dir,sub_output_dir,paste0("13_CellChat_global_communication_patterns_Show_",pattern)),
                    width = 2*(4 + length(unique(cellchat@meta[[opt$celltype]]))*0.2),
                    height = length(cellchat@netP$pathways) * 0.25,
                    res = 600,
                    formats = opt$formats
                   )
        
        # 河流图
        plot_river <- netAnalysis_river(cellchat, pattern = pattern)
        save_figure(obj = plot_river,
                    filename = file.path(opt$output_plot_dir,sub_output_dir,paste0("13_CellChat_global_communication_patterns_River_",pattern)),
                    width = 6,
                    height = 4 + length(cellchat@netP$pathways)*0.15,
                    res = 600,
                    formats = opt$formats
                   )
        
        # dot plot
        pd <- netAnalysis_dot(cellchat, pattern = pattern)
        save_figure(obj = pd,
                    filename = file.path(opt$output_plot_dir,sub_output_dir,paste0("13_CellChat_global_communication_patterns_DotPlot_",pattern)),
                    width = 4 + length(cellchat@netP$pathways)*0.15,
                    height = 4 + length(unique(cellchat@meta[[opt$celltype]]))*0.1,
                    res = 600,
                    formats = opt$formats
                   )
   }

   # Other
                                                                                       
} 

# Compare模式出图函数
compare_cellchat_plot <- function(cellchat, opt,sub_output_dir,group_levels ){
    # ANSI 颜色代码
    red <- "\033[31m"
    green <- "\033[32m"
    yellow <- "\033[33m"
    blue <- "\033[34m"
    magenta <- "\033[35m"
    cyan <- "\033[36m"
    reset <- "\033[0m"
    
    # 根据信号通路显著性和基因表达过滤信号通路
    net.pct <- netMappingDEG(object = cellchat,features.name = paste0(group_levels[2],"_vs_",group_levels[1]),variable.all = TRUE,thresh = opt$plot_Net_netP_Pvalue)
    # 将NA替换为值
    pvalue_cols <- grep("\\.pvalues$", colnames(net.pct), value = TRUE)
    logFC_cols  <- grep("\\.logFC$",   colnames(net.pct), value = TRUE)
    pct_cols    <- grep("\\.pct\\.",   colnames(net.pct), value = TRUE)
    
    # 替换NA
    net.pct[pvalue_cols] <- lapply(net.pct[pvalue_cols],\(x) replace(x, is.na(x), 1))
    net.pct[logFC_cols] <- lapply(net.pct[logFC_cols], \(x) replace(x, is.na(x), 0))
    net.pct[pct_cols] <- lapply(net.pct[pct_cols],\(x) replace(x, is.na(x), 0))
    # net
    df_net <- subsetCommunication(cellchat,
                                  net = net.pct,
                                  sources.use = NULL,
                                  slot.name = "net",
                                  thresh = opt$plot_Net_netP_Pvalue,   
                                  ligand.pvalues = opt$plot_geneExp_Pvalue,
                                  ligand.logFC = NULL,
                                  ligand.pct.1 = opt$plot_geneExp_minpct*100,
                                  ligand.pct.2 = opt$plot_geneExp_minpct*100,
                                  receptor.pvalues = opt$plot_geneExp_Pvalue,
                                  receptor.logFC = NULL,
                                  receptor.pct.1 = opt$plot_geneExp_minpct*100,
                                  receptor.pct.2 = opt$plot_geneExp_minpct*100
                                 )
    df_net <- df_net %>% dplyr::filter(abs(ligand.logFC) >= opt$plot_geneExp_logFC,abs(receptor.logFC) >= opt$plot_geneExp_logFC)
    df_net <- df_net %>% dplyr::filter(prob > quantile(df_net$prob,probs = opt$plot_prob_Quantile,na.rm = TRUE))
    pairLR_use <- data.frame(interaction_name = unique(df_net$interaction_name))

    
    # 整体信号柱子图
    gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), color.use = c("#00BFC4","#F8766D"),size.text = 14)
    gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight", color.use = c("#00BFC4","#F8766D"),size.text = 14)

    save_figure(obj = gg1,
            filename = file.path(opt$output_plot_dir,sub_output_dir,paste0("01_CellChat_Bar_interaction_",opt$group_col)),
            width = 3,
            height = 4,
            res = 600,
            formats = opt$formats
           )
    save_figure(obj = gg2,
            filename = file.path(opt$output_plot_dir,sub_output_dir,paste0("01_CellChat_Bar_weight_",opt$group_col)),
            width = 3,
            height = 4,
            res = 600,
            formats = opt$formats
           )
    
    p <- gg1 + gg2
    save_figure(obj = p,
                filename = file.path(opt$output_plot_dir,sub_output_dir,paste0("01_CellChat_Bar_interaction_weight_",opt$group_col)),
                width = 6,
                height = 4,
                res = 600,
                formats = opt$formats
               )
    


    # 配体受体整体堆叠柱形图

    df_rank <- rankNet(cellchat,mode = "comparison",measure = "weight",sources.use = NULL,targets.use = NULL,stacked = TRUE,do.stat = TRUE,cutoff.pvalue = opt$plot_Net_netP_Pvalue,thresh = opt$plot_Net_netP_Pvalue,return.data = TRUE)
    # 保存原始信息流数据
    write.csv(df_rank$signaling.contribution,file = file.path(opt$output_plot_dir,sub_output_dir,"03_CellChat_Bar_weight.csv"))
    
    gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = TRUE, do.stat = TRUE,cutoff.pvalue = opt$plot_Net_netP_Pvalue,thresh = opt$plot_Net_netP_Pvalue)
    save_figure(obj = gg1,
                filename = file.path(opt$output_plot_dir,sub_output_dir,paste0("03_CellChat_Bar_weight_Stack")),
                width = 6,
                height = 4 + length(unique(as.character(gg1@data$name))) * 0.1,
                res = 600,
                formats = opt$formats
               )
    
    gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = FALSE, do.stat = TRUE,cutoff.pvalue = opt$plot_Net_netP_Pvalue,thresh = opt$plot_Net_netP_Pvalue)
    save_figure(obj = gg2,
                filename = file.path(opt$output_plot_dir,sub_output_dir,paste0("03_CellChat_Bar_weight_NoStack")),
                width = 6,
                height = 4 + length(unique(as.character(gg2@data$name))) * 0.1,
                res = 600,
                formats = opt$formats
               ) 

    if(length(opt$source_cell) != 0 & length(opt$target_cell) != 0 ){
        cat("\n\t",yellow," Plot Compare:  Source & Target Stack",reset,"\n")
        flush.console() 
        # 设置了源和受细胞
        source_cells <- match(opt$source_cell,levels(cellchat@idents$joint))
        target_cells <- match(opt$target_cell,levels(cellchat@idents$joint))
        # 拆分成数字向量
        source_cells_list <- as.list(source_cells)
        source_cells_list[[length(source_cells_list) + 1]] <- source_cells
    
        target_cells_list <- as.list(target_cells)
        target_cells_list[[length(target_cells_list) + 1]] <- target_cells
        
        cell_pairs <- expand.grid(source = source_cells_list,target = target_cells_list,stringsAsFactors = FALSE)
        lapply(1:nrow(cell_pairs),function(i){
            source_cell_use <- cell_pairs$source[[i]]
            target_cell_use <- cell_pairs$target[[i]]
            tryCatch({
                gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = source_cell_use, targets.use = target_cell_use,
                               stacked = T, do.stat = TRUE,cutoff.pvalue = opt$plot_Net_netP_Pvalue,thresh = opt$plot_Net_netP_Pvalue
                              )
                save_figure(obj = gg1,
                            filename = file.path(opt$output_plot_dir,sub_output_dir,paste0("03_CellChat_Bar_",
                                                                                           paste(levels(cellchat@idents$joint)[source_cell_use],collapse = "_"),
                                                                                           "_to_",
                                                                                           paste(levels(cellchat@idents$joint)[target_cell_use],collapse = "_"),
                                                                                           "_weight_Stack")
                                                ),
                            width = 6,
                            height = 4 + length(unique(as.character(gg1@data$name))) * 0.1,
                            res = 600,
                            formats = opt$formats
                           )
            }, error = function(e) {NULL})
            
            tryCatch({
                gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = source_cell_use, targets.use = target_cell_use,
                               stacked = F, do.stat = TRUE,cutoff.pvalue = opt$plot_Net_netP_Pvalue,thresh = opt$plot_Net_netP_Pvalue
                              )
                save_figure(obj = gg2,
                            filename = file.path(opt$output_plot_dir,sub_output_dir,paste0("03_CellChat_Bar_",
                                                                                           paste(levels(cellchat@idents$joint)[source_cell_use],collapse = "_"),
                                                                                           "_to_",
                                                                                           paste(levels(cellchat@idents$joint)[target_cell_use],collapse = "_"),
                                                                                           "_weight_NoStack")
                                                ),
                            width = 6,
                            height = 4 + length(unique(as.character(gg2@data$name))) * 0.1,
                            res = 600,
                            formats = opt$formats
                           )
            }, error = function(e) {NULL})
        })
        
    }
    
    if( opt$plot_allcirc_allheat_allbub_singcirc_singbub[1] ){
        cat("\n\t",yellow," Plot Compare: allCircle",reset,"\n")
        flush.console()

        # ==================== 保存差异通讯矩阵 ====================
        # Control vs NAPS，与 comparison = c(1,2) 一致
        group1 <- names(cellchat@net)[1]
        group2 <- names(cellchat@net)[2]
        
        for (type in opt$plot_count_weight) {
            # 提取两组通讯矩阵
            mat1 <- cellchat@net[[group1]][[type]]
            mat2 <- cellchat@net[[group2]][[type]]
            
            # 计算差异：NAPS - Control
            mat_diff <- mat2 - mat1
            
            # 保存CSV
            write.csv(mat_diff,file = file.path(opt$output_plot_dir, sub_output_dir,paste0("02_CellChat_Circle_", type, "_",opt$group_col, "_",group2, "_vs_", group1,"_diff.csv")))
        }
        
        lapply(opt$plot_top_num,function(i){
            if( "count" %in% opt$plot_count_weight ){
                p <- function(){
                    par(mar = rep(1, 4),cex = 1.5)
                    netVisual_diffInteraction(cellchat, weight.scale = T, title.name = "Differential number of interactions",margin = c(0.2,0.2,0.4,0.2), comparison = c(1,2),
                                    arrow.size = 0.8,top = i,vertex.label.cex = 0.8)
                }
                save_figure(obj = p,
                            filename = file.path(opt$output_plot_dir,sub_output_dir,paste0("02_CellChat_Circle_interaction_",opt$group_col,"_",group_levels[2],"_vs_",group_levels[1],"_top_",i)),
                            width = 6.5,
                            height = 6.5,
                            res = 600,
                            formats = opt$formats
                           ) 
                
            }
            if( "weight" %in% opt$plot_count_weight ){
                p <- function(){
                    par(mar = rep(1, 4),cex = 1.5)
                    netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight",title.name = "Differential interaction strength",margin = c(0.2,0.2,0.4,0.2),comparison = c(1,2),
                                    arrow.size = 0.8,top = i,vertex.label.cex = 0.8 )
                }
                save_figure(obj = p,
                            filename = file.path(opt$output_plot_dir,sub_output_dir,paste0("02_CellChat_Circle_weight_",opt$group_col,"_",group_levels[2],"_vs_",group_levels[1],"_top_",i)),
                            width = 6.5,
                            height = 6.5,
                            res = 600,
                            formats = opt$formats
                           ) 
                
            }
            if( all(c("count","weight") %in% opt$plot_count_weight) ){
                p <- function(){
                    par(mar = rep(1, 4),cex = 1.5)
                    cowplot::plot_grid(netVisual_diffInteraction(cellchat, weight.scale = T, title.name = "Differential number of interactions",margin = c(0.2,0.2,0.4,0.2), comparison = c(1,2),
                                                                       arrow.size = 0.8,top = i,vertex.label.cex = 0.8),
                                       netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight",title.name = "Differential interaction strength",margin = c(0.2,0.2,0.4,0.2),comparison = c(1,2),
                                                                       arrow.size = 0.8,top = i,vertex.label.cex = 0.8)
                                      )
                }
                save_figure(obj = p,
                            filename = file.path(opt$output_plot_dir,sub_output_dir,paste0("02_CellChat_Circle_interaction_weight_",opt$group_col,"_",group_levels[2],"_vs_",group_levels[1],"_top_",i)),
                            width = 6.5,
                            height = 6.5,
                            res = 600,
                            formats = opt$formats
                           ) 
                
            }
        })
 
    }

    if( opt$plot_allcirc_allheat_allbub_singcirc_singbub[2] ){
        cat("\n\t",yellow," Plot Compare: allHeatmap",reset,"\n")
        flush.console()
        # 整体热图
        gg1 <- netVisual_heatmap(cellchat,measure = "count",comparison = c(1,2),font.size = 14,font.size.title = 16)
        #> Do heatmap based on a merged object
        gg2 <- netVisual_heatmap(cellchat, measure = "weight",comparison = c(1,2),font.size = 14,font.size.title = 16)

        #> Do heatmap based on a merged object
        save_figure(obj = gg1,
                    filename = file.path(opt$output_plot_dir,sub_output_dir,paste0("02_CellChat_Heatmap_interaction_",opt$group_col,"_",group_levels[2],"_vs_",group_levels[1])),
                    width = 8.5,
                    height = 8,
                    res = 600,
                    formats = opt$formats
                   ) 
        save_figure(obj = gg2,
                    filename = file.path(opt$output_plot_dir,sub_output_dir,paste0("02_CellChat_Heatmap_weight_",opt$group_col,"_",group_levels[2],"_vs_",group_levels[1])), 
                    width = 8.5,
                    height = 8,
                    res = 600,
                    formats = opt$formats
                   ) 
    }
    
     if( opt$plot_allcirc_allheat_allbub_singcirc_singbub[3] ){
        cat("\n\t",yellow," Plot Compare: allBuble",reset,"\n")
        flush.console()
        if(nrow(pairLR_use) != 0){
            tryCatch({
                # ALL
                p1 <- netVisual_bubble(cellchat,  comparison = c(1, 2), angle.x = 45,pairLR.use = pairLR_use) +
                    theme(plot.margin = margin(t = 5,r = 5,b = 5,l = 50,unit = "pt"))      
                # 保存数据
                write.csv(p1@data,file = file.path(opt$output_plot_dir,sub_output_dir,"05_CellChat_Buble_L-R_pairs_all.csv"))
                
                width = 4 + length(unique(as.character(p1@data$source.target))) * 0.2
                height = 5 + length(unique(as.character(p1@data$interaction_name_2))) * 0.1         
                save_figure(obj = p1,
                            filename = file.path(opt$output_plot_dir,sub_output_dir,paste0("05_CellChat_Buble_L-R_pairs_all")), 
                            width = min(50,width),
                            height = min(50,height),
                            res = 300,
                            formats = opt$formats
                           )    
            }, error = function(e) {NULL})
        }else{
            cat("\n\t",red,"In Compare mode, no signaling pathways meet the current threshold criteria.",reset,"\n")
            flush.console()
        }
     }

    
    if( opt$plot_allcirc_allheat_allbub_singcirc_singbub[4] ){
        cat("\n\t",yellow," Plot Compare: singCircle",reset,"\n")
        flush.console()
        # 单个圆环图
        #celltype_colname <- colnames(cellchat@meta)[grep("^celltype$",tolower(colnames(cellchat@meta)))]
    	celltype_colname <- opt$celltype
        all_celltypes <- unique(as.character(cellchat@meta[[celltype_colname]]))
        
        if( "source" %in%  opt$plot_source_target ){
            # 单一 圆环图
            if(is.null(opt$key_cell)){
                # 全部出图
                
                lapply(all_celltypes,function(cell_select){
                    par(mar = c(2, 2, 2, 2))   
                    if(  "count" %in% opt$plot_count_weight ){
                        p <- function(){
                            par(mar = rep(1, 4),cex = 1.5)
                            netVisual_diffInteraction(cellchat, weight.scale = T,measure = "count", title.name = "Differential number of interactions",margin = c(0.2,0.2,0.4,0.2), comparison = c(1,2),
                                                      sources.use = cell_select,arrow.size = 0.8,vertex.label.cex = 0.8)
                        }
                        save_figure(obj = p,
                                    filename = file.path(opt$output_plot_dir,sub_output_dir,paste0("02_CellChat_Circle_interaction_",opt$group_col,"_",group_levels[2],"_vs_",group_levels[1],"_top_",1,"_source_",cell_select)), 
                                    width = 6.5,
                                    height = 6.5,
                                    res = 600,
                                    formats = opt$formats
                                   )          
                    }
                    
                    if(  "weight" %in% opt$plot_count_weight ){
                        p <- function(){
                            par(mar = rep(1, 4),cex = 1.5)
                            netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight",title.name = "Differential interaction strength",margin = c(0.2,0.2,0.4,0.2),comparison = c(1,2),
                                                      sources.use = cell_select,arrow.size = 0.8,vertex.label.cex = 0.8)
                        }
                        save_figure(obj = p,
                                    filename = file.path(opt$output_plot_dir,sub_output_dir,paste0("02_CellChat_Circle_weight_",opt$group_col,"_",group_levels[2],"_vs_",group_levels[1],"_top_",1,"_source_",cell_select)),
                                    width = 6.5,
                                    height = 6.5,
                                    res = 600,
                                    formats = opt$formats
                                   )                                                  
                        
                    }   
                    
                    if( all(c("count","weight") %in% opt$plot_count_weight) ){
                        p <- function(){
                            par(mar = rep(1, 4),cex = 1.5)
                            cowplot::plot_grid(netVisual_diffInteraction(cellchat, weight.scale = T,measure = "count", title.name = "Differential number of interactions",margin = c(0.2,0.2,0.4,0.2), comparison = c(1,2),
                                                                         sources.use = cell_select,arrow.size = 0.8,vertex.label.cex = 0.8),
                                               netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight",title.name = "Differential interaction strength",margin = c(0.2,0.2,0.4,0.2),comparison = c(1,2),
                                                                         sources.use = cell_select,arrow.size = 0.8,vertex.label.cex = 0.8)
                                              )
                        }
                        
                        save_figure(obj = p,
                                    filename = file.path(opt$output_plot_dir,sub_output_dir,paste0("02_CellChat_Circle_interaction_weight_",opt$group_col,"_",group_levels[2],"_vs_",group_levels[1],"_top_",1,"_source_",cell_select)),
                                    width = 6.5,
                                    height = 6.5,
                                    res = 600,
                                    formats = opt$formats
                                   )
                        
                    }
                })
    
            }else{
                # 只出 指定的
                if(all(opt$key_cell %in% all_celltypes )){
                    lapply(opt$key_cell,function(cell_select){
                        if(  "count" %in% opt$plot_count_weight ){
                            p <- function(){
                                par(mar = rep(1, 4),cex = 1.5)
                                netVisual_diffInteraction(cellchat, weight.scale = T,measure = "count", title.name = "Differential number of interactions",margin = c(0.2,0.2,0.4,0.2), comparison = c(1,2),
                                                          sources.use = cell_select,arrow.size = 0.8, vertex.label.cex = 0.8)
                            }
                            save_figure(obj = p,
                                        filename = file.path(opt$output_plot_dir,sub_output_dir,paste0("02_CellChat_Circle_interaction_",opt$group_col,"_",group_levels[2],"_vs_",group_levels[1],"_top_",1,"_source_",cell_select)),
                                        width = 6.5,
                                        height = 6.5,
                                        res = 600,
                                        formats = opt$formats
                                       )       
                        }
                        
                        if(  "weight" %in% opt$plot_count_weight ){
                            p <- function(){
                                par(mar = rep(1, 4),cex = 1.5)
                                netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight",title.name = "Differential interaction strength",margin = c(0.2,0.2,0.4,0.2),comparison = c(1,2),
                                                          sources.use = cell_select,arrow.size = 0.8,vertex.label.cex = 0.8)
                            }
                            save_figure(obj = p,
                                        filename = file.path(opt$output_plot_dir,sub_output_dir,paste0("02_CellChat_Circle_weight_",opt$group_col,"_",group_levels[2],"_vs_",group_levels[1],"_top_",1,"_source_",cell_select)),
                                        width = 6.5,
                                        height = 6.5,
                                        res = 600,
                                        formats = opt$formats
                                       )                                              
                            
                        }   
                        
                        if( all(c("count","weight") %in% opt$plot_count_weight) ){
                            p <- function(){
                                par(mar = rep(1, 4),cex = 1.5)
                                cowplot::plot_grid(netVisual_diffInteraction(cellchat, weight.scale = T,measure = "count", title.name = "Differential number of interactions",margin = c(0.2,0.2,0.4,0.2), comparison = c(1,2),
                                                                             sources.use = cell_select,arrow.size = 0.8,vertex.label.cex = 0.8 ),
                                                   netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight",title.name = "Differential interaction strength",margin = c(0.2,0.2,0.4,0.2),comparison = c(1,2),
                                                                             sources.use = cell_select,arrow.size = 0.8,vertex.label.cex = 0.8)
                                                  )
                            }
                            
                            save_figure(obj = p,
                                        filename = file.path(opt$output_plot_dir,sub_output_dir,paste0("02_CellChat_Circle_interaction_weight_",opt$group_col,"_",group_levels[2],"_vs_",group_levels[1],"_top_",1,"_source_",cell_select)),
                                        width = 6.5,
                                        height = 6.5,
                                        res = 600,
                                        formats = opt$formats
                                       )  
                            
                        }
                    })                 
                             
                }else{
                    cat("\n",blue,setdiff(opt$key_cell,all_celltypes),"Not in Celltype",reset,"\n")
                    quit(save = "no")
                }
            }
        }
        
        if(  "target" %in%  opt$plot_source_target ){
            if(is.null(opt$key_cell)){
                lapply(all_celltypes,function(cell_select){
                    
                    if(  "count" %in% opt$plot_count_weight ){
                        p <- function(){
                            par(mar = rep(1, 4),cex = 1.5)
                            netVisual_diffInteraction(cellchat, weight.scale = T,measure = "count", title.name = "Differential number of interactions",margin = c(0.2,0.2,0.4,0.2), comparison = c(1,2),
                                                      targets.use = cell_select,arrow.size = 0.8,vertex.label.cex = 0.8)
                        }
                        save_figure(obj = p,
                                    filename = file.path(opt$output_plot_dir,sub_output_dir,paste0("02_CellChat_Circle_interaction_",opt$group_col,"_",group_levels[2],"_vs_",group_levels[1],"_top_",1,"_targets_",cell_select)),
                                    width = 6.5,
                                    height = 6.5,
                                    res = 600,
                                    formats = opt$formats
                                   )  
                    }
        
                    if(  "weight" %in% opt$plot_count_weight ){
                        
                        p <- function(){
                            par(mar = rep(1, 4),cex = 1.5)
                            netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight",title.name = "Differential interaction strength",margin = c(0.2,0.2,0.4,0.2),comparison = c(1,2),
                                                      targets.use = cell_select,arrow.size = 0.8,vertex.label.cex = 0.8)
                        }
                        save_figure(obj = p,
                                    filename = file.path(opt$output_plot_dir,sub_output_dir,paste0("02_CellChat_Circle_weight_",opt$group_col,"_",group_levels[2],"_vs_",group_levels[1],"_top_",1,"_targets_",cell_select)),
                                    width = 6.5,
                                    height = 6.5,
                                    res = 600,
                                    formats = opt$formats
                                   )  
                    }
        
                    if( all(c("count","weight") %in% opt$plot_count_weight) ){       
                        p <- function(){
                            par(mar = rep(1, 4),cex = 1.5)
                            cowplot::plot_grid(netVisual_diffInteraction(cellchat, weight.scale = T,measure = "count", title.name = "Differential number of interactions",margin = c(0.2,0.2,0.4,0.2), comparison = c(1,2),
                                                                         targets.use = cell_select,arrow.size = 0.8,vertex.label.cex = 0.8),
                                               netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight",title.name = "Differential interaction strength",margin = c(0.2,0.2,0.4,0.2),comparison = c(1,2),
                                                                         targets.use = cell_select,arrow.size = 0.8,vertex.label.cex = 0.8)
                                              )
                        }
                        save_figure(obj = p,
                                    filename = file.path(opt$output_plot_dir,sub_output_dir,paste0("02_CellChat_Circle_interaction_weight_",opt$group_col,"_",group_levels[2],"_vs_",group_levels[1],"_top_",1,"_targets_",cell_select)),
                                    width = 6.5,
                                    height = 6.5,
                                    res = 600,
                                    formats = opt$formats
                                   )
                        
                    } 
                })
            }else{
                #
                 if(all(opt$key_cell %in% all_celltypes )){
                       lapply(opt$key_cell,function(cell_select){
                        
                        if(  "count" %in% opt$plot_count_weight ){
                            p <- function(){
                                par(mar = rep(1, 4),cex = 1.5)
                                netVisual_diffInteraction(cellchat, weight.scale = T,measure = "count", title.name = "Differential number of interactions",margin = c(0.2,0.2,0.4,0.2), comparison = c(1,2),
                                                          targets.use = cell_select,arrow.size = 0.8,vertex.label.cex = 0.8)
                            }
                            save_figure(obj = p,
                                        filename = file.path(opt$output_plot_dir,sub_output_dir,paste0("02_CellChat_Circle_interaction_",opt$group_col,"_",group_levels[2],"_vs_",group_levels[1],"_top_",1,"_targets_",cell_select)),
                                        width = 6.5,
                                        height = 6.5,
                                        res = 600,
                                        formats = opt$formats
                                       )
                        }
            
                        if(  "weight" %in% opt$plot_count_weight ){
                            p <- function(){
                                par(mar = rep(1, 4),cex = 1.5)
                                netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight",title.name = "Differential interaction strength",margin = c(0.2,0.2,0.4,0.2),comparison = c(1,2),
                                                          targets.use = cell_select,arrow.size = 0.8,vertex.label.cex = 0.8)
                            }
                            save_figure(obj = p,
                                        filename = file.path(opt$output_plot_dir,sub_output_dir,paste0("02_CellChat_Circle_weight_",opt$group_col,"_",group_levels[2],"_vs_",group_levels[1],"_top_",1,"_targets_",cell_select)),
                                        width = 6.5,
                                        height = 6.5,
                                        res = 600,
                                        formats = opt$formats
                                       )
                        }
            
                        if( all(c("count","weight") %in% opt$plot_count_weight) ){
                            p <- function(){
                                par(mar = rep(1, 4),cex = 1.5)
                                cowplot::plot_grid(netVisual_diffInteraction(cellchat, weight.scale = T,measure = "count", title.name = "Differential number of interactions",margin = c(0.2,0.2,0.4,0.2), comparison = c(1,2),
                                                                             targets.use = cell_select,arrow.size = 0.8,vertex.label.cex = 0.8),
                                                    netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight",title.name = "Differential interaction strength",margin = c(0.2,0.2,0.4,0.2),comparison = c(1,2),
                                                                              targets.use = cell_select,arrow.size = 0.8,vertex.label.cex = 0.8)
                                                  )
                            }
                            save_figure(obj = p,
                                        filename = file.path(opt$output_plot_dir,sub_output_dir,paste0("02_CellChat_Circle_interaction_weight_",opt$group_col,"_",group_levels[2],"_vs_",group_levels[1],"_top_",1,"_targets_",cell_select)),
                                        width = 6.5,
                                        height = 6.5,
                                        res = 600,
                                        formats = opt$formats
                                       )
                        } 
                    })                  
                    
                }else{
                    cat("\n",blue,setdiff(opt$key_cell,all_celltypes),"Not in Celltype",reset,"\n")
                    quit(save = "no")                  
                }
            }
            
        }
        
    }


    if( opt$plot_allcirc_allheat_allbub_singcirc_singbub[5] ){
        cat("\n\t",yellow," Plot Compare: singBuble",reset,"\n")
        flush.console()
        # 源
         if(opt$key_cell %in%  levels(cellchat@meta[[opt$celltype]]) %>% all() ){
             if( "source" %in% opt$plot_source_target ){
                sources_uses <- match(opt$key_cell, levels(cellchat@idents$joint))
                df_net <- subsetCommunication(cellchat,
                                              net = net.pct,
                                              sources.use = sources_uses,
                                              slot.name = "net",
                                              thresh = opt$plot_Net_netP_Pvalue,   
                                              ligand.pvalues = opt$plot_geneExp_Pvalue,
                                              ligand.logFC = NULL,
                                              ligand.pct.1 = opt$plot_geneExp_minpct*100,
                                              ligand.pct.2 = opt$plot_geneExp_minpct*100,
                                              receptor.pvalues = opt$plot_geneExp_Pvalue,
                                              receptor.logFC = NULL,
                                              receptor.pct.1 = opt$plot_geneExp_minpct*100,
                                              receptor.pct.2 = opt$plot_geneExp_minpct*100
                                             )
                df_net <- df_net %>% dplyr::filter(abs(ligand.logFC) >= opt$plot_geneExp_logFC,abs(receptor.logFC) >= opt$plot_geneExp_logFC)
                df_net <- df_net %>% dplyr::filter(prob > quantile(df_net$prob,probs = opt$plot_prob_Quantile,na.rm = TRUE))
                pairLR_use <- data.frame(interaction_name = unique(df_net$interaction_name))

                if(nrow(pairLR_use) != 0){
                    tryCatch({
                        p1 <- netVisual_bubble(cellchat, sources.use = sources_uses,pairLR.use = pairLR_use, comparison = c(1, 2), angle.x = 45)+
                            theme(plot.margin = margin(t = 5,r = 5,b = 5,l = 50,unit = "pt"))
                        
                        # 保存气泡图数据
                        write.csv(p1@data,file = file.path(opt$output_plot_dir,sub_output_dir,paste0("05_CellChat_Buble_L-R_pairs_source_",paste(opt$key_cell,collapse = "_"),".csv")))
                        
                        width = 4 + length(unique(as.character(p1@data$source.target))) * 0.2
                        height = 5 + length(unique(as.character(p1@data$interaction_name_2))) * 0.1                   
                        save_figure(obj = p1,
                                    filename = file.path(opt$output_plot_dir,sub_output_dir,paste0("05_CellChat_Buble_L-R_pairs_source_",paste(opt$key_cell,collapse = "_"))),
                                    width = min(50,width),
                                    height = min(50,height),
                                    res = 300,
                                    formats = opt$formats
                                   )        
                    }, error = function(e) {NULL})
                }else{
                    cat("\n\t",red,"In Compare mode, no signaling pathways meet the current threshold criteria.",reset,"\n")
                    flush.console()
                }
                 # 单个
                lapply(sources_uses,function(sources_use){
                    df_net <- subsetCommunication(cellchat,
                                                  net = net.pct,
                                                  sources.use = sources_use,
                                                  slot.name = "net",
                                                  thresh = opt$plot_Net_netP_Pvalue,   
                                                  ligand.pvalues = opt$plot_geneExp_Pvalue,
                                                  ligand.logFC = NULL,
                                                  ligand.pct.1 = opt$plot_geneExp_minpct*100,
                                                  ligand.pct.2 = opt$plot_geneExp_minpct*100,
                                                  receptor.pvalues = opt$plot_geneExp_Pvalue,
                                                  receptor.logFC = NULL,
                                                  receptor.pct.1 = opt$plot_geneExp_minpct*100,
                                                  receptor.pct.2 = opt$plot_geneExp_minpct*100
                                                 )
                    df_net <- df_net %>% dplyr::filter(abs(ligand.logFC) >= opt$plot_geneExp_logFC,abs(receptor.logFC) >= opt$plot_geneExp_logFC)
                    df_net <- df_net %>% dplyr::filter(prob > quantile(df_net$prob,probs = opt$plot_prob_Quantile,na.rm = TRUE))
                    pairLR_use <- data.frame(interaction_name = unique(df_net$interaction_name))

                    
                    celltype_name <-  levels(cellchat@idents$joint)[sources_use]
                    if(nrow(pairLR_use) != 0){
                        tryCatch({
                            gg1 <- netVisual_bubble(cellchat, sources.use = sources_use,  comparison = c(1, 2), max.dataset = 2,pairLR.use = pairLR_use, title.name = paste0("Increased signaling in ",group_levels[2]), angle.x = 45, remove.isolate = T)+
                                theme(plot.margin = margin(t = 5,r = 5,b = 5,l = 50,unit = "pt"))
                            # 保存气泡图数据
                            write.csv(gg1@data,file = file.path(opt$output_plot_dir,sub_output_dir,paste0("05_CellChat_Buble_L-R_pairs_Increased_source_",paste(celltype_name,collapse = "_"),"_in_",group_levels[2],".csv")))
                            
                            width = 4 + length(unique(as.character(gg1@data$source.target))) * 0.2
                            height = 5 + length(unique(as.character(gg1@data$interaction_name_2))) * 0.1     
                            save_figure(obj = gg1,
                                        filename = file.path(opt$output_plot_dir,sub_output_dir,paste0("05_CellChat_Buble_L-R_pairs_Increased_source_",paste(celltype_name,collapse = "_"),"_in_",group_levels[2])),
                                        width = min(50,width),
                                        height = min(50,height),
                                        res = 300,
                                        formats = opt$formats
                                       )      
                        }, error = function(e) {NULL})
                      
                        tryCatch({
                            gg2 <- netVisual_bubble(cellchat, sources.use = sources_use,  comparison = c(1, 2), max.dataset = 1,pairLR.use = pairLR_use, title.name = paste0("Decreased signaling in ",group_levels[2]), angle.x = 45, remove.isolate = T)+
                                theme(plot.margin = margin(t = 5,r = 5,b = 5,l = 50,unit = "pt"))
                            # 保存气泡图数据
                            write.csv(gg2@data,file = file.path(opt$output_plot_dir,sub_output_dir,paste0("05_CellChat_Buble_L-R_pairs_Decreased_source_",paste(celltype_name,collapse = "_"),"_in_",group_levels[2],".csv")))
                            
                            width = 4 + length(unique(as.character(gg2@data$source.target))) * 0.2
                            height = 5 + length(unique(as.character(gg2@data$interaction_name_2))) * 0.1          
                            save_figure(obj = gg2,
                                        filename = file.path(opt$output_plot_dir,sub_output_dir,paste0("05_CellChat_Buble_L-R_pairs_Decreased_source_",paste(celltype_name,collapse = "_"),"_in_",group_levels[2])),
                                        width = min(50,width),
                                        height = min(50,height),
                                        res = 300,
                                        formats = opt$formats
                                       )     
                        }, error = function(e) {NULL})                      

                    }else{
                        cat("\n\t",red,"In Compare mode, no signaling pathways meet the current threshold criteria.",reset,"\n")
                        flush.console()                        
                    }    
                })

             }
             
             if("target" %in% opt$plot_source_target ){
                targets_uses <- match(opt$key_cell, levels(cellchat@idents$joint))
                df_net <- subsetCommunication(cellchat,
                                              net = net.pct,
                                              targets.use = targets_uses,
                                              slot.name = "net",
                                              thresh = opt$plot_Net_netP_Pvalue,   
                                              ligand.pvalues = opt$plot_geneExp_Pvalue,
                                              ligand.logFC = NULL,
                                              ligand.pct.1 = opt$plot_geneExp_minpct*100,
                                              ligand.pct.2 = opt$plot_geneExp_minpct*100,
                                              receptor.pvalues = opt$plot_geneExp_Pvalue,
                                              receptor.logFC = NULL,
                                              receptor.pct.1 = opt$plot_geneExp_minpct*100,
                                              receptor.pct.2 = opt$plot_geneExp_minpct*100
                                             )
                df_net <- df_net %>% dplyr::filter(abs(ligand.logFC) >= opt$plot_geneExp_logFC,abs(receptor.logFC) >= opt$plot_geneExp_logFC)
                df_net <- df_net %>% dplyr::filter(prob > quantile(df_net$prob,probs = opt$plot_prob_Quantile,na.rm = TRUE))
                pairLR_use <- data.frame(interaction_name = unique(df_net$interaction_name))

                if(nrow(pairLR_use) != 0){ 
                    tryCatch({
                        p1 <- netVisual_bubble(cellchat, targets.use = targets_uses, pairLR.use = pairLR_use,comparison = c(1, 2), angle.x = 45)+
                            theme(plot.margin = margin(t = 5,r = 5,b = 5,l = 50,unit = "pt"))    
                        # 保存气泡图数据
                        write.csv(p1@data,file = file.path(opt$output_plot_dir,sub_output_dir,paste0("05_CellChat_Buble_L-R_pairs_target_",paste(opt$key_cell,collapse = "_"),".csv")))
                        
                        width = 4 + length(unique(as.character(p1@data$source.target))) * 0.2
                        height = 5 + length(unique(as.character(p1@data$interaction_name_2))) * 0.1  
                        save_figure(obj = p1,
                                    filename = file.path(opt$output_plot_dir,sub_output_dir,paste0("05_CellChat_Buble_L-R_pairs_target_",paste(opt$key_cell,collapse = "_"))),
                                    width = min(50,width),
                                    height = min(50,height),
                                    res = 300,
                                    formats = opt$formats
                                   )      
                    }, error = function(e) {NULL})
                   
                }else{
                    cat("\n\t",red,"In Compare mode, no signaling pathways meet the current threshold criteria.",reset,"\n")
                    flush.console()
                }

                lapply(targets_uses,function(targets_use){ 
                    df_net <- subsetCommunication(cellchat,
                                                  net = net.pct,
                                                  targets.use = targets_use,
                                                  slot.name = "net",
                                                  thresh = opt$plot_Net_netP_Pvalue,   
                                                  ligand.pvalues = opt$plot_geneExp_Pvalue,
                                                  ligand.logFC = NULL,
                                                  ligand.pct.1 = opt$plot_geneExp_minpct*100,
                                                  ligand.pct.2 = opt$plot_geneExp_minpct*100,
                                                  receptor.pvalues = opt$plot_geneExp_Pvalue,
                                                  receptor.logFC = NULL,
                                                  receptor.pct.1 = opt$plot_geneExp_minpct*100,
                                                  receptor.pct.2 = opt$plot_geneExp_minpct*100
                                                 )
                    df_net <- df_net %>% dplyr::filter(abs(ligand.logFC) >= opt$plot_geneExp_logFC,abs(receptor.logFC) >= opt$plot_geneExp_logFC)
                    df_net <- df_net %>% dplyr::filter(prob > quantile(df_net$prob,probs = opt$plot_prob_Quantile,na.rm = TRUE))
                    pairLR_use <- data.frame(interaction_name = unique(df_net$interaction_name))
                    
                    celltype_name <-  levels(cellchat@idents$joint)[targets_use]
                    if(nrow(pairLR_use) != 0){
                        tryCatch({
                            gg1 <- netVisual_bubble(cellchat, targets.use = targets_use,  comparison = c(1, 2), max.dataset = 2,pairLR.use = pairLR_use, title.name = paste0("Increased signaling in ",group_levels[2]), angle.x = 45, remove.isolate = T)+
                                theme(plot.margin = margin(t = 5,r = 5,b = 5,l = 50,unit = "pt"))
                            # 保存气泡图数据
                            write.csv(gg1@data,file = file.path(opt$output_plot_dir,sub_output_dir,paste0("05_CellChat_Buble_L-R_pairs_Increased_target_",paste(celltype_name,collapse = "_"),"_in_",group_levels[2],".csv")))
                            
                            width = 4 + length(unique(as.character(gg1@data$source.target))) * 0.2
                            height = 5 + length(unique(as.character(gg1@data$interaction_name_2))) * 0.1
                            save_figure(obj = gg1,
                                        filename = file.path(opt$output_plot_dir,sub_output_dir,paste0("05_CellChat_Buble_L-R_pairs_Increased_target_",paste(celltype_name,collapse = "_"),"_in_",group_levels[2])),
                                        width = min(50,width),
                                        height = min(50,height),
                                        res = 300,
                                        formats = opt$formats
                                       )     
                        }, error = function(e) {NULL})
                                               
                        tryCatch({
                            gg2 <- netVisual_bubble(cellchat, targets.use = targets_use,  comparison = c(1, 2), max.dataset = 1,pairLR.use = pairLR_use, title.name = paste0("Decreased signaling in ",group_levels[2]), angle.x = 45, remove.isolate = T)+
                                theme(plot.margin = margin(t = 5,r = 5,b = 5,l = 50,unit = "pt"))
                            # 保存气泡图数据
                            write.csv(gg2@data,file = file.path(opt$output_plot_dir,sub_output_dir,paste0("05_CellChat_Buble_L-R_pairs_Decreased_target_",paste(celltype_name,collapse = "_"),"_in_",group_levels[2],".csv")))
                            
                            width = 4 + length(unique(as.character(gg2@data$source.target))) * 0.2
                            height = 5 + length(unique(as.character(gg2@data$interaction_name_2))) * 0.1 
                            save_figure(obj = gg2,
                                        filename = file.path(opt$output_plot_dir,sub_output_dir,paste0("05_CellChat_Buble_L-R_pairs_Decreased_target_",paste(celltype_name,collapse = "_"),"_in_",group_levels[2])),
                                        width = min(50,width),
                                        height = min(50,height),
                                        res = 300,
                                        formats = opt$formats
                                       )      
                        }, error = function(e) {NULL})   

                    }else{
                        cat("\n\t",red,"In Compare mode, no signaling pathways meet the current threshold criteria.",reset,"\n")
                        flush.console()
                    }
                })
             }       
         }else{
            cat("\n",blue,setdiff(opt$key_cell,levels(cellchat@meta[[opt$celltype]])),"Not in Celltype",reset,"\n")
            quit(save = "no")
             
         } 
        
         # 有源和受细胞
         if(length(opt$source_cell) != 0 & length(opt$target_cell) != 0){
            cat("\n\t",yellow," Plot Compare:  Source & Target Buble",reset,"\n")
            flush.console() 
            # 设置了源和受细胞
            source_cells <- match(opt$source_cell,levels(cellchat@idents$joint))
            target_cells <- match(opt$target_cell,levels(cellchat@idents$joint))
            # 拆分成数字向量
            source_cells_list <- as.list(source_cells)
            source_cells_list[[length(source_cells_list) + 1]] <- source_cells
            
            target_cells_list <- as.list(target_cells)
            target_cells_list[[length(target_cells_list) + 1]] <- target_cells
            
            cell_pairs <- expand.grid(source = source_cells_list,target = target_cells_list,stringsAsFactors = FALSE)
             
            lapply(1:nrow(cell_pairs),function(i){
                source_cell_use <- cell_pairs$source[[i]]
                target_cell_use <- cell_pairs$target[[i]] 

                # net
                df_net <- subsetCommunication(cellchat,
                                              net = net.pct,
                                              sources.use = source_cell_use,
                                              targets.use = target_cell_use,
                                              slot.name = "net",
                                              thresh = opt$plot_Net_netP_Pvalue,   
                                              ligand.pvalues = opt$plot_geneExp_Pvalue,
                                              ligand.logFC = NULL,
                                              ligand.pct.1 = opt$plot_geneExp_minpct*100,
                                              ligand.pct.2 = opt$plot_geneExp_minpct*100,
                                              receptor.pvalues = opt$plot_geneExp_Pvalue,
                                              receptor.logFC = NULL,
                                              receptor.pct.1 = opt$plot_geneExp_minpct*100,
                                              receptor.pct.2 = opt$plot_geneExp_minpct*100
                                             )
                df_net <- df_net %>% dplyr::filter(abs(ligand.logFC) >= opt$plot_geneExp_logFC,abs(receptor.logFC) >= opt$plot_geneExp_logFC)
                df_net <- df_net %>% dplyr::filter(prob > quantile(df_net$prob,probs = opt$plot_prob_Quantile,na.rm = TRUE))
                pairLR_use <- data.frame(interaction_name = unique(df_net$interaction_name))
                if(nrow(pairLR_use) != 0){
                    tryCatch({
                        gg1 <- netVisual_bubble(cellchat,sources.use = source_cell_use, targets.use = target_cell_use,comparison = c(1, 2), max.dataset = 2,pairLR.use = pairLR_use, title.name = paste0("Increased signaling in ",group_levels[2]), angle.x = 45, remove.isolate = T)+
                            theme(plot.margin = margin(t = 5,r = 5,b = 5,l = 50,unit = "pt"))
                        # 保存气泡图数据
                        write.csv(gg1@data,file = file.path(opt$output_plot_dir,sub_output_dir,paste0("05_CellChat_Buble_L-R_pairs_Increased_",
                                                                                                   paste(levels(cellchat@idents$joint)[source_cell_use],collapse = "_"),
                                                                                                   "_to_",
                                                                                                   paste(levels(cellchat@idents$joint)[target_cell_use],collapse = "_"),
                                                                                                   "_in_",group_levels[2],".csv")
                                                          )
                                 )
                        
                        width = 4 + length(unique(as.character(gg1@data$source.target))) * 0.2
                        height = 5 + length(unique(as.character(gg1@data$interaction_name_2))) * 0.1
                        save_figure(obj = gg1,
                                    filename = file.path(opt$output_plot_dir,sub_output_dir,paste0("05_CellChat_Buble_L-R_pairs_Increased_",
                                                                                                   paste(levels(cellchat@idents$joint)[source_cell_use],collapse = "_"),
                                                                                                   "_to_",
                                                                                                   paste(levels(cellchat@idents$joint)[target_cell_use],collapse = "_"),
                                                                                                   "_in_",group_levels[2])),
                                    width = min(50,width),
                                    height = min(50,height),
                                    res = 300,
                                    formats = opt$formats
                                   )
                    }, error = function(e) {NULL})
                    
                    tryCatch({
                        gg2 <- netVisual_bubble(cellchat,sources.use = source_cell_use, targets.use = target_cell_use,comparison = c(1, 2), max.dataset = 1,pairLR.use = pairLR_use, title.name = paste0("Decreased signaling in ",group_levels[2]), angle.x = 45, remove.isolate = T)+
                            theme(plot.margin = margin(t = 5,r = 5,b = 5,l = 50,unit = "pt"))
                        # 保存气泡图数据
                        write.csv(gg2@data,file = file.path(opt$output_plot_dir,sub_output_dir,paste0("05_CellChat_Buble_L-R_pairs_Decreased_",
                                                                                                   paste(levels(cellchat@idents$joint)[source_cell_use],collapse = "_"),
                                                                                                   "_to_",
                                                                                                   paste(levels(cellchat@idents$joint)[target_cell_use],collapse = "_"),
                                                                                                   "_in_",group_levels[2],".csv")
                                                          )
                                 )
                        
                        width = 4 + length(unique(as.character(gg2@data$source.target))) * 0.2
                        height = 5 + length(unique(as.character(gg2@data$interaction_name_2))) * 0.1 
                        save_figure(obj = gg2,
                                    filename = file.path(opt$output_plot_dir,sub_output_dir,paste0("05_CellChat_Buble_L-R_pairs_Decreased_",
                                                                                                   paste(levels(cellchat@idents$joint)[source_cell_use],collapse = "_"),
                                                                                                   "_to_",
                                                                                                   paste(levels(cellchat@idents$joint)[target_cell_use],collapse = "_"),
                                                                                                   "_in_",group_levels[2])),
                                    width = min(50,width),
                                    height = min(50,height),
                                    res = 300,
                                    formats = opt$formats
                                   )
                    }, error = function(e) {NULL})
                }else{
                    cat("\n\t",red,"In Compare mode, no signaling pathways meet the current threshold criteria.",reset,"\n")
                    flush.console()
                }
                
            }) 
         }
    }
    
}


#####################################################################################################################
############################################  Main Judgment and Dirs ################################################
#####################################################################################################################
# 判断是佛已有分析结果
if(opt$analysis_mode == "compare"){
    rds_file <- list.files(opt$output_rds_dir,pattern = "05_Cellchat_Compare_obj.rds",full.names = T)
    if(length(rds_file) != 0){
        # 已经有 RDS
        run_CellChat_boole = FALSE
    }else{
        # 无 RDS，进行分析
        run_CellChat_boole = TRUE
    }
}

if(opt$analysis_mode == "single"){
    rds_file <- list.files(opt$output_rds_dir,pattern = "05_Cellchat_Single_obj.rds",full.names = T)
    if(length(rds_file) != 0){
        # 已经有 RDS
        run_CellChat_boole = FALSE
    }else{
        # 无 RDS，进行分析
        run_CellChat_boole = TRUE
    }
}

if(opt$analysis_mode == "all"){
    rds_file_single <- list.files(opt$output_rds_dir,pattern = "05_Cellchat_Single_obj.rds",full.names = T)
    rds_file_compare <- list.files(opt$output_rds_dir,pattern = "05_Cellchat_Single_obj.rds",full.names = T)
    
    if(length(rds_file_single) != 0 & length(rds_file_compare) != 0 ){
        # 已经有 RDS
        run_CellChat_boole = FALSE
    }else{
        # 无 RDS，进行分析
        run_CellChat_boole = TRUE
    }
}

# 创建文件夹
seurat_all <- seurat

# 检查参数
if(!all(opt$key_cell %in% unique(seurat_all@meta.data[[opt$celltype]]))){
    cat("\n",blue,opt$key_cell,"Not all in Celltype",reset,"\n")
    quit(save = "no", status = 0, runLast = FALSE)
}
if(!all(opt$hierarchy_top %in% unique(seurat_all@meta.data[[opt$celltype]]))){
    cat("\n",blue,opt$hierarchy_top,"Not all in Celltype",reset,"\n")
    quit(save = "no", status = 0, runLast = FALSE)
}


# 分组水平、文件夹设置
sub_output_dirs <- c("01_Single","02_Compare")
if(is.null(opt$group_levels)){
    group_levels <- group_new_levels(vector = as.character(seurat_all@meta.data[[opt$group_col]]))
}else{
    if( unlist(opt$group_levels) %in% unique(as.character(seurat_all@meta.data[[opt$group_col]])) %>% all() ){
        group_levels <- opt$group_levels
    }else{
        cat("\n",blue,"The 'group_levels' Not matched Seurat!",reset,"\n")
        quit(save = "no", status = 0)
    }
}

# 文件夹设置
if( tolower(opt$analysis_mode) == "all" ){
    sub_output_dir <- c(paste0(sub_output_dirs[1],"/",group_levels),sub_output_dirs[2])
}else if( tolower(opt$analysis_mode) == "single" ){
    if(class(group_levels) != "list" ){
        sub_output_dir <- paste0(sub_output_dirs[1],"/",group_levels)  # ***************************
    }else if(class(group_levels) == "list" ){
        sub_output_dir <- paste0(sub_output_dirs[1],"/",unlist(group_levels)) 
    }
}else if(  tolower(opt$analysis_mode) == "compare" ){
    sub_output_dir <- sub_output_dirs[2]
}
if( any(opt$plot_pathway_hier_circ_heat_chord_bar_violin_singRoleheat) ){   
    lapply(file.path(opt$output_plot_dir,sub_output_dir,"01_Signa_Pathway"),function(dirname){
        if (!dir.exists(dirname)){dir.create(dirname,recursive = T)}
    })
}else{
    lapply(file.path(opt$output_plot_dir,sub_output_dir),function(dirname){
        if (!dir.exists(dirname)){dir.create(dirname,recursive = T)}
    }) 
}

#####################################################################################################################
############################################  Main CellChat Run #####################################################
#####################################################################################################################
#--运行CellChat-----------------------------------------------------------
if(run_CellChat_boole){
    cat("\n",yellow,"############################## Start Run CellChat ##############################",reset,"\n")
    flush.console()
    
    # 选择分组样本
    seurat_all <- subset(seurat_all,!!sym(opt$group_col) %in% unlist(group_levels))
    if(length(group_levels) > 1){
    	# 筛选共有Celltype
    	meta_group1 <- seurat_all@meta.data %>% dplyr::filter(!!sym(opt$group_col) == unlist(group_levels)[1])
    	meta_group2 <- seurat_all@meta.data %>% dplyr::filter(!!sym(opt$group_col) == unlist(group_levels)[2])
    	com_celltypes <- intersect(unique(meta_group1[[opt$celltype]]),unique(meta_group2[[opt$celltype]]))
    
    	seurat_all <- subset(seurat_all, !!sym(opt$celltype) %in% com_celltypes)
    
    	# Levels
    	seurat_all@meta.data[[opt$group_col]] <- factor(seurat_all@meta.data[[opt$group_col]],levels = unlist(group_levels))
    }
    options(future.globals.maxSize = 100 * 1024^3) 
    future::plan("multisession", workers = opt$ncore)
    
    # 判断 分析模式  all/single/compare
    if(tolower(opt$analysis_mode) == "compare" ){
        cat("\n",yellow," Start Run CellChat Compare ~",reset,"\n")
        flush.console()
        # 进行差异Cellchat
        compare_obj <- run_compare_analysis(seurat_all, group_levels, opt)
    
    }else if(tolower(opt$analysis_mode) == "single" ){
        cat("\n",yellow," Start Run CellChat Single ~",reset,"\n")
        flush.console()
        # 进行单独 Cellchat
        single_obj <- run_single_analysis(seurat_all, group_levels, opt)
        
    }else if(tolower(opt$analysis_mode) == "all" ){
        cat("\n",yellow," Start Run CellChat Compare & Single ~",reset,"\n")
        flush.console()
        # 进行 single + compare
        compare_obj <- run_compare_analysis(seurat_all, group_levels, opt)
        single_obj  <- run_single_analysis(seurat_all, group_levels, opt)
        
    }
    plan("sequential")
    
    # Save
    if(exists("compare_obj")){
        saveRDS(compare_obj,file = file.path(opt$output_rds_dir,"05_Cellchat_Compare_obj.rds"))
    }
    if(exists("single_obj")){
        saveRDS(single_obj,file = file.path(opt$output_rds_dir,"05_Cellchat_Single_obj.rds"))
    }
}

#####################################################################################################################
############################################  Main CellChat Plot ####################################################
#####################################################################################################################
# -- 出图---------------------------------------------------------
cat("\n\n",yellow,"############################## Start Plot CellChat ##############################",reset,"\n\n")
flush.console()

if(opt$analysis_mode == "all"){
    cat("\n",yellow," Start Plot CellChat Compare & Single ~",reset,"\n")
    flush.console()    
    # All 模式出图
    # read RDS
    cellchat_list <- lapply(c("Single","Compare"),function(type){
        readRDS(file.path(opt$output_rds_dir,paste0("05_Cellchat_",type,"_obj.rds")))
    })
    cellchat_list <- setNames(cellchat_list,c("Single","Compare"))

    # Single
    if(length(group_levels) > 1){
        invisible(suppressMessages(suppressWarnings({
            cat("\n\n",yellow," Start Plot CellChat Single: ",reset,group_levels[1],"\n")
            flush.console()  
            single_cellchat_plot(cellchat = cellchat_list$Single[[group_levels[1]]], opt = opt, sub_output_dir = sub_output_dir[1])
            cat("\n\n",yellow," Start Plot CellChat Single: ",reset,group_levels[2],"\n")
            flush.console()  
            single_cellchat_plot(cellchat = cellchat_list$Single[[group_levels[2]]], opt = opt, sub_output_dir = sub_output_dir[2])
        })))
    }else{
        invisible(suppressMessages(suppressWarnings({
            cat("\n\n",yellow," Start Plot CellChat Single: ",reset,group_levels[1],"\n")
            flush.console() 
            single_cellchat_plot(cellchat = cellchat_list$Single[[group_levels[1]]], opt = opt, sub_output_dir = sub_output_dir[1])
        })))
    } 

    # Compare
    invisible(suppressMessages(suppressWarnings({
        cat("\n\n",yellow," Start Plot CellChat Compare: ",reset,group_levels[2]," vs ",group_levels[1],"\n")
        flush.console() 
        compare_cellchat_plot(cellchat = cellchat_list$Compare, opt = opt, sub_output_dir = sub_output_dir[3], group_levels = group_levels)
    })))

}else if(opt$analysis_mode == "single"){
    cat("\n\n",yellow," Start Plot CellChat Single ~",reset,"\n")
    flush.console()   
    # Single 模式出图
    cellchat <- readRDS(file.path(opt$output_rds_dir,"05_Cellchat_Single_obj.rds"))
    # single
    invisible(suppressMessages(suppressWarnings({
        for(i in 1:length(group_levels)){
            cat("\n\n",yellow," Start Plot CellChat Single: ",reset,group_levels[[i]],"\n")
            flush.console()
            single_cellchat_plot(cellchat = cellchat[[group_levels[[i]]]], opt = opt, sub_output_dir = sub_output_dir[[i]])
        }
    })))
}else if(opt$analysis_mode == "compare"){
    cat("\n",yellow," Start Plot CellChat Compare ~",reset,"\n")
    flush.console()   
    # Compare 模式出图
    cellchat <- readRDS(file.path(opt$output_rds_dir,"05_Cellchat_Compare_obj.rds"))
    
    invisible(suppressMessages(suppressWarnings({
        cat("\n\n",yellow," Start Plot CellChat Compare: ",reset,group_levels[2]," vs ",group_levels[1],"\n")
        flush.console() 
        compare_cellchat_plot(cellchat = cellchat, opt = opt, sub_output_dir = sub_output_dir, group_levels = group_levels)
    })))
}
}
