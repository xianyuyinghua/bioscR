#' Run inferCNV analysis and plot cell-type CNV profiles
#'
#' Runs inferCNV on raw counts from a Seurat object, caches the analysis, and
#' saves a chromosome heatmap and cell-type CNV score boxplot.
#'
#' @param seurat A Seurat object with gene symbols as feature names, a counts
#' layer in the selected assay, and the metadata columns specified below.
#' @param celltype_col Name of a factor metadata column containing cell-type
#' annotations. Remove unused factor levels before calling. Defaults to
#' "Celltype".
#' @param sample_col Name of the metadata column identifying samples. Used
#' together with celltype_col for heatmap downsampling. Defaults to "sample".
#' @param species Species passed to AnnoProbe::annoGene(), for example
#' "human" or "mouse". Defaults to "human". Chromosome ordering and
#' annotation are fixed to chr1 through chr22 in the current implementation.
#' @param ref_name A single reference cell-type name present in celltype_col.
#' Required despite the NULL default. Use a name without hyphens, because
#' metadata labels are converted from hyphens to underscores internally,
#' while ref_name is not converted.
#' @param assay Assay containing the counts layer. Defaults to "RNA".
#' @param heatmap_downsample Optional positive integer giving the maximum
#' number of cells per sample and cell-type combination in the heatmap. NULL
#' uses all cells. Analysis and CNV scoring always use all retained cells.
#' @param heatmap_cnv_probs Two quantile probabilities used to set the
#' heatmap colour range symmetrically around 1. Defaults to c(0.025, 0.975).
#' @param celltype_col_colors Vector of plotting colours, with enough colours
#' for all cell-type factor levels. Defaults to basicR::get_colors(number =
#' 50.1).
#' @param chr_colors Colours repeated across the 22 chromosome annotation
#' levels. Defaults to c("#E64B35", "#4DBBD5").
#' @param ncor Number of threads passed to infercnv::run(). Defaults to 10.
#' @param output_dir Existing writable directory for input tables, RDS
#' caches, CSV tables, and figures. Defaults to the working directory. The
#' function does not create this directory. The inferCNV run directory is
#' separately fixed to infercnv_output in the working directory.
#' @param formats Character vector of figure formats passed to
#' basicR::save_figure(). Defaults to c("pdf", "png").
#'
#' @details
#' Gene positions are obtained with AnnoProbe::annoGene(). Mitochondrial and
#' sex chromosomes (chrM, chrX, and chrY) are excluded; duplicated gene
#' symbols are removed and genes are ordered by chromosome and start
#' position. Hyphens in cell names, sample labels, and cell-type labels are
#' replaced by underscores.
#'
#' The inferCNV run uses cutoff = 0.1, cluster_by_groups = FALSE,
#' hclust_method = "ward.D2", analysis_mode = "subclusters", denoise = TRUE,
#' HMM = FALSE, and leiden_resolution = "auto". Intermediate and final
#' inferCNV plotting are disabled. The function sets options(scipen = 100)
#' when running inferCNV and does not restore the previous value.
#'
#' Input files are expFile.txt, geneFile.txt, groupFiles.txt, and
#' cell_meta.rds. The initial and analysed inferCNV objects are saved as
#' infercnv_obj.rds and infercnv_obj_analysis.rds. Existing cell_meta.rds and
#' infercnv_obj_analysis.rds files are reused without checking the input or
#' parameters. For a fresh analysis, use a new output directory and a working
#' directory without a previous infercnv_output run. A partial cache
#' containing cell_meta.rds but no analysed object cannot be resumed reliably
#' because input file path variables are only assigned during input
#' preparation.
#'
#' The heatmap matrix is saved to 01_inferCNV_Heatmap.csv with cells in rows
#' and genes in columns. Heatmap figures use the same filename stem. Rows are
#' split by cell type and columns by chromosome, without clustering.
#' Downsampling sets the random seed to 123 and affects only this heatmap and
#' its CSV table.
#'
#' CNV score thresholds are derived from the reference expression matrix: the
#' lower and upper bounds are the mean of gene-wise means minus or plus twice
#' the mean of gene-wise standard deviations. Their difference defines one
#' copy interval. Each gene contributes 0 within the central interval, 1 in
#' either adjacent interval, and 2 beyond them. Scores are summed per cell.
#' Metadata and scores are saved to 02_Celltype_TotalCNV.csv; the boxplot
#' uses the same filename stem and includes a dashed line at the reference
#' median score.
#'
#' The function loads dplyr, ggplot2, Seurat, infercnv, AnnoProbe,
#' ComplexHeatmap, circlize, and grid, and also uses tibble and basicR. These
#' packages must be installed. If ref_name is NULL or absent from the
#' original cell-type metadata, the current implementation calls q("no"),
#' which terminates the R session.
#'
#' @return The return value of the final basicR::save_figure() call. Analysis
#' objects, plotting data, and figures are saved to disk; the function does
#' not return a modified Seurat object.
#'
#' @examples
#' \dontrun{
#' seu <- readRDS("seurat.rds")
#' seu$Celltype <- droplevels(factor(seu$Celltype))
#' # Replace T_cells with a reference level present in seu$Celltype.
#' dir.create("cnv_results", recursive = TRUE, showWarnings = FALSE)
#' inferCNV_analysis(
#'   seurat = seu,
#'   ref_name = "T_cells",
#'   heatmap_downsample = 200,
#'   output_dir = "cnv_results"
#' )
#' }
#' @export
inferCNV_analysis <- function(seurat,
                              celltype_col = "Celltype",
                              sample_col = "sample",
                              species = "human",
                              ref_name = NULL,
                              assay = "RNA",
                              heatmap_downsample = NULL,
                              heatmap_cnv_probs = c(0.025, 0.975),
                              celltype_col_colors = basicR::get_colors(number = 50.1),
                              chr_colors = c("#E64B35", "#4DBBD5"),
                              ncor = 10,
                              output_dir = ".",
                              formats = c("pdf", "png")
                             ){

    suppressPackageStartupMessages({
        library(dplyr)
        library(ggplot2)
        library(Seurat)
        library(infercnv)
        library(AnnoProbe)
        library(ComplexHeatmap)
        library(circlize)
        library(grid)
    })

    # ============================================================
    # 0. 参数检查与初始化
    # ============================================================

    red <- "\033[31m"
    green <- "\033[32m"
    yellow <- "\033[33m"
    magenta <- "\033[35m"
    reset <- "\033[0m"

    if (is.null(ref_name) || length(ref_name) == 0) {
        stop("Please set 'ref_name'.")
    }

    if (!all(c(celltype_col, sample_col) %in% colnames(seurat@meta.data))) {
        stop("Required metadata columns not found.")
    }

    if (!all(ref_name %in%  as.character(seurat@meta.data[[celltype_col]]))) {
        stop("Reference cell types not found: ",paste(setdiff(ref_name,unique(as.character(seurat@meta.data[[celltype_col]]))), collapse = ", "))
    }

    if (!is.null(heatmap_downsample)) {
        if (length(heatmap_downsample) != 1 ||
            is.na(heatmap_downsample) || 
            !is.numeric(heatmap_downsample) || 
            heatmap_downsample < 1 ||
            heatmap_downsample != floor(heatmap_downsample)) {
            stop("heatmap_downsample must be NULL or a positive integer.")
        }
    }

    if (length(heatmap_cnv_probs) != 2 ||
        anyNA(heatmap_cnv_probs) ||
        any(heatmap_cnv_probs < 0 | heatmap_cnv_probs > 1) ||
        heatmap_cnv_probs[1] >= heatmap_cnv_probs[2]) {
        stop("Invalid heatmap_cnv_probs.")
    }

    dir.create(output_dir,recursive = TRUE,showWarnings = FALSE)

    # 统一处理参考细胞名称
    ref_name <- unique(gsub("-", "_", ref_name))

    # ============================================================
    # 1. 统一定义输入文件路径
    # ============================================================

    file_meta <- file.path(output_dir, "cell_meta.rds")
    file_exp <- file.path(output_dir, "expFile.txt")
    file_gene <- file.path(output_dir, "geneFile.txt")
    file_group <- file.path(output_dir, "groupFiles.txt")
    file_infercnv_obj <- file.path(output_dir, "infercnv_obj.rds")
    file_infercnv_obj2 <- file.path(output_dir, "infercnv_obj_analysis.rds")
    infercnv_run_dir <- file.path(output_dir, "infercnv_output")

    # ============================================================
    # 2. 准备当前细胞元数据
    # ============================================================

    df_meta_current <- seurat@meta.data %>% dplyr::select(dplyr::all_of(c(sample_col, celltype_col)) ) %>%
        tibble::rownames_to_column("Cell")

    # 保持与 inferCNV 输入一致的细胞名称
    df_meta_current$Cell <- gsub("-", "_", df_meta_current$Cell)
    df_meta_current[[sample_col]] <- gsub("-","_",as.character(df_meta_current[[sample_col]]))
    df_meta_current[[celltype_col]] <- gsub("-","_",as.character(df_meta_current[[celltype_col]]))

    if (anyDuplicated(df_meta_current$Cell)) {
        stop("Duplicated cell names after replacing '-' with '_'.")
    }

    if (anyNA(df_meta_current[[sample_col]]) || anyNA(df_meta_current[[celltype_col]])) {
        stop("Sample or cell-type metadata contains NA.")
    }

    original_levels <- if (  is.factor(seurat@meta.data[[celltype_col]])) {
        levels(seurat@meta.data[[celltype_col]])
    }else{
        unique(as.character(seurat@meta.data[[celltype_col]]))
    }

    original_levels <- unique(gsub("-", "_", original_levels))

    df_meta_current[[celltype_col]] <- factor(df_meta_current[[celltype_col]],levels = original_levels)

    # ============================================================
    # 3. 检查输入缓存
    # ============================================================

    input_files <- c(file_meta,file_exp,file_gene,file_group)

    if (!all(file.exists(input_files))) {
        cat("\n", yellow,"Preparing inferCNV input files...",reset, "\n")
        flush.console()

        # --------------------------------------------------------
        # 3.1 获取原始 counts
        # --------------------------------------------------------

        dat <- Seurat::GetAssayData(seurat,assay = assay,layer = "counts")
        colnames(dat) <- gsub("-", "_", colnames(dat))

        # 确保表达矩阵与元数据对应
        if (!setequal(colnames(dat),df_meta_current$Cell)) {
            stop("Counts cell names do not match metadata.")
        }

        dat <- dat[,df_meta_current$Cell,drop = FALSE]

        # --------------------------------------------------------
        # 3.2 获取基因染色体位置信息
        # --------------------------------------------------------

        df_gene <- AnnoProbe::annoGene(rownames(dat),"SYMBOL",species)

        df_gene <- df_gene %>% dplyr::filter(chr %in% paste0("chr", 1:22)) %>%
            dplyr::mutate( chr = factor(chr,levels = paste0("chr", 1:22))) %>%
            dplyr::arrange(chr, start) %>% dplyr::select(SYMBOL, chr, start, end) %>%
            dplyr::distinct(SYMBOL, .keep_all = TRUE)

        df_gene <- df_gene[df_gene$SYMBOL %in% rownames(dat),,drop = FALSE]

        if (nrow(df_gene) == 0) {
            stop("No genes matched chromosome annotations.")
        }

        # --------------------------------------------------------
        # 3.3 按染色体位置排序表达矩阵
        # --------------------------------------------------------

        dat <- dat[match(df_gene$SYMBOL,rownames(dat)),, drop = FALSE]

        # --------------------------------------------------------
        # 3.4 制作细胞类型注释
        # --------------------------------------------------------

        df_meta <- df_meta_current
        df_group <- df_meta %>%  dplyr::select(Cell,dplyr::all_of(celltype_col))

        # --------------------------------------------------------
        # 3.5 保存输入文件
        # --------------------------------------------------------

        write.table(dat,file = file_exp,sep = "\t",quote = FALSE)
        write.table(df_gene,file = file_gene,sep = "\t",quote = FALSE,col.names = FALSE,row.names = FALSE)
        write.table(df_group,file = file_group,sep = "\t",quote = FALSE,col.names = FALSE,row.names = FALSE)

        saveRDS(df_meta, file_meta)

        rm(dat, df_gene, df_group)
        gc()

        cat("\n", green,"Input files prepared.",reset, "\n")

    }else{

        cat("\n", green,"Using existing inferCNV input files.",reset, "\n")

        df_meta <- readRDS(file_meta)

        # 检查缓存是否属于当前细胞数据
        if (!setequal(df_meta$Cell,df_meta_current$Cell)) {
            stop("Cached cell metadata differs from current Seurat object. ","Use a new output_dir for a fresh analysis.")
        }

        old_meta <- df_meta[match(df_meta_current$Cell,df_meta$Cell),,drop = FALSE]

        if (!identical(as.character(old_meta[[sample_col]]),as.character(df_meta_current[[sample_col]])) || !identical(
            as.character(old_meta[[celltype_col]]),as.character(df_meta_current[[celltype_col]])
        )) {
            stop("Cached sample/cell-type metadata has changed. ","Use a new output_dir.")
        }
    }

    # ============================================================
    # 4. 运行 inferCNV
    # ============================================================

    if (!file.exists(file_infercnv_obj2)) {

        cat("\n", yellow,"Start inferCNV...",reset, "\n")
        flush.console()

        infercnv_obj <- infercnv::CreateInfercnvObject(raw_counts_matrix = file_exp,
                                                       annotations_file = file_group,
                                                       delim = "\t",
                                                       gene_order_file = file_gene,
                                                       ref_group_names = ref_name
                                                      )

        saveRDS(infercnv_obj,file_infercnv_obj)

        options(scipen = 100)

        infercnv_obj2 <- infercnv::run(infercnv_obj,
                                        cutoff = 0.1,
                                        out_dir = infercnv_run_dir,
                                        cluster_by_groups = FALSE,
                                        hclust_method = "ward.D2",
                                        analysis_mode = "subclusters",
                                        denoise = TRUE,
                                        HMM = FALSE,
                                        plot_steps = FALSE,
                                        inspect_subclusters = FALSE,
                                        no_plot = TRUE,
                                        leiden_resolution = "auto",
                                        num_threads = ncor
                                      )

        saveRDS(infercnv_obj2,file_infercnv_obj2)

        rm(infercnv_obj)
        gc()
        cat("\n", green,"inferCNV analysis completed.",reset, "\n")
    } else {
        cat("\n", green,"Using cached inferCNV analysis.",reset, "\n")
        infercnv_obj2 <- readRDS(file_infercnv_obj2)
    }

    # ============================================================
    # 5. 检查参考细胞与缓存
    # ============================================================

    cached_ref_names <- names(infercnv_obj2@reference_grouped_cell_indices)

    if (!setequal(cached_ref_names,ref_name)) {
        stop(
            "Reference cell types differ from cached inferCNV results.\n",
            "Current references: ",
            paste(ref_name, collapse = ", "), "\n",
            "Cached references: ",
            paste(cached_ref_names, collapse = ", "), "\n",
            "Please use a new output_dir."
        )
    }

    df_meta <- readRDS(file_meta)

    cnv_mat <- infercnv_obj2@expr.data
    gene_order <- infercnv_obj2@gene_order

    # ============================================================
    # 6. 热图细胞降采样
    # ============================================================

    cat("\n", yellow,"Start Heatmap...",reset, "\n")
    flush.console()

    df_cell_group <- df_meta %>% dplyr::filter( Cell %in% colnames(cnv_mat) ) %>%
        dplyr::arrange(.data[[celltype_col]])

    if (!is.null(heatmap_downsample)) {
        set.seed(123)
        group_cols <- unique(c(sample_col, celltype_col))

        df_cell_group_downsample <- df_cell_group %>%
            dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
            dplyr::group_modify(
                ~ {
                    dplyr::slice_sample(.x,n = min(nrow(.x),heatmap_downsample))
                }
            ) %>%
            dplyr::ungroup() %>%
            dplyr::arrange(.data[[celltype_col]]) %>%
            dplyr::select(Cell,dplyr::all_of(group_cols))

        cat("\n\t", magenta,"Heatmap using downsample ",nrow(df_cell_group_downsample)," cells.",reset, "\n")

    } else {

        df_cell_group_downsample <- df_cell_group

        cat("\n\t", magenta,"Heatmap using all ",nrow(df_cell_group_downsample)," cells.",reset, "\n")

        if (nrow(df_cell_group_downsample) > 30000) {
            warning(
                "Heatmap uses more than 30000 cells. ",
                "Consider heatmap_downsample to reduce memory usage."
            )
        }
    }

    # ============================================================
    # 7. 构建热图矩阵
    # ============================================================

    mat_plot <- t( cnv_mat[,df_cell_group_downsample$Cell, drop = FALSE])

    write.csv(mat_plot,file = file.path(output_dir,"01_inferCNV_Heatmap.csv"))

    # ============================================================
    # 8. 热图颜色
    # ============================================================

    q <- quantile(as.numeric(mat_plot),probs = heatmap_cnv_probs,na.rm = TRUE)

    lim <- max( 1 - q[1], q[2] - 1)

    # 防止颜色断点重复
    if (!is.finite(lim) || lim <= 0) { lim <- 0.01 }

    breaks <- c(1 - lim,1,1 + lim)

    col_fun <- circlize::colorRamp2(breaks,c("#2166AC", "white", "#B2182B"))

    # ============================================================
    # 9. 细胞类型颜色注释
    # ============================================================

    group_factor <- df_cell_group_downsample[[celltype_col]]

    if (!is.factor(group_factor)) {
        group_factor <- factor(group_factor)
    }

    group_factor <- droplevels(group_factor)

    n_celltypes <- nlevels(group_factor)

    if (length(celltype_col_colors) < n_celltypes) {
        stop("Not enough celltype_col_colors.")
    }

    group_color <- setNames(celltype_col_colors[seq_len(n_celltypes)], levels(group_factor))

    ha_left <- ComplexHeatmap::rowAnnotation(
        Celltype = group_factor,
        col = list(Celltype = group_color),
        show_annotation_name = FALSE,
        show_legend = FALSE,
        simple_anno_size = grid::unit(2, "mm")
    )

    # ============================================================
    # 10. 染色体注释
    # ============================================================

    chr_factor <- factor(gsub("chr","",as.character(gene_order$chr)),levels = as.character(1:22))

    chr_color_map <- setNames(rep(chr_colors,length.out = 22),levels(chr_factor))

    ha_bottom <- ComplexHeatmap::HeatmapAnnotation(
        Chromosome = chr_factor,
        col = list(
            Chromosome = chr_color_map
        ),
        show_annotation_name = FALSE,
        show_legend = FALSE,
        simple_anno_size = grid::unit(2, "mm")
    )

    # ============================================================
    # 11. 绘制 inferCNV Heatmap
    # ============================================================

    ht <- ComplexHeatmap::Heatmap(
        mat_plot,
        name = "CNV",
        col = col_fun,

        cluster_rows = FALSE,
        cluster_columns = FALSE,

        show_row_names = FALSE,
        show_column_names = FALSE,

        row_split = group_factor,
        column_split = chr_factor,

        row_title = "%s",
        row_title_side = "left",
        row_title_rot = 0,
        row_title_gp = grid::gpar(fontsize = 11, fontface = "bold"),

        left_annotation = ha_left,

        column_title = "%s",
        column_title_side = "bottom",
        column_title_gp = grid::gpar(fontsize = 11,fontface = "bold" ),

        bottom_annotation = ha_bottom,

        row_gap = grid::unit(2, "mm"),
        column_gap = grid::unit(2, "mm"),

        use_raster = TRUE,
        raster_quality = 10
    )

    ht_drawn <- ComplexHeatmap::draw(ht,padding = grid::unit(c(10, 10, 10, 5),"mm"))

    basicR::save_figure(obj = ht_drawn,
                        filename = file.path(output_dir,"01_inferCNV_Heatmap"),
                        width = 13,
                        height = min(8,max(6,nrow(mat_plot) * 0.0005)),
                        res = 600,
                        formats = formats
                       )

    rm(ht, ht_drawn, mat_plot)
    gc()

    # ============================================================
    # 12. 获取多个参考细胞的 CNV 表达矩阵
    # ============================================================

    cat("\n", yellow,"Calculating CNV scores...",reset, "\n")

    # 支持 ref_name = c("Epi1", "Epi2")
    ref_cells <- unlist(
        infercnv_obj2@reference_grouped_cell_indices[
            ref_name
        ],
        use.names = FALSE
    )

    ref_cells <- unique(as.integer(ref_cells))

    if (length(ref_cells) == 0 ||
        anyNA(ref_cells) ||
        any(ref_cells < 1) ||
        any(ref_cells > ncol(cnv_mat))) {
        stop("Invalid reference cell indices.")
    }

    tmp <- cnv_mat[ , ref_cells, drop = FALSE]

    if (ncol(tmp) < 2) {
        stop("At least two reference cells are required.")
    }

    # ============================================================
    # 13. 计算参考细胞 CNV 范围
    # ============================================================

    ref_gene_mean <- rowMeans(tmp)

    ref_gene_sd <- apply(tmp,1,sd)

    down <- mean(ref_gene_mean) - 2 * mean(ref_gene_sd)

    up <- mean(ref_gene_mean) + 2 * mean(ref_gene_sd)

    oneCopy <- up - down

    a2 <- down - oneCopy
    a3 <- up + oneCopy

    rm(tmp, ref_gene_mean, ref_gene_sd)
    gc()

    # ============================================================
    # 14. 计算全部细胞的 CNV score
    # ============================================================

    # 保持原函数的 CNV score 计算方式
    # 分块计算以降低大矩阵的临时内存占用

    cell_scores <- numeric(ncol(cnv_mat))

    names(cell_scores) <- colnames(cnv_mat)

    block_size <- 1000

    for (start_idx in seq(1,ncol(cnv_mat),by = block_size)) {

        end_idx <- min( start_idx + block_size - 1,ncol(cnv_mat))

        idx <- start_idx:end_idx

        block <- cnv_mat[,idx,drop = FALSE]

        # 轻度 CNV 异常计 1 分
        score1 <- colSums((block >= a2 & block < down) |(block >= up & block <= a3))

        # 较强 CNV 异常计 2 分
        score2 <- 2 * colSums(block < a2 | block > a3)

        cell_scores[idx] <- score1 + score2
    }

    cell_scores_CNV <- data.frame(Cell = names(cell_scores),
                                  cnv_score = as.numeric(cell_scores),
                                  stringsAsFactors = FALSE
                                 )

    df_cell_total_cnvScore <- merge(df_meta,cell_scores_CNV,by = "Cell")

    write.csv(df_cell_total_cnvScore,file = file.path(output_dir,"02_Celltype_TotalCNV.csv"),row.names = FALSE)

    # ============================================================
    # 15. 绘制 CNV score 箱线图
    # ============================================================

    # 多个参考细胞类型：分别计算各自中位数，再取最小的那个
    ref_score_medians <- sapply(ref_name, function(x) {
        median(
            df_cell_total_cnvScore$cnv_score[
                as.character(df_cell_total_cnvScore[[celltype_col]]) == x
            ],
            na.rm = TRUE
        )
    })
    
    ref_score_median <- min(ref_score_medians, na.rm = TRUE)

    p <- ggplot(df_cell_total_cnvScore,aes(x = .data[[celltype_col]],y = cnv_score,color = .data[[celltype_col]])) +
        geom_boxplot() +
        scale_color_manual(values = group_color) +
        geom_hline(yintercept = ref_score_median,linetype = "dashed",color = "grey70",linewidth = 0.6) +
        labs(x = NULL,y = "Total CNV score") +
        theme_classic( base_size = 16) +
        theme(
            axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face = "bold"),
            axis.text.y = element_text(face = "bold"),
            axis.title = element_text(size = 17,face = "bold"),
            legend.position = "none"
        )

    basicR::save_figure(obj = p,
                        filename = file.path(output_dir,"02_Celltype_TotalCNV"),
                        width = 6,
                        height = 6,
                        res = 600,
                        formats = formats
                       )

    cat("\n", green,"inferCNV analysis and plotting completed.",reset, "\n")

    # 与原函数一致：结果主要保存到文件
    # 额外以 invisible 方式返回 CNV score 数据
    invisible(df_cell_total_cnvScore)
}
