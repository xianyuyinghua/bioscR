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
                              species = "human", # huamn/mouse
                              ref_name = NULL,
                              assay = "RNA",
                              heatmap_downsample = NULL,
                              heatmap_cnv_probs = c(0.025,0.975),
                              celltype_col_colors = basicR::get_colors(number = 50.1),
                              chr_colors = c("#E64B35", "#4DBBD5"),
                              ncor = 10,
                              output_dir = ".",
                              formats = c("pdf","png")
                             ){

    # 
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

    # colors
    red <- "\033[31m"
    green <- "\033[32m"
    yellow <- "\033[33m"
    blue <- "\033[34m"
    magenta <- "\033[35m"
    cyan <- "\033[36m"
    reset <- "\033[0m"

    
    # 参数检查
    if(is.null(ref_name)){
        cat("\n",yellow,"Please set the 'ref_name'.",reset,"\n")
        flush.console()
        q("no")
    }else if(!(ref_name %in% seurat@meta.data[[celltype_col]])){
        cat("\n",yellow,"Please Verify the 'ref_name'.",reset,"\n")
        flush.console()
        q("no")        
    }

    if(is.null(heatmap_downsample)){
        cat("\n",magenta,"Analysis and Heatmap using all ",dim(seurat)[2]," cells.",reset,"\n")
        flush.console()
    }
    
    file_meta <- file.path(output_dir,'cell_meta.rds')
    if(!file.exists(file_meta)){
        # 表达量矩阵
        dat <- GetAssayData(seurat,layer = 'counts',assay = assay)
        colnames(dat) <- gsub("-", "_", colnames(dat)) # infercnv 不识别特殊符号
    
        # 基因信息
        df_gene <- annoGene(rownames(dat),"SYMBOL",species)
        df_gene <- df_gene[!df_gene$chr %in% c("chrM", "chrX", "chrY"),]
        df_gene$chr <- factor(df_gene$chr,levels = paste0("chr",1:22))
        df_gene <- df_gene %>% arrange(chr,start) %>% dplyr::select("SYMBOL","chr","start","end")
        df_gene <- df_gene %>% distinct(SYMBOL, .keep_all = TRUE)
    
        # 表达数据筛选和排序
        dat <- dat[rownames(dat) %in% df_gene$SYMBOL,]
        dat <- dat[match(df_gene$SYMBOL,rownames(dat)),]
        
        # 文件制作2：样本的描述
        df_meta <- seurat@meta.data %>% dplyr::select(sample_col,celltype_col) %>% tibble::rownames_to_column("Cell")
        df_meta$Cell <- gsub("-", "_", df_meta$Cell)
        df_meta[[celltype_col]] <- gsub("-", "_", df_meta[[celltype_col]])
        df_meta[[sample_col]] <- gsub("-", "_", df_meta[[sample_col]])
        df_meta[[celltype_col]] <- factor(df_meta[[celltype_col]],gsub("-", "_", levels(seurat@meta.data[[celltype_col]])))
        df_group <- df_meta %>% dplyr::select("Cell",celltype_col)
    
        # 保存 表达、基因、样本分组 数据
        file_exp <- file.path(output_dir,'expFile.txt')
        file_gene <- file.path(output_dir,'geneFile.txt')
        file_group <- file.path(output_dir,'groupFiles.txt')
        
        write.table(dat, file = file_exp, sep = '\t', quote = F)
        write.table(df_gene, file = file_gene, sep = '\t',quote = F, col.names = F, row.names = F)
        write.table(df_group,file = file_group, sep = '\t',quote = F, col.names = F, row.names = F)
        saveRDS(df_meta,file = file_meta)        
    }
    
    ###########################################   Main   ###############################################
    file_infercnv_obj2 <- file.path(output_dir,"infercnv_obj_analysis.rds")
    if(!file.exists(file_infercnv_obj2)){
        cat("\n",yellow,"Start inferCNV...",reset,"\n")
        flush.console()
        # 创建 inferCNV 对象
        infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = file_exp,
                                             annotations_file = file_group,
                                             delim = "\t",
                                             gene_order_file = file_gene,
                                             ref_group_names = ref_name  # 参考
                                            )
        saveRDS(infercnv_obj,file = file.path(output_dir,"infercnv_obj.rds"))
    
        # Run inferCNV
        options(scipen = 100)
        infercnv_obj2 <- infercnv::run(infercnv_obj,
                                       cutoff =  0.1, # smart-seq选择1,10X选择0.1
                                       out_dir = "infercnv_output", # dir is auto
                                       cluster_by_groups =  F, #是否根据患者类型(由细胞注释文件中定义)
                                       hclust_method = "ward.D2",# ward.D2 方法进行层次聚类
                                       analysis_mode = "subclusters", # 默认是samples，推荐是subclusters
                                       denoise = TRUE, # 去噪音
                                       HMM = F,  ##特别耗时间,是否要去背景噪音
                                       plot_steps = F, #不在每个步骤后生成图形
                                       inspect_subclusters = FALSE, # 不画 subcluster 检查热图
                                       no_plot = TRUE, # 不画最终大热图
                                       leiden_resolution = "auto", #可以手动调参数
                                       num_threads = ncor # 4线程工作，加快速度
                                       )
        saveRDS(infercnv_obj2,file = file_infercnv_obj2)
    }

    
    #########################################   Plot Heatmap  ###############################################
    cat("\n",yellow,"Start Heatmap...",reset,"\n")
    flush.console()
    df_meta <- readRDS(file_meta)
    infercnv_obj2 <- readRDS(file_infercnv_obj2)
    
    cnv_mat <- infercnv_obj2@expr.data
    gene_order <- infercnv_obj2@gene_order

    df_cell_group <- df_meta %>% dplyr::filter(Cell %in% colnames(cnv_mat)) %>% arrange(!!sym(celltype_col))

    # 降采样出热图
    if(!is.null(heatmap_downsample)){
        # 降采样
        set.seed(123)
        group_cols <- unique(c(sample_col, celltype_col))
        df_cell_group_downsample <- df_cell_group %>%
            dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
            dplyr::group_modify(~ {dplyr::slice_sample(.x, n = min(nrow(.x), heatmap_downsample))}) %>%
            dplyr::ungroup() %>% dplyr::arrange(.data[[celltype_col]]) %>%
            dplyr::select("Cell",dplyr::all_of(group_cols))
        
        cat("\n\t",magenta,"Heatmap using downsample ",nrow(df_cell_group_downsample)," cells.",reset,"\n")
        flush.console()
    }else{
        # 全采样
        df_cell_group_downsample <- df_cell_group
        cat("\n\t",magenta,"Heatmap using all ",nrow(df_cell_group_downsample)," cells.",reset,"\n")
        flush.console()
    }
    
    # 降采样CNV矩阵
    mat_plot <- cnv_mat[, df_cell_group_downsample$Cell, drop = FALSE] %>% t()
    write.csv(mat_plot,file = file.path(output_dir,"01_inferCNV_Heatmap.csv"))
    
    # 设置颜色范围
    q <- quantile(as.numeric(mat_plot), probs = heatmap_cnv_probs, na.rm = TRUE)
    lim <- max(1 - q[1],q[2] - 1)
    breaks <- c(1 - lim,1,1 + lim)
    col_fun <- colorRamp2(breaks,c("#2166AC", "white", "#B2182B"))


    # 分组颜色
    group_factor <- df_cell_group_downsample[[celltype_col]]
    group_color <- setNames(celltype_col_colors[1:length(levels(group_factor))],
                          levels(group_factor)
                         )
    # 左侧细胞分组颜色条
    ha_left <- rowAnnotation(Celltype = group_factor,
                             col = list(Celltype = group_color),
                             show_annotation_name = FALSE,
                             show_legend = FALSE,
                             simple_anno_size = unit(2, "mm")
                            )

    # 染色体颜色
    # 只显示染色体颜色数字编号
    chr_factor <- factor(gsub("chr", "", as.character(gene_order$chr)),levels = as.character(1:22))
    chr_colors <- setNames(rep(chr_colors,length.out = 22),levels(chr_factor))
    # 底部染色体颜色条
    ha_bottom <- HeatmapAnnotation(Chromosome = chr_factor,
                                   col = list(Chromosome = chr_colors),
                                   show_annotation_name = FALSE,
                                   show_legend = FALSE,
                                   simple_anno_size = unit(2, "mm")
                                  )
    # Heatmap
    ht <- Heatmap(mat_plot,
                  name = "CNV",
                  col = col_fun,
                  cluster_rows = FALSE,
                  cluster_columns = FALSE,
                  show_row_names = FALSE,
                  show_column_names = FALSE,
                  row_split = group_factor, # 按 Epi 分块
                  column_split = chr_factor, # 按染色体分块
                  row_title = "%s", # 每个 Epi slice 只显示一次名称
                  row_title_side = "left",
                  row_title_rot = 0,
                  row_title_gp = gpar(fontsize = 11,fontface = "bold"),
                  left_annotation = ha_left, # 左侧细胞分类颜色条
                  column_title = "%s", # 染色体编号
                  column_title_side = "bottom",
                  column_title_gp = gpar(fontsize = 11,fontface = "bold"),
                  bottom_annotation = ha_bottom, # 底部染色体颜色条
                  row_gap = unit(2, "mm"), # 行间距
                  column_gap = unit(2, "mm"), # 列间距
                  use_raster = TRUE, # 栅格化处理
                  raster_quality = 10 # 栅格化清晰度
                 )
    
    basicR::save_figure(obj = draw(ht,padding = unit(c(10, 10, 10, 5), "mm")),
                    filename =  file.path(output_dir, "01_inferCNV_Heatmap"),
                    width = 13,
                    height = min(8, max(6,nrow(mat_plot) * 0.0005)),
                    res = 600,
                    formats = formats  # "pdf","tiff","png","svg"
                   )

    ############################################  Plot box  ###########################################
    cat("\n",yellow,"Start Boxplot...",reset,"\n")
    flush.console()
    tmp <- cnv_mat[,infercnv_obj2@reference_grouped_cell_indices[[ref_name]], drop = FALSE]

    down <- mean(rowMeans(tmp)) - 2 * mean(apply(tmp, 1, sd))
    up   <- mean(rowMeans(tmp)) + 2 * mean(apply(tmp, 1, sd))
    
    oneCopy <- up - down
    
    a2 <- down - oneCopy
    a3 <- up + oneCopy

    cell_scores_CNV <- data.frame(
        cnv_score = colSums( (cnv_mat >= a2 & cnv_mat < down) | (cnv_mat >= up & cnv_mat <= a3)) +
        2 * colSums(cnv_mat < a2 | cnv_mat > a3)
    ) %>% tibble::rownames_to_column("Cell")

    df_cell_total_cnvScore <- merge(df_meta,cell_scores_CNV,by = "Cell")
    write.csv(df_cell_total_cnvScore,file = file.path(output_dir,"02_Celltype_TotalCNV.csv"))
    
    ggplot(df_cell_total_cnvScore, aes(x = !!sym(celltype_col),y = cnv_score,color = !!sym(celltype_col))) +
        geom_boxplot() +
        scale_color_manual(values = celltype_col_colors) +
        geom_hline(yintercept = median(df_cell_total_cnvScore$cnv_score[df_cell_total_cnvScore[[celltype_col]] == ref_name]), linetype = "dashed", color = "grey70", size = 0.6) +
        labs(y = "Total CNV score") +
        theme_classic(base_size = 16) +
        theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face = "bold"),
              axis.text.y = element_text(,face = "bold"),
              axis.title = element_text(size = 17,face = "bold"),
              legend.position = "none"
             ) -> p

    basicR::save_figure(obj = p,
                        filename =  file.path(output_dir, "02_Celltype_TotalCNV"),
                        width = 6,
                        height = 6,
                        res = 600,
                        formats = formats  # "pdf","tiff","png","svg"
                       )
}