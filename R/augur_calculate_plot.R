#' Calculate and plot pairwise Augur cell-type prioritisation
#'
#' Runs Augur for every pair of group factor levels in a Seurat object, caches the results, and saves AUC tables and cell-type lollipop plots.
#'
#' @param seurat A Seurat object with expression data and metadata columns
#'   specified by cell_type_col and group_col. A factor metadata column named
#'   sample is also required.
#' @param cell_type_col Name of the metadata column containing cell-type
#'   annotations.
#' @param group_col Name of a factor metadata column defining the groups.
#'   Supply at least two observed levels and remove unused levels before
#'   calling.
#' @param colors Vector of plotting colours, with enough colours for the cell
#'   types. Defaults to basicR::get_colors(number = 50.1).
#' @param formats Character vector of output formats passed to
#'   basicR::save_figure(), by default PDF and PNG.
#' @param output_dir Existing writable directory for RDS caches, CSV tables,
#'   and figures. The function does not create this directory before saving RDS
#'   files.
#'
#' @details
#' Each pair is analysed with Augur::calculate_auc(), using min_cells = 10,
#' subsample_size = 10, and n_threads = 1. Comparisons run through
#' parallel::mclapply() with one core per pair; multiple pairs require a
#' platform supporting forked parallelism. No random seed is set internally.
#'
#' Files use the prefix 06_Augur_AUC_Celltype_by_ followed by the two group
#' levels joined with an underscore. Existing RDS files are reused without
#' checking whether the input data have changed. CSV tables and figures are
#' regenerated from these caches. AUC tables are sorted in ascending order;
#' plots include a reference line at AUC 0.5 and are saved at 6 by 6 inches and
#' 600 dpi.
#'
#' The implementation uses Augur, Seurat, dplyr, ggplot2, and basicR.
#' Unqualified helpers, including sym(), the pipe operator, and plotting
#' functions, must be available in the evaluation environment.
#'
#' @return A list with one element per group pair, containing the return values
#'   of basicR::save_figure(). The Augur results themselves are saved as RDS
#'   files; sorted AUC tables are saved as CSV files.
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(dplyr)
#' library(ggplot2)
#' seu <- readRDS("seurat.rds")
#' seu$group <- droplevels(factor(seu$group))
#' seu$sample <- factor(seu$sample)
#' dir.create("augur_output", recursive = TRUE, showWarnings = FALSE)
#' augur_calculate_plot(seurat = seu, output_dir = "augur_output")
#' }
#' @export
augur_calculate_plot <- function(seurat,
                                 cell_type_col = "Celltype",
                                 group_col = "group",
                                 colors = basicR::get_colors(number = 50.1),
                                 formats = c("pdf", "png"),
                                 output_dir = "03_Celltype_Marker"
                                ){
    # 分组因子判断
    
    # 分组设置
    group_combinations <- combn(levels(seurat@meta.data[[group_col]]),2,simplify = FALSE)

    parallel::mclapply(1:length(group_combinations),function(i){
        rds_file <- file.path(output_dir,paste0("06_Augur_AUC_Celltype_by_",paste(group_combinations[[i]],collapse = "_"),".rds"))
        if(!file.exists(rds_file)){
            # 筛选分组
            subseurat <- subset(seurat,subset = !!sym(group_col) %in% group_combinations[[i]])
            subseurat@meta.data[[group_col]] <- droplevels(subseurat@meta.data[[group_col]])
            subseurat@meta.data$sample <- droplevels(subseurat@meta.data$sample)
            
            augur <- Augur::calculate_auc(subseurat,
                                          cell_type_col = cell_type_col,
                                          label_col = group_col,
                                          min_cells = 10,
                                          subsample_size = 10,
                                          n_threads = 1  # 不能设置多线程，否则超级慢
                                         )
            saveRDS(augur,file = file.path(output_dir,paste0("06_Augur_AUC_Celltype_by_",paste(group_combinations[[i]],collapse = "_"),".rds")))
        }
    }, mc.cores = length(group_combinations))
        
    # Plot    
    lapply(1:length(group_combinations),function(i){
        rds_file <- file.path(output_dir,paste0("06_Augur_AUC_Celltype_by_",paste(group_combinations[[i]],collapse = "_"),".rds"))
        augur <- readRDS(rds_file)
        
        # AUC表（按AUC从低到高）
        df <- augur$AUC %>% arrange(auc) %>% dplyr::rename("Celltype" = "cell_type")
        df$Celltype <- factor(df$Celltype,levels = unique(df$Celltype))
        write.csv(df,file = file.path(output_dir,paste0("06_Augur_AUC_Celltype_by_",paste(group_combinations[[i]],collapse = "_"),".csv")))
        
        # 绘制棒棒糖
        p <- ggplot(df, aes(x = reorder(Celltype, auc), y = auc)) +
          # 棒棒糖的竖线（用相同的颜色映射）
          geom_segment(aes(x = Celltype, xend = Celltype, 
                           y = 0, yend = auc, color = Celltype),
                       size = 0.8) +
          # 顶端的点
          geom_point(aes(color = Celltype), size = 5, alpha = 1) +
          geom_hline(yintercept = 0.5,color = "grey50",linetype = "dashed",linewidth = 0.6,alpha = 0.5) +
          coord_flip() +
          scale_color_manual(values = colors) +
          labs(title = paste0(group_combinations[[i]][2]," vs ",group_combinations[[i]][1] ), x = NULL, y = "AUC(Augur)", color = "Celltype") +
          theme_classic(base_size = 16) +
          theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
                axis.text = element_text(size = 16, face = "bold"),
                axis.title.x = element_text(size = 16, face = "bold"),
                legend.position = "none",
                axis.line = element_line(linewidth = 0.8),
                plot.margin = margin(b = 10,l = 10,t = 10,r = 20)
          )
        
        basicR::save_figure(obj = p,
                            filename = file.path(output_dir,paste0("06_Augur_AUC_Celltype_by_",paste(group_combinations[[i]],collapse = "_"))),
                            width = 6,
                            height = 6,
                            res = 600,
                            formats = formats
                           )   
    })
    
}