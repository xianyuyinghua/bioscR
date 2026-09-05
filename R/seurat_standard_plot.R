#' Save standard Seurat embedding and cell-composition plots
#'
#' Creates embedding plots for combinations of grouping and splitting metadata columns, together with cell-count and cell-proportion plots for each split.
#'
#' @param seurat A Seurat object with the requested reduction and metadata
#'   columns. Supply this argument explicitly. Grouping columns should be
#'   factors with the desired cell-type order.
#' @param reduction Name of an existing dimensional reduction, by default
#'   "umap.harmony".
#' @param group_by Character vector of metadata column prefixes used for
#'   colouring cells. Prefixes are combined into a regular expression anchored
#'   at the start of the column name.
#' @param split_by Character vector of metadata column prefixes used for
#'   splitting plots, matched in the same way as group_by. An unsplit embedding
#'   is also generated.
#' @param colors Vector of colours passed to the embedding and composition
#'   plots. Defaults to basicR::get_colors(100.1).
#' @param formats Character vector of output formats passed to
#'   basicR::save_figure(), by default PDF and PNG.
#' @param output_dir Directory passed to basicR::save_figure() for saving
#'   plots. Defaults to Output/03_Celltype_Marker. Create a writable directory
#'   before use.
#'
#' @details
#' For every matched grouping column, the function creates split and unsplit
#' Seurat::DimPlot() figures, each with either a legend or cell-group labels.
#' Point size is calculated from the total cell count and clipped to [0.01, 8].
#' Facet layout and figure dimensions are selected automatically.
#'
#' For each splitting column, cell-count and within-group cell-proportion
#' figures include alluvial flows. Composition file names are shared between
#' the legend and label iterations, so these files are saved again during the
#' loop. Embedding file names include a grouping index, grouping column,
#' reduction (dots replaced by underscores), Legend or Label, and splitting
#' column (NULL for unsplit plots). Composition file names contain Cell_Number
#' or Cell_Proportion. Figures are saved at 600 dpi.
#'
#' The function attaches dplyr, Seurat, ggplot2, and basicR, and additionally
#' uses stringr, tibble, tidyr, ggalluvial, and scales. It uses
#' group_new_levels() for non-factor splitting columns.
#'
#' @return A list with one element per grouping, splitting, and legend
#'   combination. Elements are the return values of the final composition save
#'   for split combinations, and NULL for unsplit combinations. The primary
#'   output is the set of saved figure files.
#'
#' @examples
#' \dontrun{
#' seu <- readRDS("seurat.rds")
#' seu$Celltype <- factor(seu$Celltype)
#' dir.create("standard_plots", recursive = TRUE, showWarnings = FALSE)
#' seurat_standard_plot(
#'   seurat = seu, reduction = "umap",
#'   group_by = "Celltype", split_by = c("sample", "group"),
#'   output_dir = "standard_plots"
#' )
#' }
seurat_standard_plot <- function(seurat = seurat,
                                 reduction = "umap.harmony",
                                 group_by = c("Cluster","Celltype"),
                                 split_by = c("sample","group"),
                                 colors = basicR::get_colors(100.1),
                                 formats = c("pdf","png"),
                                 output_dir = file.path("Output","03_Celltype_Marker")
                                ){

    # ANSI 颜色代码
    red <- "\033[31m"
    green <- "\033[32m"
    yellow <- "\033[33m"
    blue <- "\033[34m"
    magenta <- "\033[35m"
    cyan <- "\033[36m"
    reset <- "\033[0m"
    
    suppressPackageStartupMessages({
        library(dplyr)
        library(Seurat)
        library(ggplot2)
        library(basicR)
    })

    
    ########################################## Function ######################################################
    
    calculate_cell_point_size <- function(cell_count) {
        # 单细胞定义计算点大小的函数
        a <- exp(3.84342)  # 截距的指数
        b <- -0.57370      # 回归系数
        point_size <- a * cell_count^b
        point_size <- pmin(pmax(point_size, 0.01), 8)
        return(point_size)  
    }
    
    facet_nrow_ncol <- function(facet_number){
        # 根据分面数量计算行列数
        nrow = round(sqrt(facet_number))
        ncol = ceiling(sqrt(facet_number))
        if(facet_number == 3){
            nrow = 1
            ncol = 3
        }
        return(c(nrow,ncol))
    }
  
    ######################################## Caculate ###############################################################
    reduction_name <- gsub("[.]","_",reduction)
    
    # Dimplot
    ## 设置点大小，细胞数，图片尺寸
    pt_size <- calculate_cell_point_size(cell_count = ncol(seurat) )  # 细胞点大小
    raster <- ifelse(ncol(seurat) < 10,TRUE,FALSE)  # TRUE：点以矢量方式绘制，适合于少量细胞数；FALSE：点先转换成光栅图（bitmap），然后再嵌入最终图形，适合大量细胞

    ## 设置分面的出图组合
    group_by_colnames <- grep(paste0("^(", paste(group_by, collapse = "|"), ")"),colnames(seurat@meta.data),value = TRUE)
    split_by_colnames <- c(grep(paste0("^(", paste(split_by, collapse = "|"), ")"),colnames(seurat@meta.data),value = TRUE), NA)
    legend_or_not     <- c(FALSE, TRUE)
    df_combinations <- expand.grid(group_by = group_by_colnames,
                                   split_by = split_by_colnames,
                                   legend   = legend_or_not,
                                   stringsAsFactors = FALSE
                                  ) %>% dplyr::arrange(group_by,split_by,legend)

    df_combinations <- df_combinations %>% dplyr::mutate(number = sprintf("%02d", as.integer(factor(group_by, levels = unique(group_by))))) 

    ####################################### Main Plot #################################################################
    
    lapply(seq_len(nrow(df_combinations)),function(i){

        ## Dimplot
        # 计算分面数
        if( length(unique(seurat@meta.data[[df_combinations$split_by[i]]])) == 0 ){
            facet_nums <- c(1,1)
        } else {
            facet_nums <- facet_nrow_ncol(facet_number = seurat@meta.data[[df_combinations$split_by[i]]] %>% unique() %>% length())
        }
        # 字体缩放大小
        text_scale_size <- max(facet_nums)/min(facet_nums)/10 + 1

        
        # 判断是否分面
        if (is.na(df_combinations$split_by[i])) {
            split_by_colname <- NULL
            split_by_name <- "NULL"
        } else {
            split_by_colname <- df_combinations$split_by[i]
            split_by_name <- split_by_colname
        }

        str_length_max <- stringr::str_length(unique(seurat@meta.data[[df_combinations$group_by[i]]])) %>% max()
        ncol_legend <- ceiling(length(unique(seurat@meta.data[[df_combinations$group_by[i]]]))/20)
        
        # Label and p size
        if( df_combinations$legend[i] ){
            # 有图例
            label_name <- "Legend"

            height_size <- 6 * facet_nums[1]
            width_size  <- 6 * facet_nums[2] + str_length_max * 0.1 * ncol_legend
        }else{
            # 无图例
            label_name <- "Label"

            height_size <- 6 * facet_nums[1]
            width_size  <- 6 * facet_nums[2]
        }
        
        # picture size
        scale_factor <- (6 + facet_nums[1]*1) / height_size
        height_size <- height_size * scale_factor  # 结果就是 8
        width_size  <- width_size * scale_factor

        
        cat("\n",yellow,"Plot: ",reset,df_combinations$group_by[i],"&",df_combinations$split_by[i],"&",label_name,"\n")
        flush.console()
        
        # Plot
        DimPlot(seurat,
                reduction = reduction,
                cols = colors,
                group.by = df_combinations$group_by[i],
                split.by = split_by_colname,
                ncol = facet_nums[2],
                label = !(df_combinations$legend[i]),
                label.size = 5 * text_scale_size,
                repel = T,
                pt.size = pt_size,
                raster.dpi = c(800, 800),
                raster = raster
               ) &
        guides(color = guide_legend(override.aes = list(size = 5), ncol = ncol_legend)) &
        theme(axis.title = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank(),
              plot.title = element_text(size = 20 * text_scale_size, face = "bold"),
              strip.text = element_text(size = 16 * text_scale_size, face = "bold"),
              legend.title = element_text(size = 17 * text_scale_size, face = "bold"), 
              legend.text = element_text(size = 15 * text_scale_size),
              legend.key.size = unit(0.8, "cm") 
             )-> p1
        
        if(!df_combinations$legend[i]){
            # 不要图例
            p1 <- p1 & NoLegend() 
        }
        basicR::save_figure(obj = p1,
                    filename =  file.path(output_dir, paste0(df_combinations$number[i],"_",df_combinations$group_by[i],"_",reduction_name,"_",label_name,"_by_",split_by_name)),
                    width = width_size,
                    height = height_size,
                    res = 600,
                    formats = formats # "pdf","tiff","png","svg"
                   )   

        ## Proportion
        if(!is.na(df_combinations$split_by[i])){
            group_type_length <- seurat@meta.data[[df_combinations$split_by[i]]] %>% as.character() %>% unique() %>% length()
            my_theme <- theme_classic(base_size = 18) +
                          theme(
                            plot.title = element_text(size = 18, hjust = 0.5),
                            panel.grid.minor = element_line(colour = NA, linetype = "blank"),
                            panel.background = element_rect(fill = NA),
                            plot.background = element_rect(colour = NA),
                            axis.title = element_text(face = "bold"),
                            axis.text = element_text(face = "bold"),
                            axis.title.y = element_text(size = 17, margin = ggplot2::margin(r = 10)),
                            axis.line.x = element_line(color = "black", linewidth = 0.7),
                            axis.line.y = element_line(color = "black", linewidth = 0.7),
                            legend.position = "right",
                            legend.title = element_blank(),  
                            legend.text = element_text(size = 16),
                            plot.margin = ggplot2::margin(l = 15, t = 30, r = 10, b = 10),
                            legend.key.size = unit(0.65, "cm")
                          )
    
                
                tab <- table(seurat@meta.data[[df_combinations$split_by[i]]],seurat@meta.data[[df_combinations$group_by[i]]])
                tab_df <- as.data.frame.matrix(tab)
                tab_df <- tibble::rownames_to_column(tab_df, var = "Var1")
                plot_df <- tab_df %>% tidyr::pivot_longer(cols = -Var1,names_to = "Var2",values_to = "value",names_transform = list(Var2 = as.character))
                    
                colnames(plot_df) <- c("Group","CellType","Number")
                split_var <- seurat@meta.data[[df_combinations$split_by[i]]]
                if(is.factor(split_var)){
                    grouplevels <- levels(split_var)
                }else{
                    grouplevels <- group_new_levels(
                        vector = plot_df$Group %>% as.character()
                    )
                }
                plot_df$Group <- factor(plot_df$Group,levels = grouplevels)  # 设置水平
                plot_df$CellType <- factor(plot_df$CellType,levels = levels(seurat@meta.data[[df_combinations$group_by[i]]]))  # 设置水平
                plot_df = plot_df %>% group_by(Group) %>% dplyr::mutate(percent = Number/sum(Number))
                
                ggplot(plot_df, aes(x = Group, y = Number, fill = CellType,stratum=CellType, alluvium=CellType)) +
                    geom_bar(stat = "identity", width=0.8,aes(group=CellType),position="fill") +
                    scale_y_continuous(labels = waiver())+ ####用来将y轴移动位置  
                    geom_col(width = 0.5, color=NA) + 
                    ggalluvial::geom_flow(width= 0.5,alpha=0.4, knot.pos=0.5)+
                    labs(y = "Cell Number", x = NULL) +
                    scale_fill_manual(values = colors)+
                    guides(fill = guide_legend(ncol = ncol_legend)) + 
                    my_theme -> p1
                if(group_type_length > 2){
                    p1 <- p1 + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"))
                }
            
                basicR::save_figure(obj = p1,
                                    filename =  file.path(output_dir, paste0(df_combinations$number[i],"_",df_combinations$group_by[i],"_Cell_Number_by_",split_by_name)),
                                    width = 4 + 0.3 * group_type_length + str_length_max * 0.1 * ncol_legend,
                                    height = 6,
                                    res = 600,
                                    formats = formats
                                   )

            
                ggplot(plot_df, aes(x = Group, y = percent, fill = CellType,stratum=CellType, alluvium=CellType)) +
                    geom_bar(stat = "identity", width=0.8,aes(group=CellType),position="fill") +
                    scale_y_continuous(labels = scales::percent_format())+ ####用来将y轴移动位置  
                    geom_col(width = 0.5, color=NA) + 
                    ggalluvial::geom_flow(width= 0.5,alpha=0.4, knot.pos=0.5)+
                    labs(y = "Cell Proportion", x = NULL) +
                    scale_fill_manual(values = colors)+
                    guides(fill = guide_legend(ncol = ncol_legend)) + 
                    my_theme -> p2
                if(group_type_length > 2){
                    p2 <- p2 + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"))
                }
                basicR::save_figure(obj = p2,
                                    filename =  file.path(output_dir, paste0(df_combinations$number[i],"_",df_combinations$group_by[i],"_Cell_Proportion_by_",split_by_name)),
                                    width = 3.5 + 0.3 * group_type_length + str_length_max * 0.1 * ncol_legend,
                                    height = 6,
                                    res = 600,
                                    formats = formats
                                   )   
            
        }

        # Other
        
    }) 
   
}