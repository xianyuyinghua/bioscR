#' Calculate and plot pairwise Milo neighbourhood differential abundance
#'
#' Performs neighbourhood differential abundance analysis for every pair of
#' group factor levels in a Seurat object and saves neighbourhood graphs, cell-
#' type beeswarm plots with boxplots, and annotated result tables.
#'
#' @param seurat A Seurat object with expression data suitable for conversion
#'   by Seurat::as.SingleCellExperiment(), metadata columns for groups and cell
#'   types, and a factor column named sample identifying biological samples.
#'   Each sample must belong to one group; supply biological replication for
#'   the differential abundance model.
#' @param cell_type_col Name of a factor metadata column containing cell-type
#'   annotations. Defaults to "Celltype". The current plotting code also
#'   accesses the annotated result column Celltype explicitly, so use
#'   "Celltype" for the supported workflow.
#' @param group_col Name of a factor metadata column defining comparison
#'   groups, with at least two observed levels. Remove unused levels before
#'   calling. The name must be usable directly in an R formula because the
#'   model is constructed as ~ group_col without quoting the column name.
#' @param reduction Name of the embedding used for neighbourhood graph
#'   visualisation. Its uppercase form is passed as the layout name in the
#'   converted SingleCellExperiment object. Defaults to "umap.harmony". This
#'   argument does not select the dimensional reduction used to construct the
#'   nearest-neighbour graph.
#' @param celltype_colors Vector of colours for the cell-type boxplots, with
#'   enough colours for the annotated cell types. Defaults to
#'   basicR::get_colors(number = 50.1).
#' @param continuous_colors Three colours for negative, central, and positive
#'   log fold changes, respectively.
#' @param formats Character vector of figure formats passed to
#'   basicR::save_figure(). Defaults to PDF and PNG.
#' @param output_dir Writable directory for figures and CSV tables. Defaults to
#'   "03_Celltype_Marker". Create the directory before calling; the function
#'   contains no explicit directory-creation step.
#'
#' @details
#' Comparisons follow the order returned by combn() on the levels of the group
#' factor and run sequentially. Each pair is subsetted, unused group, sample,
#' and cell-type levels are dropped, and the subset is converted to a
#' SingleCellExperiment and then a Milo object.
#'
#' The nearest-neighbour graph uses k = 30 and d = 30. Neighbourhood sampling
#' uses prop = 0.2, k = 30, d = 30, and refined = TRUE. Neighbourhood distances
#' also use d = 30. These calls use the default reduced-dimension selection of
#' miloR; provide input suitable for its 30-dimensional workflow. The reduction
#' argument is used only for graph visualisation. No random seed is set
#' internally; set one before calling if needed.
#'
#' Cells are counted by the fixed metadata column sample. The design contains
#' one row per sample and uses only the grouping variable, with no additional
#' covariates or explicit contrast supplied to testNhoods(). Under standard
#' treatment contrasts, the second retained group level is compared with the
#' first. Neighbourhoods are annotated using cell_type_col, and the plotted
#' Celltype order is reversed relative to the input factor levels.
#'
#' Both plotNhoodGraphDA() and plotDAbeeswarm() receive alpha = 1. The plots
#' therefore do not apply a conventional 0.05 significance threshold. The graph
#' fill scale is limited to [-5, 5], with outlying colours squished to the
#' endpoints; the beeswarm plot includes reference lines at log fold changes -1
#' and 1.
#'
#' For each group pair, the function saves a neighbourhood graph with prefix
#' 07_Celltype_Cell_OR_milo_by_group_ and a beeswarm/boxplot figure and CSV
#' table with prefix 07_Celltype_Cell_Beeswarm_and_Boxplots_by_group_. Each
#' prefix is followed by the two group levels joined by an underscore. Despite
#' the OR text in the graph filename, the displayed values are log fold
#' changes, not odds ratios. Graph figures are 6.5 by 5.5 inches and beeswarm
#' figures are 6.5 by 6 inches, all at 600 dpi. CSV files contain the annotated
#' testNhoods() results. There is no result cache or explicit save of the Milo
#' object.
#'
#' The function attaches dplyr, Seurat, miloR, scran, scater, SeuratWrappers,
#' SingleCellExperiment, patchwork, ggbeeswarm, scales, forcats, data.table,
#' stringr, and basicR. It also calls ggnewscale::new_scale_color() and uses
#' unqualified ggplot2 plotting functions; these must be available in the
#' evaluation environment.
#'
#' @return A list with one element per group pair, containing the return values
#'   of the final basicR::save_figure() call for the beeswarm/boxplot figure.
#'   The primary outputs are the saved figures and CSV tables; the function
#'   does not directly return the Milo objects or differential abundance
#'   tables.
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(ggplot2)
#' seu <- readRDS("seurat.rds")
#' seu$group <- droplevels(factor(seu$group))
#' seu$sample <- droplevels(factor(seu$sample))
#' seu$Celltype <- droplevels(factor(seu$Celltype))
#' dir.create("milo_output", recursive = TRUE, showWarnings = FALSE)
#'
#' # Supply an object suitable for the 30-dimensional graph workflow
#' # and containing the requested embedding.
#' set.seed(123)
#' miloR_calculate_plot(
#'   seurat = seu,
#'   cell_type_col = "Celltype",
#'   group_col = "group",
#'   reduction = "umap.harmony",
#'   output_dir = "milo_output"
#' )
#' }
#' @importFrom dplyr distinct sym
#' @importFrom Seurat as.SingleCellExperiment
#' @importFrom SeuratObject "Idents<-"
#' @importFrom miloR makeNhoods countCells calcNhoodDistance testNhoods buildNhoodGraph plotNhoodGraphDA annotateNhoods plotDAbeeswarm
#' @importFrom SummarizedExperiment colData
#' @importFrom ggplot2 aes element_blank element_text geom_boxplot geom_hline labs margin scale_color_gradient2 scale_color_manual scale_fill_gradient2 theme theme_classic
#' @importFrom scales squish
#' @importFrom stats as.formula
#' @importFrom utils combn write.csv
#' @export
miloR_calculate_plot <- function(seurat,
                                 cell_type_col = "Celltype",
                                 group_col = "group",
                                 reduction = "umap.harmony", 
                                 celltype_colors = basicR::get_colors(number = 50.1),
                                 continuous_colors = c("#070091","lightgrey","#910000"),
                                 formats = c("pdf","png"),
                                 output_dir = "03_Celltype_Marker"
                                ){
    
    suppressPackageStartupMessages({
        library(dplyr)
        library(Seurat)
        library(miloR)
        library(scran)
        library(scater)
        library(SeuratWrappers)
        library(SingleCellExperiment)
        library(patchwork)
        library(ggbeeswarm)
        library(scales)
        library(forcats)
        library(data.table)
        library(stringr)
        library(basicR)
    })

    # 分组 
    group_combinations <- combn(levels(seurat@meta.data[[group_col]]),2,simplify = FALSE)

    # Run
    lapply(1:length(group_combinations),function(i){
        subseurat <- subset(seurat,subset = !!sym(group_col) %in% group_combinations[[i]])
        subseurat@meta.data[[group_col]] <- droplevels(subseurat@meta.data[[group_col]])
        subseurat@meta.data$sample <- droplevels(subseurat@meta.data$sample)
        subseurat@meta.data[[cell_type_col]] <- droplevels(subseurat@meta.data[[cell_type_col]])
        
        # 1. 构建milo对象
        Idents(subseurat) <- subseurat@meta.data[[group_col]]
        
        #miloR输入对象是SingleCellExperiment，所以我们是常用的seurat对象的话，转化一下
        seurat_sc <- as.SingleCellExperiment(subseurat)
        scmilo <- miloR::Milo(seurat_sc)
    
        # 2. 构建KNN-Graph
        scmilo <- miloR::buildGraph(scmilo,k = 30, d = 30)
    
        # 3. 在图上定义邻域
        scmilo <- makeNhoods(scmilo, 
                             prop = 0.2, #定义要随机抽样的图顶点的比例，通常为0.1-0.2
                             k = 30, #建议使用与buildGraph一样的k值
                             d = 30, #KNN降维数，建议使用与buildGraph一样的d值
                             refined = TRUE
                            )
        
        # 4. 量化每个邻域的细胞数量
        scmilo <- countCells(scmilo,meta.data = data.frame(colData(scmilo)), sample="sample")
    
        # 5. 构建分组矩阵
        traj_design <- data.frame(colData(scmilo))[,c("sample", group_col)]#分别是重复样本ID和分组
        traj_design$sample <- as.factor(traj_design$sample)
        traj_design <- distinct(traj_design)
        rownames(traj_design) <- traj_design$sample
    
        # 6. 计算细胞距离
        scmilo <- calcNhoodDistance(scmilo, d = 30)
    
        # 7.两个分组之间的差异分析
        design_formula <- as.formula(paste0("~ ", group_col))
        da_results <- testNhoods(scmilo,
                                 design = design_formula, # 可以增加协变量
                                 design.df = traj_design
                                )
    
        # 为可视化建立一个抽象的邻域图
        scmilo <- buildNhoodGraph(scmilo)
        plotNhoodGraphDA(scmilo, layout  = toupper(reduction), da_results, alpha=1)+
        scale_fill_gradient2(low = continuous_colors[1],mid = continuous_colors[2],high = continuous_colors[3],name = "log2FC",limits = c(-5,5),oob = squish)+
        labs(title = paste0(group_combinations[[i]][2]," vs ",group_combinations[[i]][1])) +
        theme_classic(base_size = 16) +
        theme(axis.text = element_blank(),
              axis.title = element_blank(),
              axis.line = element_blank(),
              axis.ticks = element_blank(),
              legend.title = element_text(size = 14,face = "bold"),
              legend.text = element_text(size = 13),
              plot.title = element_text(face = "bold",hjust = 0.5),
              plot.margin = margin(t = 10,r = 10,b = 5,l = 5,unit = "mm")
             ) -> p2
        basicR::save_figure(obj = p2,
                            filename =  file.path(output_dir,paste0("07_Celltype_Cell_OR_milo_by_group_",paste(group_combinations[[i]],collapse = "_"))),
                            width = 6.5,
                            height = 5.5,
                            res = 600,
                            formats = formats
                           )
    
    
        # 蜂群图展示celltype logFC变化
        da_results <- annotateNhoods(scmilo, da_results, coldata_col = cell_type_col)
        da_results$Celltype <- factor(da_results$Celltype,levels = rev(levels(subseurat@meta.data[[cell_type_col]])) )
        write.csv(da_results,file = file.path(output_dir,paste0("07_Celltype_Cell_Beeswarm_and_Boxplots_by_group_",paste(group_combinations[[i]],collapse = "_"),".csv")))
    
        plotDAbeeswarm(da_results, group.by = "Celltype", alpha = 1) +
            scale_color_gradient2(midpoint = 0, low = continuous_colors[1], mid = continuous_colors[2], high = continuous_colors[3], space = "Lab") +
            ggnewscale::new_scale_color() +
            geom_boxplot(aes(group = Celltype,colour = Celltype), outlier.shape = NA,show.legend = FALSE) +  # 显式指定 group
            scale_color_manual(values = celltype_colors) +
            geom_hline(yintercept = c(-1,1), linetype = "dashed",color = "grey50",linewidth = 0.8,alpha = 0.8)+
            labs(y = paste0(group_combinations[[i]][1],"  ←  Log(FC)  →  ",group_combinations[[i]][2]),x = "")+
            theme_classic(base_size = 18)+
            theme(axis.text = element_text(face = "bold"),
                  axis.text.y = element_text(,face = "bold"),
                  axis.title = element_text(face ="bold"),
                  plot.margin = margin(t = 10,r = 20,b = 5,l = 5,unit = "mm")
                 ) -> p2
        p2$layers[[1]]$aes_params$size <- 0.7  # 修改点大小
        basicR::save_figure(obj = p2,
                            filename = file.path(output_dir,paste0("07_Celltype_Cell_Beeswarm_and_Boxplots_by_group_",paste(group_combinations[[i]],collapse = "_"))),
                            width = 6.5,
                            height = 6,
                            res = 600,
                            formats = formats
                           )
        
    })
    
}