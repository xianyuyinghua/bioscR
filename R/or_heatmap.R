#' Plot cell-type enrichment odds ratios across groups
#'
#' Calculates Fisher exact-test odds ratios for cell-type and group combinations and returns an annotated heatmap with its underlying odds ratios and significance symbols.
#'
#' @param seurat A Seurat object containing the group and cell-type metadata
#'   columns. Supply this argument explicitly.
#' @param group Name of the metadata column defining groups.
#' @param celltype Name of the metadata column containing cell-type
#'   annotations.
#' @param row_col With "celltype_group", cell types form rows and groups form
#'   columns. Any other value transposes the displayed matrix and returned
#'   data.
#' @param colors Three colours for the lower endpoint, log2(OR + 1) = 1, and
#'   upper endpoint of the heatmap scale.
#'
#' @details
#' The function builds a cell-type by group count table and replaces zero
#' counts with one before testing. For each combination, a two-sided Fisher
#' exact test compares that combination with the remaining cells. P values are
#' adjusted together using the Benjamini-Hochberg method. Symbols are shown
#' only when OR > 1.5 or OR < 0.5: * for adjusted P < 0.05, ** for adjusted P <
#' 0.01, and *** for adjusted P < 0.001.
#'
#' Heatmap colours encode log2(OR + 1), with endpoints at its minimum and 99th
#' percentile and a middle colour at 1. Rows and columns are not clustered.
#' Factor levels determine their order; non-factor cell types use first-
#' occurrence order, while non-factor groups are ordered by group_new_levels().
#' Use at least two observed cell types and groups, without unused factor
#' levels.
#'
#' The implementation uses data.table, plyr, tidyr, circlize, ComplexHeatmap,
#' grid, gridExtra, and cowplot. Unqualified grid and heatmap drawing functions
#' must be available; group_new_levels() is required when group is not a
#' factor. This is a pooled cell-level comparison and does not model biological
#' replicates.
#'
#' @return A list with plot, a cowplot/ggplot heatmap, and data, a list
#'   containing or_df (untransformed odds ratios) and p_df (a matrix of
#'   significance strings, not numeric P values). The orientation of both data
#'   components follows row_col.
#'
#' @examples
#' \dontrun{
#' library(ComplexHeatmap)
#' library(grid)
#' seu <- readRDS("seurat.rds")
#' seu$group <- droplevels(factor(seu$group))
#' seu$Celltype <- droplevels(factor(seu$Celltype))
#' result <- or_heatmap(seurat = seu)
#' result$plot
#' result$data$or_df
#' result$data$p_df
#' }

or_heatmap <- function(seurat = seurat,       # seurat 对象
                       group = "group",       # 样本分组或其他分组
                       celltype = "Celltype",    # 细胞类型列
                       row_col = "celltype_group",   # 图和数据的行列设置
                       colors = c("#0073C2FF","white","#CD534CFF")
                      ){
    
    # 计算单细胞数据组间的各细胞类型的富集比和差异并出热图

    # 检查 group 和 celltype 是否都是因子型并转化
    if(is.null(levels(seurat@meta.data[[celltype]]))){
        seurat@meta.data[[celltype]] <- factor(seurat@meta.data[[celltype]],levels = unique(seurat@meta.data[[celltype]]))
    }
    if(is.null(levels(seurat@meta.data[[group]]))){
        seurat@meta.data[[group]] <- factor(seurat@meta.data[[group]],levels = group_new_levels(vector =seurat@meta.data[[group]]) )
    }
    
    do.tissueDist <- function(cellInfo.tb = cellInfo.tb,
                              meta.cluster = cellInfo.tb$meta.cluster,
                              colname.patient ="patient",
                              loc = cellInfo.tb$loc,
                              pdf.width=3,
                              pdf.height=5,
                              verbose=0
                             ){
        
            # 不同组间细胞亚群分布差异
        
            ##input data 
            library(data.table)
            cellInfo.tb = data.table(cellInfo.tb)
            cellInfo.tb$meta.cluster =as.character(meta.cluster)
            
            if(is.factor(loc)){
                cellInfo.tb$loc = loc
            }else{
                cellInfo.tb$loc =as.factor(loc)
            }
            
            loc.avai.vec <- levels(cellInfo.tb[["loc"]])
            count.dist <-unclass(cellInfo.tb[,table(meta.cluster,loc)])[,loc.avai.vec]
            freq.dist <- sweep(count.dist,1,rowSums(count.dist),"/")
            freq.dist.bin <-floor(freq.dist *100/10)
            
            {
            count.dist[count.dist == 0] <- 1  # 把为0的变为1
            count.dist.melt.ext.tb <- test.dist.table(count.dist)
            p.dist.tb <- dcast(count.dist.melt.ext.tb,rid~cid,value.var="p.value")
            OR.dist.tb <- dcast(count.dist.melt.ext.tb,rid~cid,value.var="OR")
            OR.dist.mtx <-as.matrix(OR.dist.tb[,-1])
            rownames(OR.dist.mtx)<- OR.dist.tb[[1]]
            }
            
            if(verbose == 1){
                return(list("count.dist.melt.ext.tb"=count.dist.melt.ext.tb,
                            "p.dist.tb"=p.dist.tb,
                            "OR.dist.tb"=OR.dist.tb,
                            "OR.dist.mtx"=OR.dist.mtx))
            }else{
                return(OR.dist.mtx)
            }
    }
    
    
    test.dist.table <- function(count.dist,min.rowSum=0){
                count.dist <- count.dist[rowSums(count.dist)>=min.rowSum,,drop=F]
                sum.col <- colSums(count.dist)
                sum.row <- rowSums(count.dist)
                count.dist.tb <-as.data.frame(count.dist)
                setDT(count.dist.tb,keep.rownames=T)
                count.dist.melt.tb <- melt(count.dist.tb,id.vars="rn")
                colnames(count.dist.melt.tb)<-c("rid","cid","count")
                count.dist.melt.ext.tb <- as.data.table(plyr::ldply(seq_len(nrow(count.dist.melt.tb)),function(i){
                            this.row <- count.dist.melt.tb$rid[i]
                            this.col <- count.dist.melt.tb$cid[i]
                            this.c <- count.dist.melt.tb$count[i]
                                other.col.c <- sum.col[this.col]-this.c
                            this.m <- matrix(c(this.c,
                                                   sum.row[this.row]-this.c,
                                                   other.col.c,
                            sum(sum.col)-sum.row[this.row]-other.col.c),
                                                 ncol=2)
                                res.test <- fisher.test(this.m)
                                data.frame(rid=this.row,
                                           cid=this.col,
                                           p.value=res.test$p.value,
                                           OR=res.test$estimate)
                        }))
                count.dist.melt.ext.tb <- merge(count.dist.melt.tb,count.dist.melt.ext.tb, by=c("rid","cid"))
                count.dist.melt.ext.tb$adj.p.value <- p.adjust(count.dist.melt.ext.tb$p.value, method = "BH")
                return(count.dist.melt.ext.tb)
    }

    # 该分析只需要 分组信息 和 cluster/celltype结果 ，也就是meta.data 中的两列信息。
    # 修改我们输入矩阵的名字来适配函数 
    meta <- seurat@meta.data
    meta$group <- meta[[group]]
    meta$celltype <- meta[[celltype]]

    #主分析函数
    OR.immune.list <- do.tissueDist(cellInfo.tb = meta,
                                    meta.cluster = meta$celltype,
                                    loc = meta$group,
                                    pdf.width = 4,
                                    pdf.height = 8,
                                    verbose=1
                                    )
    #可视化
    # a 存OR值结果
    a=OR.immune.list[["OR.dist.tb"]]
    a <-as.data.frame(a,check.names = F)
    rownames(a) <- a$rid
    a <- a[,-1]
    a <- na.omit(a)
    
    # b存P值结果
    b <- OR.immune.list$count.dist.melt.ext.tb[,c(1,2,6)]
    b <- tidyr::spread(b,key ="cid", value ="adj.p.value")
    b <- data.frame(b[,-1],row.names = b$rid,check.names = F)
    b <- b[rownames(a),]
    
    #考虑到OR值在文献中定义的0.5 和 1.5 值，这里设置bk参数
    p_df = ifelse(b >= 0.05&(a > 1.5|a < 0.5),"",
               #ifelse(b<0.0001&(a>1.5|a<0.5),"****",
                      ifelse(b < 0.001 & (a > 1.5| a < 0.5),"***",
                             ifelse(b < 0.01 & (a > 1.5| a < 0.5),"**",
                                    ifelse(b < 0.05 & (a > 1.5 | a < 0.5),"*","")
                                   )
                            )
                     #)
                 )
    # 排序
    a <- a[levels(meta$celltype),levels(meta$group)] 
    p_df <- p_df[levels(meta$celltype),levels(meta$group)]

    if(row_col != "celltype_group"){
        a <- t(a); p_df <- t(p_df)
    }

    a_log2 <- log2(a + 1) 
    
    # 定义颜色映射
    value_top <- quantile(as.vector(as.matrix(a_log2)), 0.99, na.rm = TRUE)
    value_bottom <- quantile(as.vector(as.matrix(a_log2)), 0, na.rm = TRUE)
    col_fun <- circlize::colorRamp2(c(value_bottom, 1, value_top),colors)
    
    legend <- ComplexHeatmap::Legend(
        col_fun = col_fun,  # 颜色映射
        title = "log2(OR+1)",  # 图例标题
        title_gp = gpar(fontsize = 14, col = "black"),
        at = pretty(c(value_bottom, value_top), n = 4),  # 图例刻度位置
        labels_gp = gpar(fontsize = 12),  # 图例字体大小
        legend_height = unit(5, "cm"),  # 图例高度
        grid_width = unit(0.6, "cm"),  # 图例宽度
        direction = "vertical",  # 设置图例方向为垂直
        title_position = "topcenter",
        title_gap = unit(0.5, "cm")
        )

    ComplexHeatmap::Heatmap(
          a_log2, 
          col = col_fun,  
          cluster_columns = FALSE,  
          cluster_rows = FALSE,     
          show_heatmap_legend = FALSE, 
          row_names_side = "left",  
          column_names_side = "bottom",
          column_names_rot = 0,
          column_names_centered = TRUE,
          row_names_gp = gpar(fontsize = 12,fontface = "bold"),  
          column_names_gp = gpar(fontsize = 13,fontface = "bold"),
          cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(p_df[i, j], x, y, gp = gpar(fontsize = 10, col = "black"))
          }
    ) -> p

    # 合并为一个图
    # 抓取主图和图例的 grob（主图支持 padding）
    p_grob <- grid::grid.grabExpr(draw(p, padding = unit(c(5, 5, 5, 0), "mm")))
    
    # 图例对象不能加 padding 参数，去掉即可
    legend_grob <- grid::grid.grabExpr(draw(legend, y = unit(0.8, "npc"),just = c("center", "top")))
    
    # 包装到 frameGrob 里，缩小内部尺寸实现边距
    legend_padded <- grid::frameGrob()
    legend_padded <- grid::packGrob(
      legend_padded, legend_grob,
      side = "center",
      width = unit(1, "npc") - unit(10, "mm"),
      height = unit(1, "npc") - unit(10, "mm")
    )
    
    # 合并 grob：主图 + 图例（水平布局）
    combined_grob <- gridExtra::arrangeGrob(
      grobs = list(p_grob, legend_padded),
      ncol = 2,
      widths = c(0.8, 0.2)
    )

    ggplot_like <- cowplot::ggdraw() + cowplot::draw_grob(combined_grob)

    return(list(plot = ggplot_like,data = list(or_df = a,p_df = p_df)))
}