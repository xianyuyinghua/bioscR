#' Calculate CytoTRACE2 scores and plot an existing cell embedding
#'
#' Runs CytoTRACE2 on a Seurat object, caches the scored object, exports
#' cell metadata and scores, and plots CytoTRACE2_Relative on an existing
#' two-dimensional embedding.
#'
#' @param seurat A Seurat object with counts and feature identifiers accepted by
#' CytoTRACE2::cytotrace2(), metadata columns sample, group, and celltype_col,
#' and the requested reduction. Required, with no default. No assay is
#' explicitly selected by this wrapper; counts access follows CytoTRACE2.
#' Ignored when output_rds already exists.
#' @param species A single species name passed to CytoTRACE2::cytotrace2(),
#' "human" or "mouse". Defaults to "human". Not validated by this wrapper
#' and only used when computing a new result.
#' @param celltype_col Name of the cell-type metadata column, as a single
#' character string. Defaults to "Celltype". The column is checked and
#' exported during a fresh calculation; it does not define analysis groups
#' or plot colours. No factor requirement is imposed.
#' @param reduction Name of an existing reduction, as a single character string.
#' Defaults to "umap.harmony". Its cell.embeddings matrix must have exactly
#' two columns, because the code renames all embedding columns to dim1 and
#' dim2 after adding cell IDs. This is required in the cached object too.
#' @param colors Vector of valid plotting colours for the continuous score
#' gradient, ordered from low to high values. Defaults to
#' basicR::get_colors(number = 10.1, package = "viridis", name = "D").
#' Passed to scale_color_gradientn(); it is not a cell-type palette.
#' @param formats Character vector of figure formats accepted by
#' basicR::save_figure(), such as "pdf", "png", "tiff", and "svg".
#' Defaults to c("pdf", "png"). Appropriate graphics devices and fonts
#' must be available for the selected formats.
#' @param ncores Positive integer number of cores passed to
#' CytoTRACE2::cytotrace2(). Defaults to 5. Not validated by this wrapper;
#' only used for a fresh calculation. The environment must support the
#' parallel execution requested by CytoTRACE2.
#' @param output_dir Existing writable directory for the score CSV and figures.
#' Defaults to "./04_Pseudotime". The function does not create it. Relative
#' paths are resolved from the working directory.
#' @param output_rds Single path to the scored Seurat RDS cache. Defaults to
#' "./Rdata/04_Seurat_Celltype_CytoTRACE2.rds". For a new calculation its
#' parent directory must already exist and be writable. An existing file
#' must be a readable, compatible scored Seurat object; it is reused
#' without checking the supplied object or calculation parameters. This
#' path is independent of output_dir.
#'
#' @details
#' When output_rds does not exist, the function checks for the fixed metadata
#' columns sample and group and the column named by celltype_col. Missing
#' columns cause a message followed by q(save = "no", status = 0,
#' runLast = FALSE), which terminates the R session with status zero rather
#' than signalling an ordinary R error. Other inputs are not explicitly
#' validated.
#'
#' A fresh calculation calls CytoTRACE2::cytotrace2() with is_seurat = TRUE,
#' slot_type = "counts", batch_size = 10000, smooth_batch_size = 1000,
#' parallelize_models = TRUE, parallelize_smoothing = TRUE, and seed = 123.
#' These settings are fixed, including batch_size for objects with more than
#' 10000 cells. Any additional preprocessing or filtering is delegated to
#' CytoTRACE2; the wrapper performs no explicit gene or cell filtering before
#' the call. It does not split by sample, group, or cell type, perform pairwise
#' comparisons, define a comparison direction, or apply significance cutoffs.
#' It neither constructs a trajectory nor computes an embedding.
#'
#' After calculation, the CSV 01_CytoTRACE2_Score_in_Celltype.csv is written
#' under output_dir, followed by the scored object at output_rds. CSV rows
#' are cells, with row names included by write.csv(); columns are sample,
#' group, celltype_col, and every metadata column whose name contains the
#' case-sensitive substring "Cyto". This selection can include pre-existing
#' unrelated columns and is not restricted to a fixed set of score fields.
#' The RDS contains the complete object returned by CytoTRACE2, including
#' its metadata and reductions; this wrapper does not define additional
#' result components.
#'
#' If output_rds exists, readRDS() replaces the supplied seurat object.
#' Calculation, metadata checks, and CSV export are skipped, so a missing CSV
#' is not regenerated. Use a new cache path to recompute for changed inputs
#' or parameters. In both branches, the object must contain a numeric
#' CytoTRACE2_Relative metadata column and the requested reduction.
#'
#' For plotting, metadata and embedding rows are joined by cell ID using
#' merge(), retaining only matching IDs. No explicit score or missing-value
#' filter is applied. The plot uses CytoTRACE2_Relative unchanged with a
#' continuous colour gradient; it does not reverse the score or interpret it
#' as pseudotime. Point size is computed from the total object cell count by
#' bioscR::calculate_cell_point_size(), not the number of matched plot rows.
#' Figures use the stem 01_CytoTRACE2_Score_in_Celltype in output_dir, with
#' width = 6, height = 6, and res = 600 passed to basicR::save_figure().
#' Existing figure files, and the CSV on a fresh calculation, may be overwritten.
#'
#' The workflow uses CytoTRACE2, Seurat/SeuratObject, dplyr, tibble, ggplot2,
#' basicR, and utils, plus bioscR::calculate_cell_point_size(). The default
#' palette uses the viridis support in basicR. The function attaches no
#' packages. In the current package namespace, the unqualified calls to
#' %>%, select(), rownames_to_column(), ggplot(), geom_point(),
#' scale_color_gradientn(), and theme_minimal() lack explicit imports.
#' Attach dplyr, tibble, and ggplot2 before calling and avoid masking select().
#' CytoTRACE2 and tibble are not declared in DESCRIPTION; viridis is also
#' undeclared there and is an indirect dependency of the default palette.
#' The nested plotting helper has a viridis() default, but this is not
#' evaluated because the outer function always supplies colors explicitly.
#' No external script is sourced or invoked by this wrapper. Apart from an
#' optional existing RDS cache, input is supplied as an in-memory object;
#' any further installation or runtime requirements belong to CytoTRACE2
#' and the plotting dependencies.
#'
#' @return The return value of the final basicR::save_figure() call; its type and
#' visibility depend on the installed basicR implementation. No structured
#' analysis result, modified Seurat object, or plot is explicitly returned.
#' The primary results are the scored Seurat RDS, the cell metadata/score CSV
#' (on a fresh calculation only), and the embedding figures described above.
#' Read output_rds to retrieve the scored object.
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(dplyr)
#' library(tibble)
#' library(ggplot2)
#' library(bioscR)
#'
#' # Requires CytoTRACE2 and its runtime dependencies.
#' seu <- readRDS("seurat.rds")
#' # Supply counts, sample/group/Celltype metadata, and a two-column embedding.
#' stopifnot(all(c("sample", "group", "Celltype") %in% colnames(seu@meta.data)))
#' stopifnot(ncol(seu@reductions[["umap.harmony"]]@cell.embeddings) == 2L)
#' dir.create("cytotrace2_results", recursive = TRUE, showWarnings = FALSE)
#' dir.create("Rdata", recursive = TRUE, showWarnings = FALSE)
#' # Choose an unused cache path for a fresh calculation.
#' cache <- "Rdata/Seurat_CytoTRACE2.rds"
#' Trajectory_Cytotrace2(
#'   seurat = seu,
#'   species = "human",
#'   reduction = "umap.harmony",
#'   output_dir = "cytotrace2_results",
#'   output_rds = cache
#' )
#' scored_seurat <- readRDS(cache)
#' }
#' @export
Trajectory_Cytotrace2 <- function(seurat,
                                  species = "human",  #  mouse / human
                                  celltype_col = "Celltype",
                                  reduction = "umap.harmony",
                                  colors = basicR::get_colors(number = 10.1,package = "viridis",name = "D"), 
                                  formats = c("pdf","png"),
                                  ncores = 5,
                                  output_dir = "./04_Pseudotime",
                                  output_rds = "./Rdata/04_Seurat_Celltype_CytoTRACE2.rds"
                                 ){
    
    # ANSI 颜色代码
    red <- "\033[31m"
    green <- "\033[32m"
    yellow <- "\033[33m"
    blue <- "\033[34m"
    magenta <- "\033[35m"
    cyan <- "\033[36m"
    reset <- "\033[0m"

    if(!file.exists(output_rds)){
        # 检查参数和数据
        diff_colnames <- setdiff(c("sample","group",celltype_col),colnames(seurat@meta.data))
        if(length(diff_colnames) > 0){
            cat("\n",red,diff_colnames," Not in Seurat Cols.",reset,"\n")
            q(save = "no", status = 0, runLast = FALSE)
        }
        
        # cytotrace2 Calculation
        seurat <- CytoTRACE2::cytotrace2(seurat, 
                                         is_seurat = TRUE,   # 默认为非Seurat对象
                                         slot_type = "counts", # 当提供 Seurat对象时，默认是counts
                                         species = species,   # 默认是 mouse，human
                                         batch_size = 10000, # 进行样本抽样的细胞数，默认为10000，当细胞数超过10000时候，建议为NULL
                                         smooth_batch_size = 1000, # 进一步样本抽样的细胞数量。默认值为 1000，推荐输入数据大小 > 1K 个细胞时使用
                                         parallelize_models = TRUE, # 是否在多线程上并行运行模型预测函数（默认值为 TRUE）
                                         parallelize_smoothing = TRUE,
                                         ncores = ncores,
                                         seed = 123
                                        )
        
        df_score <-  seurat@meta.data %>% dplyr::select("sample","group",celltype_col,grep("Cyto",colnames(.)))
        write.csv(df_score,file = file.path(output_dir,paste0("01_CytoTRACE2_Score_in_Celltype",".csv")))
        saveRDS(seurat,file = output_rds)
        
    }else{
        seurat <- readRDS(output_rds)
    }
    
    Score_plot <- function(seurat,
                           reduction = "umap",
                           col = "CytoTRACE2_Relative",
                           colors = viridis(10),
                           point_size = bioscR::calculate_cell_point_size(cell_count = dim(seurat)[2])
                          ){
        
        dim_df <-  seurat@reductions[[reduction]]@cell.embeddings %>% as.data.frame() %>% rownames_to_column(var = "cell_id")
        colnames(dim_df) <- c('cell_id','dim1','dim2')
        color_df <- seurat@meta.data %>% select(col) %>% as.data.frame() %>% rownames_to_column(var = "cell_id")
        colnames(color_df) <- c("cell_id","color")
        df <- merge(color_df,dim_df,by = "cell_id")
        ggplot(df,aes(dim1,dim2,color = color)) +
          geom_point(size = point_size) +
          scale_color_gradientn(colors = colors) +
          labs(title = "",x = "",y = "",color = "CytoTRACE2\nScore") +
        theme_minimal()+
        theme(panel.grid = element_blank(),
              legend.title = element_text(size = 15, face = "bold", margin = margin(b = 20),hjust = 0.5),
              legend.text = element_text(size = 14),
              axis.text = element_blank(),
              axis.title.y = element_text(hjust = 0, vjust = 2)
             )
    }

    p <- Score_plot(seurat,reduction = reduction,col = "CytoTRACE2_Relative",colors = colors,point_size = bioscR::calculate_cell_point_size(cell_count = dim(seurat)[2]))
    basicR::save_figure(obj = p,
                        filename =  file.path(output_dir,"01_CytoTRACE2_Score_in_Celltype"),
                        width = 6,
                        height = 6,
                        res = 600,
                        formats = formats # "pdf","tiff","png","svg"
                       ) 
}