#' Infer a Monocle3 trajectory using an existing Seurat embedding
#'
#' Builds and caches a Monocle3 cell_data_set from RNA counts and an existing
#' Seurat embedding, orders cells from selected roots, tests graph-associated
#' genes, and exports pseudotime tables and trajectory figures.
#'
#' @param seurat A Seurat object, required with no default. A fresh analysis uses
#' SeuratObject::LayerData() with assay = "RNA" and layer = "counts",
#' cell metadata with matching cell names, and the requested embedding.
#' Supply a compatible counts layer and enough cells and features for the
#' fixed 30-dimensional preprocessing. All counts cells must occur in the
#' metadata and embedding. Ignored when output_rds exists.
#' @param reduction Single name of a Seurat reduction. Defaults to "umap.harmony".
#' It must contain at least two columns; only the first two are used,
#' reordered to the CDS cell order and stored under the name "UMAP".
#' No new UMAP is computed. Only used for a fresh analysis.
#' @param celltype_col Single name of a character or factor cell-type metadata
#' column. Defaults to "Celltype". Used for candidate-root annotation and
#' cell-type plotting, including plotting a cached object. With
#' point_start_use = FALSE, actual root-cell selection always uses the
#' literal Celltype column, irrespective of this argument. Use "Celltype"
#' for that workflow, or ensure both columns have consistent labels.
#' @param celltype_start Single starting cell-type label. Defaults to NULL,
#' which triggers readline() during a fresh analysis. Provide an observed,
#' non-missing label. It identifies candidate leaf nodes for display; with
#' point_start_use = FALSE, all cells with this label in the literal
#' Celltype column are passed as root_cells. No explicit label validation
#' is performed. Ignored when output_rds exists.
#' @param point_start_use Single non-missing logical value. Defaults to TRUE.
#' For a fresh analysis, TRUE always prompts for comma-separated principal
#' node IDs, for example "Y_1,Y_2", and passes them as root_pr_nodes.
#' There is no argument for supplying these IDs programmatically. Input
#' is split on commas without trimming whitespace or checking candidates.
#' FALSE selects root cells by celltype_start without a node prompt.
#' For non-interactive use, set FALSE and supply celltype_start explicitly.
#' Ignored when output_rds exists.
#' @param celltype_colors Vector of valid colours with enough entries for the
#' cell-type categories in the final labelled trajectory plot. Defaults to
#' basicR::get_colors(number = 50.1). Named colours can match category
#' labels. This palette is not applied to the candidate-root figure or
#' the pseudotime-coloured figures.
#' @param output_dir Existing writable directory for CSV tables and figures.
#' Defaults to "./04_Pseudotime". The function does not create it. Relative
#' paths use the working directory. On successful completion, files with
#' names containing "02_Monocle3_RootNodes" are deleted from this directory,
#' including pre-existing matching files.
#' @param output_rds Single path to the CDS cache, independent of output_dir.
#' Defaults to "./Rdata/04_Monocle3.cds". Despite the extension, this is an
#' RDS file written with saveRDS(). Its parent directory must exist and be
#' writable for a fresh analysis. An existing file must contain a readable,
#' compatible ordered cell_data_set with UMAP coordinates, graph information,
#' and the requested cell-type metadata; it is reused without validation
#' against the supplied Seurat object or trajectory parameters.
#' @param formats Character vector of formats accepted by basicR::save_figure(),
#' such as "pdf", "png", "tiff", or "svg". Defaults to c("pdf", "png").
#' Used for all figures, including the temporary root-selection figure.
#' The corresponding graphics devices and fonts must be available.
#'
#' @details
#' For a fresh analysis, counts come from the fixed RNA assay and counts
#' layer. Metadata are reordered to counts columns, and gene_short_name is
#' set to the counts row names. The wrapper performs no explicit gene or
#' cell filtering before constructing the CDS. It calls preprocess_cds()
#' with num_dim = 30 and norm_method = "log", copies the first two embedding
#' coordinates into UMAP, calls cluster_cells() with reduction_method =
#' "UMAP" and resolution = 0.1, then learn_graph() with use_partition = TRUE,
#' close_loop = FALSE, and learn_graph_control = list(minimal_branch_len = 30,
#' prune_graph = TRUE). There is no explicit random seed in this wrapper.
#'
#' Root candidates are degree-one principal-graph nodes with at least one
#' assigned cell of celltype_start, using the closest-vertex mapping and the
#' column named by celltype_col. Counts and within-node proportions are
#' computed by node and cell type. No minimum proportion, majority rule,
#' or highest-proportion selection is applied. The temporary figure
#' 02_Monocle3_RootNodes highlights these candidates. Manually entered nodes
#' are not restricted to the displayed candidates by this wrapper. With
#' point_start_use = FALSE, root cells are selected from the literal Celltype
#' column across the whole CDS, not only from candidate leaf nodes.
#'
#' Pseudotime is obtained from order_cells() using the selected roots; the
#' wrapper does not reverse it, select an endpoint, or verify root coverage
#' of every partition. Missing or non-finite pseudotime values are not
#' explicitly filtered before CSV export. There are no sample/group
#' comparisons, treatment contrasts, or fold-change directions. Clustering
#' and graph learning use Monocle3 partitions; cell-type labels guide root
#' selection and colouring rather than separate analyses.
#'
#' The ordered CDS is saved before gene testing. If output_rds already
#' exists, construction, preprocessing, graph learning, root prompts, and
#' root selection are skipped. The cached object replaces the supplied
#' input without checking parameter consistency. Gene testing, CSV export,
#' final plotting, and root-figure cleanup still run on every call. Use an
#' unused cache path to change inputs or roots. A failure after caching can
#' leave a usable cache alongside incomplete outputs.
#'
#' All CDS genes are split in row order into batches of at most 5000 genes,
#' not 1000 as stated in an internal comment. Each batch is passed to
#' graph_test() with neighbor_graph = "principal_graph" and cores = 1.
#' The batch results are row-bound; q_value values are taken from each
#' separate call without a new correction across all batches. The exported
#' subset requires status == "OK", q_value < 0.05, and morans_I > 0.1,
#' and is sorted by decreasing morans_I. Missing filter values are dropped
#' by dplyr::filter(). These are graph-associated genes, not a pairwise
#' differential-expression contrast or a fitted pseudotime regression.
#'
#' The file 02_Monocle3_Deg_PseudoTime_sig.csv contains a gene column derived
#' from the combined result row names and all graph_test() result columns,
#' including status, q_value, and morans_I, for retained genes. The full
#' unfiltered result is neither saved separately nor added to the cached CDS.
#' The file 02_Monocle3_pseudotime.csv contains one column named
#' monocle3_pseudotime with cell IDs as row names. Both CSV files include
#' row names through write.csv(). The cache retains counts, cell and gene
#' metadata, preprocessing results, UMAP coordinates, clusters/partitions,
#' and principal-graph and ordering information in the cell_data_set.
#' Retrieve it with readRDS(output_rds).
#'
#' Four final figure stems are written under output_dir:
#' 03_Monocle3_pseudotime_noNode_noTrajectory (pseudotime colours, no graph
#' or node labels); 03_Monocle3_pseudotime_noNode_Trajectory (pseudotime
#' colours and graph, no node labels); 03_Monocle3_pseudotime_Node_Trajectory
#' (pseudotime colours, graph, and root labels only); and
#' 03_Monocle3_pseudotime_Label_Trajectory (cell-type colours and graph,
#' with principal-node IDs for branch, leaf, and root nodes, not cell-type
#' text labels or every principal node). All figures, including the root
#' candidate figure, pass width = 7, height = 6, and res = 600 to
#' basicR::save_figure(). Existing output tables and figures may be overwritten.
#'
#' Plotting uses a locally modified plot_cells() helper, not a direct call
#' to monocle3::plot_cells(). The numeric colour scale uses viridis option
#' "C"; the helper adds overlapping point layers. It reads Monocle3 internal
#' slots, so compatibility with the installed Monocle3 version is required.
#' On successful completion, list.files() with the unanchored pattern
#' "02_Monocle3_RootNodes" identifies files to remove. This also runs with
#' a cached CDS and can remove unrelated files containing that substring.
#' Earlier errors prevent this final cleanup.
#'
#' The function attaches dplyr, basicR, monocle3, Seurat, ggplot2, and
#' SingleCellExperiment even when reading a cache. Executed paths also use
#' SeuratObject, SummarizedExperiment, assertthat, methods, igraph, tibble,
#' viridis, ggrepel, stats, and utils. Packages must be installed with their
#' runtime dependencies. The local plotting helper also references Matrix,
#' reshape2, and ggrastr in gene-expression or rasterisation branches, but
#' those branches are not used by this wrapper. DESCRIPTION currently omits
#' monocle3, assertthat, methods, igraph, tibble, viridis, and ggrepel, as well
#' as the unused-branch references Matrix, reshape2, and ggrastr. No external
#' script is sourced or launched. Input is the in-memory Seurat object or
#' an existing RDS cache. A fresh run with default root-selection settings
#' requires an interactive console and valid root input; invalid or empty
#' input is not handled by this wrapper.
#'
#' @return An invisible logical vector returned by file.remove() if matching root
#' figure files are found: one element per file, TRUE for successful removal
#' and FALSE for failure. If no matching files exist, the final conditional
#' returns invisible NULL. The function does not return the CDS, gene-test
#' table, pseudotime vector, or plots. Its analysis results are the RDS cache,
#' two CSV tables, and figures described above.
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(bioscR)
#' seu <- readRDS("seurat.rds")
#' # Requires RNA counts, at least two embedding columns, and Celltype metadata.
#' # Replace Progenitor with an observed starting cell-type label.
#' stopifnot("Celltype" %in% colnames(seu@meta.data))
#' stopifnot("Progenitor" %in% seu$Celltype)
#' stopifnot(ncol(Embeddings(seu, "umap.harmony")) >= 2L)
#' dir.create("monocle3_results", recursive = TRUE, showWarnings = FALSE)
#' dir.create("Rdata", recursive = TRUE, showWarnings = FALSE)
#' cache <- "Rdata/Monocle3.cds"
#' # Use an unused cache path for a fresh trajectory with these roots.
#' set.seed(123)
#' Trajectory_Monocle3(
#'   seurat = seu,
#'   celltype_col = "Celltype",
#'   celltype_start = "Progenitor",
#'   point_start_use = FALSE,
#'   output_dir = "monocle3_results",
#'   output_rds = cache
#' )
#' cds <- readRDS(cache)
#' # For manual node selection in an interactive console, use
#' # point_start_use = TRUE and a new cache path; inspect the root figure
#' # when prompted and enter valid node IDs separated by commas, without spaces.
#' }
#' @export
Trajectory_Monocle3 <- function(seurat,
                                reduction = "umap.harmony",
                                celltype_col = "Celltype",
                                celltype_start = NULL,
                                point_start_use = TRUE, # TRUE/FALSE
                                celltype_colors = basicR::get_colors(number = 50.1),
                                output_dir = "./04_Pseudotime",
                                output_rds = "./Rdata/04_Monocle3.cds",
                                formats = c("pdf","png")
                               ){
    
    suppressPackageStartupMessages({
        library(dplyr)
        library(basicR)
        library(monocle3)
        library(Seurat)
        library(ggplot2)
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

    ############################ Function ##################################
    # 修改后，点的边框色和填充色一致（第360行）
    # 修改过 label_cell_groups 的标签文本（原始脚本是以cell_group分组的）
    plot_cells <- function(cds,
                           x=1,
                           y=2,
                           reduction_method = c("UMAP", "tSNE", "PCA", "LSI", "Aligned"),
                           color_cells_by="cluster",
                           group_cells_by="cluster",
                           genes=NULL,
                           show_trajectory_graph=TRUE,
                           trajectory_graph_color="grey28",
                           trajectory_graph_segment_size=0.75,
                           norm_method = c("log", "size_only"),
                           label_cell_groups = TRUE,
                           label_groups_by_cluster=TRUE,
                           group_label_size=2,
                           labels_per_group=1,
                           label_branch_points=TRUE,
                           label_roots=TRUE,
                           label_leaves=TRUE,
                           graph_label_size=2,
                           cell_size=0.35,
                           cell_stroke= I(cell_size / 2),
                           alpha = 1,
                           min_expr=0.1,
                           rasterize=FALSE,
                           scale_to_range=TRUE,
                           label_principal_points = FALSE
                          ){
      feature_label <- min_val_for_feature <- max_val_for_feature <- NULL # no visible binding
      cell_group <- cell_color <- cells_in_cluster <- per <- NULL # no visible binding
      reduction_method <- match.arg(reduction_method)
      assertthat::assert_that(methods::is(cds, "cell_data_set"))
      assertthat::assert_that(!is.null(SingleCellExperiment::reducedDims(cds)[[reduction_method]]),
                              msg = paste("No dimensionality reduction for",
                                          reduction_method, "calculated.",
                                          "Please run reduce_dimension with",
                                          "reduction_method =", reduction_method,
                                          "before attempting to plot."))
      low_dim_coords <- SingleCellExperiment::reducedDims(cds)[[reduction_method]]
      assertthat::assert_that(ncol(low_dim_coords) >=max(x,y),
                              msg = paste("x and/or y is too large. x and y must",
                                          "be dimensions in reduced dimension",
                                          "space."))
      if(!is.null(color_cells_by)) {
        assertthat::assert_that(color_cells_by %in% c("cluster", "partition",
                                                      "pseudotime") |
                                  color_cells_by %in% names(colData(cds)),
                                msg = paste("color_cells_by must one of",
                                            "'cluster', 'partition', 'pseudotime,",
                                            "or a column in the colData table."))
    
        if(color_cells_by == "pseudotime") {
          tryCatch({pseudotime(cds, reduction_method = reduction_method)},
                   error = function(x) {
                     stop("No pseudotime for ", reduction_method,
                                " calculated. Please run order_cells with ",
                                "reduction_method = ", reduction_method,
                                " before attempting to color by pseudotime.")})
    
        }
      }
      assertthat::assert_that(!is.null(color_cells_by) || !is.null(markers),
                              msg = paste("Either color_cells_by or markers must",
                                          "be NULL, cannot color by both!"))
    
      norm_method = match.arg(norm_method)
    
      assertthat::assert_that((group_cells_by %in% c("cluster", "partition")) ||
                              (group_cells_by %in% names(SummarizedExperiment::colData(cds))),
                              msg = paste("group_cells_by must be one of",
                              "'cluster', 'partition', or a column in the colData table."))
    
      assertthat::assert_that(!is.null(color_cells_by) || !is.null(genes),
                              msg = paste("Either color_cells_by or genes must be",
                                          "NULL, cannot color by both!"))
    
      if (show_trajectory_graph &&
          is.null(principal_graph(cds)[[reduction_method]])) {
        message("No trajectory to plot. Has learn_graph() been called yet?")
        show_trajectory_graph = FALSE
      }
      if (label_principal_points &&
          is.null(principal_graph(cds)[[reduction_method]])) {
        message("Cannot label principal points when no trajectory to plot. Has learn_graph() been called yet?")
        label_principal_points = FALSE
      }
    
      if (label_principal_points) {
        label_branch_points <- FALSE
        label_leaves <- FALSE
        label_roots <- FALSE
      }
    
    
      gene_short_name <- NA
      sample_name <- NA
      #sample_state <- colData(cds)$State
      data_dim_1 <- NA
      data_dim_2 <- NA
      if (rasterize){
        plotting_func <- ggrastr::geom_point_rast
      }else{
        plotting_func <- ggplot2::geom_point
      }
    
      S_matrix <- SingleCellExperiment::reducedDims(cds)[[reduction_method]]
      data_df <- data.frame(S_matrix[,c(x,y)])
    
      colnames(data_df) <- c("data_dim_1", "data_dim_2")
      data_df$sample_name <- row.names(data_df)
    
      data_df <- as.data.frame(cbind(data_df, colData(cds)))
      if (group_cells_by == "cluster") {
        data_df$cell_group <-
          tryCatch({clusters(cds,
                             reduction_method = reduction_method)[
                               data_df$sample_name]},
                   error = function(e) {NULL})
      } else if (group_cells_by == "partition") {
        data_df$cell_group <-
          tryCatch({partitions(cds,
                               reduction_method = reduction_method)[
                                 data_df$sample_name]},
                   error = function(e) {NULL})
      } else{
        data_df$cell_group <- SummarizedExperiment::colData(cds)[data_df$sample_name, group_cells_by]
      }
    
      if (color_cells_by == "cluster"){   #---------------------------------------------------------------
        data_df$cell_color <-
          tryCatch({clusters(cds,
                             reduction_method = reduction_method)[
                               data_df$sample_name]},
                   error = function(e) {NULL})
      } else if (color_cells_by == "partition") {
        data_df$cell_color <-
          tryCatch({partitions(cds,
                               reduction_method = reduction_method)[
                                 data_df$sample_name]},
                   error = function(e) {NULL})
      } else if (color_cells_by == "pseudotime") {
        data_df$cell_color <-
          tryCatch({pseudotime(cds,
                               reduction_method = reduction_method)[
                                 data_df$sample_name]}, error = function(e) {NULL})
      } else{
        data_df$cell_color <- colData(cds)[data_df$sample_name, color_cells_by]
      }
    
      ## Graph info
      if (show_trajectory_graph) {
    
        ica_space_df <- t(cds@principal_graph_aux[[reduction_method]]$dp_mst) %>%
          as.data.frame() %>%
          dplyr::select(prin_graph_dim_1 = {{x}}, prin_graph_dim_2 = {{y}}) %>%
          dplyr::mutate(sample_name = rownames(.),
                        sample_state = rownames(.))
    
        dp_mst <- cds@principal_graph[[reduction_method]]
    
        edge_df <- dp_mst %>%
          igraph::as_data_frame() %>%
          dplyr::select(source = "from", target = "to") %>%
          dplyr::left_join(ica_space_df %>%
                             dplyr::select(
                               source="sample_name",
                               source_prin_graph_dim_1="prin_graph_dim_1",
                               source_prin_graph_dim_2="prin_graph_dim_2"),
                           by = "source") %>%
          dplyr::left_join(ica_space_df %>%
                             dplyr::select(
                               target="sample_name",
                               target_prin_graph_dim_1="prin_graph_dim_1",
                               target_prin_graph_dim_2="prin_graph_dim_2"),
                           by = "target")
      }
    
      ## Marker genes
      markers_exprs <- NULL
      expression_legend_label <- NULL
      if (!is.null(genes)) {
        if (!is.null(dim(genes)) && dim(genes)[[2]] >= 2){
          markers = unlist(genes[,1], use.names=FALSE)
        } else {
          markers = genes
        }
    	markers_rowData <- rowData(cds)[(rowData(cds)$gene_short_name %in% markers) |
    							        (row.names(rowData(cds)) %in% markers),,drop=FALSE]
    	markers_rowData <- as.data.frame(markers_rowData)
        if (nrow(markers_rowData) == 0) {
          stop("None of the provided genes were found in the cds")
        }
        if (nrow(markers_rowData) >= 1) {
          cds_exprs <- SingleCellExperiment::counts(cds)[row.names(markers_rowData), ,drop=FALSE]
          # assertthat::assert_that(!is.null(size_factors(cds_exprs)))
          cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/size_factors(cds))
    
          if (!is.null(dim(genes)) && dim(genes)[[2]] >= 2){
            #genes = as.data.frame(genes)
            #row.names(genes) = genes[,1]
            #genes = genes[row.names(cds_exprs),]
    
            agg_mat = as.matrix(aggregate_gene_expression(cds, genes, norm_method=norm_method, gene_agg_fun="mean", scale_agg_values=FALSE))
            markers_exprs = agg_mat
            markers_exprs <- reshape2::melt(markers_exprs)
            colnames(markers_exprs)[1:2] <- c('feature_id','cell_id')
            if (is.factor(genes[,2]))
              markers_exprs$feature_id = factor(markers_exprs$feature_id,
                                                levels=levels(genes[,2]))
    
            markers_exprs$feature_label <- markers_exprs$feature_id
            norm_method = "size_only"
            expression_legend_label = "Expression score"
          } else {
            cds_exprs@x = round(10000*cds_exprs@x)/10000
            markers_exprs = matrix(cds_exprs, nrow=nrow(markers_rowData))
            colnames(markers_exprs) = colnames(SingleCellExperiment::counts(cds))
            row.names(markers_exprs) = row.names(markers_rowData)
            markers_exprs <- reshape2::melt(markers_exprs)
            colnames(markers_exprs)[1:2] <- c('feature_id','cell_id')
            markers_exprs <- merge(markers_exprs, markers_rowData,
                                   by.x = "feature_id", by.y="row.names")
            if (is.null(markers_exprs$gene_short_name)) {
              markers_exprs$feature_label <-
                as.character(markers_exprs$feature_id)
            } else {
              markers_exprs$feature_label <-
                as.character(markers_exprs$gene_short_name)
            }
    
            markers_exprs$feature_label <- ifelse(is.na(markers_exprs$feature_label) | !as.character(markers_exprs$feature_label) %in% markers,
                                                  as.character(markers_exprs$feature_id),
                                                  as.character(markers_exprs$feature_label))
    
            markers_exprs$feature_label <- factor(markers_exprs$feature_label,
                                                  levels = markers)
            if (norm_method == "size_only")
              expression_legend_label = "Expression"
            else
              expression_legend_label = "log10(Expression)"
          }
    
          if (scale_to_range){
    				markers_exprs = dplyr::group_by(markers_exprs, feature_label) %>%
    				  dplyr::mutate(max_val_for_feature = max(value),
    								min_val_for_feature = min(value)) %>%
    				  dplyr::mutate(value = 100 * (value - min_val_for_feature) / (max_val_for_feature - min_val_for_feature))
    				expression_legend_label = "% Max"
    			  }
    			}
    		  }
    
          if (label_cell_groups && is.null(color_cells_by) == FALSE){
    			if (is.null(data_df$cell_color)){
    			  if (is.null(genes)){
    				message(color_cells_by, " not found in colData(cds), cells will ",
    							  "not be colored")
    			  }
    			  text_df = NULL
    			  label_cell_groups = FALSE
    			}else{
    			  if(is.character(data_df$cell_color) || is.factor(data_df$cell_color)) {
    
    				if (label_groups_by_cluster && is.null(data_df$cell_color) == FALSE){
    				  text_df = data_df %>%
    					dplyr::group_by(cell_color) %>%
    					dplyr::mutate(cells_in_cluster= dplyr::n()) %>%
    					dplyr::group_by(cell_color, .add=TRUE) %>%
    					dplyr::mutate(per=dplyr::n()/cells_in_cluster)
    				  median_coord_df = text_df %>%
    					dplyr::summarize(fraction_of_group = dplyr::n(),
    									 text_x = stats::median(x = data_dim_1),
    									 text_y = stats::median(x = data_dim_2))
    				  text_df = suppressMessages(text_df %>% dplyr::select(per) %>%
    											   dplyr::distinct())
    				  text_df = suppressMessages(dplyr::inner_join(text_df,
    															   median_coord_df))
    				  text_df = text_df %>% dplyr::group_by(cell_color) %>%
    					dplyr::top_n(labels_per_group, per)
    				} else {
    				  text_df = data_df %>% dplyr::group_by(cell_color) %>%
    					dplyr::mutate(per=1)
    				  median_coord_df = text_df %>%
    					dplyr::summarize(fraction_of_group = dplyr::n(),
    									 text_x = stats::median(x = data_dim_1),
    									 text_y = stats::median(x = data_dim_2))
    				  text_df = suppressMessages(text_df %>% dplyr::select(per) %>%
    											   dplyr::distinct())
    				  text_df = suppressMessages(dplyr::inner_join(text_df,
    															   median_coord_df))
    				  text_df = text_df %>% dplyr::group_by(cell_color) %>%
    					dplyr::top_n(labels_per_group, per)
    				}
    
    				text_df$label = as.character(text_df %>% dplyr::pull(cell_color))
    				# I feel like there's probably a good reason for the bit below, but I
    				# hate it and I'm killing it for now.
    				# text_df$label <- paste0(1:nrow(text_df))
    				# text_df$process_label <- paste0(1:nrow(text_df), '_',
    				# as.character(as.matrix(text_df[, 1])))
    				# process_label <- text_df$process_label
    				# names(process_label) <- as.character(as.matrix(text_df[, 1]))
    				# data_df[, group_by] <-
    				#  process_label[as.character(data_df[, group_by])]
    				# text_df$label = process_label
    			  } else {
    				message("Cells aren't colored in a way that allows them to ",
                            "be grouped.")
    				text_df = NULL
    				label_cell_groups = FALSE
    			  }
    			}
    		  }
    
    		  if (!is.null(markers_exprs) && nrow(markers_exprs) > 0){
    			data_df <- merge(data_df, markers_exprs, by.x="sample_name",
    							 by.y="cell_id")
    			data_df$value <- with(data_df, ifelse(value >= min_expr, value, NA))
    			ya_sub <- data_df[!is.na(data_df$value),]
    			na_sub <- data_df[is.na(data_df$value),]
    			if(norm_method == "size_only"){
    			  g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2)) +
    				plotting_func(aes(data_dim_1, data_dim_2), size=I(cell_size),
    							  stroke = I(cell_stroke), color = "grey80", alpha = alpha,
    							  data = na_sub) +
    				plotting_func(aes(color=value), size=I(cell_size),
    							  stroke = I(cell_stroke),
    							  data = ya_sub[order(ya_sub$value),]) +
    				viridis::scale_color_viridis(option = "viridis",
    											 name = expression_legend_label,
    											 na.value = NA, end = 0.8,
    											 alpha = alpha) +
    				guides(alpha = 'none') + facet_wrap(~feature_label)
        } else {
          g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2)) +
            plotting_func(aes(data_dim_1, data_dim_2), size=I(cell_size),
                          stroke = I(cell_stroke), color = "grey80",
                          data = na_sub, alpha = alpha)
          if (scale_to_range){
            g <- g + plotting_func(aes(color=value),
                                   size=I(cell_size), stroke = I(cell_stroke),
                                   data = ya_sub[order(ya_sub$value),],
                                   alpha = alpha)
          }else{
            g <- g + plotting_func(aes(color=log10(value+min_expr)),
                                   size=I(cell_size), stroke = I(cell_stroke),
                                   data = ya_sub[order(ya_sub$value),],
                                   alpha = alpha)
          }
    
          g <- g +
            viridis::scale_color_viridis(option = "viridis",
                                         name = expression_legend_label,
                                         na.value = NA, end = 0.8,
                                         alpha = alpha) +
            guides(alpha = 'none') + facet_wrap(~feature_label)
        }
      } else {
        g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2))
    
        g <- g + geom_point(aes(color = cell_color), size=1.5*cell_size,    #  修改的地方，边框颜色和填充色一致----------------------------------------
                            stroke = I(cell_stroke), na.rm = TRUE,
                            alpha = I(alpha))
    
        # We don't want to force users to call order_cells before even being able
        # to look at the trajectory, so check whether it's null and if so, just
        # don't color the cells
        if(color_cells_by %in% c("cluster", "partition")){
          if (is.null(data_df$cell_color)){
            g <- g + geom_point(aes(color = cell_color), size=I(cell_size),  # 修改的地方，边框颜色和填充色一致----------------------------------------
                                stroke = I(cell_stroke), na.rm = TRUE,
                                alpha = I(alpha))
            message("cluster_cells() has not been called yet, can't ",
                          "color cells by cluster")
          } else{
            g <- g + geom_point(aes(color = cell_color), size=I(cell_size),
                                stroke = I(cell_stroke), na.rm = TRUE,
                                alpha = alpha)
          }
          g <- g + guides(color = guide_legend(title = color_cells_by,
                                               override.aes = list(size = 4)))
        } else if (methods::is(data_df$cell_color, "numeric")) {
          g <- g + geom_point(aes(color = cell_color), size=I(cell_size),
                              stroke = I(cell_stroke), na.rm = TRUE, alpha = alpha)
          g <- g + viridis::scale_color_viridis(name = color_cells_by, option="C")
        } else {
          g <- g + geom_point(aes(color = cell_color), size=I(cell_size),
                              stroke = I(cell_stroke), na.rm = TRUE, alpha = alpha)
          g <- g + guides(color = guide_legend(title = color_cells_by,
                                               override.aes = list(size = 4)))
        }
    
      }
      if (show_trajectory_graph){
        g <- g + geom_segment(aes_string(x="source_prin_graph_dim_1",
                                         y="source_prin_graph_dim_2",
                                         xend="target_prin_graph_dim_1",
                                         yend="target_prin_graph_dim_2"),
                              size=trajectory_graph_segment_size,
                              color=I(trajectory_graph_color),
                              linetype="solid",
                              na.rm=TRUE,
                              data=edge_df)
    
    
        if (label_principal_points) {
          mst_branch_nodes <- branch_nodes(cds, reduction_method)
          mst_leaf_nodes <- leaf_nodes(cds, reduction_method)
          mst_root_nodes <- root_nodes(cds, reduction_method)
          pps <- c(mst_branch_nodes, mst_leaf_nodes, mst_root_nodes)
          princ_point_df <- ica_space_df %>%
            dplyr::slice(match(names(pps), sample_name))
    
          g <- g +
            geom_point(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2"),
                       shape = 21, stroke=I(trajectory_graph_segment_size),
                       color="white",
                       fill="black",
                       size=I(graph_label_size * 1.5),
                       na.rm=TRUE, princ_point_df) +
            ggrepel::geom_text_repel(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2",
                                 label="sample_name"),
                      size=I(graph_label_size * 1.5), color="Black", na.rm=TRUE,
                      princ_point_df)
        }
        if (label_branch_points){
          mst_branch_nodes <- branch_nodes(cds, reduction_method)
          branch_point_df <- ica_space_df %>%
            dplyr::slice(match(names(mst_branch_nodes), sample_name)) %>%
            dplyr::mutate(branch_point_idx = seq_len(dplyr::n()))
    
          g <- g +
            geom_point(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2"),
                       shape = 21, stroke=I(trajectory_graph_segment_size),
                       color="white",
                       fill="black",
                       size=I(graph_label_size * 1.5),
                       na.rm=TRUE, branch_point_df) +
            geom_text(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2",
                                 label="branch_point_idx"),
                      size=I(graph_label_size), color="white", na.rm=TRUE,
                      branch_point_df)
        }
    
        if (label_leaves){
          mst_leaf_nodes <- leaf_nodes(cds, reduction_method)
          leaf_df <- ica_space_df %>%
            dplyr::slice(match(names(mst_leaf_nodes), sample_name)) %>%
            dplyr::mutate(leaf_idx = seq_len(dplyr::n()))
    
          g <- g +
            geom_point(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2"),
                       shape = 21, stroke=I(trajectory_graph_segment_size),
                       color="black",
                       fill="lightgray",
                       size=I(graph_label_size * 1.5),
                       na.rm=TRUE,
                       leaf_df) +
            geom_text(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2",
                                 label="leaf_idx"),
                      size=I(graph_label_size), color="black", na.rm=TRUE, leaf_df)
        }
    
        if (label_roots){
          mst_root_nodes <- root_nodes(cds, reduction_method)
          root_df <- ica_space_df %>%
            dplyr::slice(match(names(mst_root_nodes), sample_name)) %>%
            dplyr::mutate(root_idx = seq_len(dplyr::n()))
    
          g <- g +
            geom_point(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2"),
                       shape = 21, stroke=I(trajectory_graph_segment_size),
                       color="black",
                       fill="white",
                       size=I(graph_label_size * 1.5),
                       na.rm=TRUE,
                       root_df) +
            geom_text(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2",
                                 label="root_idx"),
                      size=I(graph_label_size), color="black", na.rm=TRUE, root_df)
        }
      }
    
      if(label_cell_groups) {
        g <- g + ggrepel::geom_text_repel(data = text_df,
                                          mapping = aes_string(x = "text_x",
                                                               y = "text_y",
                                                               label = "label"),
                                          size=I(group_label_size))
        # If we're coloring by gene expression, don't hide the legend
        if (is.null(markers_exprs))
          g <- g + theme(legend.position="none")
      }
    
      g <- g +
        #scale_color_brewer(palette="Set1") +
        monocle_theme_opts() +
        xlab(paste(reduction_method, x)) +
        ylab(paste(reduction_method, y)) +
        #guides(color = guide_legend(label.position = "top")) +
        theme(legend.key = element_blank()) +
        theme(panel.background = element_rect(fill='white'))
      g
    }
    
    #--root_nodes------------------------------------------------------------------------------------
    root_nodes <- function(cds, reduction_method="UMAP"){
      g = principal_graph(cds)[[reduction_method]]
      root_pr_nodes <- which(names(igraph::V(g)) %in%
                        cds@principal_graph_aux[[reduction_method]]$root_pr_nodes)
      names(root_pr_nodes) <-
        cds@principal_graph_aux[[reduction_method]]$root_pr_nodes
      return(root_pr_nodes)
    }
    
    #--monocle_theme_opts----------------------------------------------------------------------------
    monocle_theme_opts <- function(){
      theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
        theme(panel.border = element_blank()) +
        theme(axis.line.x = element_line(linewidth =0.25, color="black")) +
        theme(axis.line.y = element_line(linewidth=0.25, color="black")) +
        theme(panel.grid.minor.x = element_blank(),
              panel.grid.minor.y = element_blank()) +
        theme(panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_blank()) +
        theme(panel.background = element_rect(fill='white')) +
        theme(legend.key=element_blank())
    }
    
    #--branch_nodes---------------------------------------------------------------------------------
    branch_nodes <- function(cds,reduction_method="UMAP"){
      g = principal_graph(cds)[[reduction_method]]
      branch_points <- which(igraph::degree(g) > 2)
      branch_points = branch_points[branch_points %in% root_nodes(cds, reduction_method) == FALSE]
      return(branch_points)
    }
    
    #--leaf_nodes-----------------------------------------------------------------------------------
    leaf_nodes <- function(cds,reduction_method="UMAP"){
      g = principal_graph(cds)[[reduction_method]]
      leaves <- which(igraph::degree(g) == 1)
      leaves = leaves[leaves %in% root_nodes(cds, reduction_method) == FALSE]
      return(leaves)
    }

    #########################################################################################
    #########################################   Main    #####################################
    #########################################################################################
    
    if(!file.exists(output_rds)){
        # 无运行好的cds
        # ==================== 1. 构建 CDS ====================
        counts <- SeuratObject::LayerData(seurat,assay = "RNA",layer = "counts")
        cell_metadata <- seurat@meta.data
        gene_metadata <- data.frame(gene_short_name = rownames(counts),row.names = rownames(counts))
    
        # 确保细胞顺序一致
        stopifnot(all(colnames(counts) %in% rownames(cell_metadata)))
        cell_metadata <- cell_metadata[colnames(counts), , drop = FALSE]
        cds <- new_cell_data_set(expression_data = counts,
                                 cell_metadata = cell_metadata,
                                 gene_metadata = gene_metadata
                                )
        
        # ==================== 2. 预处理 ====================
        cds <- preprocess_cds(cds,num_dim = 30,norm_method = "log")
    
        # ==================== 3. 导入 Seurat UMAP ====================
        umap <- Embeddings(seurat,reduction = reduction)
        stopifnot(all(colnames(cds) %in% rownames(umap)))
        umap <- umap[colnames(cds), 1:2, drop = FALSE]
        colnames(umap) <- c("UMAP_1", "UMAP_2")
        SingleCellExperiment::reducedDim(cds, "UMAP") <- umap
    
        # ==================== 4. 聚类与 Partition ====================
        cds <- cluster_cells(cds,reduction_method = "UMAP",resolution = 0.1)
        
        # ==================== 5. 轨迹推断 ====================
        cds <- learn_graph(cds,
                           use_partition = TRUE,
                           close_loop = FALSE,
                           learn_graph_control = list(minimal_branch_len = 30,prune_graph = TRUE)
                          )
    
        #############################  出节点图  #######################################
        g <- principal_graph(cds)[["UMAP"]]
        
        # 获取 principal graph 的所有叶节点
        leaf_node_ids <- names(which(igraph::degree(g) == 1))
    
        # 获取节点信息
        closest_vertex <- principal_graph_aux(cds)[["UMAP"]]$pr_graph_cell_proj_closest_vertex
        closest_vertex <- as.matrix(closest_vertex)
        
        df_root <- data.frame(cell = rownames(closest_vertex),
                              node = paste0("Y_", closest_vertex[, 1]),
                              Celltype = colData(cds)[rownames(closest_vertex), celltype_col]
                             )
        df_leaf <- df_root %>% filter(node %in% leaf_node_ids)
        leaf_stat <- df_leaf %>%
            dplyr::count(node, Celltype) %>% dplyr::group_by(node) %>%
            dplyr::mutate( total = sum(n), prop = n / total) %>%
            dplyr::arrange(node, desc(prop))
    
        if(is.null(celltype_start)){
            cat("\n",yellow,"请输入起始细胞类型: ",reset,"\n")
            flush.console()
            celltype_start <- readline(prompt = "start celltype: ")
        }
        
        candidate_root <- leaf_stat %>% filter(Celltype == celltype_start)
        node_coords <- t(principal_graph_aux(cds)[["UMAP"]]$dp_mst)
        
        node_coords <- as.data.frame(node_coords)
        colnames(node_coords) <- c("UMAP_1", "UMAP_2")
        node_coords$node <- rownames(node_coords)
        
        candidate_coords <- node_coords %>% filter(node %in% candidate_root$node)
        
        p <- plot_cells(cds,
                        color_cells_by = celltype_col,
                        label_cell_groups = FALSE,
                        label_leaves = FALSE,
                        label_branch_points = FALSE,
                        show_trajectory_graph = TRUE
                    ) + 
            geom_point( data = candidate_coords,aes(x = UMAP_1, y = UMAP_2),inherit.aes = FALSE,size = 4,shape = 21,fill = "white") +
            geom_text(data = candidate_coords,
                      aes(x = UMAP_1,y = UMAP_2,label = node),
                        inherit.aes = FALSE,
                        size = 4,
                        nudge_y = 0.15
                    )
        basicR::save_figure(obj = p,
                        filename =  file.path(output_dir, paste0("02_Monocle3_RootNodes")),
                        width = 7,
                        height = 6,
                        res = 600,
                        formats = formats
                       )
    
        if(point_start_use){
            cat("\n",yellow,"请输入起始细胞节点: ",reset,"\n")
            cat("\n",yellow,"Eg: Y_794,Y_795",reset,"\n")
            flush.console()
            point_start_character <- readline(prompt = "start point: ")
            point_start = point_start_character %>% strsplit(",") %>% unlist()
        }
    
        # 设置起点拟时序
        if(point_start_use){
            # 设置拟时序起点
            cds <- order_cells(cds,root_pr_nodes = point_start)
            
        }else{
            # 提取 Macro_CIITA 细胞
            root_cells <- colnames(cds)[colData(cds)$Celltype == celltype_start]
            # 设置拟时序起点
            cds <- order_cells(cds,root_cells = root_cells)  
            
        }
    
        saveRDS(cds,file = output_rds)
        
    }else{
        # 已有运行好的cds
        cds <- readRDS(output_rds)
        
    }
        
    ####################################### Pseudotime DEGs ##########################################
    cat("\n",yellow,"Pseudotime DEGs Start ...",reset,"\n")
    # 拟时序差异基因，每批 1000 个基因运行,控制内存占用
    genes <- rownames(cds)
    gene_batches <- split(genes,ceiling(seq_along(genes) / 5000))
    length(gene_batches)
    
    deg_list <- vector("list", length(gene_batches))
    
    for (i in seq_along(gene_batches)) {
        cat("Batch:", i, "/", length(gene_batches)," genes:", length(gene_batches[[i]]), "\n")
        deg_list[[i]] <- graph_test(cds[gene_batches[[i]], ],neighbor_graph = "principal_graph",cores = 1)
        gc()
    }
    deg_pt <- do.call(rbind,deg_list)
    
    deg_pt_sig <- deg_pt %>%
        tibble::rownames_to_column("gene") %>%
        dplyr::filter(status == "OK",q_value < 0.05,morans_I > 0.1) %>%
        dplyr::arrange(desc(morans_I))
    
    write.csv(deg_pt_sig,file = file.path(output_dir,"02_Monocle3_Deg_PseudoTime_sig.csv"))

    df_pt <- pseudotime(cds) %>% as.data.frame() %>% setNames("monocle3_pseudotime")
    write.csv(df_pt,file = file.path(output_dir,"02_Monocle3_pseudotime.csv") )

    ######################################## Plot #################################################
    # 无节点无轨迹拟时序
    plot_cells(cds, 
               show_trajectory_graph  = FALSE,
               reduction_method = "UMAP",
               cell_size = 1,
               color_cells_by = "pseudotime", 
               label_cell_groups = F, 
               label_leaves = F,  
               label_branch_points = F,
               label_roots = T,
               labels_per_group = F) + 
        theme_minimal(base_size = 16)+
        labs(x= "",y = "")+
        theme(panel.grid = element_blank(),
              legend.title = element_text(size = 10, face = "bold", margin = margin(b = 20),hjust = 0.5),
              axis.text = element_blank()
             ) -> p
    basicR::save_figure(obj = p,
                        filename =  file.path(output_dir, paste0("03_Monocle3_pseudotime_noNode_noTrajectory")),
                        width = 7,
                        height = 6,
                        res = 600,
                        formats = formats
                       )
    
    # 无节点有轨迹拟时序
    plot_cells(cds, 
               reduction_method = "UMAP",
               cell_size = 1,
               color_cells_by = "pseudotime", 
               label_cell_groups = F, 
               label_leaves = F,  
               label_branch_points = F,
               label_roots = F,
               labels_per_group = F
              ) + 
        theme_minimal(base_size = 16)+
        labs(x= "",y = "")+
        theme(panel.grid = element_blank(),
              legend.title = element_text(size = 18, face = "bold", margin = margin(b = 20),hjust = 0.5),
              axis.text = element_blank()
             ) -> p
    basicR::save_figure(obj = p,
                        filename =  file.path(output_dir, paste0("03_Monocle3_pseudotime_noNode_Trajectory")),
                        width = 7,
                        height = 6,
                        res = 600,
                        formats = formats
                       )
    
    # 有节点有轨迹拟时序
    plot_cells(cds, 
               reduction_method = "UMAP",
               cell_size = 1,
               color_cells_by = "pseudotime", 
               label_cell_groups = F, 
               label_leaves = F,  
               label_branch_points = F,
               label_roots = T,
               labels_per_group = F
              ) + 
        theme_minimal(base_size = 16)+
        labs(x= "",y = "")+
        theme(panel.grid = element_blank(),
              legend.title = element_text(size = 18, face = "bold", margin = margin(b = 20),hjust = 0.5),
              axis.text = element_blank()
             ) -> p
    basicR::save_figure(obj = p,
                        filename =  file.path(output_dir, paste0("03_Monocle3_pseudotime_Node_Trajectory")),
                        width = 7,
                        height = 6,
                        res = 600,
                        formats = formats
                       )
    
    # 有标签有轨迹拟时序
    plot_cells(cds,
               reduction_method="UMAP",
               color_cells_by= celltype_col,
               label_principal_points = T,
               label_branch_points = T,
               label_roots =T,
               label_cell_groups = F,
               cell_size = 1
              ) + 
    scale_color_manual(values = celltype_colors) + 
    theme_minimal(base_size = 16)+
    labs(x= "",y = "")+
    theme(panel.grid = element_blank(),
          legend.title = element_text(size = 18, face = "bold", margin = margin(b = 20),hjust = 0.5),
          axis.text = element_blank()
         ) -> p
    basicR::save_figure(obj = p,
                        filename =  file.path(output_dir, paste0("03_Monocle3_pseudotime_Label_Trajectory")),
                        width = 7,
                        height = 6,
                        res = 600,
                        formats = formats
                       )
    # 删除节点辅助图
    RootNodes_files <- list.files(output_dir,pattern = "02_Monocle3_RootNodes",full.names = T)
    if(length(RootNodes_files) > 0){
        invisible(file.remove(RootNodes_files))
    }
    
}