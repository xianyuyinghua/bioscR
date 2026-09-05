#' Calculate plotting point size from the number of cells
#'
#' Calculates an empirical point size for single-cell scatter plots, decreasing the size as the number of cells increases.
#'
#' @param cell_count A positive numeric cell count or a numeric vector of cell
#'   counts.
#'
#' @details
#' The size is calculated as exp(3.84342) * cell_count^(-0.57370), then clipped
#' to the interval [0.01, 8]. The calculation is vectorised. Inputs are not
#' validated; missing values propagate to the result.
#'
#' @return A numeric vector of point sizes with the same length as cell_count.
#'
#' @examples
#' calculate_cell_point_size(10000)
#' calculate_cell_point_size(c(100, 1000, 10000))
calculate_cell_point_size <- function(cell_count) {
    # 单细胞定义计算点大小的函数
    a <- exp(3.84342)  # 截距的指数
    b <- -0.57370      # 回归系数
    point_size <- a * cell_count^b
    point_size <- pmin(pmax(point_size, 0.01), 8)
    return(point_size)  
}