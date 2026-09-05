% Documentation for R/calculate_cell_point_size.R; keep in sync with its opening help comments.
\name{calculate_cell_point_size}
\alias{calculate_cell_point_size}
\title{
Calculate plotting point size from the number of cells
}
\description{
Calculates an empirical point size for single-cell scatter plots, decreasing the size as the number of cells increases.
}
\usage{
calculate_cell_point_size(cell_count)
}
\arguments{
\item{cell_count}{A positive numeric cell count or a numeric vector of cell counts.}
}
\details{
The size is calculated as exp(3.84342) * cell_count^(-0.57370), then clipped to the interval [0.01, 8]. The calculation is vectorised. Inputs are not validated; missing values propagate to the result.
}
\value{
A numeric vector of point sizes with the same length as cell_count.
}
\examples{
calculate_cell_point_size(10000)
calculate_cell_point_size(c(100, 1000, 10000))
}
