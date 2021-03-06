% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coloc_heatmap_fun.R
\name{get_colocalization_heatmap}
\alias{get_colocalization_heatmap}
\title{This function takes in deconvolution matrix and returns a network igraph object}
\usage{
get_colocalization_heatmap(decon_mtrx)
}
\arguments{
\item{decon_mtrx}{Object of class matrix, SPOTxCELL-TYPE, with the cell-type proportions per spot.}
}
\value{
This function returns a heatmap representing how the cell type colocalize on the tissue
}
\description{
This function takes in deconvolution matrix and returns a network igraph object
}
\examples{

}
