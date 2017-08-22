### iterative quantile bin plotting functions

#--------------------------------------
#' Plotting 2-dimensional iterative quantile bins
#'
#' @description used for binning the numeric values by empirical quantile
#'
#' @param iq_def iterative quantile definition: either from iq 
#' @param nbin An integer defining the number of bins to partion the xs into
#' @param output Output Structure: "data" for just the binned data,"definition" for a list of bin centers and boundaries, or "both" for list containing both data and definition
#' @param jit non-negative value to specify a random uniform jitter to the observed values prior to partitioning by quantiles.
#'
#' @return output as specified
#' @examples

# iqbin_plot_2d <- function(iq_def, ){
#   
# }