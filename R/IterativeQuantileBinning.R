### iterative quantile binning functions
# first source in all helper functions
# source("IterativeQuantileBinningSupportFunctions.R")
# Note: a #!# tag will be added on items that need improved efficiency/clarity


# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#--------------------------------------
#' Iterative Quantile Binning
#'
#' @description used for binning the numeric values by empirical quantile
#'
#' @param data data frame to be binned (will coerce matrix or tibble to simple data frame)
#' @param bin_cols vector of column names of variables to iteratively bin, ordered first to last
#' @param nbins vector of number of bins per step of iterative binning, ordered first to last
#' @param jit vector of margins for uniform jitter to each dimension to create seperability of tied obs due to finite precision
#' @param output Output Structure: "data" for just the binned data,"definition" for a list of bin centers and boundaries, or "both" for list containing both data and definition
#'
#' @return output list containing values specified by user
#' @family iterative quantile binning functions
#' @export
#' @examples
#' iqbin(data=iris, bin_cols=c("Sepal.Length","Sepal.Width","Petal.Width"),
#'                     nbins=c(3,5,2), output="both",jit=rep(0.001,3))
#' iqbin(data=iris, bin_cols=c("Sepal.Length","Sepal.Width","Petal.Width"),
#'                     nbins=c(3,5,2), output="data")
iqbin <- function(data, bin_cols, nbins, jit = rep(0,length(bin_cols)), output="data"){
  data <- as.data.frame(data) 
  bin_dim <- length(bin_cols)
  # Initialize with first binning step
  bin_jit <- NULL
  step_bin_info <- quant_bin_1d(data[,bin_cols[1]], nbins[1],output="both",jit[1])
  # if(sum(jit>0)>0) bin_jit <- data.frame("V1"=step_bin_info$jit_values) #!# drop for timing study
  bound_list <- list(NULL) 
  bound_list[[nbins[1]+1]] <- step_bin_info$bin_bounds
  bin_bounds <- matrix(c(step_bin_info$bin_bounds[1:nbins[1]],
                         step_bin_info$bin_bounds[2:(nbins[1]+1)]),
                       nrow=nbins[1],byrow=FALSE )
  bin_indeces <- matrix(step_bin_info$bin_number, ncol=1)
  # Loop over remaining variables to use quantile binning WITHIN each of previous state bins
  #!# At some stage need to think about restructuring how binning definition is structured, more included with nested lists perhaps or indexed dataframe with list columns of bin attributes
  for(d in 2:bin_dim){
    stack_size <- nrow(bin_bounds)
    stack_matrix <- make_stack_matrix(stack_size,nbins[d])
    bin_bounds <- cbind(stack_matrix %*% bin_bounds,matrix(rep(NA,2*stack_size*nbins[d]),ncol=2))
    bin_indeces <- cbind(bin_indeces,NA)
    unique_bin_indeces <- unique(bin_indeces)
    # iterate through unique bins from prior step which are the {1,1+nbins[d],1+2*nbins[d],...} rows of the bin matrices
    for(b in 1:prod(nbins[1:(d-1)]) ){
      in_bin_b <- apply(bin_indeces,1,identical,y=unique_bin_indeces[b,])
      step_bin_info <- quant_bin_1d(data[in_bin_b,bin_cols[d]], nbins[d],output="both",jit[d])
      bin_bounds[(b-1)*nbins[d]+1:nbins[d],c(2*d-1,2*d)] <- matrix(c(step_bin_info$bin_bounds[1:nbins[d]],
                                                                     step_bin_info$bin_bounds[2:(nbins[d]+1)]),
                                                                   nrow=nbins[d],byrow=FALSE)
      bin_indeces[in_bin_b,d] <- step_bin_info$bin_number
      # if(sum(jit>0)>0) bin_jit[in_bin_b,d] <- step_bin_info$jit_values #!# drop for timing study
    }
  }
  # add bin index column to data with collapse indeces in for each bin into single 
  data$bin_index <- sapply(1:nrow(data), function(x) index_collapser(nbins=nbins, indeces=bin_indeces[x,]))
  
  # Create list tree structure for fast queries
  bin_list <- make_bin_list(bin_bounds,nbins)
  if(output=="data") iqbin_obj <- list(data=data,bin_data=bin_data,bin_jit=bin_jit)
  if(output=="definition") {
    iqbin_obj <- list(bin_bounds=bin_bounds, bin_cols=bin_cols, nbins=nbins, jit=jit, bin_list=bin_list)
  }
  if(output=="both"){
    iqbin_obj <- list(bin_data=list(data=data,bin_jit=bin_jit),
                      bin_def=list(bin_bounds=bin_bounds, bin_cols=bin_cols, nbins=nbins, jit=jit, bin_list=bin_list))
    attributes(iqbin_obj)$iq_obj_type <- "iqbin"
  }
  return(iqbin_obj)
}
  
  

#--------------------------------------
#' Stretching to Add Tolerance Buffer to outermost Iterative Quantile Bins
#'
#' @description We may wish to include values just outside the constructed bins in future applications. \code{iqbin_stretch} redefines the outermost bin in each dimension from a definition created by \code{\link{iqbin}}
#'
#' @param iq_def  iterative quantile binning definition list from \code{\link{iqbin}} function
#' @param tol vector of tolerance values to stretch each dimension for future binning
#'
#' @return updated binning definition list with outermost bin boundaries extended by tolerance values
#' @family iterative quantile binning functions
#' @export
#' @examples
#' iq_def <- iqbin(data=iris, bin_cols=c("Sepal.Length","Sepal.Width","Petal.Width"),
#'                              nbins=c(3,2,2), output="both")
#' stretch_iq_def <- iqbin_stretch(iq_def, tol = c(1,1,1))
#' iq_def$bin_def$bin_bounds
#' stretch_iq_def$bin_def$bin_bounds

iqbin_stretch <- function(iq_def, tol){
  b = nrow(iq_def$bin_def$bin_bounds)
  p = length(iq_def$bin_def$nbins)
  for (d in 1:p){
    blocks <- prod(iq_def$bin_def$nbins[1:d-1])
    blocks_n <- b/blocks
    subblocks <- prod(iq_def$bin_def$nbins[1:d])
    subblocks_n <- b/subblocks
    # stretch the bins
    for(block in 1:blocks){
      lowest_in_block <- seq(1+((block-1)*blocks_n),subblocks_n+((block-1)*blocks_n),by=1)
      highest_in_block <- seq(blocks_n-subblocks_n+1+((block-1)*blocks_n), blocks_n+((block-1)*blocks_n),by=1)
      iq_def$bin_def$bin_bounds[lowest_in_block,2*d-1] <- iq_def$bin_def$bin_bounds[lowest_in_block,2*d-1] - tol[d]
      iq_def$bin_def$bin_bounds[highest_in_block,2*d] <- iq_def$bin_def$bin_bounds[highest_in_block,2*d] + tol[d]
    }
  }
  iq_def$bin_def$bin_list <- make_bin_list(iq_def$bin_def$bin_bounds,iq_def$bin_def$nbins)
  return(iq_def)
}

#--------------------------------------
#' Assigning New observations to Existing Iterative Quantile Bins
#'
#' @description Assigning New observations to Existing Iterative Quantile Bins by using binning definition created by \code{\link{iqbin}} function 
#'
#' @param iq_def Iterative quantile binning definition list as output from \code{\link{iqbin}} function
#' @param new_data Data frame with column names matching the binned columns from bin-training data
#' @param output Matches format of \code{\link{iqbin}} and inherets properties from \code{\link{iqnn}} if applicable {"data","both"}
#' @param strict TRUE/FALSE: If TRUE Observations must fall within existing bins to be assigned; if FALSE the outer bins in each dimension are unbounded to allow outlying values to be assigned.
#'
#' @return list of new_data, binned data and indecies. Optionally the bin definition may also be included
#' @family iterative quantile binning functions
#' @export
#' @examples
#' withhold_index <- c(1,2,51,52,101,102)
#' iq_def <- iqbin(data=iris[-withhold_index,], bin_cols=c("Sepal.Length","Sepal.Width","Petal.Width"),
#'                               nbins=c(3,2,2), output="definition")
#' iqbin_assign(bin_def=iq_def, new_data=iris[withhold_index,], output="data")

iqbin_assign <- function(bin_def, new_data, output="data", strict=FALSE){
  # new_data <- as.data.frame(new_data)
  #!# need to introduce similar jitter to new data as in definition so "boundary" points allocated randomly
  #
  # loop over each obs in new data, identify the bin indeces then return bin centers for associated bins
  new_data$bin_index <- sapply(1:nrow(new_data), function(i){
    val <- bin_index_finder_nest(new_data[i,bin_def$bin_cols],bin_def, strict=strict)
    if(is.null(val)) val = NA
    return(val)
  })
  
  if(output=="bin_index") return(new_data$bin_index)
  if(output=="data") return(new_data)
  if(output=="both"){
    return(list(bin_data=list(data=new_data,bin_data=bin_def$bin_centers[bin_indeces,],bin_indeces=bin_indeces),
                bin_def=iq_def))
  }
}



