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
#' @return output as specified
#' @examples
#' iqbin(data=iris, bin_cols=c("Sepal.Length","Sepal.Width","Petal.Width"),
#'                     nbins=c(3,5,2), output="both",jit=rep(0.001,3))
#' iqbin(data=iris, bin_cols=c("Sepal.Length","Sepal.Width","Petal.Width"),
#'                     nbins=c(3,5,2), output="data")
iqbin <- function(data, bin_cols, nbins, jit = rep(0,length(bin_cols)), output="data"){
  data <- as.data.frame(data)
  # row.names(data) <- 1:nrow(data)
  bin_dim <- length(bin_cols)
  bin_data <- matrix(NA,nrow=nrow(data),ncol=bin_dim, dimnames=list(row.names(data),paste(bin_cols,"binned",sep="_")))
  # Initialize with first binning step
  bin_jit <- NULL
  step_bin_info <- quant_bin_1d(data[,bin_cols[1]], nbins[1],output="both",jit[1])
  if(sum(jit>0)>0) bin_jit <- data.frame("V1"=step_bin_info$jit_values)
  bin_bounds <- matrix(c(step_bin_info$bin_bounds[1:nbins[1]],
                         step_bin_info$bin_bounds[2:(nbins[1]+1)]),
                       nrow=nbins[1],byrow=FALSE )
  bin_centers <- matrix(step_bin_info$bin_centers, nrow=nbins[1])
  bin_data[,1] <- step_bin_info$bin_data
  # Loop over remaining variables to use quantile binning WITHIN each of previous state bins
  #!# At some stage need to think about restructuring how binning definition is structured, more included with nested lists perhaps or indexed dataframe with list columns of bin attributes
  for(d in 2:bin_dim){
    stack_size <- nrow(bin_centers)
    stack_matrix <- make_stack_matrix(stack_size,nbins[d])
    bin_centers <- cbind(stack_matrix %*% bin_centers,matrix(rep(NA,stack_size*nbins[d]),ncol=1))
    bin_bounds <- cbind(stack_matrix %*% bin_bounds,matrix(rep(NA,2*stack_size*nbins[d]),ncol=2))
    # iterate through unique bins from prior step which are the {1,1+nbins[d],1+2*nbins[d],...} rows of the bin matrices
    for(b in seq(1,1+(stack_size-1)*nbins[d],by=nbins[d]) ){
      in_bin_b <- apply(matrix(bin_data[,1:(d-1)],ncol=(d-1)),1,identical,y=bin_centers[b,-d])
      step_bin_info <- quant_bin_1d(data[in_bin_b,bin_cols[d]], nbins[d],output="both",jit[d])
      bin_bounds[b:(b+nbins[d]-1),c(2*d-1,2*d)] <- matrix(c(step_bin_info$bin_bounds[1:nbins[d]],
                                                            step_bin_info$bin_bounds[2:(nbins[d]+1)]),
                                                          nrow=nbins[d],byrow=FALSE)
      bin_centers[b:(b+nbins[d]-1),d] <- matrix(step_bin_info$bin_centers, nrow=nbins[d])
      bin_data[in_bin_b,d] <- step_bin_info$bin_data
      if(sum(jit>0)>0) bin_jit[,d] <- step_bin_info$jit_values
    }
  }
  # add bin index column to bin_data
  bin_centers_idx <- as.data.frame(bin_centers)
  bin_data_idx <- as.data.frame(bin_data)
  names(bin_centers_idx) <- colnames(bin_data)
  bin_centers_idx$bin_index <- 1:nrow(bin_centers_idx)
  bin_data_idx$order <- 1:nrow(bin_data_idx)
  bin_data <- merge(round_df(bin_data_idx,10),round_df(bin_centers_idx,10), all.x=TRUE)
  bin_data <- bin_data[order(bin_data$order),]
  row.names(bin_data) <- 1:nrow(bin_data)
  #
  bin_list <- make_bin_list(bin_bounds,nbins)
  if(output=="data") iqbin_obj <- list(data=data,bin_data=bin_data,bin_jit=bin_jit)
  if(output=="definition") iqbin_obj <- list(bin_centers=bin_centers, bin_bounds=bin_bounds, bin_cols=bin_cols, nbins=nbins, jit=jit, bin_list=bin_list)
  if(output=="both"){
    iqbin_obj <- list(bin_data=list(data=data,bin_data=bin_data,bin_jit=bin_jit),
                bin_def=list(bin_centers=bin_centers, bin_bounds=bin_bounds, bin_cols=bin_cols, nbins=nbins, jit=jit, bin_list=bin_list))
  }
  class(iqbin_obj) <- "iqbin"
  return(iqbin_obj)
}

#--------------------------------------
#' Adding Tolerance Buffer to outermost bins from Iterative Quantile Binning
#'
#' @description New observations selected from the same population as the data used to build bin definitions may fall just outside the bins. If we wish to include nearby values we can either allow outer bins to be extended (this function) or to leave the outer bins unbounded.
#'
#' @param iq_def  iterative quantile binning definition list
#' @param tol vector of tolerance values to stretch each dimension for future binning
#'
#' @return updated binning definition with bins extended by tolerance values
#' @examples
# iq_def <- iqbin(data=iris, bin_cols=c("Sepal.Length","Sepal.Width","Petal.Width"),
#                               nbins=c(3,2,2), output="both")
# stretch_iq_def <- iqbin_stretch(iq_def, tol = c(1,1,1))
# iq_def$bin_def$bin_bounds
# stretch_iq_def$bin_def$bin_bounds

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
#' Iterative Quantile Binning New Data from defined bins
#'
#' @description New observations selected from the same population as the data used to build bin definitions may fall just outside the bins. If we wish to include nearby values we can either allow outer bins to be extended (this function) or to leave the outer bins unbounded.
#'
#' @param iq_def Iterative quantile binning definition list
#' @param new_data Data frame with column names matching the binned columns from bin-training data
#' @param output Matches format of iqbin and inherets properties from iqnn if applicable {"data","both"}
#' @param strict TRUE/FALSE: If TRUE Observations must fall within existing bins to be assigned; if FALSE the outer bins in each dimension are unbounded to allow outlying values to be assigned.
#'
#' @return updated binning definition with bins extended by tolerance values
#' @examples
#' withhold_index <- c(1,2,51,52,101,102)
#' iq_def <- iqbin(data=iris[-withhold_index,], bin_cols=c("Sepal.Length","Sepal.Width","Petal.Width"),
#'                               nbins=c(3,2,2), output="definition")
#' iqbin_assign(bin_def=iq_def, new_data=iris[withhold_index,], output="data")

iqbin_assign <- function(bin_def, new_data, output="data", strict=FALSE){
  new_data <- as.data.frame(new_data)
  #!# need to introduce similar jitter to new data as in definition so "boundary" points allocated randomly
  #
  # loop over each obs in new data, identify the bin indeces then return bin centers for associated bins
  bin_indeces <- sapply(1:nrow(new_data), function(i){
    # bin_index_finder(new_data[i,iq_def$bin_cols],iq_def$bin_bounds, iq_def$nbins, strict=strict)
    bin_index_finder_nest(new_data[i,bin_def$bin_cols],bin_def, strict=strict)
  })

  if(output=="data") return(list(data=new_data,bin_data=bin_def$bin_centers[bin_indeces,],bin_indeces=bin_indeces))
  if(output=="both"){
    return(list(bin_data=list(data=new_data,bin_data=bin_def$bin_centers[bin_indeces,],bin_indeces=bin_indeces),
                bin_def=iq_def))
  }
}



