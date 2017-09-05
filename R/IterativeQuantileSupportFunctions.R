### Helper functions for Iterative Quantile Binning Functions

#--------------------------------------
#' One-Dimensional Empirical Quantile Binning
#'
#' @description used for binning the numeric values by empirical quantile
#'
#' @param xs Numeric vector of values to be binned
#' @param nbin An integer defining the number of bins to partion the xs into
#' @param output Output Structure: "data" for just the binned data,"definition" for a list of bin centers and boundaries, or "both" for list containing both data and definition
#' @param jit non-negative value to specify a random uniform jitter to the observed values prior to partitioning by quantiles.
#'
#' @return output as specified
#' @export
#' @examples
#' quant_bin_1d(ggplot2::diamonds$price,4,output="data")
#' quant_bin_1d(ggplot2::diamonds$price,4,output="definition")
#' quant_bin_1d(runif(100,0,10),nbin=4,output="both")
#' quant_bin_1d(runif(100,0,10),nbin=4,output="both", jit=.1)

quant_bin_1d <- function(xs, nbin, output="data",jit=0){
  jit_values <- NULL
  if(jit > 0) {
    jit_values <- runif(length(xs),-jit,jit)
    xs <- xs + jit_values 
  }
  quants <- quantile(xs, seq(0, 1, by=1/(2*nbin)))
  bin_centers <- quants[seq(2,length(quants)-1, by=2)]
  bin_bounds <- quants[seq(1,length(quants)+1, by=2)]
  if(jit > 0) bin_bounds[c(1,length(bin_bounds))] <- bin_bounds[c(1,length(bin_bounds))]+c(-jit,jit)
  if(output=="definition") {
    return(list(bin_centers=bin_centers,bin_bounds=bin_bounds))
  } else{
    bin_data <- bin_centers[.bincode(xs,bin_bounds,T,T)]
    if(output=="data") return(bin_data)
    if(output=="both") return(list(bin_data=bin_data,bin_centers=bin_centers,bin_bounds=bin_bounds,jit_values=jit_values))
  }
}
# Speed test
# load("~/onePercentSample.Rdata")
# timer <- Sys.time()
# quant_bin_1d(onePercentSample$total_amount,100,output="data", jit=.00001)
# Sys.time()-timer
# Note: using .bincode() this take ~2 seconds instead of ~10 seconds with for loop overwrite from original

#--------------------------------------
#' Convert bin bounds data frame into R-tree structure using nested lists
#'
#' @description Convert bin bounds data frame into R-tree structure using nested lists
#'
#' @param bin_bounds Data frame of bin bounds; one row per bin, columns for lower and upper bound on each dimension
#' @param nbins number of bins in each dimension
#' 
#' @return R-tree nested list of bin boundaries
#' @export
#' @examples 
#' iq_def <- iqbin(data=iris, bin_cols=c("Sepal.Length","Sepal.Width","Petal.Width"),
#'                     nbins=c(3,2,2), output="definition",jit=rep(0.001,3))
#' iq_def$bin_bounds
#' make_bin_list(bin_bounds=iq_def$bin_bounds,nbins=c(3,2,2))

make_bin_list <- function(bin_bounds,nbins){
  bin_dim = length(nbins)
  ### build nested list version of bin_bounds to speed up future searching for bins
  lower_level_list <- list(NULL)
  for(i in 1:nrow(bin_bounds)){
    lower_level_list[[i]] <- i
  } 
  for(d in bin_dim:1){
    # for each dimension from second lowest to highest, group up observations from lower_level_list into items in upper_level_list 
    upper_level_list <- list(NULL)
    upper_blocks <- ifelse(d==1,1,prod(nbins[1:(d-1)]))
    lower_block_size <- nbins[d]
    upper_indeces <- prod(nbins[d:length(nbins)])
    lower_indeces <- ifelse(d==bin_dim,1,prod(nbins[(d+1):bin_dim]))
    
    for(ul in 1:upper_blocks){
      # create upper level groups 
      upper_level_list[[ul]] <- list(NULL)
      for(ll in 1:lower_block_size){
        upper_level_list[[ul]][[ll]] <- lower_level_list[[(ul-1)*lower_block_size+ll]]
      }
      # upper_level_list[[ul]][[lower_block_size+1]] <- bin_bounds[(ul-1)*upper_indeces + 1:lower_block_size*lower_indeces,(d-1)*2+1:2]
      upper_level_list[[ul]][[lower_block_size+1]] <- unique(as.vector(bin_bounds[(ul-1)*upper_indeces + 1:lower_block_size*lower_indeces,(d-1)*2+1:2]))
    }
    lower_level_list <- upper_level_list
  }
  bin_list <- lower_level_list
  return(bin_list)
}

#--------------------------------------
#' Find bin index from R-tree structure
#'
#' @description Use R-tree structure using from make_bin_list() function to find bin index for new observation
#' @param x vector of input values for each of the binned dimensions
#' @param bin_def Iterative quantile binning definition list
#' @param strict TRUE/FALSE: If TRUE Observations must fall within existing bins to be assigned; if FALSE the outer bins in each dimension are unbounded to allow outlying values to be assigned.
#' 
#' @return bin index for new observation
#' @export
#' @examples 
#' iq_def <- iqbin(data=iris[-test_index,], bin_cols=c("Sepal.Length","Sepal.Width","Petal.Width"),
#'                               nbins=c(3,2,2), output="both")
#' bin_index_finder_nest(x=c(6,3,1.5),bin_def=iq_def$bin_def, strict=TRUE)

bin_index_finder_nest <- function(x, bin_def, strict=TRUE){ 
  bin_dim = length(bin_def$nbins)
  nest_list <- bin_def$bin_list[[1]]
  x <- as.numeric(x)
    for(d in 1:bin_dim){
      nest_index <- .bincode(x[[d]], nest_list[[bin_def$nbins[d]+1]],T,T)
      if(strict == FALSE){
        if( x[[d]] < min(nest_list[[bin_def$nbins[d]+1]]) ) nest_index <- 1
        if( x[[d]] > max(nest_list[[bin_def$nbins[d]+1]]) ) nest_index <- bin_def$nbins[d]
      }
      # if(length(nest_index)==0) return(print("Observation outside of observed bins, set strict=FALSE "))
      nest_list <- nest_list[[nest_index]]
    }
    idx <- nest_list
  return(idx)
} 

#--------------------------------------
#' Stack Matrix builder
#'
#' @description function to make a J matrix for duplicating the N rows of a matrix M times each. (support funciton for iqbin function)
#' @param N number of rows
#' @param M number of duplicates
#' 
#' @return bin index for new observation
#' @export
#' @examples 
#' make_stack_matrix(3,4)

make_stack_matrix <- function(N,M){ 
  mat <- unname(model.matrix(~as.factor(rep(1:N,each=M))-(1)))
  attributes(mat)[2:3]<-NULL
  return(mat)
} 

#--------------------------------------
#' CV cohort additions
#'
#' @description Create set of indeces for tracking observations to cv folds. Creates equally balanced folds
#'
#' @param dat data frame for cross validation
#' @param cv_K number of folds needed in indexing
#' 
#' @return Vector of fold indeces
#' @export
#' @examples 
#' make_cv_cohorts(iris,cv_K=10)
make_cv_cohorts <- function(dat,cv_K){
  if(nrow(dat) %% cv_K == 0){ # if perfectly divisible
    cv_cohort <- sample(rep(1:cv_K, each=(nrow(dat)%/%cv_K)))
  } else { # if not perfectly divisible
    cv_cohort <- sample(c(rep(1:(nrow(dat) %% cv_K), each=(nrow(dat)%/%cv_K + 1)),
                              rep((nrow(dat) %% cv_K + 1):cv_K,each=(nrow(dat)%/%cv_K)) ) )
  }
  return(cv_cohort)
}

#--------------------------------------
#' Create data frame with rounded number of trailing digits
#'
#' @description Based on https://gist.github.com/avsmith/e6f4f654451da139230b to round all numeric variables
#' @param x data frame 
#' @param digits number of digits to round
#' 
#' @return data frame with rounded numeric variables
#' @export
#' @examples 
#' round_df(head(faithful),digits=1)

round_df <- function(x, digits=2) {
  numeric_columns <- sapply(x, class) == 'numeric'
  x[numeric_columns] <-  round(x[numeric_columns], digits)
  x
}

#--------------------------------------
#' Simple Majority Vote Counter
#'
#' @description Identify the maximum vote earners, then randomly pick winner if there is a tie to break
#'
#' @param votes character or factor vector
#' @export
#' @examples 
#' votes <- c("a","a","a","b","b","c")
#' majority_vote(votes)

majority_vote <- function(votes){
  top_votes <- names(which.max(table(as.character(votes)))) # collect top vote earner (ties allowed)
  return(sample(top_votes,1)) # randomly select to break any ties for best
}

#--------------------------------------
#' Function to create list of nbins vectors to put into tuning iqnn 
#'
#' @description create a list of nbins vectors, use progression that increases number of bins in each dimension while always staying balanced between dimensions. 
#'
#' @param nbin_range positive integer vector containing lower and upper bounds on number of bins in each dimension
#' @param p number of binning dimensions
#' 
#' @return list of nbins vectors
#' @export
#' @examples 
#' make_nbins_list(c(2,4),3)

make_nbins_list <- function(nbin_range, p){
  nbins_list <- list(rep(nbin_range[1],p))
  counter = 1
  for(i in 1:(nbin_range[2]-nbin_range[1])){
    for(j in 1:p){
      nbins_list[[counter+1]] <- nbins_list[[counter]]
      nbins_list[[counter+1]][j] <- nbins_list[[counter+1]][j] + 1
      counter <- counter+1
    }
  }
  return(nbins_list)
}

#--------------------------------------
### Helper function for suggesting parameters for jittering number of bins in each dimension
# based on number of ties and data resolution
#!#
# roots_for_nbins <- function(x, p, k){
#   nbin_opt <- x/k
#   rep(ceiling(nbin_opt^(1/p)),p)
# }
# # Test with goal to mimic 10-nn with p=3 dimensions
# roots_for_nbins(270,3,10)
