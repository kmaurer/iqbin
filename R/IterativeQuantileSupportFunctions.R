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
quant_bin_1d <- function(xs, nbin, output="data",jit=0){
  jit_values <- NULL
  if(jit > 0) {
    jit_values <- stats::runif(length(xs),-jit,jit)
    xs <- xs + jit_values 
  }
  quants <- stats::quantile(xs, seq(0, 1, by=1/(2*nbin)), type=7)
  bin_bounds <- quants[seq(1,length(quants)+1, by=2)]
  if(jit > 0) bin_bounds[c(1,length(bin_bounds))] <- bin_bounds[c(1,length(bin_bounds))]+c(-jit,jit)
  if(output=="definition") {
    return(list(bin_bounds=bin_bounds))
  } else{
    bin_number <- .bincode(xs,bin_bounds,T,T)
    if(output=="data") return(bin_number)
    if(output=="both") return(list(bin_number=bin_number,bin_bounds=bin_bounds,jit_values=jit_values))
  }
}

#--------------------------------------
#' Convert bin bounds data frame into R-tree structure using nested lists
#'
#' @description Convert bin bounds data frame into R-tree structure using nested lists
#'
#' @param bin_bounds Data frame of bin bounds; one row per bin, columns for lower and upper bound on each dimension
#' @param nbins number of bins in each dimension
#' 
#' @return R-tree nested list of bin boundaries
make_bin_list <- function(bin_bounds,nbins){
  bin_dim = length(nbins)
  if(!is.numeric(nbins)){
    stop("nbins needs to be an integer vector")
  } else if(bin_dim == 1){
    bin_list <- list(NULL)
    bin_list[[1]] <- list(sort(unique(as.numeric(bin_bounds))))
  } else {
    ### build nested list version of bin_bounds to speed up future searching for bins
    upper_blocks <- prod(nbins[1:(bin_dim-1)])
    lower_level_list <-lapply(1:upper_blocks, function(x){
      list(unique(as.vector(bin_bounds[(x-1)*nbins[bin_dim] + 1:nbins[bin_dim],(bin_dim-1)*2+1:2])))
    })
    
    # at each subsequent level higher in tree, bundle all lower_level_list objects associated with p-th dimensional interval l
    for(p in (bin_dim-1):1){
      upper_level_list <- NULL
      
      upper_blocks <- ifelse(p==1,1,prod(nbins[1:(p-1)]))
      lower_block_size <- nbins[p]
      upper_indeces <- prod(nbins[p:length(nbins)])
      lower_indeces <- ifelse(p==bin_dim,1,prod(nbins[(p+1):bin_dim]))
      
      for(ul in 1:upper_blocks){
        # create upper level groups 
        upper_level_list[[ul]] <- list(NULL)
        for(ll in 1:lower_block_size){
          upper_level_list[[ul]][[ll]] <- lower_level_list[[(ul-1)*lower_block_size+ll]]
        }
        # upper_level_list[[ul]][[lower_block_size+1]] <- bin_bounds[(ul-1)*upper_indeces + 1:lower_block_size*lower_indeces,(d-1)*2+1:2]
        upper_level_list[[ul]][[lower_block_size+1]] <- unique(as.vector(bin_bounds[(ul-1)*upper_indeces + 1:lower_block_size*lower_indeces,(p-1)*2+1:2]))
      }
      lower_level_list <- upper_level_list
    }
    bin_list <- lower_level_list
  }
  return(bin_list)
}

#--------------------------------------
#' Function to turn indeces from p-dimensions into one set of unique indeces
#'
#' @description take p-dim indeces and convert fast/slow running indeces
#'
#' @param nbins number of bins in each dimension
#' @param indeces p-dimensional indeces for bin
#' 
#' @return integer index unique to p-dim bin
index_collapser <- function(nbins, indeces){
  p = length(nbins)
  bin_index <- indeces[p]
  for(j in 1:(p-1)){
    bin_index <- bin_index + (indeces[j]-1)*prod(nbins[(j+1):p]) 
  }
  return(bin_index)
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
bin_index_finder_nest <- function(x, bin_def, strict=TRUE){ 
  bin_dim <- length(bin_def$nbins)
  nest_list <- bin_def$bin_list[[1]]
  nest_index <- rep(NA,bin_dim)
  x <- as.numeric(x) 
  for(d in 1:bin_dim){
    bound_vec_idx <- ifelse(d==bin_dim,1,bin_def$nbins[d]+1)
    nest_index[d] <- .bincode(x[[d]], nest_list[[bound_vec_idx]],T,T)
    if(strict == FALSE){
      if( x[[d]] < min(nest_list[[bound_vec_idx]]) ) nest_index[d] <- 1
      if( x[[d]] > max(nest_list[[bound_vec_idx]]) ) nest_index[d] <- bin_def$nbins[d]
    }
    if(length(nest_index[d])==0 |  is.na(nest_index[d])) return(print("Observation outside of observed bins, set strict=FALSE "))
    if(d<bin_dim) nest_list <- nest_list[[nest_index[d]]]
  }
  
  if(bin_dim==1){
    idx <- nest_index
  } else{
    idx <- index_collapser(bin_def$nbins,nest_index)
  }
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
make_stack_matrix <- function(N,M){ 
  mat <- unname(stats::model.matrix(~as.factor(rep(1:N,each=M))-(1)))
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
