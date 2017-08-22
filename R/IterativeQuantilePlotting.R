### iterative quantile bin plotting functions

#--------------------------------------
#' Plotting 2-dimensional iterative quantile bins
#'
#' @description used for binning the numeric values by empirical quantile
#'
#' @param iq_obj iterative quantile object: either from iqbin or iqnn function with data and definition
#'
#' @return create ggplot2 object with iq bins displayed
#' @examples

iqbin_plot_2d <- function(iq_obj){
  if(class(iq_obj)=="iqbin"){
    bounds <- as.data.frame(iq_obj$bin_def$bin_bounds)
    training <- iq_obj$bin_data$data
  }
  if(class(iq_obj)=="iqnn") bounds <- as.data.frame(iq_obj$bin_bounds)
  
  p1 <- ggplot() +
    geom_rect(aes(xmin=V1, xmax=V2, ymin=V3, ymax=V4),color="black",fill=NA, data=bounds)+
    theme_bw()
  if(class(iq_obj)=="iqbin") p1 <- p1 + geom_point(aes_string(x=iq_obj$bin_def$bin_cols[1],y=iq_obj$bin_def$bin_cols[2]), data=iq_obj$bin_data$data)
  p1
}

iq_obj <- iqbin(data=iris, bin_cols=c("Sepal.Length","Sepal.Width"),
                    nbins=c(5,3), output="both",jit=rep(0.001,2))
str(iq_obj)
class(iq_obj)

iq_obj <- iqnn(data=iris, y="Species", mod_type="class", bin_cols=c("Sepal.Length","Sepal.Width"),
                 nbins=c(5,3), jit=rep(0.001,2))
str(iq_obj)
class(iq_obj)
