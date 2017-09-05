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
#' 2+2
iqbin_plot_2d <- function(iq_obj){
  if(class(iq_obj)=="iqbin"){
    bounds <- as.data.frame(iq_obj$bin_def$bin_bounds)
    train <- iq_obj$bin_data$data[,iq_obj$bin_def$bin_cols]
    bin_jit <- iq_obj$bin_data$bin_jit
    train_jit <- train + bin_jit
  }
  if(class(iq_obj)=="iqnn") bounds <- data.frame(as.data.frame(iq_obj$bin_bounds),iq_obj$bin_stats)
  
  p1 <- ggplot() +
    theme_bw()
  if(class(iq_obj)=="iqbin"){
    p1 <- p1 + geom_rect(aes(xmin=V1, xmax=V2, ymin=V3, ymax=V4),color="black",fill=NA, data=bounds)+
      geom_point(aes_string(x=iq_obj$bin_def$bin_cols[1],
                                     y=iq_obj$bin_def$bin_cols[2]), 
                          data=train_jit)
  } 
  if(class(iq_obj)=="iqnn"){
    p1 <- p1 + geom_rect(aes(xmin=V1, xmax=V2, ymin=V3, ymax=V4,fill=pred),color="black", data=bounds)+
      labs(x=iq_obj$bin_cols[1],y=iq_obj$bin_cols[2])
    if(iq_obj$mod_type=="reg") p1 <- p1 + scale_fill_continuous(paste0("Predicted \n",iq_obj$y),low="gray90",high="gray20")
  } 
  p1
}
# iq_obj1 <- iqbin(data=iris, bin_cols=c("Sepal.Length","Sepal.Width"),
# nbins=c(5,3), output="both",jit=rep(0.1,2))
# str(iq_obj1)
# class(iq_obj1)
# 
# iq_obj2 <- iqnn(data=iris, y="Species", mod_type="class", bin_cols=c("Sepal.Length","Sepal.Width"),
#                nbins=c(5,3), jit=rep(0.001,2))
# str(iq_obj2)
# class(iq_obj2)
# 
# iq_obj <- iqnn(data=iris, y="Petal.Length", mod_type="reg", bin_cols=c("Sepal.Length","Sepal.Width"),
#                 nbins=c(5,3), jit=rep(0.001,2))
# str(iq_obj)
# class(iq_obj)

# iq_obj <- iqbin(data=iris, bin_cols=c("Sepal.Length","Sepal.Width","Petal.Width"),
# nbins=c(2,5,3), output="both",jit=rep(0.1,3))
# str(iq_obj)
# class(iq_obj)
# 
# iq_obj <- iqnn(data=iris, y="Species", mod_type="class", bin_cols=c("Sepal.Length","Sepal.Width","Petal.Width"),
#                nbins=c(2,5,3), jit=rep(0.001,3))
# str(iq_obj)
# class(iq_obj)
# 
# iq_obj <- iqnn(data=iris, y="Petal.Length", mod_type="reg", bin_cols=c("Sepal.Length","Sepal.Width","Petal.Width"),
#                nbins=c(2,5,3), jit=rep(0.001,3))
# str(iq_obj)
# class(iq_obj)

iqbin_plot_3d <- function(iq_obj, round_digits=2){
  if(class(iq_obj)=="iqbin"){
    bounds <- as.data.frame(iq_obj$bin_def$bin_bounds)
    bounds$x1 <- paste0(iq_obj$bin_cols[1]," in (",round(bounds$V1,round_digits),",",round(bounds$V2,round_digits),"]")
    bounds$x1 <- factor(bounds$x1, levels=unique(bounds$x1))
    train <- iq_obj$bin_data$data[,iq_obj$bin_def$bin_cols]
    bin_jit <- iq_obj$bin_data$bin_jit
    train_jit <- train + bin_jit
  }
  if(class(iq_obj)=="iqnn") {
    bounds <- data.frame(as.data.frame(iq_obj$bin_bounds),iq_obj$bin_stats)
    bounds$x1 <- paste0(iq_obj$bin_cols[1]," in (",round(bounds$V1,round_digits),",",round(bounds$V2,round_digits),"]")
    bounds$x1 <- factor(bounds$x1, levels=unique(bounds$x1))
  }
  
  p1 <- ggplot() +
    theme_bw()
  if(class(iq_obj)=="iqbin"){
    p1 <- p1 + geom_rect(aes(xmin=V1, xmax=V2, ymin=V3, ymax=V4),color="black",fill=NA, data=bounds)+
      geom_point(aes_string(x=iq_obj$bin_def$bin_cols[1],
                            y=iq_obj$bin_def$bin_cols[2]), 
                 data=train_jit)
  } 
  if(class(iq_obj)=="iqnn"){
    p1 <- p1 + geom_rect(aes(xmin=V3, xmax=V4, ymin=V5, ymax=V6,fill=pred),color="black", data=bounds)+
      facet_grid(.~x1)+
      labs(x=iq_obj$bin_cols[2],y=iq_obj$bin_cols[3])
    if(iq_obj$mod_type=="reg") p1 <- p1 + scale_fill_continuous(paste0("Predicted \n",iq_obj$y),low="gray90",high="gray20")
  } 
  p1
}

