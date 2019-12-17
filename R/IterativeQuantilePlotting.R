### iterative quantile bin plotting functions

#--------------------------------------
#' Plotting 2-dimensional iterative quantile bins
#'
#' @description used for binning the numeric values by empirical quantile
#'
#' @param iq_obj iterative quantile object: either from iqbin or iqnn function with data and definition
#'
#' @return create ggplot2 object with iq bins displayed
#' @export
#' @examples
#' iq_obj <- iqbin(data=iris, bin_cols=c("Sepal.Length","Sepal.Width"), 
#'                 nbins=c(5,3), output="both",jit=rep(0.1,2))
#' iqbin_plot_2d(iq_obj)
#' 
#' iq_obj <- iqbin(data=iris, bin_cols=c("Sepal.Length","Sepal.Width"), 
#'                 nbins=c(5,3), output="both")
#' iqbin_plot_2d(iq_obj)
#' 
#' iq_obj <- iqnn(data=iris, y="Species", mod_type="class", 
#'                bin_cols=c("Sepal.Length","Sepal.Width"), 
#'                nbins=c(5,3), jit=rep(0.001,2))
#' iqbin_plot_2d(iq_obj)
#' 
#' iq_obj <- iqnn(data=iris, y="Petal.Length", mod_type="reg", 
#'                bin_cols=c("Sepal.Length","Sepal.Width"), 
#'                nbins=c(5,3), jit=rep(0.001,2))
#' iqbin_plot_2d(iq_obj)

iqbin_plot_2d <- function(iq_obj){
  if(attributes(iq_obj)$iq_obj_type == "iqbin"){
    bounds <- as.data.frame(iq_obj$bin_def$bin_bounds)
    train_points <- iq_obj$bin_data$data[,iq_obj$bin_def$bin_cols] 
    if(!is.null(iq_obj$bin_data$bin_jit)) train_points <- train_points + iq_obj$bin_data$bin_jit
    bin_cols <- iq_obj$bin_def$bin_cols
  }
  if(attributes(iq_obj)$iq_obj_type == "iqnn"){
    bounds <- data.frame(as.data.frame(iq_obj$bin_bounds),iq_obj$bin_stats)
    bin_cols <- iq_obj$bin_cols
  }
  
  p1 <- ggplot2::ggplot() +
    ggplot2::theme_bw()
  if(attributes(iq_obj)$iq_obj_type == "iqbin"){
    p1 <- p1 + ggplot2::geom_rect(ggplot2::aes(xmin=train_points$V1, xmax=train_points$V2,
                                               ymin=train_points$V3, ymax=train_points$V4), 
                                  data=bounds,color="black",fill=NA) +
      ggplot2::geom_point(ggplot2::aes_string(x=bin_cols[1],y=bin_cols[2]),data=train_points)
  } 
  if(attributes(iq_obj)$iq_obj_type == "iqnn"){
    p1 <- p1 + ggplot2::geom_rect(ggplot2::aes(xmin=train_points$V1, xmax=train_points$V2,
                                               ymin=train_points$V3, ymax=train_points$V4,
                                               fill=train_points$pred),
                                  data=bounds,color="black")+
      ggplot2::labs(x=bin_cols[1],y=bin_cols[2])
    if(iq_obj$mod_type=="class") p1 <- p1 + ggplot2::scale_fill_discrete(paste0("Predicted \n",iq_obj$y))
    if(iq_obj$mod_type=="reg") p1 <- p1 + ggplot2::scale_fill_continuous(paste0("Predicted \n",iq_obj$y),low="gray90",high="gray20")
  } 
  p1
}

#--------------------------------------
#' Plotting 3-dimensional iterative quantile bins
#'
#' @description Plotting 3-dimensional iterative quantile bins using faceting
#'
#' @param iq_obj iterative quantile object: either from iqbin or iqnn function with data and definition
#' @param round_digits decimal places to round the data used in structure of plot
#'
#' @return create ggplot2 object with iq bins displayed
#' @export
#' @examples
#' iq_obj <- iqbin(data=iris, bin_cols=c("Sepal.Length","Sepal.Width","Petal.Width"), 
#'                 nbins=c(2,5,3), output="both",jit=rep(0.1,3))
#' iqbin_plot_3d(iq_obj)
#' 
#' iq_obj <- iqnn(data=iris, y="Species", mod_type="class", 
#'                bin_cols=c("Sepal.Length","Sepal.Width","Petal.Width"), 
#'                nbins=c(2,5,3), jit=rep(0.001,3))
#' iqbin_plot_3d(iq_obj)
#'
#' iq_obj <- iqnn(data=iris, y="Petal.Length", mod_type="reg", 
#'                bin_cols=c("Sepal.Length","Sepal.Width","Petal.Width"), 
#'                nbins=c(2,5,3), jit=rep(0.001,3))
#' iqbin_plot_3d(iq_obj)

iqbin_plot_3d <- function(iq_obj, round_digits=2){
  if(attributes(iq_obj)$iq_obj_type == "iqbin"){
    bin_cols <- iq_obj$bin_def$bin_cols
    nbins <- iq_obj$bin_def$nbins
    bounds <- as.data.frame(iq_obj$bin_def$bin_bounds)
    train_points <- iq_obj$bin_data$data[,iq_obj$bin_def$bin_cols] 
    if(!is.null(iq_obj$bin_data$bin_jit)) train_points <- train_points + iq_obj$bin_data$bin_jit
    train_points$bin_index <- iq_obj$bin_data$data[,"bin_index"]
    closed_lower <- ifelse(bounds$V1[train_points$bin_index]==min(bounds$V1[train_points$bin_index]),"[","(")
    train_points$x1 <- paste0(bin_cols[1]," in ",closed_lower,round(bounds$V1,round_digits)[train_points$bin_index],",",round(bounds$V2,round_digits)[train_points$bin_index],"]")
  }
  if(attributes(iq_obj)$iq_obj_type == "iqnn") {
    bin_cols <- iq_obj$bin_cols
    nbins <- iq_obj$nbins
    bounds <- data.frame(as.data.frame(iq_obj$bin_bounds),iq_obj$bin_stats)
  }
  lower_bounds <- c(rep("[",sum(bounds$V1==min(bounds$V1))),rep("(",nrow(bounds)-sum(bounds$V1==min(bounds$V1))))
  bounds$x1 <- paste0(bin_cols[1]," in ",lower_bounds,round(bounds$V1,round_digits),",",round(bounds$V2,round_digits),"]")
  bounds$x1 <- factor(bounds$x1, levels=unique(bounds$x1))
  if(attributes(iq_obj)$iq_obj_type == "iqbin") train_points$x1 <- factor(train_points$x1, levels=levels(bounds$x1))
  
  p1 <- ggplot2::ggplot() +
    ggplot2::theme_bw()
  if(attributes(iq_obj)$iq_obj_type == "iqbin"){
    p1 <- p1 + ggplot2::geom_rect(ggplot2::aes(xmin=train_points$V3, xmax=train_points$V4,
                                               ymin=train_points$V5, ymax=train_points$V6),
                                  color="black",fill=NA, data=bounds)+
      ggplot2::facet_grid(.~x1)+
      ggplot2::geom_point(ggplot2::aes_string(x=bin_cols[2],y=bin_cols[3]), 
                 data=train_points)+
      ggplot2::labs(title="Iterative Quantile Bins Defined by Training Data",
                    subtitle=paste(paste(bin_cols,collapse=" X ")," (",paste(nbins,collapse="X"),")",sep="") )
  } 
  if(attributes(iq_obj)$iq_obj_type == "iqnn"){
    p1 <- p1 + ggplot2::geom_rect(ggplot2::aes(xmin=train_points$V3, xmax=train_points$V4, 
                                               ymin=train_points$V5, ymax=train_points$V6,
                                               fill=train_points$pred),color="black", data=bounds)+
      ggplot2::facet_grid(.~x1) + 
      ggplot2::labs(x=iq_obj$bin_cols[2],y=iq_obj$bin_cols[3],
                    title="Iterative Quantile Nearest Neighbor Predicted Responses",
                    subtitle=paste(paste(bin_cols,collapse=" X ")," (",paste(nbins,collapse="X"),")",sep="") )+
      ggplot2::theme(legend.position="bottom")
    if(iq_obj$mod_type=="class") p1 <- p1 +       ggplot2::scale_fill_discrete(paste0("Predicted \n",iq_obj$y))
    if(iq_obj$mod_type=="reg") p1 <- p1 + ggplot2::scale_fill_continuous(paste0("Predicted \n",iq_obj$y),low="gray90",high="gray20")
  } 
  p1
}

