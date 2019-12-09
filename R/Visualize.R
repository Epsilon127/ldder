#' Matrix to data frame
#'
#' Transforms the data matrix into a data frame usealbe by ggplot 2
#'
#' @param data A matrix where each column represents a data point in R^d
#' @return A data frame where each row represents a data point and each column
#'   contains one component of each data point. The colums are named with X1, X2,...
MatToDf <- function(data){
  df <- data.frame(t(data))
}

#' Plot raw data
#'
#' Plots the raw data projected in 2 dimensions.
#'
#' @export
#' @param data A matrix containing in each column a data point or a Dataframe
#' @return Plots the raw data.
PlotData <- function(data, ...){
  if(!is.data.frame(data)){
    data <- as.matrix(data)
    d    <- dim(data)[1]
    df   <- MatToDf(data)
    names(df) <- c("X1", "X2")
  }
  if(d == 1){
    warning("Plots for univariate data aren't yet available...")
  }
  if(d == 2){
    # plot(df[,1],df[,2])
        p <- ggplot(df, aes(X1, X2)) + geom_point(...)
  }
  if(d > 2){
    warning("Plots for more then two dimensions aren't yet available...")
  }
  return(p)
}
#' Colorplot of the estimator
#'
#' gives a colorplot of the estimator. This means 6 plots (d = 2), each
#' representing one component of the estimator and the color codes the value of
#' component
#'
#' @export
#' @param est The estimator as a dataframe. This is the first element returned
#'   by LogDensity2.
#' @param data The data matrix
#' @param max The maximum value of the theta. All values with absolut value
#'   bigger then \code{max} will be set to \code{max}. That makes it easier to
#'   see the different in the color where theta is small.
#' @return The plot as a ggplot2-object.
PlotEst <- function(est, max = NULL, color = NULL, data = NULL, ...){
  ds_long <- gather(est, key = theta, value = value, theta1:theta6)
  if(!is.null(max)){
    ds_long$value.max <- pmin(ds_long$value, max)
    ds_long$value.max <- pmax(ds_long$value.max, -max)
  } else{
    ds_long$value.max <- ds_long$value
  }
  p <- ggplot(ds_long, aes(X1, X2)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    facet_wrap(~theta, nrow = 2, dir = "v") +
    geom_tile(aes(fill = value.max))
  if(is.null(color)){
    p <- p + scale_fill_gradient2(low = "darkblue", high = "darkred", na.value = "grey50")
  } else{
    p <- p + scale_fill_gradient2(low = color[1], mid = color[2], high = color[3], na.value = "grey50")
  }
  if(!is.null(data)){
    arg <- list(...)
    p <- p + geom_point(data = MatToDf(data), aes(X1, X2), size = 0.5, alpha = 0.6)
    if(length(arg) > 0){
      p <- p + geom_point(data = MatToDf(data), aes(X1, X2), ...)
    }
  }
  return(p)
}

geom_circle_ridge <- function(r = 1, sd = 0.1, n_points = 360){
  r_ridge <- circle_radius(r, sd)
  geom_circle(aes(x0 = 0, y0 = 0, r = r_ridge), n = n_points)
}
