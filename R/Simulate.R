#' Simulated circle data
#'
#' Gives simulated 2-dimensional data which is uniformally distributed on a
#' circle with Gaussian noise.
#'
#' @param n Sample size.
#' @param r Radius of the circle.
#' @param sd Standard deviation of the Gaussian noise.
#' @examples
#' data <- SimCircle(200)
#' PlotData(data)
#' @export
SimCircle <- function(n,r = 1, sd = 0.1){
  unif <- runif(n)
  x    <- cos(2 * pi * unif) + rnorm(n, sd = sd)
  y   <- sin(2 * pi * unif) + rnorm(n, sd = sd)
  data <- rbind(X1 = x, X2 = y)

  return(data = data)
}
