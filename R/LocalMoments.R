#' Sample local 0-moment
#'
#' \code{LM0} returns the value of the sample local 0-moment of the data
#' points.
#'
#' @param x A vector or matrix containing the points at which the moment will be
#'   calculated.
#' @param h A number specifying the bandwidth.
#' @param type A character string specifying the used kernel function.
#' @param data A data matrix containing in each column a data point.
#' @return If \code{x} is a vector returns the value \eqn{s_n^0(x)}. If the
#'   input is a matrix a vector \eqn{s_n^0(x_1), \ldots, s_n^0(x_m)} is
#'   returned, where \eqn{x = [x_1, \ldots, x_m]}.
#' @examples
#' data <- matrix(rnorm(20),2,10)
#' x <- matrix(c(seq(-2, 2, length.out = 9), rep(0, 9)), 2, 9)
#' LM0(x, 0.5, "Gaussian", data)
LM0 <- function(x, h, type = "Gaussian", data){
  x <- as.matrix(x)
  if(dim(x)[1] != dim(data)[1]){
    stop("The dimensions of the data and points do not match.")
  }
  d <- dim(x)[1]
  n <- dim(data)[2]
  m <- dim(x)[2]
  s0 <- apply(x, 2, function(x) mean(K((data - x) / h, type = type)) / h^d)
  return(s0 = s0)
}

#' Sample local 1-moments
#'
#' \code{LM1} returns a vector with the 1st order sample local moments.
#'
#' @inheritParams LM0
#' @return If \code{x} is a vector returns the vector \eqn{s_n^{\otimes 1}(x)}. If the
#'   input is a matrix the matrix \eqn{s_n^{\otimes 1}(x_1), \ldots, s_n^{\otimes 1}(x_m)} is
#'   returned, where \eqn{x = [x_1, \ldots, x_m]}.
#' @examples
#' data <- matrix(rnorm(20),2,10)
#' x <- matrix(c(seq(-2, 2, length.out = 9), rep(0, 9)), 2, 9)
#' LM1(x, 0.5, "Gaussian", data)
LM1 <- function(x, h, type = "Gaussian", data){
  x <- as.matrix(x)
  if(dim(x)[1] != dim(data)[1]){
    stop("The dimensions of the data and points do not match.")
  }
  d  <- dim(x)[1]
  n  <- dim(data)[2]
  m  <- dim(x)[2]
  s1 <- apply(x, 2, function(x) ((data - x) / h)
                   %*% matrix(K((data - x) / h, type = type), ncol = 1)) / n / h^d
  return(s1 = s1)
}

#' Sample local 2-moments
#'
#' \code{LM2} returns a vector with the 2nd order sample local moments.
#'
#' @inheritParams LM0
#' @return If \code{x} is a vector returns the vector \eqn{s_n^{\otimes 2}(x)}. If the
#'   input is a matrix the matrix \eqn{s_n^{\otimes 2}(x_1), \ldots, s_n^{\otimes 2}(x_m)} is
#'   returned, where \eqn{x = [x_1, \ldots, x_m]}.
#' @examples
#' data <- matrix(rnorm(20),2,10)
#' x <- matrix(c(seq(-2, 2, length.out = 9), rep(0, 9)), 2, 9)
#' LM2(x, 0.5, "Gaussian", data)
LM2 <- function(x, h, type = "Gaussian", data){
  x <- as.matrix(x)
  if(dim(x)[1] != dim(data)[1]){
    stop("The dimensions of the data and points do not match.")
  }
  d <- dim(x)[1]
  n <- dim(data)[2]
  m <- dim(x)[2]
  s2 <- matrix(0, d^2, m)
  for(j in 1:m){
    data_local <- (data - x[ , j]) / h
    s2[ , j]   <- rowMeans(apply(data_local, 2, function(x) kronecker(x, x) *
                           K(x, type = type))) / h^d
  }
  return(s2 = s2)
}

#' Sample local 0-moment of \eqn{DK}
#'
#' Returns the sample local 0-moment of the 1st order derivative of the kernel
#' function.
#'
#' @inheritParams LM0
#' @return If \code{x} is a vector it returns the vector \eqn{Q_n^{\otimes
#'   0}(x)}. If \code{x} is a matrix it returns the matrix \eqn{Q_n^{\otimes
#'   0}(x_1), \ldots, Q_n^{\otimes 0}(x_m)}, where \eqn{x = [x_1, \ldots, x_m]}.
LM0_DK <- function(x, h, type = "Gaussian", data){
  x  <- as.matrix(x)
  d  <- dim(x)[1]
  q0 <- apply(x, 2, function(x) rowMeans(DK((data - x) / h, type = type))) / h^(d + 1)
  return(q0 = q0)
}

#' Sample local 1-moment of \eqn{DK}
#'
#' Returns the sample local 1-moment of the 1st order derivatives of the kernel
#' function.
#'
#' @inheritParams LM0
#' @return If \code{x} is a vector it returns the vector \eqn{Q_n^{\otimes
#'   1}(x)}. If \code{x} is a matrix it returns the matrix \eqn{Q_n^{\otimes
#'   1}(x_1), \ldots, Q_n^{\otimes 1}(x_m)}, where \eqn{x = [x_1, \ldots, x_m]}.
LM1_DK <- function(x, h, type = "Gaussian", data){
  x  <- as.matrix(x)
  d  <- dim(x)[1]
  q1 <- apply(x, 2, function(x) rowMeans(apply(data, 2, function(y) kronecker(DK((y - x) / h, type = type),
                                                   (y - x) / h)))) / h^(d + 1)
  return(q1 = q1)
}

#' Sample local 0-moment of \eqn{D^{\otimes 2}K}
#'
#' Returns the sample local 0-moment of the 2nd order derivatives of the kernel
#' function.
#'
#' @inheritParams LM0
#' @return If \code{x} is a vector it returns the vector \eqn{\frac{1}{nh^2}
#'   \sum_{i=1}^n (D^{\otimes 2}K)_h (X_i - x)}. If \code{x} is a matrix it
#'   returns the matrix \eqn{\frac{1}{nh^2} \sum_{i=1}^n (D^{\otimes 2}K)_h (X_i
#'   - x_1), \ldots, \frac{1}{nh^2} \sum_{i=1}^n (D^{\otimes 2}K)_h (X_i -
#'   x_m)}, where \eqn{x = [x_1, \ldots, x_m]}.
LM0_D2K <- function(x, h, type = "Gaussian", data){
  x  <- as.matrix(x)
  d  <- dim(x)[1]
  Q0 <- apply(x, 2, function(x) rowMeans(D2K((data - x) / h, type = type))) / h^(d + 2)
  return(Q0 = Q0)
}
