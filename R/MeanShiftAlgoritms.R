#' Mode Seeking Algorithm
#'
#' Iterative mean shift algorithm detecting the mode of the underlaying pdf.
#'
#' The algorithm uses LHSE, LMME or KDE for finding the modes of the underlaying
#' pdf.
#'
#' @param data Data matrix [d, n].
#' @param h Kernel bandwidth parameter
#' @param method Used method for the log-density estimation. Either \code{LHSE},
#'   \code{LMME} or \code{KDE}.
#' @param type Used kernel. Either \code{Gaussian} or \code{Epan}-(echikov).
#' @param startpoints Matrix [d, m] containing in each raw a starting point for
#'   the algorithm. If \code{NULL} the data points will be used.
#' @param delta Stopping criteria.
#' @param iter_max Maximal number of iterations before the procedure stops. Will
#'   print a warning if it did not converge.
#' @return Fixpoints of the iteration as matrix [d, m]
ModeSeekingVersion1 <- function(data, h, method = "LHSE", type = "Gaussian", startpoints = NULL,
                        delta = 1e-10, iter_max = Inf){
  if (is.null(startpoints)) startpoints <- data
  d <- dim(data)[1]
  range <- max(range(data[1, ]), range(data[2, ]))
  startpoints <- as.matrix(startpoints)
  m <- dim(startpoints)[2]
  iter <- 1
  points <- as.matrix(startpoints)
  index <- rep(TRUE, m)
  while ((sum(index) > 0) & (iter <= iter_max)){
    theta1 <- as.matrix(LogDensity(points[ , index], h, type, method, data)[1:d, ])
    ddelta <- apply(theta1, 2, function(x) sqrt(sum(x^2)))
    # naive stepsize correction
    theta1 <- theta1 * pmin(range / ddelta, 1)
    points[ , index] <- as.matrix(points[ , index] + theta1)
    index[index] <- (ddelta > delta)
    print(paste("Iteration:" , iter, "Points not yet converge:", sum(index)))
    p <- PlotData(data) + geom_point(aes(x, y), data = data.frame(x = points[1, ],
                                                                 y = points[2, ]), color = "red")
    print(p)
    iter <- iter + 1
  }
  return(points)
}

#' Mode Seeking Algorithm
#'
#' Gradient acent with step size correction according to the direction of the gradient
#'
#' @param x Starting point [d]
#' @param h Kernel bandwidth parameter
#' @param method Used method for the log-density estimation. Either \code{LHSE},
#'   \code{LMME} or \code{KDE}.
#' @param type Used kernel. Either \code{Gaussian} or \code{Epan}-(echikov).
#' @param data Data matrix [d, n].
#' @param delta Stopping criteria.
#' @return Fixpoints of the iteration [d]
#' @export
ModeSeeking <- function(x, h, type = "Gaussian", method = "LHSE", data, delta = 1e-10){
  d <- length(x)
  step_size <- 2 * delta
  while (step_size > delta){
    step <- LogDensity(x, h, type, method, data)[1:d]
    x_new <- x + step
    step2 <- LogDensity(x_new, h, type, method, data)[1:d]
    direction <- sign(sum(step * step2))
    while (direction < 0){
      x_new <- (x + x_new) / 2
      step2 <- LogDensity(x_new, h, type, method, data)[1:d]
      direction <- sign(sum(step * step2))
    }
    step_size <- sqrt(sum((x_new - x)^2))
    x <- x_new
  }
  return(x)
}

#' Ridge Seeking Algorithm
#'
#' Gradient acent with step size correction according to the direction of the gradient
#'
#' @param x Starting point [d]
#' @param h Kernel bandwidth parameter
#' @param method Used method for the log-density estimation. Either \code{LHSE},
#'   \code{LMME} or \code{KDE}.
#' @param type Used kernel. Either \code{Gaussian} or \code{Epan}-(echikov).
#' @param data Data matrix [d, n].
#' @param delta Stopping criteria.
#' @param iter_max Maximal number of iterations before the procedure stops. Will
#'   print a warning if it did not converge.
#' @param path \code{TRUE} if the whole path of the iteration should be returned.
#' @return Fixpoints of the iteration as vector [d]. If \code{path = TRUE} it returns a list with
#' \describe{\item{x}{Fixpoint [d] of the iteration.}
#' \item{path}{Path of the iteration [d, iter], where iter is the number of steps.}}
#' @export
RidgeSeeking_V1 <- function(x, h, type = "Gaussian", method = "LHSE", data, delta = 1e-6, iter_max = Inf, path = FALSE){
  path_values <- x
  iter <- 1
  d <- length(x)
  step_size <- 2 * delta
  while (step_size > delta | iter <= iter_max){
    tmp <- LogDensity(x, h, type, method, data)
    step <- tmp[1:d]
    x_new <- x + step
    step2 <- LogDensity(x_new, h, type, method, data)[1:d]
    direction <- sign(sum(step * step2))
    while(direction < 0){
      x_new <- (x + x_new) / 2
      step2 <- LogDensity(x_new, h, type, method, data)[1:d]
      direction <- sign(sum(step * step2))
    }
    step_projected <- SubspaceProjection(x_new - x, tmp[(d+1):length(tmp)])
    step_size <- sqrt(sum(step_projected^2))
    x <- x + step_projected
    path_values <- cbind(path_values, x)
    iter <- iter + 1
  }
  if (path == TRUE){
    return(list(x = x, path = path_values))
  } else{
    return(x)
  }
}

#' Ridge Seeking Algorithm
#'
#' Gradient acent with step size correction according to the direction of the gradient
#'
#' @param x Starting point [d]
#' @param h Kernel bandwidth parameter
#' @param method Used method for the log-density estimation. Either \code{LHSE},
#'   \code{LMME} or \code{KDE}.
#' @param type Used kernel. Either \code{Gaussian} or \code{Epan}-(echikov).
#' @param data Data matrix [d, n].
#' @param delta Stopping criteria.
#' @param iter_max Maximal number of iterations before the procedure stops. Will
#'   print a warning if it did not converge.
#' @param path \code{TRUE} if the whole path of the iteration should be returned.
#' @return Fixpoints of the iteration as vector [d]. If \code{path = TRUE} it returns a list with
#' \describe{\item{x}{Fixpoint [d] of the iteration.}
#' \item{path}{Path of the iteration [d, iter], where iter is the number of steps.}}
#' @export
RidgeSeeking <- function(x, h, type = "Gaussian", method = "LHSE", data, delta = 1e-6, iter_max = Inf, path = FALSE){
  path_values <- x
  iter <- 1
  d <- length(x)
  step_size <- 2 * delta
  while (step_size > delta){
    tmp <- LogDensity(x, h, type, method, data)
    step <- tmp[1:d]
    step_projected <- SubspaceProjection(x + step, tmp[(d+1):length(tmp)])
    v <- step_projected / sqrt(sum(step_projected^2))
    a  <- 0
    b  <- 1
    tmp <- LogDensity(x + b * step_projected, h, type, method, data)
    Dv_b <- sum(v * tmp[1:d])
    while (Dv_b > 0 & !is.na(Dv_b)){
      a <- b
      b <- 2 * b
      tmp <- LogDensity(x + b * step_projected, h, type, method, data)
      Dv_b <- sum(v * tmp[1:d])
    }
    while (b - a > delta / 10){
      m <- mean(c(a, b))
      tmp <- LogDensity(x + m * step_projected, h, type, method, data)
      Dv_m <- sum(v * tmp[1:d])
      if (Dv_m > 0 & !is.na(Dv_m)){
        a <- m
      } else{
        b <- m
      }
    }
    m <- mean(c(a, b))
    x <- x + m * step_projected
    path_values <- cbind(path_values, x)
    step_size <- m * sqrt(sum(step_projected^2))
    if (iter >= iter_max){
      warning(paste0("RidgeSeeking stopped after ", iter, "iterations!"))
      return(list(x = x, path = path_values))
    }
    iter <- iter + 1
  }
  if (path == TRUE){
    return(list(x = x, path = path_values))
  } else{
    return(x)
  }
}

#' Subspace Projection
#'
#' @param x Vector [d]
#' @param theta2 The vectorized Hessian of the log-density [d^2]
#' @return
#' Projects a vector on the d-1 dimensional subspace spanned by the d-1
#' eigenvectors corresponding to the d-1 smallest eigenvalues of the Hessian of
#' the log-density.
#' @export
SubspaceProjection <- function(x, theta2){
  d <- length(x)
  V <- eigen(matrix(theta2, d), symmetric = TRUE)$vectors[ , 2:d, drop = FALSE]
  projection <- V %*% t(V) %*% as.matrix(x)
  return(as.vector(projection))
}
