#' Kernel function
#'
#' \code{K} returns the value of the kernel function at its argument(s).
#'
#' @param z A vector or matrix containing where the kernel is evaluated. If
#'   it is a matrix each column is interpreted as a vector.
#' @param type A character string specify the used kernel function. Current
#'   options are "Gaussian" and...
#' @return If the input is a vector the value \eqn{K(z)} will be returned. If
#'   the input is a matrix a vector \eqn{K(z_1), \ldots, K(z_m)} will be
#'   returned, where \eqn{z = [z_1, \ldots, z_m]}.
#' @examples
#' z <- rnorm(3)
#' K(z, "Gaussian")
#'
#' Z <- matrix(rnorm(12), 3, 4)
#' K(Z, "Gaussian")
#' @export
K <- function(z, type = "Gaussian"){
  z <- as.matrix(z)
  res <- switch(type,
                Gaussian = K.Gaussian(z),
                Epan     = K.Epan(z),
                Uniform  = K.Uniform(z),
                stop("This Kernel type is not available!
                     (Error occured in K)"))
  return(Kz = res)
}

#' Theoretical moments of the kernel functions
#'
#' List of the theoretical moments of the kernel functions used for the LMME.
#'
#' @param order A number interpreted as multi-index giving the order of the
#'   moment.
#' @param type A character string specify the used kernel function.
#' @return Returns the theoratical moment of the kernel function. Either
#'   \eqn{\mu_{4} = \int K(z) z^4 dz} or \eqn{\mu_{22} = \int K(z) z_1^2 z_2^2
#'   dz}.
KernelMoment <- function(order, type = "Gaussian"){
  all_order <- c(22, 4)
  if(!is.element(order, all_order)){
    stop("This order is not available for the kernel function")
  }
  res <- switch(type,
                Gaussian = c(1, 3),
                Epan     = c(1, 15 / 7),
                Uniform  = c(NA, NA),
                stop("This kernel is not available in the function KernelMoment"))
  return(mu = res[match(order, all_order)])
}




#' Gradient of the kernel function
#'
#' \code{DK} returns the gradient of K at its argument(s).
#'
#' @inheritParams K
#' @return If the input is avector the value \eqn{DK(z)} will be returned as a
#'   vector. If the input is a matrix a matrix \eqn{DK(z_1), \ldots, DK(z_m)}
#'   will be returned, where \eqn{z = [z_1, \ldots, z_m]}.
#' @examples
#' z <- rnorm(3)
#' DK(z, "Gaussian")
#'
#' Z <- matrix(rnorm(12), 3, 4)
#' DK(Z, "Gaussian")
DK <- function(z, type = "Gaussian"){
  z <- as.matrix(z)
  res <- switch(type,
                Gaussian = DK.Gaussian(z),
                Epan     = DK.Epan(z),
                stop("This Kernel type is not available!
                     (Error occured in DK)"))
  return(DKz = res)
}

#' 2nd order derivative of the kernel function
#'
#' \code{D2K} returns the values of the 2nd order derivative of the kernel
#' function at its argument(s).
#'
#' @inheritParams K
#' @return If the input is a vector the vector \eqn{D^{\otimes 2} K(z)} will be
#'   returned. If it is a matrix a matrix \eqn{D^{\otimes 2}K(z_1), \ldots,
#'   D^{\otimes 2}K(z_m)} will be returned, where \eqn{z = [z_1, \ldots, z_m]}.
#' @examples
#' z <- rnorm(3)
#' D2K(z, "Gaussian")
#'
#' Z <- matrix(rnorm(12), 3, 4)
#' D2K(Z, "Gaussian")
D2K <- function(z, type = "Gaussian"){
  z <- as.matrix(z)
  res <- switch(type,
                Gaussian = D2K.Gaussian(z),
                Epan     = D2K.Epan(z),
                stop("This Kernel type is not available!
                     (Error occured in D2K)"))
  return(DK2z = res)
}

#' Epanechnikov Kernel
#'
#' This function is used inside the functions /code{K, DK, D2K} if /code{type =
#' "Epan"}.
#'
#' @param z A matrix containing on which the kernel is evaluated. If
#'   it is a matrix each column is interpreted as a vector.
#' @return If the input is a vector the value \eqn{K(z)} will be returned. If
#'   the input is a matrix a vector \eqn{K(z_1), \ldots, K(z_m)} will be
#'   returned, where \eqn{z = [z_1, \ldots, z_m]}.
K.Epan <- function(z){
  d <- dim(z)[1]
  m <- dim(z)[2]
  const <- 3 * sqrt(5) / 20
  pos_ind     <- which(apply(z, 2, function(y) sum(-sqrt(5) < y & y < sqrt(5)) == d))
  Kz          <- rep(0, m)
  if(length(pos_ind) != 0){
    if(length(pos_ind) == 1){
      Kz <- prod(1 - z[ , pos_ind]^2 / 5) * const^d
    } else{
      Kz[pos_ind] <- apply(z[ , pos_ind], 2, function(y) prod(1 - y^2 / 5)) * const^d
    }
  }
  return(Kz = Kz)
}

#' 1st derivative of the Epanechnikov Kernel
#'
#' This function is used inside the functions \code{DK} if \code{type =
#' "Epan"}.
#'
#' @param z A matrix containing on which the kernel is evaluated. If
#'   it is a matrix each column is interpreted as a vector.
#' @return If the input is a vector the value \eqn{DK(z)} will be returned as a
#'   vector. If the input is a matrix a matrix \eqn{DK(z_1), \ldots, DK(z_m)}
#'   will be returned, where \eqn{z = [z_1, \ldots, z_m]}.
DK.Epan <- function(z){
  d <- dim(z)[1]
  m <- dim(z)[2]
  const <- 3 * sqrt(5) / 20
  pos_ind <- which(apply(z, 2, function(y) sum(-sqrt(5) < y & y < sqrt(5)) == d))
  DKz <- matrix(0, d, m)
  if(length(pos_ind) != 0){
    if(length(pos_ind) == 1){
      tmp <- prod(1 - z[, pos_ind]^2) * const^d
    } else{
      tmp <- apply(z[ , pos_ind], 2, function(y) prod(1 - y^2 / 5)) * const^d
    }
    DKz[ , pos_ind] <- rep(tmp, each = d) * (-2 / 5 * z[ , pos_ind] / (1 - z[ , pos_ind]^2 / 5))
  }
  return(DKz = DKz)
  }

#' 2nd derivative of the Epanechnikov Kernel
#'
#' This function is used inside the function \code{D2K} if \code{type =
#' "Epan"}.
#'
#' @param z A matrix containing on which the kernel is evaluated. If
#'   it is a matrix each column is interpreted as a vector.
#' @return If the input is a vector the vector \eqn{D^{\otimes 2} K(z)} will be
#'   returned. If it is a matrix a matrix \eqn{D^{\otimes 2}K(z_1), \ldots,
#'   D^{\otimes 2}K(z_m)} will be returned, where \eqn{z = [z_1, \ldots, z_m]}.
D2K.Epan <- function(z){
  # only works if d = 2
  d <- dim(z)[1]
  if(d != 2) stop("D2K.Epan only works for d == 2")
  m <- dim(z)[2]
  const <- 3 * sqrt(5) / 20
  pos_ind <- which(apply(z, 2, function(y) sum(-sqrt(5) < y & y < sqrt(5)) == d))
  D2Kz <- matrix(0, d^2, m)
  if(length(pos_ind) != 0){
    if(length(pos_ind) == 1){
      tmp <- prod(1 - z[ , pos_ind]^2 / 5) * const^d
    } else{
      tmp <- apply(z[ , pos_ind], 2, function(y) prod(1 - y^2 / 5)) * const^d
    }
    D2Kz[c(1, 4) , pos_ind] <- rep(tmp, each = 2) * (-2 / 5) /
                      (1 - z[ , pos_ind]^2 / 5)
    D2Kz[c(2, 3) , pos_ind] <- rep(tmp * (-2 / 5 * z[1, pos_ind] / (1 - z[1, pos_ind]^2 / 5) *
                                    -2 / 5 * z[2, pos_ind] / (1 - z[2, pos_ind]^2 / 5)), each = 2)
  }
  return(D2Kz = D2Kz)
}

#' Gaussian Kernel
#'
#' This function is used inside the functions /code{K, DK, D2K} if /code{type =
#' "Gaussian"}.
#'
#' @param z A matrix containing on which the kernel is evaluated. If
#'   it is a matrix each column is interpreted as a vector.
#' @return If the input is a vector the value \eqn{K(z)} will be returned. If
#'   the input is a matrix, a vector \eqn{K(z_1), \ldots, K(z_m)} will be
#'   returned, where \eqn{z = [z_1, \ldots, z_m]}.
K.Gaussian <- function(z){
  d <- dim(z)[1]
  return(Kz = (2 * pi)^(-d / 2) * exp(- 0.5 * colSums(z^2)))
}

#' 1st derivative of the Gaussian Kernel
#'
#' This function is used inside the functions \code{DK} if \code{type =
#' "Gaussian"}.
#'
#' @param z A matrix containing on which the kernel is evaluated. If
#'   it is a matrix each column is interpreted as a vector.
#' @return If the input is avector the value \eqn{DK(z)} will be returned as a
#'   vector. If the input is a matrix a matrix \eqn{DK(z_1), \ldots, DK(z_m)}
#'   will be returned, where \eqn{z = [z_1, \ldots, z_m]}.
DK.Gaussian <- function(z){
  d <- dim(z)[1]
  Kz = (2 * pi)^(-d / 2) * exp(- 0.5 * colSums(z^2))
  return(DKz = - z * rep(Kz, each = d))
}

#' 2nd derivative of the Gaussian Kernel
#'
#' This function is used inside the function \code{D2K} if \code{type =
#' "Gaussian"}.
#'
#' @param z A matrix containing on which the kernel is evaluated. If
#'   it is a matrix each column is interpreted as a vector.
#' @return If the input is a vector the vector \eqn{D^{\otimes 2} K(z)} will be
#'   returned. If it is a matrix a matrix \eqn{D^{\otimes 2}K(z_1), \ldots,
#'   D^{\otimes 2}K(z_m)} will be returned, where \eqn{z = [z_1, \ldots, z_m]}.
D2K.Gaussian <- function(z){
  d <- dim(z)[1]
  D2Kz <- rep(K(z, "Gaussian"), each = d^2) *
    (apply(z, 2, function(z) kronecker(z, z)) - get_vecId(d))
  return(D2Kz = D2Kz)
}

#' Uniform Kernel
#'
#' This function is used inside the functions \code{K} if \code{type =
#' "Uniform"}.
#'
#' @param z A matrix containing where the kernel is evaluated. If
#'   it is a matrix each column is interpreted as a vector.
#' @return If the input is a vector the value \eqn{K(z)} will be returned. If
#'   the input is a matrix a vector \eqn{K(z_1), \ldots, K(z_m)} will be
#'   returned, where \eqn{z = [z_1, \ldots, z_m]}.
K.Uniform <- function(z){
  d <- dim(z)[1]
  return(Kz = as.numeric(sqrt(colSums(z^2)) <= 1))
}

