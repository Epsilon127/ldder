#' Local Estimator of the log-density
#'
#' Returns the local estimator of \eqn{\hat{\boldsymbols{\theta}}} for the
#' location(s) \eqn{x}
#'
#' @export
#' @param x A vector or matrix containing the points at which the estimator will be
#'   calculated.
#' @param h A number specifying the bandwidth.
#' @param type A character string specifying the used kernel function.
#' @param method A cherecter string specifying the uesed method. The possible
#'   methods are \code{LHSE, LLLE, LMME, KDE}.
#' @param data A data matrix containing in each column a data point.
#' @return Returnes the the value \eqn{\hat{\boldsymbols{\theta}}} at \code{x}.
#'   Eihter as a vector or matrix.
LogDensity <- function(x, h, type = "Gaussian", method = "LHSE", data){
  x   <- as.matrix(x)
  res <- switch(method,
                LHSE  = LSME (x, h, type, data),
                LSME  = LSME (x, h, type, data),
                LLLE  = LLLE (x, h, type, data),
                LMME  = LMME (x, h, type, data),
                LMME2 = LMME2(x, h, type, data),
                KDE   = KDE  (x, h, type, data),
                stop("This method is not available!"))
  return(theta = res)
}

#' Local Estimator of the log-density
#'
#' Returns the local estimator of \eqn{\hat{\boldsymbols{\theta}}} for the
#' location(s) \eqn{x} and the input arguments as a list object. This may be
#' useful for further use as plotting.
#'
#' @export
#' @param x A vector or matrix containing the points at which the moment will be
#'   calculated.
#' @param h A number specifying the bandwidth.
#' @param type A character string specifying the used kernel function.
#' @param method A cherecter string specifying the uesed method. The possible
#'   methods are \code{LSME, LLLE, LMME, KDE}.
#' @param data A data matrix containing in each column a data point.
#' @return Returnes the the value \eqn{\hat{\boldsymbols{\theta}}} at \code{x}.
#'   Eihter as a vector or matrix.
LogDensity2 <- function(x, h, type = "Gaussian", method = "LSME", data){
  x   <- as.matrix(x)
  res <- switch(method,
                LSME = LSME(x, h, type, data),
                LHSE = LSME(x, h, type, data),
                LLLE = LLLE(x, h, type, data),
                LMME = LMME(x, h, type, data),
                LMME2 = LMME2(x, h, type, data),
                KDE  = KDE (x, h, type, data),
                stop("This method is not available!"))
  d <- dim(x)[1]
  rownames(res) <- paste0("theta", 1:(d + d^2))
  return(list(theta = data.frame(t(x), t(res)), method = method,
              h = h, type = type, data = data))
}

#' Local Estimator of the log-density with LSME
#'
#' Returns the local estimator of \eqn{\hat{\boldsymbols{\theta}}} for the
#' location(s) \eqn{x} calculated with the LSME
#'
#' @param x A vector or matrix containing the points at which the moment will be
#'   calculated.
#' @param h A number specifying the bandwidth.
#' @param type A character string specifying the used kernel function.
#' @param data A data matrix containing in each column a data point.
#' @return Returnes a list with
#' \describe{
#' \item{theta}{Dataframe containig the evaluated points and estimators.}
#' \item{method}{Used method.}
#' \item{h}{Used bandwidth.}
#' \item{type}{Used kernel type.}
#' \item{data}{Data [d, n].}}
LSME <- function(x, h, type = "Gaussian", data){
  tmp <- LHSE2(x, h, type, data)
  return(theta = tmp$theta)
}
#' Local Hyvärinen Score estimator used for calculating p-values
#'
#' Returns the local estimator of \eqn{\hat{\boldsymbols{\theta}}} for the
#' location(s) \eqn{x} calculated with the LSME and S.star for the caluculation
#' of the p-values.
#'
#' In case the local 0-moment at \code{x} is equal to 0, the corresponding theta
#' will be set to \code{NA}
#'
#' @param x A vector or matrix containing the points at which the moment will be
#'   calculated.
#' @param h A number specifying the bandwidth.
#' @param type A character string specifying the used kernel function.
#' @param data A data matrix containing in each column a data point.
#' @return Returnes the the value \eqn{\hat{\boldsymbols{\theta}}} at \code{x}.
#'   Eihter as a vector or matrix and s0, s1, s2 either as a vector or a matrix
#'   gathered in a list-object.
LHSE2 <- function(x, h, type = "Gaussian", data){
  if(dim(as.matrix(x))[1] != dim(data)[1]){
    stop("Dimensions of x and data do not match in LSME")
  }
  d  <-     dim(data)[1]
  p  <-     dim(x)[2]
  s0 <-     LM0(x, h, type, data)
  s1.tmp <- LM1(x, h, type, data)
  s2.tmp <- LM2(x, h, type, data)
  ind <- which(s0 >= 1e-10)
  if (length(ind) == 0){
    return(list(theta = matrix(NA, d + d^2, p), s0 = s0, s1 = s1.tmp, s2 = s2.tmp))
  }
  s0 <- s0[ind]
  s1.tmp <- s1.tmp[ , ind, drop = FALSE]
  s2.tmp <- s2.tmp[ , ind, drop = FALSE]
  x  <- x[ , ind, drop = FALSE]
  # Normalized local moments
  s1 <- s1.tmp / rep(s0, each = d)
  s2 <- s2.tmp / rep(s0, each = d^2)
  # Normalized local derivatve moments
  q0 <- LM0_DK(x, h, type, data) / rep(s0, each = d)
  q1 <- LM1_DK(x, h, type, data) / rep(s0, each = d^2)
  S_tilde_as_vector <- s2 - apply(s1, 2, function(y) kronecker(y, y))
  S_tilde_inverse_as_vector <- apply(S_tilde_as_vector, 2,
                                     function(y) as.vector(qr.solve(matrix(y, d, d))))
  s1_times_S_tilde_inverse <- apply(as.matrix(s1[rep(1:d, d), ] * S_tilde_inverse_as_vector), 2,
                                    function(y) colSums(matrix(y, d, d)))
  s1_times_S_tilde_inverse_times_s1 <- colSums(s1_times_S_tilde_inverse * s1)
  theta_1 <- - rep(1 + s1_times_S_tilde_inverse_times_s1, each = d) * q0 +
             apply(as.matrix(s1_times_S_tilde_inverse[rep(1:d, each = d), ] * q1), 2,
                   function(y) rowSums(matrix(y, d, d))) +
             s1_times_S_tilde_inverse / h
  index1  <- rep(1:d, each = d)
  index2  <- rep(1:d, d)
  index3  <- rep(1:d^2, each = d)
  index4  <- d * rep(seq(0, d-1, 1), each = d^2) + rep(seq(1, d, 1), d^2)
  theta_2 <- s1_times_S_tilde_inverse [index1, ] * q0[index2, ] / h -
             apply(as.matrix(S_tilde_inverse_as_vector[index3, ] * q1[index4, ]), 2,
                   function(y) rowSums(matrix(y, d^2))) / h -
             S_tilde_inverse_as_vector / h^2
  # Add the 0-columns correspondig to s0 == 0
  theta <- matrix(NA, d + d^2, p)
  theta[ , ind] <- rbind(theta_1, theta_2)
  return(list(theta = theta, s0 = s0, s1 = s1.tmp, s2 = s2.tmp))
}


#' Local Estimator of the log-density with LLLE
#'
#' Returns the local estimator of \eqn{\hat{\boldsymbols{\theta}}} for the
#' location(s) \eqn{x} calculated with the LLLE
#'
#' @param x A vector or matrix containing the points at which the moment will be
#'   calculated.
#' @param h A number specifying the bandwidth.
#' @param type A character string specifying the used kernel function.
#' @param data A data matrix containing in each column a data point.
#' @return Returnes the the value \eqn{\hat{\boldsymbols{\theta}}} at \code{x}.
#'   Eihter as a vector or matrix.
LLLE <- function(x, h, type = "Gaussian", data){
  if(dim(as.matrix(x))[1] != dim(data)[1]){
    stop("Dimensions of x and data do not match in LLLE")
  }
  return("TO DO")
}

#' Local Estimator of the log-density with LMME
#'
#' Returns the local estimator of \eqn{\hat{\boldsymbols{\theta}}} for the
#' location(s) \eqn{x} calculated with the LMME
#'
#' @param x A vector or matrix containing the points at which the moment will be
#'   calculated.
#' @param h A number specifying the bandwidth.
#' @param type A character string specifying the used kernel function.
#' @param data A data matrix containing in each column a data point.
#' @return Returnes the the value \eqn{\hat{\boldsymbols{\theta}}} at \code{x}.
#'   Eihter as a vector or matrix.
LMME <- function(x, h, type = "Gaussian", data){
  if(dim(as.matrix(x))[1] != dim(data)[1]){
    stop("Dimensions of x and data do not match in LMME")
  }
  S_overline_2 <- rbind(LM0(x, h, type = type, data),
                        LM1(x, h, type = type, data),
                        LM2(x, h, type = type, data))
  d     <- dim(data)[1]
  mu_22 <- KernelMoment(22, type)
  mu_4  <- KernelMoment(4,  type)
  vecId <- get_vecId(d)
  B     <- rbind(t(vecId), matrix(0, d, d^2))
  M_33_inv <- qr.solve(kronecker(vecId,t(vecId)) * (mu_22 - 1) +
                       diag(rep(1, d^2)) * 2 * mu_22 +
                       diag(vecId) * (mu_4 - 3 * mu_22))
  M_inv <- cbind(rbind(diag(rep(1, d + 1)) + B %*% M_33_inv %*% t(B),
                       - M_33_inv %*% t(B)),
                 rbind(- B %*% M_33_inv,
                       M_33_inv))
  tmp <- rep(c(1, 1 / h, 2 / h^2), c(1, d, d^2)) * (M_inv %*% S_overline_2)
  return(theta = LogEstimation(tmp))
}

#' Local Estimator of the log-density with LMME2
#'
#' Returns the local estimator of \eqn{\hat{\boldsymbols{\theta}}} for the
#' location(s) \eqn{x} calculated with the LMME2. This method does not take into
#' account the second order Taylor expansion for calculating \eqn{\hat{f}(x)}.
#'
#' @param x A vector or matrix containing the points at which the moment will be
#'   calculated.
#' @param h A number specifying the bandwidth.
#' @param type A character string specifying the used kernel function.
#' @param data A data matrix containing in each column a data point.
#' @return Returnes the the value \eqn{\hat{\boldsymbols{\theta}}} at \code{x}.
#'   Eihter as a vector or matrix.
LMME2 <- function(x, h, type = "Gaussian", data){
  if(dim(as.matrix(x))[1] != dim(data)[1]){
    stop("Dimensions of x and data do not match in LMME")
  }
  S_overline_2 <- rbind(LM0(x, h, type = type, data),
                        LM1(x, h, type = type, data),
                        LM2(x, h, type = type, data))
  d     <- dim(data)[1]
  mu_22 <- KernelMoment(22, type)
  mu_4  <- KernelMoment(4,  type)
  vecId <- get_vecId(d)
  C     <- cbind(vecId, matrix(0, d^2, d))
  m_inv <- qr.solve(kronecker(vecId,t(vecId)) * mu_22 +
                       diag(rep(1, d^2)) * 2 * mu_22 +
                       diag(vecId) * (mu_4 - 3 * mu_22))
  M_inv <- cbind(rbind(diag(rep(1, d + 1)), - m_inv %*% C),
                 rbind(matrix(0, d + 1, d^2), m_inv))
  tmp <- rep(c(1, 1 / h, 2 / h^2), c(1, d, d^2)) * (M_inv %*% S_overline_2)
  return(theta = LogEstimation(tmp))
}


#' Local Estimator of the log-density with KDE
#'
#' Returns the local estimator of \eqn{\hat{\boldsymbols{\theta}}} for the
#' location(s) \eqn{x} calculated with the KDE
#'
#' @param x A vector or matrix containing the points at which the moment will be
#'   calculated.
#' @param h A number specifying the bandwidth.
#' @param type A character string specifying the used kernel function.
#' @param data A data matrix containing in each column a data point.
#' @return Returnes the the value \eqn{\hat{\boldsymbols{\theta}}} at \code{x}.
#'   Eihter as a vector or matrix.
KDE <- function(x, h, type = "Gaussian", data){
  if(dim(as.matrix(x))[1] != dim(data)[1]){
    stop("Dimensions of x and data do not match in KDE")
  }
  est <- rbind(LM0    (x, h, type = type, data = data),
               - LM0_DK (x, h, type = type, data = data),
               LM0_D2K(x, h, type = type, data = data))
  return(theta = LogEstimation(est))
}

#' Log-derivative transformation
#'
#' Transforms the vector of density derivatives up to second order into a vector
#' of log-density derivatives of 1st and 2nd order.
#'
#' @param est Estimator of the density derivatives in R^{1 + d + d^2}.
#' @return  Estimator of the log-density derivatives in R^{d + d^2}.
LogEstimation <- function(est){
  est <- as.matrix(est)
  p   <- dim(est)[1]
  d   <- - 0.5 + sqrt(-3 + 4*p) / 2
  if(abs(d - round(d)) > sqrt(.Machine$double.eps)){
    stop(paste("The length p of the vector does not fit any original dimension d. \n",
               "p =", p, ", calculated d = ",d))
  }
  xi1 <- est[1, ]
  xi2 <- est[2:(d + 1), ]
  xi3 <- est[-(1:(d + 1)), ]
  theta1 <- as.matrix(xi2 / rep(xi1, each = d))
  theta2 <- as.matrix(xi3 / rep(xi1, each = d^2) - apply(theta1, 2, function(y) kronecker(y, y)))
  return(theta = as.matrix(rbind(theta1, theta2)))
}

#' LHSE for Gaussian kernel
#'
#' Calculates the LHSE in case of the Gaussian kernel with a simpler formulea.
#'
LHSE.Gaussian <- function(x, h, data){
  if(dim(as.matrix(x))[1] != dim(data)[1]){
    stop("Dimensions of x and data do not match in LSME")
  }
  d  <-     dim(data)[1]
  s0 <-     LM0(x, h, type, data)
  s1.tmp <- LM1(x, h, type, data)
  s2.tmp <- LM2(x, h, type, data)
  # Normalized local moments
  s1 <- s1.tmp / rep(s0, each = d)
  s2 <- s2.tmp / rep(s0, each = d^2)
  S_tilde_as_vector <- s2 - apply(s1, 2, function(y) kronecker(y, y))
  S_tilde_inverse_as_vector <- apply(S_tilde_as_vector, 2,
                                     function(y) as.vector(qr.solve(matrix(y, d, d))))
  s1_times_S_tilde_inverse <- apply(as.matrix(s1[rep(1:d, d), ] * S_tilde_inverse_as_vector), 2,
                                    function(y) colSums(matrix(y, d, d)))
  theta_1 <- s1_times_S_tilde_inverse / h
  theta_2 <- (get_vecId(d) - S_tilde_inverse_as_vector) / h^2
    return(list(theta = rbind(theta_1, theta_2), s0 = s0, s1 = s1.tmp, s2 = s2.tmp))
}
