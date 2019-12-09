#' Test of ridge point
#'
#' \code{test.ridge} returns TRUE if with confidence 1-alpha there is a ridge
#' point at x.
#'
#' @export
#' @param est Vector with the estimator of the log-density derivatives.
#' @param x A vector or matrix containing the points at which the test will be
#'   executed.
#' @param h A number specifying the bandwidth.
#' @param n Sample size
#' @param type A character string specifying the used kernel function.
#' @param alpha Confidence level
#' @param S.tilde Hessian matrix of the sample score
#' @return Returns a logical wheter or not \code{x} is a ridge point or not.
test.ridge <- function(est, x, n, h, type, alpha = 0.05, S.tilde){
  if(length(x) != 2){
    stop("Dimension has to be 2!")
  }
  # Calculation of the critical value
  kappa.sq <- critical.value(alpha, type) / n / h^4
  # Calculate the eigenvalues and eigenvectors of \hat{A}
  b      <- est[1:2]
  A      <- matrix(est[3:6], 2)
  # calculate the gradient of the function at x
  # b      <- b + A %*% x
  tmp    <- eigen(A)
  lambda <- tmp$values[2:1]
  V2     <- tmp$vectors[ , 1]
  V1     <- c(V2[2], -V2[1])
  V      <- cbind(V1, V2)
  # Modify the gradient s.t. x is on the ridgeline
  b.tilde <- sum(b * V2) * V2
  # Modify A s.t. x is on the ridgeline
  # Angel between V2 and b
  gamma2 <- acos(sum(b * V2) / sqrt(sum(b^2)))
  # Angel between V1 and b
  gamma1 <- acos(sum(b * V1) / sqrt(sum(b^2)))
  # minimal Angel between V2 and +/- b
  if(0 < gamma1 & gamma1 <= pi / 2){
    if(0 < gamma2 & gamma2 <= pi / 2){
      gamma <- -gamma2
    } else{# pi / 2 < gamma2 <= pi
      gamma <- pi - gamma2
    }
  } else{# pi / 2 < gamma1 <= pi
    if(0 < gamma2 & gamma2 <= pi / 2){
      gamma <- gamma2
    } else{# pi /2 < gamma1 <= pi
      gamma <- -pi + gamma2
      }
  }
  if(-pi / 4 < gamma & gamma < pi / 4){
    lambda.tilde <- lambda * cos(gamma)^2 + lambda[2:1] * sin(gamma)^2
  } else{
    lambda.tilde <- rep(sum(lambda) / 2, 2)
  }
  lambda.tilde[1] <- min(0, lambda.tilde[1])
  V2.tilde <- c(cos(gamma) * V2[1] - sin(gamma) * V2[2],
                sin(gamma) * V2[1] + cos(gamma) * V2[2])
  V1.tilde <- c(cos(gamma) * V1[1] - sin(gamma) * V1[2],
                sin(gamma) * V1[1] + cos(gamma) * V1[2])
  V.tilde  <- cbind(V1.tilde, V2.tilde)
  A.tilde  <- V.tilde %*% diag(lambda.tilde) %*% t(V.tilde)
  # Calculating the distances
  eigen.S.tilde <- eigen(S.tilde)
  sqrt.S.tilde  <- eigen.S.tilde$vectors %*%
                   diag(sqrt(eigen.S.tilde$values)) %*% t(eigen.S.tilde$vectors)
  # Modifiy b and lambda[1] (if needed)
  A.mod  <- V %*% diag(c(min(0, lambda[1]), lambda[2])) %*% t(V)
  b.diff <- sqrt.S.tilde %*% c(b - b.tilde, A[1, 1] - A.mod[1, 1],
                               A[1, 2] - A.mod[1, 2], A[2, 2] - A.mod[2, 2])
  b.dist.sq <- sum(b.diff^2)
  # Modify A
  A.diff    <- A - V.tilde %*% diag(lambda.tilde) %*% t(V.tilde)
  A.dist.sq <- sum((sqrt.S.tilde[3:5, 3:5] %*% c(A.diff[1, 1], A.diff[1, 2], A.diff[2, 2]))^2)
  # check if the modified vector is contained in the confidence region
  res     <- !(min(b.dist.sq, A.dist.sq) < kappa.sq)
  weights <- c(sqrt(2) / 2 / pi, rep(1 / 2 / pi, 2), sqrt(2) / 4 / pi, sqrt(2) / 16 / pi)
  pv      <- 1 - pchisqsum(n * h^4 * min(b.dist.sq, A.dist.sq), rep(1, 5), weights,
                                     TRUE, "sat")
  return(list(res = res, pv = pv, b.dist.sq = b.dist.sq, A.dist.sq = A.dist.sq, kappa.sq = kappa.sq))
}

# Calculation of the critical value of the weighted least square distribution.
critical.value <- function(alpha = 0.05, type = "Gaussian", d = 2, tol = 10^(-10)){
  if(type != "Gaussian"){
    stop("The kernel has to be Gaussian!")
  }
  if(d != 2){
    stop("Dimension has to be 2!")
  }
  # Eigenvalues of Sigma^(1/2) Gamma^(-1) Sigma^(1/2)
  lambda <- c(sqrt(2) / 2 / pi, rep(1 / 2 / pi, 2), sqrt(2) / 4 / pi, sqrt(2) / 16 / pi)
  # using bisection to find the alpha-quantile
  return(qchisqsum(1 - alpha, rep(1, 5), lambda, tol))
}


# Quantil function of a weighted chi-square distribution
qchisqsum <- function(alpha, df, a, tol = 10^(-10)){
  kappa.left  <- 0
  kappa.right <- 1
  q.right <- pchisqsum(kappa.right, df, a, TRUE, "sat")
  q.left  <- 0
  while(q.right < alpha){
    kappa.left  <- kappa.right
    q.left      <- q.right
    kappa.right <- 2 * kappa.left
    q.right     <- pchisqsum(kappa.right, df, a, TRUE, "sat")
  }
  while(abs(q.right - alpha) > tol){
    kappa.tmp <- (kappa.left + kappa.right) / 2
    q.tmp     <- pchisqsum(kappa.tmp, df, a, TRUE, "sat")
    if(q.tmp > alpha){
      kappa.right <- kappa.tmp
      q.right     <- q.tmp
    } else{
      kappa.left <- kappa.tmp
      q.left     <- q.tmp
    }
  }
  return(kappa.right)
}



