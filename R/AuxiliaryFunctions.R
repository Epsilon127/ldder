#' vec(I_d)
#'
#' Returning the vectorized identity matrix of dimension d.
#'
#' @param d Dimension of the vectorized identity matrix.
#' @return A vector of dimension d^2.
get_vecId <- function(d){
  v <- rep(0,d^2)
  v[(0:(d - 1)) * (d + 1) + 1] <- 1
  return(v)
}

#' Stochastic Integration
#'
#' Estimates the integral \eqn{\int K(x) g(x) dx} if K is a density function (or
#' kernel function). \code{espilon, alpha} are control parameters.
#'
#' @param g A function we want to integrate. Takes a vector as input and returns
#'   a vector.
#' @param d Dimension of the input argument of \code{g}.
#' @param type the used kernel
#' @param eps \code{eps} is the we allow our estimator to vary compared to
#'   the true value with certainty \code{1 - alpha}.
#' @param alpha \code{alpha} is the certainty with which our estimator is in a
#'   certained interval around the true value.
#' @param min The minimal number of simulations.
#' @param max The maximal number of simulations.
#' @return An vector estimating the integral.
Integrate <- function(g, d, type = "Gaussian", eps = 1e-4, alpha = 0.01, min = 100, max = Inf){
  m   <- min # Inital number of simulated points
  sim <- switch(type,
                Gaussian = matrix(rnorm(d * m), d, m),
                stop("This kernel is not available in the function Integrate"))
  g_sim <- apply(sim, 2, g)
  p     <- dim(g_sim)[1]
  mean  <- mean_old <- rowMeans(g_sim)
  var   <- var_old  <- matrix(rowSums(apply(g_sim, 2,
                 function(x) kronecker(x - mean, x - mean))) / (m - 1), p, p)
  eigen <- eigen(var, TRUE, TRUE)$values
  cdf <- hbe(eigen, m * eps^2)
  while(cdf < 1 - alpha & m < max){
    sim <- switch(type,
                  Gaussian = rnorm(d))
    sim <- matrix(sim, d, 1)
    g_sim <- g(sim)
    mean  <- (m / (m + 1)) * mean_old + g_sim / (m + 1)
    var   <- ((m - 1) / m) * var_old +
             (m + 1) * matrix(kronecker(mean - mean_old, mean - mean_old), p, p)
    mean_old <- mean
    eigen <- eigen(var, TRUE, TRUE)$values
    m <- m + 1
    cdf <- hbe(eigen, m * eps^2)
  }
  return(list(int = mean, var = var, m = m, eigen = eigen, cdf = cdf))
}
#' Grid in 2D
#'
#' Creat a matrix which corresponds to a grid of points on the plain
#'
#' @param x The interval of the x coordinates
#' @param y The interval of the y coordinates
#' @param lenght Horizontal and vertical distance between two points on the grid.
#'   If only one number is given it will be used for the horizontal and vertical
#'   distance. Ignored if \code{N} is specified
#' @param N Number of gridpoints.
#' @return Returns a 2 x N matrix which contains all the gridpoints.
#' @export
Grid2D <- function(x, y, length = NULL, N = NULL){
  if(!is.null(N)){
    ratio <- (x[2] - x[1]) / (y[2] - y[1])
    tmp   <- floor(sqrt(N / ratio))
    length <- c((x[2] - x[1]) / (tmp * ratio), (y[2] - y[1]) / tmp)
  } else{
    if(is.null(length)) stop("Neither length or N is specified.")
    if(length(length) == 1) length <- rep(length, 2)
  }
  x_grid <- seq(x[1], x[2], length[1])
  y_grid <- seq(y[1], y[2], length[2])
  x_length <- length(x_grid)
  y_length <- length(y_grid)
  res <- rbind(rep(x_grid, y_length),
               rep(y_grid, each = x_length))
  return(grid = res)
}

#' Test.KDE
#'
#' Testfunction comparing LHSE and KDE with Gauss data. The output is the
#' adjusted LHSE (with the bias of the KDE), which should now fit the KDE up to
#' order h^2.
#'
#' @param tehta.LHSE Dataframe with the estimators of the LHSE
#' @param h bandwith used for the estimators
#' @param var1 Variance on X1
#' @param var2 Variance on X2
Test.KDE <- function(theta.LHSE, h, var1, var2){
  theta <- theta.LHSE
  B1 <- theta$X1^2 / var1^2 - 1 / var2
  B2 <- theta$X1 * theta$X2 / sqrt(var1) / sqrt(var2)
  B3 <- theta$X2^2 / var2^2 - 1 / var1
  theta[, 3] <- theta[, 3] + h^2 * (B1 * theta.LHSE[, 3] + B2 * theta.LHSE[, 4])
  theta[, 4] <- theta[, 4] + h^2 * (B2 * theta.LHSE[, 3] + B3 * theta.LHSE[, 4])
  theta[, 5] <- theta[, 5] + h^2 * (B1 * theta.LHSE[, 5] + B2 * theta.LHSE[, 6])
  theta[, 6] <- theta[, 6] + h^2 * (B2 * theta.LHSE[, 5] + B3 * theta.LHSE[, 6])
  theta[, 7] <- theta[, 6]
  theta[, 8] <- theta[, 8] + h^2 * (B2 * theta.LHSE[, 7] + B3 * theta.LHSE[, 8])
  return(theta = theta)
}

#' S_tilde
#'
#' Returns the matrix S_tile as dataframe
#'
S_tilde <- function(x, h, type = "Gaussian", data){
  if(dim(as.matrix(x))[1] != dim(data)[1]){
    stop("Dimensions of x and data do not match in LSME")
  }
  d  <- dim(data)[1]
  s0 <- LM0(x, h, type, data)
  # Normalized local moments
  s1 <- LM1(x, h, type, data)    / rep(s0, each = d)
  s2 <- LM2(x, h, type, data)    / rep(s0, each = d^2)
  S_tilde_as_vector <- s2 - apply(s1, 2, function(y) kronecker(y, y))
  return(list(S_tilde = data.frame(X1 = Grid[1, ],
                              X2 = Grid[2, ],
                              S11 = S_tilde_as_vector[1, ],
                              S12 = S_tilde_as_vector[2, ],
                              S21 = S_tilde_as_vector[3, ],
                              S22 = S_tilde_as_vector[4, ]),
         h = h,
         type = type,
         data = data))
}

PlotS_tilde <- function(S_tilde){
  df <- gather(S_tilde, key = S_tilde, value = value, S11:S22)
  p <- ggplot(df, aes(X1, X2)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    facet_wrap(~S_tilde, nrow = 2) +
    geom_tile(aes(fill = value))
  return(p)
}

#' Bivariate Gaussian
#'
#' Calculates the true values of a bivarate Gaussian random variable at a grid
#' of given points
#'
#' @param x A vector or matrix with the values where the density derivatives will be calculated.
#' @param sigma1 The variance of X1
#' @param sigma2 The variance of X2
#' @return Returnes the  value \eqn{\boldsymbols{\theta}} at \code{x}.
#'   Eihter as a vector or matrix.
BiGaussian <- function(x, sigma1, sigma2){
  x <- as.matrix(x)
  m <- dim(x)[2]
  theta1 <- -x / rep(c(sigma1, sigma2), m)
  theta <- rbind(theta1, - 1 / sigma1, 0, 0, - 1 / sigma2)
  rownames(theta) <- paste0("theta", 1:6)
  return(theta = data.frame(t(x), t(theta)))
}

#' Ridge Radius Circle Data
#'
#' Calculates the radius of the true ridge underlying the circle data.
#'
#' The calculation uses a bisection algorithm.
#'
#' @param r Radius used in the generation of the circle data.
#' @param sd Standard deviation used in the generation of the circle data.
#' @return Returns the radius of the true ridge underlying the circle data
#'   generated with radius \code{r} and standard deviation \code{sd}.
#' @export
circle_radius <- function(r = 1, sd = 0.1){
  ratio <- r / sd
  threshold <- ratio^2
  alpha <- ratio / sd
  if (ratio <= sqrt(2)){
    return(0)
  }
  nu <- function(t){
    alpha * t * besselI(alpha * t, 0) / besselI(alpha * t, 1)
  }
  a <- 0
  b <- 1
  nu_a <- 2
  nu_b <- nu(b)
  while(nu_b < threshold){
    # Increase a and b till a < threshold < b.
    a <- b
    nu_a <- nu_b
    b <- 2 * b
    nu_b <- nu(b)
  }
  while (b - a > 0.0001){
    tmp <- (b + a) / 2
    nu_tmp <- nu(tmp)
    if (nu_tmp > threshold){
      b <- tmp
      nu_b <- nu_tmp
    } else{
      a <- tmp
      nu_a <- nu_tmp
    }
  }
  return((b + a) / 2)
}






