library(ldder)

context("Estimators")

# Simulate data
set.seed(10)
data <- SimCircle(100, 1, 0.1)
z    <- matrix(c(0,  0,
                -1,  0,
                 1,  0,
                 0, -1,
                 0,  1), nrow = 2)
h    <- 0.5
d    <- 2

test_that("general estimator", {
  expect_error(LogDensity(z, h, method = "NA",   data = data))
  expect_equal(LogDensity(z, h, method = "LSME", data = data),
               LSME(z, h , data = data))
  expect_equal(LogDensity(z, h, method = "LLLE", data = data),
               LLLE(z, h , data = data))
  expect_equal(LogDensity(z, h, method = "LMME", data = data),
               LMME(z, h , data = data))
  expect_equal(LogDensity(z, h, method = "KDE", data = data),
               KDE(z, h , data = data))
})

# Naive caluclation of the LSME
m <- dim(z)[2]
s0 <- LM0(z, h, data = data)
s1 <- LM1(z, h, data = data)
s2 <- LM2(z, h, data = data)
q0 <- LM0_DK(z, h, data = data)
q1 <- LM1_DK(z, h, data = data) + rep(s0, each = d^2) * matrix(get_vecId(d), d^2, m) / h
v  <- rbind(q0, q1)
theta <- matrix(NA, d + d^2, m)
for(j in 1:m){
  S_star_11 <- diag(rep(1,d)) * s0[j]
  S_star_12 <- kronecker(diag(rep(1, d)), t(s1[ , j]))
  S_star_22 <- kronecker(diag(rep(1, d)), matrix(s2[ , j], d, d))
  S_star    <- rbind(cbind(S_star_11, S_star_12),
                     cbind(t(S_star_12), S_star_22))
  S_star_inv  <- qr.solve(S_star)
  theta[ , j] <- - S_star_inv %*% v[ , j]
  theta[-(1:d), j] <- theta[-(1:d), j] / h
  }


test_that("LSME", {
  expect_equal(dim(LSME(z, h, data = data)), c(6, 5))
  expect_length(LSME(z[ , 1], h ,data = data), 6)
  expect_equal(as.vector(LSME(z, h, data = data)[ , 2]),
               as.vector(LSME(z[ , 2], h, data = data)))
  expect_equal(LSME(z, h, data = data), theta)
})

# Naive calculation of the LMME
vecId <- get_vecId(d)
mu22  <- KernelMoment(22)
mu4   <- KernelMoment(4)
m     <- kronecker(vecId, t(vecId)) * mu22 + diag(rep(2 * mu22, d^2)) +
         diag(vecId) * (mu4 - 3 * mu22)
B     <- rbind(vecId, matrix(0, d, d^2))
dimnames(B) <- NULL
M     <- rbind(cbind(diag(rep(1, d + 1)), B),
               cbind(t(B), m))
M_inv <- qr.solve(M)
S_overline_2 <- rbind(LM0(z, h, data = data),
                      LM1(z, h, data = data),
                      LM2(z, h, data = data))
est <- rep(c(1,1 / h, 2 / h^2), c(1, d, d^2)) * (M_inv %*% S_overline_2)
theta <- LogEstimation(est)

test_that("LMME", {
  expect_equal(dim(LMME(z, h, data = data)), c(6, 5))
  expect_length(LMME(z[ , 1], h ,data = data), 6)
  expect_equal(as.vector(LMME(z, h, data = data)[ , 2]),
               as.vector(LMME(z[ , 2], h, data = data)))
  expect_equal(LMME(z, h, data = data), theta)
})

test_that("KDE", {
  expect_equal(dim(KDE(z, h, data = data)), c(6, 5))
  expect_length(KDE(z[ , 1], h ,data = data), 6)
  expect_equal(as.vector(KDE(z, h, data = data)[ , 2]),
               as.vector(KDE(z[ , 2], h, data = data)))
})

# Compare the KDE and LSME with the transformation in case of
# Gaussian kernel.

theta_KDE <- KDE(z, h, data = data)
S0 <- LM0(z, h, data = data)
S1 <- LM1(z, h, data = data)
S2 <- LM2(z, h, data = data)
S_tilde_vec <- S2 / rep(S0, each = d^2) -
               apply(S1, 2, function(x) kronecker(x, x)) / rep(S0^2, each = d^2)
theta_LSME <- LSME(z, h, data = data)

LSME_trans <- matrix(NA, d + d^2, dim(z)[2])
for(j in 1:dim(z)[2]){
  LSME_trans[ ,j] <- kronecker(diag(rep(1, d + 1)),
                               matrix(S_tilde_vec[ , j], d, d)) %*% theta_LSME[ , j]
}
test_that("Compare KDE vs LSME",{
  expect_equal(LSME_trans, theta_KDE)
})
