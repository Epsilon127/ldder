library(ldder)
context("Epanechnikov Kernel")

# Data
z <- matrix(c(0,0,
              1,1,
              -1,-1,
              -sqrt(5),sqrt(5),
              4,4),nrow = 2)
const <- 3 * sqrt(5) / 20

test_that("K.Epan",{
  expect_equal(K.Epan(z), K(z, type = "Epan"))
  expect_length(K.Epan(z), 5)
  expect_equal(K.Epan(z), c(3 * sqrt(5) / 20,
                           rep(3 * sqrt(5) / 20 * 4 / 5, 2),
                           0,0)^2)
  })

test_that("DK.Epan",{
  expect_equal(DK.Epan(z), DK(z, type = "Epan"))
  expect_equal(dim(DK.Epan(z)), c(2, 5))
  expect_equal(DK.Epan(z), matrix(c(0, 0,
                                   c(-1, -1, 1, 1) * (3 * sqrt(5) / 20)^2 * 4 / 5 * 2 / 5,
                                   rep(0, 4)), nrow = 2))
  expect_length(DK(z[ , 1], type = "Epan"), 2)
})

test_that("D2K.Epan",{
  expect_equal(D2K.Epan(z), D2K(z, type = "Epan"))
  expect_equal(dim(D2K.Epan(z)), c(4, 5))
  expect_equal(D2K.Epan(z), matrix(const^2 *
                                  c(-2 / 5, 0, 0, -2 / 5,
                                    rep(c(-2 / 5 * 4 / 5, (2 / 5)^2, (2 / 5)^2, -2 / 5 * 4 / 5), 2),
                                    rep(0, 8)), nrow = 4))
  expect_length(D2K(z[ , 1], type = "Epan"), 4)
})
