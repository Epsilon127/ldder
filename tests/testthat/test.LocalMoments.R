library(ldder)
context("Local Moments")

# Simulate data
data <- SimCircle(100, 1, 0.1)
z    <- matrix(c(0,  0,
                -1,  0,
                 1,  0,
                 0, -1,
                 0,  1), nrow = 2)

test_that("Dimension problem in local moments", {
  expect_length(    LM0    (z, 0.5, data = data),       5)
  expect_equal (dim(LM1    (z, 0.5, data = data)), c(2, 5))
  expect_equal (dim(LM2    (z, 0.5, data = data)), c(4, 5))
  expect_equal (dim(LM0_DK (z, 0.5, data = data)), c(2, 5))
  expect_equal (dim(LM0_D2K(z, 0.5, data = data)), c(4, 5))
  expect_equal (dim(LM1_DK (z, 0.5, data = data)), c(4, 5))
})
