library(ldder)
context("Kernel functions")

# Simulate data
data <- SimCircle(100, 1, 0.1)
z    <- matrix(c(0,  0,
                -1,  0,
                 1,  0,
                 0, -1,
                 0,  1), nrow = 2)


test_that("KernelMoment is working",{
  expect_equal(KernelMoment(4, "Gaussian"), 3)
  expect_error(KernelMoment(0.5))
  expect_error(KernelMoment(22, "NA"))
})

test_that("K is working",{
  expect_length(K(z), 5)
  expect_length(K(z[ , 1]), 1)
  expect_error(K(z, "NA"))
})

test_that("DK is working",{
  expect_equal(dim(DK(z)), c(2, 5))
  expect_length(DK(z[ , 1]), 2)
  expect_error(DK(z, "NA"))
})

test_that("D2K is working",{
  expect_equal(dim(D2K(z)), c(4, 5))
  expect_length(D2K(z[ , 1]), 4)
  expect_error(D2K(z, "NA"))
})
