# The curve-evaluation helpers live in validation/curve-study.R, which is
# excluded from the package build by .Rbuildignore. On CRAN (and in any
# installed copy of the package) that file does not exist, so this test is
# skipped there instead of failing. It still runs in your local source tree.
curve_file <- testthat::test_path("..", "..", "validation", "curve-study.R")
if (file.exists(curve_file)) {
  source(curve_file, local = TRUE)
}

test_that("curve evaluation preserves inputs and agrees with the probability formula", {
  skip_if_not(file.exists(curve_file),
              "validation/curve-study.R is not part of the built package")
  a <- array(0, c(3, 2, 4), dimnames = list(iteration = as.character(1:3),
    chain = as.character(1:2), variable = c("delta0", "delta1", "beta_Intercept", "beta_x")))
  a[, , "delta0"] <- .1
  a[, , "delta1"] <- .85
  a[, , "beta_Intercept"] <- -.5
  a[, , "beta_x"] <- 1.2
  original <- a
  grid <- seq(-4, 4, length.out = 31)
  curves <- curve_draws(posterior::as_draws_array(a), grid)
  expect_identical(a, original)
  expect_identical(dim(curves), c(3L, 2L, 31L))
  expected <- .1 + .75 * plogis(-.5 + 1.2 * grid)
  for (i in seq_along(grid)) {
    expect_equal(as.numeric(curves[, , i]), rep(expected[i], 6), tolerance = 1e-12)
  }
})
