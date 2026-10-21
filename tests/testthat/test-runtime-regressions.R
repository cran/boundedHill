# Guards against bugs that occurred earlier (regression tests).
#  1) Posterior draw matrices keep their matrix shape for 1 draw / 1 chain / 1 row.
#  2) BayesLogit::rpg() returns positive, finite Polya-Gamma draws.

test_that("draw matrices are ordinary matrices and rowwise probabilities agree", {
  for (ni in c(1L, 3L)) for (nc in c(1L, 2L)) for (nx in c(1L, 4L)) {
    a <- array(0, c(ni, nc, 4), dimnames = list(
      iteration = as.character(seq_len(ni)), chain = as.character(seq_len(nc)),
      variable = c("delta0", "delta1", "beta_Intercept", "beta_x")))
    a[, , "delta0"] <- .2; a[, , "delta1"] <- .8
    a[, , "beta_Intercept"] <- -.3; a[, , "beta_x"] <- 1.2
    x <- seq(-2, 2, length.out = nx)
    object <- structure(list(draws = posterior::as_draws_array(a),
      beta_names = c("beta_Intercept", "beta_x"), X = cbind(1, x)), class = "bounded_reg")
    m <- bh_posterior_draws(object, "matrix")
    expect_identical(class(m), c("matrix", "array"))
    expect_null(dim(m[, "delta0"]))
    actual <- .posterior_probability_matrix(object)
    expected <- matrix(rep(.2 + .6 * plogis(-.3 + 1.2 * x), each = ni * nc),
                       nrow = ni * nc)
    expect_identical(dim(actual), dim(expected))
    expect_equal(as.numeric(actual), as.numeric(expected))
  }
})

test_that("PG draws used for precision are positive with double arguments", {
  set.seed(816)
  counts <- c(1L, 1L, 10L, 10L)
  z <- c(-4, 0, 2, 5)
  omega <- BayesLogit::rpg(num = length(counts), h = as.double(counts), z = as.double(z))
  expect_length(omega, length(counts))
  expect_true(all(is.finite(omega) & omega > 0))
  X <- cbind(1, z)
  precision <- crossprod(X, sweep(X, 1L, as.numeric(omega), `*`)) + diag(.16, 2)
  expect_true(all(is.finite(chol(precision))))
})
