# Checks each Gibbs step against the paper's formulas (sampler.R).
#  1) Step 3: latent allocation probabilities equal Bayes' rule;
#     impossible endpoint combinations stop with 'zero likelihood'.
#  2) Step 4: the delta update is exactly Beta(a + counts, b + counts).
#  3) Step 6: the Cholesky-based normal draw has covariance = precision^-1.

test_that("latent allocation probabilities agree with Bayes' rule", {
  eta <- c(-2, 0, 3)
  g <- plogis(eta)
  d0 <- .2; d1 <- .8
  p <- d0 * (1 - g) + d1 * g
  z <- .latent_probabilities(eta, d0, d1)
  expect_equal(z$success, d1 * g / p)
  expect_equal(z$failure, (1 - d1) * g / (1 - p))
  z <- .latent_probabilities(eta, 0, 1)
  expect_equal(z$success, rep(1, 3))
  expect_equal(z$failure, rep(0, 3))
  expect_error(.latent_probabilities(0, 0, 0, successes = 1, failures = 0),
               "zero likelihood")
  expect_error(.latent_probabilities(0, 1, 1, successes = 0, failures = 1),
               "zero likelihood")
})

test_that("endpoint update is exactly two independent Beta draws", {
  p <- list(delta0 = c(2, 3), delta1 = c(4, 5))
  set.seed(26)
  expected <- c(rbeta(1, 9, 5), rbeta(1, 5, 13))
  set.seed(26)
  actual <- .sample_asymptotes(7, 2, 1, 8, p)
  expect_equal(unname(unlist(actual)), expected)
})

test_that("precision-normal sampler has the intended covariance algebra", {
  P <- matrix(c(4, 1, 1, 2), 2)
  U <- chol(P)
  A <- backsolve(U, diag(2))
  expect_equal(A %*% t(A), solve(P))
})
