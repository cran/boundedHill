# Checks input validation for priors and MCMC settings (prior.R, control.R).
#  1) bh_prior(): Beta shapes and fixed endpoints are accepted;
#     values outside [0, 1] or zero Beta shapes are rejected.
#  2) bh_control(): warmup >= iter, thin = 0, Inf or NA settings are rejected.

test_that("endpoint priors support fixed and random values", {
  p <- bh_prior(delta0 = c(1, 4), delta1 = c(4, 1))
  expect_s3_class(p, "bh_prior")
  expect_equal(p$delta0, c(1, 4))
  expect_equal(p$delta1, c(4, 1))
  for (pair in list(c(.1, .9), c(.5, .5), c(0, 1), c(0, 0), c(1, 1))) {
    expect_s3_class(bh_prior(pair[1], pair[2]), "bh_prior")
  }
  expect_error(bh_prior(-.1, 1))
  expect_error(bh_prior(0, 1.1))
  expect_error(bh_prior(c(0, 1), 1))
  expect_error(bh_prior(beta_cov = numeric()))
})

test_that("control validates settings and retains thinning convention", {
  expect_error(bh_control(iter = 100, warmup = 100))
  expect_error(bh_control(thin = 0))
  expect_error(bh_control(iter = Inf))
  expect_error(bh_control(iter = 1e20))
  expect_error(bh_control(chains = NA_real_))
  expect_s3_class(bh_control(iter = 100, warmup = 50, chains = 1), "bh_control")
})
