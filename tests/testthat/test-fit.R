# End-to-end checks of bh_fit() and the methods (fit.R, methods.R).
#  1) A short fit returns finite draws, deltas in [0, 1], predictions in [0, 1].
#  2) Fixed deltas stay fixed in every draw; the same seed gives the same draws.
#  3) Grouped data: predictive draws, log-likelihood and error messages are correct.

test_that("fit keeps all unordered posterior draws and predicts probabilities", {
  set.seed(11)
  dat <- data.frame(x = seq(-2, 2, length.out = 40))
  dat$y <- rbounded_binom(-.2 + dat$x, delta0 = .1, delta1 = .85)
  fit <- bh_fit(y ~ x, dat, prior = bh_prior(c(2, 8), c(8, 2)),
                control = bh_control(iter = 100, warmup = 50, chains = 2), seed = 9)
  expect_s3_class(fit, "bounded_reg")
  draws <- bh_posterior_draws(fit, "matrix")
  expect_equal(nrow(draws), 100L)
  expect_true(all(is.finite(draws)))
  expect_true(all(draws[, c("delta0", "delta1")] >= 0 &
                  draws[, c("delta0", "delta1")] <= 1))
  pred <- predict(fit)
  expect_true(all(pred$estimate >= 0 & pred$estimate <= 1))
  expect_s3_class(summary(fit), "summary.bounded_reg")
})

test_that("fixed endpoints are preserved in every retained draw", {
  dat <- data.frame(x = seq(-2, 2, length.out = 20), y = rep(c(0, 1), 10))
  ctl <- bh_control(iter = 30, warmup = 10, thin = 3, chains = 2)
  a <- bh_fit(y ~ x, dat, bh_prior(.1, .9), ctl, seed = 7)
  b <- bh_fit(y ~ x, dat, bh_prior(.1, .9), ctl, seed = 7)
  d <- bh_posterior_draws(a, "matrix")
  expect_equal(nrow(d), 14L)
  expect_equal(unname(d[, "delta0"]), rep(.1, 14))
  expect_equal(unname(d[, "delta1"]), rep(.9, 14))
  expect_equal(a$draws, b$draws)
  expect_warning(diag <- bh_diagnose(a), NA)
  expect_setequal(diag$fixed_parameters, c("delta0", "delta1"))
  expect_false(any(c("delta0", "delta1") %in% diag$invalid_flag))
})

test_that("grouped data and predictive/log-likelihood methods have correct rows", {
  dat <- data.frame(x = seq(-2, 2, length.out = 5),
                    success = c(1, 3, 6, 8, 9), failure = c(9, 7, 4, 2, 1))
  fit <- bh_fit(cbind(success, failure) ~ x, dat, bh_prior(0, 1),
                bh_control(iter = 40, warmup = 20, chains = 1), seed = 3)
  expect_equal(sum(fit$trials), 50)
  expect_equal(nrow(predict(fit)), 5L)
  pp <- bh_posterior_predict(fit, ndraws = 10, seed = 4)
  expect_equal(dim(pp), c(10L, 5L))
  expect_true(all(pp >= 0 & pp <= 10))
  prob <- predict(fit, type = "draws")
  ll <- bh_log_lik(fit)
  expect_equal(unname(ll[, 1]), unname(dbinom(1, 10, prob[, 1], log = TRUE)))
  expect_equal(bh_log_lik(fit, dat), ll)
  bad <- dat; bad$success[2] <- NA
  expect_error(bh_log_lik(fit, bad), "Missing responses")
  bad <- dat; bad$x[2] <- NA
  expect_error(predict(fit, bad), "Missing")
  expect_error(bh_posterior_predict(fit, trials = c(2, 3)))
  expect_error(bh_posterior_predict(fit, ndraws = NA_real_))
})
