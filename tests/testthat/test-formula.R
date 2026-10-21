# Checks formula handling (utils.R): subset, factors, missing rows,
# unseen factor levels, unsupported offset(), invalid seed, parameter names.

test_that("subset, factors and missing training rows are handled", {
  dat <- data.frame(x = 1:12, g = factor(rep(c("a", "b"), 6)), y = rep(c(1, 0, 0), 4))
  dat$x[3] <- NA
  fit <- bh_fit(y ~ x + g, dat, bh_prior(0, 1),
                bh_control(20, 10, chains = 1), subset = x < 11,
                na.action = na.omit, seed = 2)
  expect_equal(nrow(fit$X), 9L)
  expect_equal(nrow(predict(fit, data.frame(x = c(1, 2), g = c("a", "b")))), 2L)
  expect_error(predict(fit, data.frame(x = 1, g = "c")))
  expect_error(bh_fit(y ~ x + offset(x), dat), "Offsets")
  expect_error(bh_fit(y ~ x, dat, seed = 1e20), "seed")
})

test_that("parameter names cannot collide with the special intercept name", {
  expect_equal(anyDuplicated(.parameter_names(c("(Intercept)", "Intercept"))), 0L)
})
