# Checks the curve formula itself (probability.R).
#  1) bounded_probability(): eta -> -Inf gives delta0, eta -> +Inf gives delta1,
#     delta0 = delta1 gives a flat line, (0, 1) gives plain plogis();
#     invalid endpoints or lengths must stop with an error.
#  2) rbounded_binom(): simulated counts are integers between 0 and size.

test_that("probabilities allow increasing and constant endpoints", {
  eta <- c(-100, 0, 100)
  expect_equal(bounded_probability(eta, .2, .8), c(.2, .5, .8), tolerance = 1e-12)
  expect_equal(bounded_probability(eta, .3, .3), rep(.3, 3))
  expect_equal(bounded_probability(eta, 0, 1), plogis(eta))
  expect_error(bounded_probability(0, -.1, .8))
  expect_error(bounded_probability(0, .1, 1.1))
  expect_error(bounded_probability(numeric(), .1, .9))
  expect_error(bounded_probability(1:3, c(.1, .2), .9))
})

test_that("binomial simulation respects trial counts and input lengths", {
  set.seed(1)
  y <- rbounded_binom(rep(0, 20), size = 3, delta0 = .2, delta1 = .8)
  expect_true(all(y >= 0 & y <= 3 & y == floor(y)))
  expect_error(rbounded_binom(1:3, size = c(1, 2)))
  expect_error(rbounded_binom(0, size = -1))
})
