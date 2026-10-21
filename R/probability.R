#' Bounded logistic probabilities
#'
#' Evaluate \code{delta0 + (delta1 - delta0) * plogis(eta)}. No ordering of
#' the two endpoints is required.
#'
#' @param eta Finite numeric linear predictor.
#' @param delta0 Endpoint at eta tending to minus infinity, in [0, 1].
#' @param delta1 Endpoint at eta tending to plus infinity, in [0, 1].
#' @return A numeric vector of probabilities. Inputs must be scalar or have a
#'   common length.
#' @examples
#' bounded_probability(c(-2, 0, 2), delta0 = 0.1, delta1 = 0.8)
#' @export
bounded_probability <- function(eta, delta0, delta1) {
  args <- list(eta = eta, delta0 = delta0, delta1 = delta1)
  for (nm in names(args)) {
    z <- args[[nm]]
    if (!is.numeric(z) || !length(z) || any(!is.finite(z))) {
      stop(sprintf("`%s` must be nonempty, numeric, and finite.", nm), call. = FALSE)
    }
  }
  n <- max(lengths(args))
  if (any(!lengths(args) %in% c(1L, n))) {
    stop("Inputs must be scalar or have a common length.", call. = FALSE)
  }
  if (any(delta0 < 0 | delta0 > 1) || any(delta1 < 0 | delta1 > 1)) {
    stop("Both endpoints must lie in [0, 1].", call. = FALSE)
  }
  delta0 + (delta1 - delta0) * stats::plogis(eta)
}

#' Generate bounded binomial outcomes
#'
#' Simulate binomial counts whose success probability follows the bounded
#' logistic curve \code{\link{bounded_probability}}.
#'
#' @param eta Finite numeric linear predictor.
#' @param size Nonnegative trial counts, scalar or matching probability length.
#' @param delta0 Endpoint at eta tending to minus infinity.
#' @param delta1 Endpoint at eta tending to plus infinity.
#' @return A numeric vector of integer-valued counts.
#' @examples
#' set.seed(1)
#' rbounded_binom(c(-2, 0, 2), size = 10, delta0 = 0.1, delta1 = 0.8)
#' @export
rbounded_binom <- function(eta, size = 1L, delta0 = 0, delta1 = 1) {
  p <- bounded_probability(eta, delta0, delta1)
  .check_counts(size, "`size`")
  if (!length(size) %in% c(1L, length(p))) {
    stop("`size` must be scalar or match the probability length.", call. = FALSE)
  }
  stats::rbinom(length(p), size = rep_len(size, length(p)), prob = p)
}
