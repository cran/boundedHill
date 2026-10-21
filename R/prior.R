#' Beta shapes from a prior mean and standard deviation
#'
#' Convert a prior mean and standard deviation for an asymptote into the two
#' Beta shape parameters used by \code{\link{bh_prior}}. With mean \eqn{m} and
#' standard deviation \eqn{s}, the shapes are \eqn{a = m k} and
#' \eqn{b = (1 - m) k} with \eqn{k = m(1 - m)/s^2 - 1}.
#'
#' @param mean Prior mean, strictly between 0 and 1.
#' @param sd Prior standard deviation; must satisfy
#'   \code{sd^2 < mean * (1 - mean)}.
#'
#' @return A named numeric vector \code{c(shape1, shape2)}.
#' @examples
#' bh_beta_moments(0.10, 0.12)   # baseline response around 10 percent
#' bh_beta_moments(0.85, 0.12)   # plateau around 85 percent
#' bh_prior(delta0 = bh_beta_moments(0.2, 0.1),
#'          delta1 = bh_beta_moments(0.7, 0.1))
#' @export
bh_beta_moments <- function(mean, sd) {
  if (!is.numeric(mean) || length(mean) != 1L || !is.finite(mean) ||
      mean <= 0 || mean >= 1) {
    stop("`mean` must be a single number strictly between 0 and 1.", call. = FALSE)
  }
  if (!is.numeric(sd) || length(sd) != 1L || !is.finite(sd) || sd <= 0) {
    stop("`sd` must be a single positive number.", call. = FALSE)
  }
  max_sd <- sqrt(mean * (1 - mean))
  if (sd >= max_sd) {
    stop(sprintf("`sd` must be smaller than sqrt(mean * (1 - mean)) = %.4g.", max_sd),
         call. = FALSE)
  }
  k <- mean * (1 - mean) / sd^2 - 1
  c(shape1 = mean * k, shape2 = (1 - mean) * k)
}

#' Specify priors for bounded logistic regression
#'
#' Construct independent Beta priors for the two endpoints and a multivariate
#' normal prior for regression coefficients. No ordering constraint is imposed.
#' A scalar fixes an endpoint; two positive values specify Beta shapes, which
#' can be obtained from a prior mean and standard deviation with
#' \code{\link{bh_beta_moments}}. The defaults put the baseline response
#' (\code{delta0}) at mean 0.10 and the plateau (\code{delta1}) at mean 0.85,
#' both with standard deviation 0.12, and give each regression coefficient an
#' independent normal prior with mean 0 and standard deviation 2.
#'
#' \code{delta0} and \code{delta1} are the limits as the linear predictor
#' tends to minus and plus infinity, respectively; they need not be the lower
#' and upper probability bounds. Swapping the endpoints and negating
#' \code{beta} leaves the likelihood unchanged, so symmetric priors produce a
#' label-symmetric posterior. Asymmetric priors fix the labels.
#'
#' @param delta0 Length-two Beta shape vector, or a scalar fixed eta-minus-infinity
#'   asymptote. Use \code{bh_beta_moments(mean, sd)} to give a mean and SD.
#' @param delta1 Length-two Beta shape vector, or a scalar fixed eta-plus-infinity
#'   asymptote. Use \code{bh_beta_moments(mean, sd)} to give a mean and SD.
#' @param beta_mean Prior mean for regression coefficients. A scalar is
#'   recycled after the model matrix is known.
#' @param beta_sd Prior standard deviation for independent regression
#'   coefficients. Ignored when \code{beta_cov} is supplied.
#' @param beta_cov Optional prior covariance matrix. A scalar creates a
#'   diagonal covariance matrix with that variance.
#'
#' @return An object of class \code{bh_prior}: a list with elements
#'   \code{delta0}, \code{delta1}, \code{beta_mean}, \code{beta_sd} and
#'   \code{beta_cov}. The \code{print} method returns its input invisibly.
#' @examples
#' bh_prior()
#' # Asymptotes from mean and SD, wider N(0, 4^2) priors for the logistic part
#' bh_prior(delta0 = bh_beta_moments(0.15, 0.10),
#'          delta1 = bh_beta_moments(0.75, 0.10), beta_sd = 4)
#' # Standard logistic regression: both asymptotes fixed
#' bh_prior(delta0 = 0, delta1 = 1)
#' @export
bh_prior <- function(delta0 = bh_beta_moments(0.10, 0.12),
                     delta1 = bh_beta_moments(0.85, 0.12),
                     beta_mean = 0,
                     beta_sd = 2,
                     beta_cov = NULL) {
  delta0 <- .validate_asymptote_prior(delta0, "delta0")
  delta1 <- .validate_asymptote_prior(delta1, "delta1")

  if (!is.numeric(beta_mean) || !length(beta_mean) ||
      any(!is.finite(beta_mean))) {
    stop("`beta_mean` must be a finite numeric vector.", call. = FALSE)
  }
  if (!is.numeric(beta_sd) || !length(beta_sd) ||
      any(!is.finite(beta_sd)) || any(beta_sd <= 0)) {
    stop("`beta_sd` must contain positive finite values.", call. = FALSE)
  }
  if (!is.null(beta_cov)) {
    if (!is.numeric(beta_cov) || !length(beta_cov) || any(!is.finite(beta_cov))) {
      stop("`beta_cov` must be numeric and finite.", call. = FALSE)
    }
    if (length(beta_cov) == 1L && beta_cov <= 0) {
      stop("A scalar `beta_cov` must be a positive variance.",
           call. = FALSE)
    }
    if (length(beta_cov) > 1L && !is.matrix(beta_cov)) {
      stop("`beta_cov` must be a covariance matrix or a positive scalar.",
           call. = FALSE)
    }
  }

  structure(
    list(
      delta0 = delta0,
      delta1 = delta1,
      beta_mean = as.numeric(beta_mean),
      beta_sd = as.numeric(beta_sd),
      beta_cov = beta_cov
    ),
    class = "bh_prior"
  )
}

#' @rdname bh_prior
#' @param x A \code{bh_prior} object.
#' @param ... Further arguments; currently unused.
#' @export
print.bh_prior <- function(x, ...) {
  describe_delta <- function(z) {
    if (length(z) == 1L) {
      sprintf("fixed at %.4g", z)
    } else {
      m <- z[1] / sum(z)
      s <- sqrt(z[1] * z[2] / (sum(z)^2 * (sum(z) + 1)))
      sprintf("Beta(%.4g, %.4g)  [mean %.3g, SD %.3g]", z[1], z[2], m, s)
    }
  }

  cat("Prior for Bayesian bounded logistic regression\n")
  cat("  delta0:", describe_delta(x$delta0), "\n")
  cat("  delta1:", describe_delta(x$delta1), "\n")
  cat("  beta mean:", paste(x$beta_mean, collapse = ", "), "\n")
  if (is.null(x$beta_cov)) {
    cat("  beta SD:", paste(x$beta_sd, collapse = ", "), "\n")
  } else {
    cat("  beta covariance: supplied\n")
  }
  invisible(x)
}

.validate_asymptote_prior <- function(x, name) {
  if (!is.numeric(x) || !(length(x) %in% c(1L, 2L)) ||
      any(!is.finite(x))) {
    stop(sprintf("`%s` must be a fixed scalar or two Beta shapes.", name),
         call. = FALSE)
  }
  if (length(x) == 1L) {
    x <- unname(x)
    valid <- x >= 0 && x <= 1
    if (!valid) {
      interval <- "[0, 1]"
      stop(sprintf("Fixed `%s` must lie in %s.", name, interval),
           call. = FALSE)
    }
  } else if (any(x <= 0)) {
    stop(sprintf("The Beta shapes in `%s` must be positive.", name),
         call. = FALSE)
  }
  as.numeric(unname(x))
}

.resolve_prior <- function(prior, p, coef_names) {
  if (!inherits(prior, "bh_prior")) {
    stop("`prior` must be created by `bh_prior()`.", call. = FALSE)
  }

  beta_mean <- prior$beta_mean
  if (length(beta_mean) == 1L) {
    beta_mean <- rep(beta_mean, p)
  }
  if (length(beta_mean) != p) {
    stop(sprintf(
      "`beta_mean` has length %d, but the model matrix has %d columns.",
      length(beta_mean), p
    ), call. = FALSE)
  }

  if (is.null(prior$beta_cov)) {
    beta_sd <- prior$beta_sd
    if (length(beta_sd) == 1L) {
      beta_sd <- rep(beta_sd, p)
    }
    if (length(beta_sd) != p) {
      stop(sprintf(
        "`beta_sd` has length %d, but the model matrix has %d columns.",
        length(beta_sd), p
      ), call. = FALSE)
    }
    beta_cov <- diag(beta_sd^2, p)
  } else if (length(prior$beta_cov) == 1L) {
    beta_cov <- diag(as.numeric(prior$beta_cov), p)
  } else {
    beta_cov <- prior$beta_cov
    if (!all(dim(beta_cov) == c(p, p))) {
      stop(sprintf("`beta_cov` must be a %d by %d matrix.", p, p),
           call. = FALSE)
    }
    if (!isTRUE(all.equal(beta_cov, t(beta_cov), tolerance = 1e-10))) {
      stop("`beta_cov` must be symmetric.", call. = FALSE)
    }
  }

  chol_cov <- tryCatch(chol(beta_cov), error = function(e) NULL)
  if (is.null(chol_cov)) {
    stop("`beta_cov` must be positive definite.", call. = FALSE)
  }

  names(beta_mean) <- coef_names
  dimnames(beta_cov) <- list(coef_names, coef_names)

  list(
    delta0 = prior$delta0,
    delta1 = prior$delta1,
    beta_mean = beta_mean,
    beta_cov = beta_cov,
    beta_precision = chol2inv(chol_cov)
  )
}
