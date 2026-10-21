.sample_asymptotes <- function(success0, failure0, success1, failure1, prior) {
  d0_spec <- prior$delta0
  d1_spec <- prior$delta1
  delta0 <- if (length(d0_spec) == 1L) {
    d0_spec
  } else {
    .guard_unit(stats::rbeta(1L, d0_spec[1L] + success0, d0_spec[2L] + failure0))
  }
  delta1 <- if (length(d1_spec) == 1L) {
    d1_spec
  } else {
    .guard_unit(stats::rbeta(1L, d1_spec[1L] + success1, d1_spec[2L] + failure1))
  }
  list(delta0 = delta0, delta1 = delta1)
}

# Beta draws with very small shape parameters can be exactly 0 or 1 in double
# precision, which would make the likelihood zero. Estimated asymptotes are
# therefore kept a tiny distance away from 0 and 1 (fixed values are not).
.guard_unit <- function(x, eps = 1e-10) min(max(x, eps), 1 - eps)

.prior_mean_asymptote <- function(spec) {
  if (length(spec) == 1L) spec else spec[1L] / sum(spec)
}

# Starting values. The asymptotes start from a prior draw, relabelled to
# follow the order of their prior means. Because the likelihood is invariant
# to (delta0, delta1, beta) -> (delta1, delta0, -beta), a chain started with
# the wrong sign of beta can stay in the mirror-image mode that the prior
# disfavours. beta therefore starts near a standard logistic fit to the data,
# with its sign matched to the prior order of the asymptotes, plus a small
# random perturbation so that chains start at different points.
.initial_state <- function(prior, X = NULL, successes = NULL, failures = NULL) {
  asymptotes <- .sample_asymptotes(
    success0 = 0, failure0 = 0,
    success1 = 0, failure1 = 0,
    prior = prior
  )
  # Start away from the boundary even when the prior puts mass near 0 or 1.
  delta0 <- if (length(prior$delta0) == 2L) min(max(asymptotes$delta0, 0.001), 0.999) else asymptotes$delta0
  delta1 <- if (length(prior$delta1) == 2L) min(max(asymptotes$delta1, 0.001), 0.999) else asymptotes$delta1
  m0 <- .prior_mean_asymptote(prior$delta0)
  m1 <- .prior_mean_asymptote(prior$delta1)
  if (length(prior$delta0) == 2L && length(prior$delta1) == 2L &&
      (m0 < m1) != (delta0 < delta1)) {
    tmp <- delta0; delta0 <- delta1; delta1 <- tmp
  }
  beta <- NULL
  if (!is.null(X) && m0 != m1) {
    start <- tryCatch(
      suppressWarnings(stats::glm.fit(X, cbind(successes, failures),
                                      family = stats::binomial())$coefficients),
      error = function(e) NULL
    )
    if (!is.null(start) && all(is.finite(start))) {
      start <- pmax(pmin(start, 10), -10) * sign(m1 - m0)
      beta <- start + stats::rnorm(length(start), sd = 0.1)
    }
  }
  if (is.null(beta)) beta <- .rmvn_covariance(prior$beta_mean, prior$beta_cov)
  list(delta0 = delta0, delta1 = delta1, beta = beta)
}

.latent_probabilities <- function(eta, delta0, delta1,
                                  successes = NULL, failures = NULL) {
  if (delta0 == 0 && delta1 == 0) {
    if (any(successes > 0)) {
      stop("Both endpoints are zero but successes are observed (zero likelihood).",
           call. = FALSE)
    }
    success <- stats::plogis(eta)
  } else {
    success <- stats::plogis(eta + log(delta1) - log(delta0))
  }
  if (delta0 == 1 && delta1 == 1) {
    if (any(failures > 0)) {
      stop("Both endpoints are one but failures are observed (zero likelihood).",
           call. = FALSE)
    }
    failure <- stats::plogis(eta)
  } else {
    failure <- stats::plogis(eta + log1p(-delta1) - log1p(-delta0))
  }
  list(success = success, failure = failure)
}


.run_bh_chain <- function(X, successes, failures, trials,
                          prior, control, seed) {
  set.seed(seed)
  n <- nrow(X)
  keep_at <- seq.int(control$warmup + 1L, control$iter,
                     by = control$thin)
  keep_lookup <- integer(control$iter)
  keep_lookup[keep_at] <- seq_along(keep_at)

  state <- .initial_state(prior, X, successes, failures)
  delta0 <- state$delta0
  delta1 <- state$delta1
  beta <- state$beta

  parameter_names <- .parameter_names(colnames(X))
  kept <- matrix(
    NA_real_,
    nrow = length(keep_at),
    ncol = length(parameter_names),
    dimnames = list(NULL, parameter_names)
  )

  prior_rhs <- drop(prior$beta_precision %*% prior$beta_mean)

  for (iteration in seq_len(control$iter)) {
    # Algorithm step 3: eta and allocation probabilities use the old state.
    eta <- drop(X %*% beta)
    allocation_probability <- .latent_probabilities(
      eta, delta0, delta1, successes, failures
    )

    z1_success <- stats::rbinom(
      n, size = successes, prob = allocation_probability$success
    )
    z1_failure <- stats::rbinom(
      n, size = failures, prob = allocation_probability$failure
    )
    z1_total <- z1_success + z1_failure

    # Step 4: independent Beta full conditionals given the new allocation.
    asymptotes <- .sample_asymptotes(
      success0 = sum(successes - z1_success),
      failure0 = sum(failures - z1_failure),
      success1 = sum(z1_success),
      failure1 = sum(z1_failure),
      prior = prior
    )
    delta0 <- asymptotes$delta0
    delta1 <- asymptotes$delta1

    # Step 5: PG(m_i, eta_old); m_i=1 for Bernoulli observations.
    omega <- BayesLogit::rpg(
      num = n,
      h = as.double(trials),
      z = as.double(eta)
    )
    # Step 6: Gaussian beta update using new allocation and omega.
    omega <- as.numeric(omega)
    if (length(omega) != n || any(!is.finite(omega)) || any(omega <= 0)) {
      stop(sprintf("Invalid Polya-Gamma draws at iteration %d; check BayesLogit.",
                   iteration), call. = FALSE)
    }
    precision <- crossprod(X, sweep(X, 1L, omega, `*`)) + prior$beta_precision
    precision_chol <- chol(precision)
    rhs <- drop(crossprod(X, z1_total - trials / 2)) + prior_rhs
    beta_mean <- backsolve(
      precision_chol,
      forwardsolve(t(precision_chol), rhs)
    )
    beta <- .rmvn_precision(beta_mean, precision_chol)

    # Step 8: keep by warmup/thin only, regardless of endpoint order.
    keep_row <- keep_lookup[iteration]
    if (keep_row > 0L) {
      kept[keep_row, ] <- c(delta0, delta1, beta)
    }
  }

  list(draws = kept)
}
