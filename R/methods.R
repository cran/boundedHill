# ---------------------------------------------------------------------------
# Help page 1: "bounded_reg-methods"
# print, summary, print.summary, coef, predict, fitted, simulate share ONE
# help page. Every argument is documented exactly once, in this first block.
# The other functions only say `@rdname bounded_reg-methods`.
# (plot has its own page below because its `type` means something different.)
# ---------------------------------------------------------------------------

#' Methods for bounded logistic regression fits
#'
#' Standard S3 methods for objects of class \code{bounded_reg} returned by
#' \code{\link{bh_fit}}: printing, posterior summaries, point estimates,
#' predictions, fitted values, and posterior predictive simulation.
#'
#' @param x A \code{bounded_reg} fit (for \code{print}) or a
#'   \code{summary.bounded_reg} object (for \code{print.summary}).
#' @param object A \code{bounded_reg} fit.
#' @param digits Number of displayed digits.
#' @param statistic Either the posterior \code{"mean"} or \code{"median"}.
#' @param newdata Optional data frame of predictor values. If \code{NULL},
#'   the data used for fitting are used.
#' @param type For \code{predict}: \code{"response"} returns posterior
#'   summaries of the bounded probabilities, \code{"linear"} summarizes the
#'   linear predictor, and \code{"draws"} returns the probability draw matrix.
#' @param level Credible interval level, strictly between 0 and 1.
#' @param nsim Number of posterior predictive simulations.
#' @param seed Optional random seed.
#' @param trials Binomial trial counts for the rows of \code{newdata};
#'   defaults to one. For the original data the observed trial counts are used.
#' @param ... Further arguments; currently unused.
#'
#' @return
#' \itemize{
#'   \item \code{print}: the input object, invisibly.
#'   \item \code{summary}: an object of class \code{summary.bounded_reg}
#'     containing a table of posterior summaries and MCMC diagnostics.
#'   \item \code{coef}: a named numeric vector containing both asymptotes and
#'     the regression coefficients.
#'   \item \code{predict}: a data frame with columns \code{estimate},
#'     \code{median}, \code{lower} and \code{upper} for
#'     \code{type = "response"} or \code{"linear"}; a matrix of posterior
#'     probability draws (draws in rows, observations in columns) for
#'     \code{type = "draws"}.
#'   \item \code{fitted}: a numeric vector of posterior mean probabilities.
#'   \item \code{simulate}: a data frame with one column per simulation
#'     (\code{sim_1}, \code{sim_2}, ...) and one row per observation.
#' }
#'
#' @examples
#' set.seed(1)
#' x <- seq(-3, 3, length.out = 60)
#' y <- rbounded_binom(-0.5 + 1.2 * x, delta0 = 0.1, delta1 = 0.8)
#' fit <- bh_fit(y ~ x, data = data.frame(x = x, y = y),
#'               control = bh_control(iter = 400, warmup = 200, chains = 2),
#'               seed = 1)
#' fit
#' summary(fit)
#' coef(fit)
#' predict(fit, newdata = data.frame(x = c(-2, 0, 2)))
#' head(fitted(fit))
#' sims <- simulate(fit, nsim = 2, seed = 3)
#' @aliases bounded_reg-methods
#' @rdname bounded_reg-methods
#' @export
print.bounded_reg <- function(x, ...) {
  kept_per_chain <- dim(x$draws)[1L]
  cat("Bayesian bounded logistic regression\n")
  cat("Formula: ")
  print(x$formula)
  cat("Observations:", nrow(x$X), "\n")
  cat("Total binomial trials:", sum(x$trials), "\n")
  cat("Regression coefficients:", ncol(x$X), "\n")
  cat(
    "Posterior draws:",
    kept_per_chain * x$control$chains,
    sprintf("(%d chains)\n", x$control$chains)
  )
  invisible(x)
}

#' @rdname bounded_reg-methods
#' @export
summary.bounded_reg <- function(object, ...) {
  table <- posterior::summarise_draws(object$draws)
  structure(
    list(
      call = object$call,
      formula = object$formula,
      n = nrow(object$X),
      trials = sum(object$trials),
      table = table
    ),
    class = "summary.bounded_reg"
  )
}

#' @rdname bounded_reg-methods
#' @export
print.summary.bounded_reg <- function(x, digits = 3L, ...) {
  cat("Bayesian bounded logistic regression\n")
  cat("Formula: ")
  print(x$formula)
  cat("Observations:", x$n, "  Trials:", x$trials, "\n\n")

  wanted <- intersect(
    c("variable", "mean", "sd", "q5", "median", "q95",
      "rhat", "ess_bulk", "ess_tail"),
    names(x$table)
  )
  print.data.frame(
    as.data.frame(x$table[, wanted, drop = FALSE]),
    digits = digits,
    row.names = FALSE
  )
  invisible(x)
}

#' @rdname bounded_reg-methods
#' @export
coef.bounded_reg <- function(object,
                             statistic = c("mean", "median"),
                             ...) {
  statistic <- match.arg(statistic)
  draws <- .draw_matrix(object)
  if (statistic == "mean") {
    colMeans(draws)
  } else {
    apply(draws, 2L, stats::median)
  }
}

#' @rdname bounded_reg-methods
#' @export
predict.bounded_reg <- function(object,
                                newdata = NULL,
                                type = c("response", "linear", "draws"),
                                level = 0.95,
                                ...) {
  type <- match.arg(type)
  if (type == "linear") {
    eta <- .posterior_linear_matrix(object, newdata)
    colnames(eta) <- paste0("row", seq_len(ncol(eta)))
    return(.summarize_columns(eta, level))
  }

  probability <- .posterior_probability_matrix(object, newdata)
  colnames(probability) <- paste0("row", seq_len(ncol(probability)))
  if (type == "draws") {
    return(probability)
  }
  .summarize_columns(probability, level)
}

#' @rdname bounded_reg-methods
#' @export
fitted.bounded_reg <- function(object, ...) {
  predict(object, newdata = NULL, type = "response", ...)$estimate
}

#' @rdname bounded_reg-methods
#' @export
simulate.bounded_reg <- function(object,
                                 nsim = 1L,
                                 seed = NULL,
                                 newdata = NULL,
                                 trials = NULL,
                                 ...) {
  draws <- bh_posterior_predict(
    object,
    newdata = newdata,
    trials = trials,
    ndraws = nsim,
    seed = seed
  )
  result <- as.data.frame(t(draws))
  names(result) <- paste0("sim_", seq_len(ncol(result)))
  result
}

# ---------------------------------------------------------------------------
# Help page 2: plot (separate page; its `type` differs from predict's `type`)
# ---------------------------------------------------------------------------

#' Plot a bounded logistic regression fit
#'
#' Draw MCMC trace plots, marginal posterior densities, or fitted
#' probabilities with credible bands against the observed proportions.
#'
#' @param x A \code{bounded_reg} fit.
#' @param type Plot type: \code{"fitted"} (fitted probabilities and observed
#'   proportions), \code{"trace"} (MCMC traces by chain), or
#'   \code{"density"} (marginal posterior densities).
#' @param variables Optional character vector of parameter names for trace
#'   or density plots. Defaults to all parameters.
#' @param ... Additional graphical arguments passed to
#'   \code{\link[graphics]{plot}}.
#'
#' @return The input object, invisibly. Called for its side effect of
#'   drawing a plot.
#' @examples
#' set.seed(1)
#' x <- seq(-3, 3, length.out = 60)
#' y <- rbounded_binom(-0.5 + 1.2 * x, delta0 = 0.1, delta1 = 0.8)
#' fit <- bh_fit(y ~ x, data = data.frame(x = x, y = y),
#'               control = bh_control(iter = 400, warmup = 200, chains = 2),
#'               seed = 1)
#' plot(fit)
#' plot(fit, type = "trace", variables = c("delta0", "delta1"))
#' @export
plot.bounded_reg <- function(x,
                             type = c("fitted", "trace", "density"),
                             variables = NULL,
                             ...) {
  type <- match.arg(type)
  if (type == "fitted") {
    fit <- predict(x, type = "response")
    observed <- x$successes / x$trials
    index <- seq_along(observed)
    ylim <- range(c(0, 1, fit$lower, fit$upper), finite = TRUE)
    graphics::plot(
      index, fit$estimate,
      ylim = ylim,
      xlab = "Observation",
      ylab = "Response probability",
      type = "n",
      ...
    )
    graphics::polygon(
      c(index, rev(index)),
      c(fit$lower, rev(fit$upper)),
      col = grDevices::adjustcolor("steelblue", alpha.f = 0.25),
      border = NA
    )
    graphics::lines(index, fit$estimate, col = "steelblue", lwd = 2)
    graphics::points(index, observed, pch = 16, col = "black")
    return(invisible(x))
  }

  draws_df <- posterior::as_draws_df(x$draws)
  available <- setdiff(names(draws_df), posterior::reserved_variables())
  if (is.null(variables)) {
    variables <- available
  }
  if (!length(variables) || !all(variables %in% available)) {
    stop("Unknown parameter in `variables`.", call. = FALSE)
  }

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mfrow = grDevices::n2mfrow(length(variables)))

  if (type == "trace") {
    chain <- draws_df$.chain
    colors <- seq_len(max(chain))
    for (variable in variables) {
      values <- draws_df[[variable]]
      graphics::plot(
        draws_df$.iteration, values,
        type = "n",
        xlab = "Stored draw",
        ylab = variable,
        main = variable,
        ...
      )
      for (chain_id in unique(chain)) {
        use <- chain == chain_id
        graphics::lines(draws_df$.iteration[use], values[use], col = colors[chain_id])
      }
    }
  } else {
    for (variable in variables) {
      values <- draws_df[[variable]]
      if (length(unique(values)) == 1L) {
        graphics::plot(values[1L], 1, type = "h", xlim = values[1L] + c(-0.1, 0.1),
                       main = paste(variable, "(constant)"), xlab = variable,
                       ylab = "Point mass", ...)
      } else {
        d <- stats::density(values)
        graphics::plot(d, main = variable, xlab = variable, ...)
      }
    }
  }
  invisible(x)
}

# ---------------------------------------------------------------------------
# Separate help pages for the bh_* helper functions
# ---------------------------------------------------------------------------

#' Extract posterior draws
#'
#' @param object A \code{bounded_reg} fit.
#' @param format One of \code{"draws_array"} (a \pkg{posterior}
#'   \code{draws_array}), \code{"matrix"}, or \code{"data.frame"}.
#'
#' @return Posterior draws in the requested format, with the asymptotes
#'   \code{delta0}, \code{delta1} and the regression coefficients
#'   \code{beta_*} as variables.
#' @examples
#' set.seed(1)
#' x <- seq(-3, 3, length.out = 60)
#' y <- rbounded_binom(-0.5 + 1.2 * x, delta0 = 0.1, delta1 = 0.8)
#' fit <- bh_fit(y ~ x, data = data.frame(x = x, y = y),
#'               control = bh_control(iter = 400, warmup = 200, chains = 2),
#'               seed = 1)
#' draws <- bh_posterior_draws(fit, format = "matrix")
#' dim(draws)
#' @export
bh_posterior_draws <- function(object,
                               format = c("draws_array", "matrix",
                                          "data.frame")) {
  if (!inherits(object, "bounded_reg")) {
    stop("`object` must be a `bounded_reg` fit.", call. = FALSE)
  }
  format <- match.arg(format)
  switch(
    format,
    draws_array = object$draws,
    matrix = .draw_matrix(object),
    data.frame = as.data.frame(posterior::as_draws_df(object$draws))
  )
}

#' Pointwise log-likelihood draws
#'
#' Compute the binomial log likelihood of each observation under each
#' posterior draw, for example for use with the \pkg{loo} package.
#'
#' @param object A \code{bounded_reg} fit.
#' @param newdata Optional new data containing both the response and the
#'   predictors. If \code{NULL}, the data used for fitting are used.
#'
#' @return A numeric matrix with posterior draws in rows and observations in
#'   columns.
#' @examples
#' set.seed(1)
#' x <- seq(-3, 3, length.out = 60)
#' y <- rbounded_binom(-0.5 + 1.2 * x, delta0 = 0.1, delta1 = 0.8)
#' fit <- bh_fit(y ~ x, data = data.frame(x = x, y = y),
#'               control = bh_control(iter = 400, warmup = 200, chains = 2),
#'               seed = 1)
#' ll <- bh_log_lik(fit)
#' dim(ll)
#' @export
bh_log_lik <- function(object, newdata = NULL) {
  if (!inherits(object, "bounded_reg")) {
    stop("`object` must be a `bounded_reg` fit.", call. = FALSE)
  }
  if (is.null(newdata)) {
    probability <- .posterior_probability_matrix(object)
    successes <- object$successes
    trials <- object$trials
  } else {
    mf <- stats::model.frame(object$terms, newdata, na.action = stats::na.pass,
                             xlev = object$xlevels)
    prepared <- .prepare_model_frame(mf, object$contrasts)
    probability <- .posterior_probability_matrix(object, newdata)
    successes <- prepared$successes
    trials <- prepared$trials
  }

  out <- matrix(NA_real_, nrow(probability), ncol(probability))
  for (j in seq_len(ncol(probability))) {
    out[, j] <- stats::dbinom(
      successes[j], size = trials[j], prob = probability[, j], log = TRUE
    )
  }
  colnames(out) <- paste0("row", seq_len(ncol(out)))
  out
}

#' Draw from the posterior predictive distribution
#'
#' @param object A \code{bounded_reg} fit.
#' @param newdata Optional data frame of predictor values. If \code{NULL},
#'   the data used for fitting are used.
#' @param trials Trial counts for \code{newdata}; defaults to one. For the
#'   original data, the observed trial counts are used.
#' @param ndraws Number of posterior predictive draws, or \code{NULL} to use
#'   all stored posterior draws.
#' @param seed Optional random seed.
#'
#' @return An integer matrix with predictive draws in rows and observations
#'   in columns.
#' @examples
#' set.seed(1)
#' x <- seq(-3, 3, length.out = 60)
#' y <- rbounded_binom(-0.5 + 1.2 * x, delta0 = 0.1, delta1 = 0.8)
#' fit <- bh_fit(y ~ x, data = data.frame(x = x, y = y),
#'               control = bh_control(iter = 400, warmup = 200, chains = 2),
#'               seed = 1)
#' yrep <- bh_posterior_predict(fit, ndraws = 10, seed = 4)
#' dim(yrep)
#' @export
bh_posterior_predict <- function(object,
                                 newdata = NULL,
                                 trials = NULL,
                                 ndraws = 100L,
                                 seed = NULL) {
  if (!inherits(object, "bounded_reg")) {
    stop("`object` must be a `bounded_reg` fit.", call. = FALSE)
  }
  if (!is.null(seed)) {
    set.seed(.validate_seed(seed))
  }
  probability <- .posterior_probability_matrix(object, newdata)
  if (is.null(trials)) {
    trials <- if (is.null(newdata)) {
      object$trials
    } else {
      rep.int(1L, ncol(probability))
    }
  }
  .check_counts(trials, "`trials`")
  if (!length(trials) %in% c(1L, ncol(probability))) {
    stop("`trials` must be scalar or match prediction rows.", call. = FALSE)
  }
  trials <- rep_len(trials, ncol(probability))

  total_draws <- nrow(probability)
  if (is.null(ndraws)) {
    selected <- seq_len(total_draws)
  } else {
    if (!is.numeric(ndraws) || length(ndraws) != 1L ||
        !is.finite(ndraws) || ndraws < 1L ||
        ndraws > .Machine$integer.max || ndraws != floor(ndraws)) {
      stop("`ndraws` must be NULL or a positive integer.", call. = FALSE)
    }
    ndraws <- as.integer(ndraws)
    selected <- sample.int(
      total_draws,
      ndraws,
      replace = ndraws > total_draws
    )
  }
  probability <- probability[selected, , drop = FALSE]

  out <- matrix(NA_integer_, nrow(probability), ncol(probability))
  for (j in seq_len(ncol(probability))) {
    out[, j] <- stats::rbinom(
      nrow(probability), size = trials[j], prob = probability[, j]
    )
  }
  colnames(out) <- paste0("row", seq_len(ncol(out)))
  out
}
