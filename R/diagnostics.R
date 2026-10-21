#' Diagnose computation and data-range information
#'
#' Check MCMC convergence and whether the data identify the bounded
#' logistic curve. Diagnostics use the original, unmodified posterior chains.
#' Endpoints fixed by the prior are excluded from convergence flags. Undefined
#' diagnostics on other parameters are flagged rather than treated as
#' evidence of convergence. Switching of the endpoint labels can give poor
#' diagnostics for the raw parameters even when the fitted probabilities
#' (see \code{\link{predict.bounded_reg}}) mix well.
#'
#' @param object A \code{bounded_reg} fit.
#' @param rhat_threshold R-hat warning threshold.
#' @param ess_threshold Bulk and tail ESS warning threshold.
#' @return An object of class \code{bh_diagnostic}: a list with the MCMC
#'   summary table (\code{mcmc}), the names of fixed and flagged parameters,
#'   the design-matrix condition number, posterior quantiles of the logistic
#'   gate span and of the endpoint gap \code{delta1 - delta0}, the posterior
#'   probability that \code{delta0 > delta1}, the number of endpoint-order
#'   switches per chain, and a character vector of warning messages.
#' @examples
#' set.seed(1)
#' x <- seq(-3, 3, length.out = 60)
#' y <- rbounded_binom(-0.5 + 1.2 * x, delta0 = 0.1, delta1 = 0.8)
#' fit <- bh_fit(y ~ x, data = data.frame(x = x, y = y),
#'               control = bh_control(iter = 400, warmup = 200, chains = 2),
#'               seed = 1)
#' bh_diagnose(fit)
#' @export
bh_diagnose <- function(object, rhat_threshold = 1.01, ess_threshold = 400) {
  if (!inherits(object, "bounded_reg")) {
    stop("`object` must be a `bounded_reg` fit.", call. = FALSE)
  }
  for (z in list(rhat_threshold, ess_threshold)) {
    if (!is.numeric(z) || length(z) != 1L || !is.finite(z) || z <= 0) {
      stop("Diagnostic thresholds must be positive finite scalars.", call. = FALSE)
    }
  }
  tab <- posterior::summarise_draws(object$draws)
  fixed <- c(if (length(object$resolved_prior$delta0) == 1L) "delta0",
             if (length(object$resolved_prior$delta1) == 1L) "delta1")
  active <- !tab$variable %in% fixed
  rhat_flag <- tab$variable[active & is.finite(tab$rhat) & tab$rhat > rhat_threshold]
  low <- pmin(tab$ess_bulk, tab$ess_tail)
  ess_flag <- tab$variable[active & is.finite(low) & low < ess_threshold]
  invalid_flag <- tab$variable[active & (!is.finite(tab$rhat) | !is.finite(low))]
  draws <- .draw_matrix(object)
  gate <- stats::plogis(draws[, .beta_columns(object), drop = FALSE] %*% t(object$X))
  span <- apply(gate, 1L, function(z) max(z) - min(z))
  gap <- draws[, "delta1"] - draws[, "delta0"]
  # Zero-length chain differences are intentionally not interpreted as switches.
  a <- array(as.numeric(object$draws), dim = dim(object$draws),
             dimnames = dimnames(object$draws))
  switches <- vapply(seq_len(dim(a)[2L]), function(ch) {
    s <- a[, ch, "delta1"] > a[, ch, "delta0"]
    sum(diff(s) != 0)
  }, numeric(1))
  msgs <- c(
    if (length(rhat_flag)) "At least one free parameter exceeds the R-hat threshold.",
    if (length(ess_flag)) "At least one free parameter has low bulk or tail ESS.",
    if (length(invalid_flag)) "At least one free parameter has undefined MCMC diagnostics.",
    if (dim(a)[2L] < 2L) "Only one chain: between-chain convergence cannot be assessed.",
    if (stats::median(span) < 0.5) "Observed predictors cover little of the logistic transition.",
    if (stats::median(abs(gap)) < 0.1) "Small endpoint separation: beta may be weakly identified.",
    if (!length(fixed)) paste("The likelihood is invariant to swapping endpoints and negating beta.",
                             "Use asymmetric endpoint priors to fix the labels; fitted probabilities are unaffected.")
  )
  structure(list(mcmc = tab, fixed_parameters = fixed,
                 rhat_flag = rhat_flag, ess_flag = ess_flag,
                 invalid_flag = invalid_flag,
                 design_condition_number = kappa(object$X),
                 gate_span = stats::quantile(span, c(.05, .5, .95)),
                 delta_gap = stats::quantile(gap, c(.05, .5, .95)),
                 probability_delta0_gt_delta1 = mean(gap < 0),
                 endpoint_switches = switches, warnings = msgs),
            class = "bh_diagnostic")
}

#' @rdname bh_diagnose
#' @param x A \code{bh_diagnostic} object.
#' @param ... Further arguments; currently unused.
#' @param digits Number of displayed digits.
#' @export
print.bh_diagnostic <- function(x, digits = 3L, ...) {
  cat("boundedHill diagnostics (original posterior draws)\n")
  cat("  P(delta0 > delta1):", round(x$probability_delta0_gt_delta1, digits), "\n")
  cat("  endpoint-order switches per chain:", paste(x$endpoint_switches, collapse = ", "), "\n")
  cat("  design condition number:", round(x$design_condition_number, digits), "\n")
  cat("  gate span (5%, 50%, 95%):", paste(round(x$gate_span, digits), collapse = ", "), "\n")
  cat("  endpoint gap (5%, 50%, 95%):", paste(round(x$delta_gap, digits), collapse = ", "), "\n")
  if (length(x$warnings)) for (s in x$warnings) cat(" -", s, "\n")
  invisible(x)
}
