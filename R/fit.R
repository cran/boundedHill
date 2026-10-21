#' Fit Bayesian bounded logistic regression
#'
#' Fit a bounded logistic regression model to Bernoulli or grouped binomial
#' outcomes using a fully Gibbs posterior sampler. The right-hand side of
#' \code{formula} is processed by \code{\link[stats]{model.matrix}}, so
#' continuous predictors, factors, interactions, and standard transformations
#' are supported.
#'
#' @param formula Model formula. Use \code{y ~ x1 + x2} for Bernoulli data or
#'   \code{cbind(successes, failures) ~ x1 + x2} for grouped binomial data.
#' @param data Data frame containing variables in \code{formula}.
#' @param prior Prior specification from \code{\link{bh_prior}}.
#' @param control Sampling controls from \code{\link{bh_control}}.
#' @param seed Optional integer seed. Chain-specific seeds are generated from
#'   it, so the same \code{seed} reproduces the same draws. Fitting changes
#'   the state of the random number generator in the calling session, also
#'   when \code{seed = NULL}.
#' @param subset Optional row subset passed to \code{\link[stats]{model.frame}}.
#' @param na.action Missing-value action.
#' @param contrasts Optional contrasts passed to
#'   \code{\link[stats]{model.matrix}}.
#'
#' @references
#' Polson, N. G., Scott, J. G., and Windle, J. (2013). Bayesian inference for
#' logistic models using Polya-Gamma latent variables. \emph{Journal of the
#' American Statistical Association}, 108(504), 1339--1349.
#' \doi{10.1080/01621459.2013.829001}
#'
#' Lee, M., Lee, C. H., and Kim, R. S. A novel Bayesian phase 1/2 design
#' incorporating the Hill equation. Manuscript in preparation.
#'
#' Lee, M., Lee, C. H., and Kim, R. S. boundedHill: Bayesian logistic
#' regression with lower and upper asymptotes in R. Manuscript in preparation.
#'
#' @return An object of class \code{bounded_reg}: a list containing the
#'   call, formula, model matrix, response counts, prior, control settings,
#'   chain seeds, and the posterior draws (a \pkg{posterior}
#'   \code{draws_array} in element \code{draws}). Use the methods in
#'   \code{\link{bounded_reg-methods}}, \code{\link{plot.bounded_reg}} and
#'   \code{\link{bh_diagnose}} to work with it.
#' @examples
#' set.seed(1)
#' x <- seq(-3, 3, length.out = 60)
#' y <- rbounded_binom(-0.5 + 1.2 * x, delta0 = 0.1, delta1 = 0.8)
#' fit <- bh_fit(y ~ x, data = data.frame(x = x, y = y),
#'               control = bh_control(iter = 400, warmup = 200, chains = 2),
#'               seed = 1)
#' summary(fit)
#'
#' # Grouped binomial data
#' grouped <- data.frame(dose = c(50, 100, 150, 200, 250),
#'                       resp = c(2, 11, 14, 15, 14), n = c(15, 20, 21, 22, 20))
#' fit2 <- bh_fit(cbind(resp, n - resp) ~ log(dose / 150), data = grouped,
#'                control = bh_control(iter = 400, warmup = 200, chains = 2),
#'                seed = 2)
#' coef(fit2)
#' @export
bh_fit <- function(formula,
                   data,
                   prior = bh_prior(),
                   control = bh_control(),
                   seed = NULL,
                   subset,
                   na.action = stats::na.omit,
                   contrasts = NULL) {
  call <- match.call()
  if (!inherits(formula, "formula")) {
    stop("`formula` must be a model formula.", call. = FALSE)
  }
  if (!inherits(control, "bh_control")) {
    stop("`control` must be created by `bh_control()`.", call. = FALSE)
  }
  if (!is.null(seed)) {
    set.seed(.validate_seed(seed))
  }

  model_frame_arguments <- match(
    c("formula", "data", "subset", "na.action"),
    names(call),
    nomatch = 0L
  )
  model_frame_arguments <- model_frame_arguments[
    model_frame_arguments > 0L
  ]
  model_frame_call <- call[c(1L, model_frame_arguments)]
  model_frame_call[[1L]] <- quote(stats::model.frame)
  if (is.null(model_frame_call$na.action)) model_frame_call$na.action <- na.action
  model_frame_call$drop.unused.levels <- TRUE
  model_frame <- eval(model_frame_call, parent.frame())
  prepared <- .prepare_model_frame(model_frame, contrasts)
  colnames(prepared$X) <- prepared$coef_names
  resolved_prior <- .resolve_prior(
    prior, ncol(prepared$X), prepared$coef_names
  )

  chain_seeds <- sample.int(
    .Machine$integer.max,
    size = control$chains,
    replace = FALSE
  )
  chain_results <- lapply(chain_seeds, function(chain_seed) {
    .run_bh_chain(
      X = prepared$X,
      successes = prepared$successes,
      failures = prepared$failures,
      trials = prepared$trials,
      prior = resolved_prior,
      control = control,
      seed = chain_seed
    )
  })

  n_keep <- nrow(chain_results[[1L]]$draws)
  parameter_names <- colnames(chain_results[[1L]]$draws)
  draw_array <- array(
    NA_real_,
    dim = c(n_keep, control$chains, length(parameter_names)),
    dimnames = list(
      iteration = as.character(seq_len(n_keep)),
      chain = as.character(seq_len(control$chains)),
      variable = parameter_names
    )
  )
  for (chain in seq_len(control$chains)) {
    draw_array[, chain, ] <- chain_results[[chain]]$draws
  }
  draw_array <- posterior::as_draws_array(draw_array)

  structure(
    list(
      call = call,
      formula = formula,
      terms = prepared$terms,
      model_frame = prepared$model_frame,
      X = prepared$X,
      coef_names = prepared$coef_names,
      beta_names = .parameter_names(prepared$coef_names)[-(1:2)],
      response = prepared$response,
      successes = prepared$successes,
      failures = prepared$failures,
      trials = prepared$trials,
      prior = prior,
      resolved_prior = resolved_prior,
      control = control,
      seeds = chain_seeds,
      draws = draw_array,
      contrasts = prepared$contrasts,
      xlevels = prepared$xlevels,
      na.action = prepared$na.action
    ),
    class = "bounded_reg"
  )
}
