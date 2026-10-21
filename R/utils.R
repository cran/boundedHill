.check_counts <- function(x, name) {
  if (!is.numeric(x) || !length(x) || any(!is.finite(x)) || any(x < 0) ||
      any(x > .Machine$integer.max) ||
      any(abs(x - round(x)) > sqrt(.Machine$double.eps))) {
    stop(sprintf("%s must contain non-negative integer counts.", name),
         call. = FALSE)
  }
  invisible(TRUE)
}

.prepare_model_frame <- function(mf, contrasts) {
  trm <- attr(mf, "terms")
  if (!is.null(stats::model.offset(mf))) {
    stop("Offsets are not supported in version 0.1.0.", call. = FALSE)
  }

  if (!nrow(mf)) stop("No observations remain in the model frame.", call. = FALSE)
  response <- stats::model.response(mf)
  if (anyNA(response)) stop("Missing responses are not supported here.", call. = FALSE)
  if (is.matrix(response)) {
    if (ncol(response) != 2L) {
      stop(
        "A matrix response must be `cbind(successes, failures)`.",
        call. = FALSE
      )
    }
    successes <- response[, 1L]
    failures <- response[, 2L]
    .check_counts(successes, "Successes")
    .check_counts(failures, "Failures")
    trials <- successes + failures
    .check_counts(trials, "Total trials")
    if (any(trials == 0)) {
      stop("Every grouped-binomial row must contain at least one trial.",
           call. = FALSE)
    }
  } else {
    if (is.factor(response)) {
      if (nlevels(response) != 2L) {
        stop("A factor response must have exactly two levels.",
             call. = FALSE)
      }
      successes <- as.integer(response == levels(response)[2L])
    } else if (is.logical(response)) {
      successes <- as.integer(response)
    } else {
      if (!is.numeric(response) || any(!is.finite(response)) ||
          any(!(response %in% c(0, 1)))) {
        stop("A Bernoulli response must contain only 0 and 1.",
             call. = FALSE)
      }
      successes <- as.integer(response)
    }
    trials <- rep.int(1L, length(successes))
    failures <- trials - successes
  }

  X <- stats::model.matrix(trm, mf, contrasts.arg = contrasts)
  if (!ncol(X)) {
    stop("The formula must produce at least one regression coefficient.",
         call. = FALSE)
  }
  if (any(!is.finite(X))) {
    stop("The model matrix contains non-finite values.", call. = FALSE)
  }
  if (qr(X)$rank < ncol(X)) {
    warning(
      "The model matrix is rank deficient; coefficients may be weakly identified.",
      call. = FALSE
    )
  }

  list(
    model_frame = mf,
    terms = trm,
    X = unname(X),
    coef_names = colnames(X),
    successes = as.integer(successes),
    failures = as.integer(failures),
    trials = as.integer(trials),
    response = response,
    contrasts = attr(X, "contrasts"),
    xlevels = stats::.getXlevels(trm, mf),
    na.action = attr(mf, "na.action")
  )
}

.prepare_model <- function(formula, data, na.action, contrasts) {
  mf <- stats::model.frame(
    formula = formula,
    data = data,
    na.action = na.action,
    drop.unused.levels = TRUE
  )
  .prepare_model_frame(mf, contrasts)
}

.rmvn_precision <- function(mean, precision_chol) {
  mean + backsolve(precision_chol, stats::rnorm(length(mean)))
}

.rmvn_covariance <- function(mean, covariance) {
  mean + drop(t(chol(covariance)) %*% stats::rnorm(length(mean)))
}

.parameter_names <- function(coef_names) {
  clean_names <- make.names(coef_names, unique = TRUE)
  clean_names[coef_names == "(Intercept)"] <- "Intercept"
  c("delta0", "delta1", paste0("beta_", make.unique(clean_names)))
}

.draw_matrix <- function(object) {
  x <- posterior::as_draws_matrix(object$draws)
  matrix(as.numeric(x), nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x))
}

.beta_columns <- function(object) {
  object$beta_names
}

.new_model_matrix <- function(object, newdata = NULL) {
  if (is.null(newdata)) {
    return(object$X)
  }
  trm <- stats::delete.response(object$terms)
  mf <- stats::model.frame(trm, data = newdata, na.action = stats::na.pass,
                           xlev = object$xlevels)
  if (!nrow(mf)) stop("`newdata` must contain observations.", call. = FALSE)
  if (anyNA(mf)) stop("Missing newdata predictors are not supported.", call. = FALSE)
  X <- stats::model.matrix(trm, data = mf, contrasts.arg = object$contrasts)
  missing_columns <- setdiff(object$coef_names, colnames(X))
  extra_columns <- setdiff(colnames(X), object$coef_names)
  if (length(missing_columns) || length(extra_columns)) {
    stop("`newdata` does not produce the same model-matrix columns as the fit.",
         call. = FALSE)
  }
  X <- X[, object$coef_names, drop = FALSE]
  if (any(!is.finite(X))) {
    stop("The new model matrix contains non-finite values.", call. = FALSE)
  }
  unname(X)
}

.posterior_probability_matrix <- function(object, newdata = NULL) {
  X <- .new_model_matrix(object, newdata)
  draws <- .draw_matrix(object)
  beta <- draws[, .beta_columns(object), drop = FALSE]
  eta <- beta %*% t(X)
  d0 <- as.numeric(draws[, "delta0"])
  d1 <- as.numeric(draws[, "delta1"])
  sweep(sweep(stats::plogis(eta), 1L, d1 - d0, `*`), 1L, d0, `+`)
}

.posterior_linear_matrix <- function(object, newdata = NULL) {
  X <- .new_model_matrix(object, newdata)
  draws <- .draw_matrix(object)
  beta <- draws[, .beta_columns(object), drop = FALSE]
  beta %*% t(X)
}

.summarize_columns <- function(x, level = 0.95) {
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) {
    stop("`level` must lie strictly between 0 and 1.", call. = FALSE)
  }
  alpha <- (1 - level) / 2
  qs <- apply(x, 2L, stats::quantile,
              probs = c(alpha, 0.5, 1 - alpha), names = FALSE)
  if (is.null(dim(qs))) {
    qs <- matrix(qs, nrow = 3L)
  }
  data.frame(
    estimate = colMeans(x),
    median = qs[2L, ],
    lower = qs[1L, ],
    upper = qs[3L, ],
    row.names = colnames(x)
  )
}

.validate_seed <- function(seed) {
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed < 0 || seed > .Machine$integer.max || seed != floor(seed)) {
    stop("`seed` must be an integer from 0 to .Machine$integer.max.", call. = FALSE)
  }
  as.integer(seed)
}
