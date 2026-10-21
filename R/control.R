#' Control Gibbs sampling
#'
#' Set the number of iterations, warmup, thinning, and chains used by
#' \code{\link{bh_fit}}.
#'
#' @param iter Total iterations per chain, including warmup.
#' @param warmup Number of initial iterations discarded per chain.
#' @param thin Positive thinning interval.
#' @param chains Number of chains.
#'
#' @return An object of class \code{bh_control}: a list with integer
#'   elements \code{iter}, \code{warmup}, \code{thin} and \code{chains}.
#'   The \code{print} method returns its input invisibly.
#' @examples
#' bh_control(iter = 2000, warmup = 1000, chains = 2)
#' @export
bh_control <- function(iter = 4000L,
                       warmup = 2000L,
                       thin = 1L,
                       chains = 4L) {
  values <- list(
    iter = iter,
    warmup = warmup,
    thin = thin,
    chains = chains
  )
  for (nm in names(values)) {
    value <- values[[nm]]
    if (length(value) != 1L || !is.numeric(value) || !is.finite(value) ||
        abs(value) > .Machine$integer.max || value != floor(value)) {
      stop(sprintf("`%s` must be a finite integer.", nm), call. = FALSE)
    }
    values[[nm]] <- as.integer(value)
  }
  if (values$iter < 2L) {
    stop("`iter` must be at least 2.", call. = FALSE)
  }
  if (values$warmup < 0L || values$warmup >= values$iter) {
    stop("`warmup` must lie between 0 and `iter - 1`.", call. = FALSE)
  }
  if (values$thin < 1L || values$chains < 1L) {
    stop("`thin` and `chains` must be positive.",
         call. = FALSE)
  }

  structure(values, class = "bh_control")
}

#' @rdname bh_control
#' @param x A \code{bh_control} object.
#' @param ... Further arguments; currently unused.
#' @export
print.bh_control <- function(x, ...) {
  kept <- length(seq.int(x$warmup + 1L, x$iter, by = x$thin))
  cat("Gibbs sampling control\n")
  cat("  chains:", x$chains, "\n")
  cat("  iterations per chain:", x$iter, "\n")
  cat("  warmup:", x$warmup, "\n")
  cat("  thinning:", x$thin, "\n")
  cat("  retained draws:", kept * x$chains, "\n")
  invisible(x)
}
