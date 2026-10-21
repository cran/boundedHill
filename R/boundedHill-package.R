#' Bayesian bounded logistic regression
#'
#' The \pkg{boundedHill} package fits Bernoulli and grouped binomial
#' regression models of the form
#' \deqn{P(Y_i = 1 \mid x_i) =
#' \delta_0 + (\delta_1 - \delta_0)
#' \operatorname{logit}^{-1}(x_i^\top \beta),}
#' with \eqn{0 \leq \delta_0, \delta_1 \leq 1}.
#'
#' @keywords internal
#' @importFrom stats coef fitted predict simulate
#' @importFrom graphics plot
"_PACKAGE"
