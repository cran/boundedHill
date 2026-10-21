# boundedHill

`boundedHill` fits Bayesian bounded logistic regression models, also called
four-parameter logistic or bounded Hill models, to binary and grouped-binomial
data:

$$
\Pr(Y_i = 1 \mid x_i)
= \delta_0 + (\delta_1-\delta_0)
  \operatorname{logit}^{-1}(x_i^\top\beta),
\qquad 0 \leq \delta_0, \delta_1 \leq 1.
$$

In ordinary logistic regression the response probability always goes from 0
to 1. Here the curve instead goes from a lower level `delta0` to an upper level
`delta1`. For example, a combination therapy can give some response even at a
very low dose (`delta0 > 0`), and the response can level off below 100% at high
doses (`delta1 < 1`). Each of `delta0` and `delta1` can be estimated with a
Beta prior or fixed at a chosen value. Fixing `delta0 = 0` and `delta1 = 1`
gives ordinary logistic regression.

Estimation uses a Gibbs sampler, so no proposal distributions need tuning.
In each iteration, every observation is assigned to either the lower level or
the upper level, `delta0` and `delta1` are drawn from Beta distributions, and
the regression coefficients are drawn from a normal distribution using
Pólya-Gamma augmentation (see `?bh_fit` for references).

### A note on labels

The same curve can be written in two ways: swapping `delta0` and `delta1` and
changing the sign of `beta` gives exactly the same probabilities. The data
alone therefore cannot tell which level is which. Giving `delta0` a prior with
a low mean and `delta1` a prior with a high mean, as the default prior does,
fixes the labels. Either way, the fitted probabilities from `predict()` are
the same; only the labels of `delta0`, `delta1` and `beta` are affected.

## Installation

After the package is released on CRAN:

```r
install.packages("boundedHill")
```

For a local source tree:

```r
install.packages("boundedHill_0.1.0.tar.gz",
                 repos = NULL, type = "source")
```

## Bernoulli regression

```r
library("boundedHill")

set.seed(18)
n <- 160
dat <- data.frame(
  x = runif(n, -2.5, 2.5),
  group = factor(sample(c("A", "B"), n, replace = TRUE))
)
eta <- -0.4 + 1.2 * dat$x + 0.5 * (dat$group == "B")
dat$y <- rbounded_binom(eta, delta0 = 0.12, delta1 = 0.78)

fit <- bh_fit(
  y ~ x + group,
  data = dat,
  prior = bh_prior(
    delta0 = c(2, 10),
    delta1 = c(8, 2),
    beta_sd = 2
  ),
  control = bh_control(iter = 3000, warmup = 1500, chains = 4),
  seed = 2026
)

print(fit)
summary(fit)
bh_diagnose(fit)
plot(fit, type = "trace")
predict(fit, newdata = data.frame(
  x = c(-2, 0, 2),
  group = factor("A", levels = levels(dat$group))
))
```

## Grouped binomial regression

Grouped responses are supplied without expanding individual trials:

```r
grouped <- data.frame(
  x = seq(-2, 2, length.out = 7),
  success = c(3, 5, 9, 14, 17, 18, 19),
  failure = c(17, 15, 11, 6, 3, 2, 1)
)

fit_grouped <- bh_fit(
  cbind(success, failure) ~ x,
  data = grouped,
  prior = bh_prior(delta0 = c(2, 8), delta1 = c(8, 2)),
  seed = 45
)
```

## Model family

The same engine covers several nested models:

```r
# Standard Bayesian logistic regression
prior_logistic <- bh_prior(delta0 = 0, delta1 = 1)

# Unknown lower asymptote, upper asymptote fixed at one
prior_lower <- bh_prior(delta0 = c(1, 4), delta1 = 1)

# Both asymptotes unknown (Beta shapes)
prior_bounded <- bh_prior(delta0 = c(1, 4), delta1 = c(4, 1))

# The same, specified by prior mean and SD
prior_moments <- bh_prior(delta0 = bh_beta_moments(0.2, 0.15),
                          delta1 = bh_beta_moments(0.8, 0.15))
```

The default `bh_prior()` uses `delta0` with mean 0.10 and `delta1` with mean
0.85, both with SD 0.12, and independent N(0, 2^2) priors for the regression
coefficients.

## Public interface

| Function | Main input | Output |
|---|---|---|
| `bh_fit()` | formula, data, prior, MCMC control | `bounded_reg` fit |
| `bh_prior()` | asymptote priors and coefficient prior | `bh_prior` |
| `bh_beta_moments()` | prior mean and SD of an asymptote | Beta shapes |
| `bh_control()` | iterations, warmup, thinning, chains | `bh_control` |
| `predict()` | fit and optional new predictor data | posterior probability summaries or draws |
| `bh_posterior_draws()` | fit, output format | posterior draws |
| `bh_posterior_predict()` | fit, trial counts, number of draws | posterior predictive count matrix |
| `bh_log_lik()` | fit and optional response data | pointwise log-likelihood matrix |
| `bh_diagnose()` | fit and diagnostic thresholds | MCMC and predictor-range diagnostics |
| `bounded_probability()` | linear predictor and asymptotes | bounded probabilities |
| `rbounded_binom()` | linear predictor, trials, asymptotes | simulated binomial counts |

## References

Lee, M., Lee, C. H., and Kim, R. S. A novel Bayesian phase 1/2 design
incorporating the Hill equation. Manuscript in preparation.

Lee, M., Lee, C. H., and Kim, R. S. boundedHill: Bayesian logistic
regression with lower and upper asymptotes in R. Manuscript in preparation.
