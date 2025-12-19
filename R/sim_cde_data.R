#' Simulate longitudinal data for `cde_fit`
#'
#' This function mirrors the `sim_data()` generators used in the provided LMM and GLMM
#' simulation scripts: baseline covariate `x1`, baseline binary `Treatment`, a time
#' index `Time`, and a time-dependent covariate `x2`. Missingness is induced in `x1`
#' at the subject level and in `x2` at the observation level with time-varying rates.
#'
#' The returned `data_mis` is a long-format `data.frame` that can be passed directly to
#' `cde_fit()`, e.g. `cde_fit(y ~ x1 + Treatment + Time + x2, data = data_mis, ...)`.
#'
#' @param model One of `"lmm"` (Gaussian) or `"glmm"` (Bernoulli-logit).
#' @param beta Numeric vector of length 5: `(beta0, beta_x1, beta_trt, beta_time, beta_x2)`.
#' @param m Number of subjects.
#' @param t Number of timepoints per subject.
#' @param q Time-dependent missingness slope for `x2`; at time `j` probability is `min(1, q*j)`.
#' @param p Deprecated placeholder (kept for signature compatibility with scripts). Ignored by default.
#' @param sd_alpha Random-intercept SD for `"lmm"` generation (default 0.3).
#' @param sd_epsilon Noise SD for `"lmm"` generation (default 0.1).
#' @param seed Optional integer seed for reproducibility.
#'
#' @return A list with:
#' \describe{
#'   \item{data_mis}{Long-format `data.frame` with induced missingness in `x1` and/or `x2`.}
#'   \item{data_comp}{Corresponding complete-data `data.frame` (no missing covariates).}
#'   \item{missing}{A `data.frame` of missingness indicators `I_1` (x1 missing) and `I_4` (x2 missing).}
#' }
#' @export
sim_cde_data <- function(
  model = c("lmm", "glmm"),
  beta = c(1, 0.5, 0.2, 0.2, 0.2),
  m = 400,
  t = 5,
  p = 0.1,
  q = 0.1,
  sd_alpha = 0.3,
  sd_epsilon = 0.1,
  seed = NULL
) {
  model <- match.arg(model)
  if (!is.null(seed)) set.seed(seed)

  if (!is.numeric(beta) || length(beta) != 5) {
    stop("`beta` must be a numeric vector of length 5.", call. = FALSE)
  }

  # Covariates (match scripts)
  x1 <- rep(stats::runif(m, 0, 1), each = t)
  trt_base <- stats::rbinom(m, 1, 0.5)
  Treatment <- rep(trt_base, each = t)
  Time <- rep(seq_len(t), times = m)
  x2 <- stats::rnorm(m * t, mean = 1, sd = 1)
  id <- rep(seq_len(m), each = t)

  # Outcome (generate before inducing missingness)
  if (model == "lmm") {
    alpha <- rep(stats::rnorm(m, 0, sd_alpha), each = t)
    epsilon <- stats::rnorm(m * t, 0, sd_epsilon)
    y <- beta[1] + beta[2] * x1 + beta[3] * Treatment + beta[4] * Time + beta[5] * x2 + alpha + epsilon
  } else {
    eta <- beta[1] + beta[2] * x1 + beta[3] * Treatment + beta[4] * Time + beta[5] * x2
    pr <- stats::plogis(eta)
    y <- stats::rbinom(length(pr), 1, pr)
  }

  # Missingness mechanisms (match scripts)
  # x1 missingness: subject-level, depends on Treatment via logistic link
  p_subj <- stats::plogis(3 * trt_base - 2)
  I_1_subj <- stats::rbinom(length(trt_base), 1, p_subj)
  I_1 <- rep(I_1_subj, each = t)

  # x2 missingness: observation-level, depends on time
  # At time j: min(1, q*j)
  I_4 <- as.vector(t(sapply(seq_len(t), function(j) {
    stats::rbinom(m, 1, ifelse(q * j > 1, 1, q * j))
  })))

  # Apply missingness
  x1_mis <- x1
  x1_mis[I_1 == 1] <- NA_real_

  x2_mis <- x2
  x2_mis[I_4 == 1] <- NA_real_

  data_comp <- data.frame(
    y = y,
    id = id,
    x1 = x1,
    Treatment = Treatment,
    Time = Time,
    x2 = x2
  )

  data_mis <- data.frame(
    y = y,
    id = id,
    x1 = x1_mis,
    Treatment = Treatment,
    Time = Time,
    x2 = x2_mis
  )

  missing <- data.frame(
    id = id,
    Time = Time,
    I_1 = I_1,
    I_4 = I_4
  )

  list(data_mis = data_mis, data_comp = data_comp, missing = missing)
}
