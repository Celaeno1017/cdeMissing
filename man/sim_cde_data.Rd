\name{sim_cde_data}
\alias{sim_cde_data}
\title{Simulate longitudinal data for \code{cde_fit}}
\description{
Simulates long-format longitudinal data with baseline and time-dependent covariates and
induces missingness in baseline and time-dependent covariates. The output \code{data_mis}
is compatible with \code{\link{cde_fit}}.
}
\usage{
sim_cde_data(
  model = c("lmm", "glmm"),
  beta = c(1, 0.5, 0.2, 0.2, 0.2),
  m = 400,
  t = 5,
  p = 0.1,
  q = 0.1,
  sd_alpha = 0.3,
  sd_epsilon = 0.1,
  seed = NULL
)
}
\arguments{
\item{model}{Either \code{"lmm"} for Gaussian outcomes or \code{"glmm"} for Bernoulli-logit outcomes.}

\item{beta}{Numeric vector of length 5: \code{(beta0, beta_x1, beta_trt, beta_time, beta_x2)}.}

\item{m}{Number of subjects.}

\item{t}{Number of timepoints per subject.}

\item{p}{Deprecated placeholder; included for signature compatibility.}

\item{q}{Slope for time-dependent missingness probability in \code{x2}; at time \code{j} the missingness
probability is \code{min(1, q*j)}.}

\item{sd_alpha}{Random-intercept SD used for \code{"lmm"} generation.}

\item{sd_epsilon}{Noise SD used for \code{"lmm"} generation.}

\item{seed}{Optional integer seed.}
}
\value{
A list with elements:
\itemize{
  \item \code{data_mis}: long-format \code{data.frame} with induced missingness in covariates.
  \item \code{data_comp}: complete-data version.
  \item \code{missing}: \code{data.frame} with missingness indicators \code{I_1} (baseline) and \code{I_4} (time-dependent).
}
}
\examples{
sim <- sim_cde_data(model = "lmm", m = 50, t = 5, q = 0.1, seed = 1)
head(sim$data_mis)

fit <- cde_fit(
  y ~ x1 + Treatment + Time + x2,
  data = sim$data_mis,
  id = "id",
  time = "Time",
  missing = list(baseline = "x1", time_dependent = "x2"),
  family = gaussian(),
  engine = "nlme"
)
}
\seealso{
\code{\link{cde_fit}}
}
