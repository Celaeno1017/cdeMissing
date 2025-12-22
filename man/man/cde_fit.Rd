\name{cde_fit}
\alias{cde_fit}
\title{Fit CDE-style LMM/GLMM with missing covariates}
\description{
\code{cde_fit} augments a long-format dataset with missingness indicators and timepoint dummies,
then fits a Gaussian LMM (via \code{nlme::lme}) or a GLMM (via \code{lme4::glmer}) using a
consistent interface.
}
\usage{
cde_fit(
  formula,
  data,
  id,
  time,
  random,
  correlation,
  missing = list(baseline = character(0), time_dependent = character(0)),
  family = stats::gaussian(),
  engine = c("auto", "nlme", "lme4"),
  control = list(),
  return_augmented_data = TRUE
)
}
\arguments{
\item{formula}{A model formula for the substantive fixed effects (e.g., \code{y ~ x1 + Treatment + Time + x2}).
Missingness-related terms are appended automatically.}

\item{data}{Long-format \code{data.frame}.}

\item{id}{Column name for subject ID.}

\item{time}{Column name for time variable. Required if \code{missing$time_dependent} is non-empty.}

\item{random}{A model formula for the substantive random effects. Should be in the \code{"nlme"} style.}

\item{correlation}{Additional correlation structure passed to \code{"nlme"}. Not available for GLMM.}

\item{missing}{A list with elements \code{baseline} and \code{time_dependent}, each a character vector
indicating covariates that can be missing.}

\item{family}{A GLM family. Use \code{gaussian()} for LMM; e.g., \code{binomial()} for GLMM.}

\item{engine}{One of \code{"auto"}, \code{"nlme"}, \code{"lme4"}. If \code{"auto"}, selects \code{"nlme"} for
Gaussian and \code{"lme4"} otherwise.}

\item{control}{A list of control arguments. For \code{nlme} fits, passed to \code{nlme::lmeControl()}.
For \code{lme4} fits, currently supports \code{nAGQ} (default 0).}

\item{return_augmented_data}{Logical; whether to store the augmented dataset in the returned object.}
}
\details{
Data augmentation is performed by \code{\link{cde_augment_data}} and typically includes:
\itemize{
  \item \code{mis_any}: any missingness among selected covariates in the row.
  \item \code{mis_base_*} and \code{mis_base_any}: baseline missingness indicators.
  \item \code{mis_td_*} and \code{mis_td_any}: time-dependent missingness indicators.
  \item \code{u_*}: timepoint dummies indicating when a time-dependent covariate is missing.
  \item \code{ct}: per-row observation index (used by GLMM defaults).
}
Default model components:
\itemize{
  \item Gaussian: \code{nlme::lme} with random effects \code{~ 1 + mis_base_any + mis_td_any | id} and
  group-dependent residual variance via \code{varIdent} on \code{interaction(mis_base_any, mis_td_any)}.
  \item Non-Gaussian: \code{lme4::glmer} with random effects \code{(mis_any - 1 | id) + (mis_td_any - 1 | ct)}
  and \code{nAGQ = 0} by default.
}
}
\value{
An object of class \code{"cde_fit"} with components:
\itemize{
  \item \code{fit}: the underlying fitted model object (\code{nlme::lme} or \code{lme4::glmerMod}).
  \item \code{engine}, \code{family}.
  \item \code{data_aug}: augmented data if \code{return_augmented_data = TRUE}, otherwise \code{NULL}.
  \item \code{call}: the matched call.
}
}
\examples{
## Gaussian / LMM
sim <- sim_cde_data(model = "lmm", beta = c(1,0.5,0.2,0.2,0.2))
dat <- sim$data_mis

fit <- cde_fit(
  y ~ x1 + Treatment + Time + x2,
  data = dat,
  id = "id",
  time = "Time",
  random = ~1|id,
  missing = list(baseline = "x1", time_dependent = "x2"),
  family = gaussian(),
  engine = "nlme",
  control = list(msMaxIter = 100)
)

cde_fixef(fit)

## Binary / GLMM
sim2 <- sim_cde_data(model = "glmm",beta = c(-1,0.5,0.2,0.2,0.2))
dat2 <- sim2$data_mis

fit2 <- cde_fit(
  y ~ x1 + Treatment + Time + x2,
  data = dat2,
  id = "id",
  time = "Time",
  missing = list(baseline = "x1", time_dependent = "x2"),
  family = binomial(),
  engine = "lme4",
  control = list(nAGQ = 0)
)
}
\seealso{
\code{\link{cde_augment_data}}, \code{\link{cde_fixef}}, \code{\link{cde_augmented_data}}
}
