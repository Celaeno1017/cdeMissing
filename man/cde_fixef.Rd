\name{cde_fixef}
\alias{cde_fixef}
\title{Extract fixed effects from a \code{cde_fit} object}
\description{
Returns fixed-effect estimates from the underlying fitted model.
}
\usage{
cde_fixef(object, ...)
}
\arguments{
\item{object}{A \code{cde_fit} object.}
\item{...}{Unused.}
}
\value{
A named numeric vector of fixed-effect estimates.
}
\examples{
sim <- sim_cde_data(model = "lmm", m = 30, t = 5, q = 0.1, seed = 1)
fit <- cde_fit(
  y ~ x1 + Treatment + Time + x2,
  data = sim$data_mis,
  id = "id",
  time = "Time",
  missing = list(baseline = "x1", time_dependent = "x2"),
  family = gaussian(),
  engine = "nlme"
)
cde_fixef(fit)
}
\seealso{
\code{\link{cde_fit}}
}
