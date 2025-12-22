\name{cde_augmented_data}
\alias{cde_augmented_data}
\title{Retrieve augmented data from a \code{cde_fit} object}
\description{
Returns the augmented dataset stored in a \code{cde_fit} object, if available.
}
\usage{
cde_augmented_data(object)
}
\arguments{
\item{object}{A \code{cde_fit} object.}
}
\value{
A \code{data.frame} containing augmented data, or \code{NULL} if \code{return_augmented_data = FALSE}
was used when fitting.
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
  engine = "nlme",
  return_augmented_data = TRUE
)
aug <- cde_augmented_data(fit)
head(aug)
}
\seealso{
\code{\link{cde_augment_data}}, \code{\link{cde_fit}}
}
