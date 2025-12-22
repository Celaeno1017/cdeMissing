\name{cde_augment_data}
\alias{cde_augment_data}
\title{Augment long-format data with CDE-style missingness features}
\description{
\code{cde_augment_data} creates missingness indicators and timepoint dummy variables for missing
time-dependent covariates, and replaces missing covariate values with a constant (default 0).
The response variable is not modified.
}
\usage{
cde_augment_data(
  data,
  id,
  time,
  baseline_missing = character(0),
  time_dependent_missing = character(0),
  covariates = NULL,
  replace_na_with = 0
)
}
\arguments{
\item{data}{A long-format \code{data.frame}.}

\item{id}{Column name for subject ID.}

\item{time}{Column name for time variable. Required if \code{time_dependent_missing} is non-empty.}

\item{baseline_missing}{Character vector of baseline covariates that may be missing.}

\item{time_dependent_missing}{Character vector of time-dependent covariates that may be missing.}

\item{covariates}{Optional character vector specifying the covariate set used to define \code{mis_any}.
If \code{NULL}, uses the union of \code{baseline_missing} and \code{time_dependent_missing}.}

\item{replace_na_with}{Value used to replace \code{NA} in the selected covariates. Default is 0.}
}
\details{
Adds the following columns to the returned dataset:
\itemize{
  \item \code{mis_any}: row-wise indicator for any missingness among \code{covariates}.
  \item \code{mis_base_any}: row-wise indicator for any missing baseline covariate.
  \item \code{mis_td_any}: row-wise indicator for any missing time-dependent covariate.
  \item \code{mis_base_<var>}: per-variable baseline missingness indicators.
  \item \code{mis_td_<var>}: per-variable time-dependent missingness indicators.
  \item \code{u_<var>_t<time>}: timepoint dummies indicating that variable \code{<var>} is missing at time \code{<time>}.
  \item \code{ct}: row index \code{1:nrow(data)}.
}
}
\value{
A \code{data.frame} in the same long format as input with additional columns described above.
}
\examples{
dat <- data.frame(
  id = rep(1:3, each = 3),
  Time = rep(1:3, times = 3),
  y = rnorm(9),
  x1 = c(1, NA, 3, 4, 5, 6, 7, 8, 9),
  x2 = c(NA, 2, 3, 4, 5, NA, 7, 8, 9)
)

aug <- cde_augment_data(
  dat,
  id = "id",
  time = "Time",
  baseline_missing = "x1",
  time_dependent_missing = "x2"
)
names(aug)
}
\seealso{
\code{\link{cde_fit}}
}
