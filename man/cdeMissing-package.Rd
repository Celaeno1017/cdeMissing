\name{cdeMissing-package}
\alias{cdeMissing}
\alias{cdeMissing-package}
\title{cdeMissing: random effect approach to generalized linear mixed model Analysis of Incomplete Longitudinal Data}
\description{
The \pkg{cdeMissing} package provides utilities to (i) augment long-format longitudinal data
with missingness indicators and timepoint dummy variables, and (ii) fit modified LMM/GLMM
models using \code{nlme} (Gaussian LMM) or \code{lme4} (GLMM).
}
\details{
Key functions:
\itemize{
  \item \code{\link{cde_fit}}: augment + fit LMM/GLMM.
  \item \code{\link{cde_augment_data}}: augmentation only.
  \item \code{\link{sim_cde_data}}: simulate long-format data compatible with \code{cde_fit}.
  \item \code{\link{cde_fixef}}: extract fixed-effect estimates.
  \item \code{\link{cde_augmented_data}}: retrieve augmented data stored in a \code{cde_fit} object.
}
}
\keyword{package}
