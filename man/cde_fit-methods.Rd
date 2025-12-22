\name{cde_fit-methods}
\alias{print.cde_fit}
\alias{summary.cde_fit}
\alias{coef.cde_fit}
\title{Methods for \code{cde_fit} objects}
\description{
S3 methods for \code{cde_fit} objects: printing, summaries, and coefficient extraction.
}
\usage{
\method{print}{cde_fit}(x, ...)

\method{summary}{cde_fit}(object, ...)

\method{coef}{cde_fit}(object, ...)
}
\arguments{
\item{x}{A \code{cde_fit} object.}
\item{object}{A \code{cde_fit} object.}
\item{...}{Additional arguments passed to the underlying model methods.}
}
\details{
\code{summary.cde_fit} dispatches to \code{summary()} of the underlying fitted model object.
\code{coef.cde_fit} dispatches to \code{coef()} of the underlying fitted model object.
}
\seealso{
\code{\link{cde_fit}}
}
