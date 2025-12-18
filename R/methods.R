#' @export
print.cde_fit <- function(x, ...) {
  cat("<cde_fit>\n")
  cat("  engine: ", x$engine, "\n", sep = "")
  cat("  family: ", x$family$family, "\n", sep = "")
  cat("  fixed effects (head):\n")
  print(utils::head(cde_fixef(x)))
  invisible(x)
}

#' @export
summary.cde_fit <- function(object, ...) {
  summary(object$fit, ...)
}

#' @export
coef.cde_fit <- function(object, ...) {
  stats::coef(object$fit, ...)
}
