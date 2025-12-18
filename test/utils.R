\
#' Internal helper: null-coalescing operator
#' @keywords internal
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

#' Internal helper: ensure a character scalar
#' @keywords internal
.as_chr1 <- function(x, arg) {
  if (!is.character(x) || length(x) != 1) {
    stop(sprintf("`%s` must be a single character string.", arg), call. = FALSE)
  }
  x
}

#' Internal helper: ensure character vector (possibly empty)
#' @keywords internal
.as_chrv <- function(x, arg) {
  if (is.null(x)) return(character(0))
  if (!is.character(x)) stop(sprintf("`%s` must be a character vector.", arg), call. = FALSE)
  x
}

#' Internal helper: add columns safely
#' @keywords internal
.add_col <- function(dat, name, value) {
  if (!name %in% names(dat)) dat[[name]] <- value
  dat
}
