#' Build an augmented fixed-effect RHS by appending u_* and baseline missingness terms
#' @keywords internal
.cde_fixed_rhs <- function(formula, data_aug) {
  base_terms <- attr(stats::terms(formula), "term.labels")

  # u_* dummies and baseline missingness flags
  u_terms <- grep("^u_", names(data_aug), value = TRUE)
  base_mis_terms <- grep("^mis_base_", names(data_aug), value = TRUE)

  rhs <- paste(c(base_terms, u_terms, base_mis_terms), collapse = " + ")
  rhs
}

#' Internal: extract response name from formula
#' @keywords internal
.cde_response <- function(formula) {
  all.vars(formula)[1]
}
