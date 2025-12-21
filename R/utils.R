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
.trim <- function(x) gsub("^\\s+|\\s+$", "", x)

.has_bar <- function(x) grepl("\\|", x, fixed = FALSE)



# Parse a "~ <re> | <group>" string into list(re=..., group=...)
.parse_bar <- function(s) {
  s <- gsub("\\s+", " ", s)
  s <- sub("^\\s*~\\s*", "", s)
  parts <- strsplit(s, "\\|", fixed = FALSE)[[1]]
  if (length(parts) != 2) stop("Expected a single '|' in random formula.", call. = FALSE)
  list(re = .trim(parts[1]), group = .trim(parts[2]))
}

# Extract RHS from a one-sided formula "~ 1 + x"
.formula_rhs <- function(f) {
  if (!inherits(f, "formula")) stop("Not a formula.", call. = FALSE)
  # one-sided: f[[2]] exists; deparse safely
  rhs <- paste(deparse(f[[2]]), collapse = " ")
  rhs <- gsub("\\s+", " ", rhs)
  rhs <- .trim(rhs)
  if (rhs == "") rhs <- "1"
  rhs
}

# Convert a single lme-style random spec into ONE lme4 random term string "(... | group)" or "(... || group)"
# spec can be: formula, or call like pdDiag(~...) (optionally with "| group" embedded)
.lme_random_one_to_lme4 <- function(spec, group_name = NULL) {

  # Case A: formula
  if (inherits(spec, "formula")) {
    s <- paste(deparse(spec), collapse = " ")
    s <- gsub("\\s+", " ", s)

    if (.has_bar(s)) {
      parsed <- .parse_bar(s)
      re    <- parsed$re
      group <- parsed$group
      if (!is.null(group_name) && nzchar(group_name) && group_name != group) {
        # Prefer list name (what lme uses for list format); still accept mismatch.
        group <- group_name
      }
      if (re == "") re <- "1"
      return(paste0("(", re, " | ", group, ")"))
    } else {
      # lme list-format: group comes from list name; spec is "~ 1 + x"
      if (is.null(group_name) || !nzchar(group_name)) {
        stop("Random formula has no '| group'. Provide a grouping factor via list name (e.g., list(id=~1+x)) or use ~ ... | id.",
             call. = FALSE)
      }
      re <- .formula_rhs(spec)
      return(paste0("(", re, " | ", group_name, ")"))
    }
  }

  # Case B: pdMat call, e.g. pdDiag(~ 1 + Time) OR pdDiag(~ 1 + Time) | id
  if (is.call(spec)) {
    s <- paste(deparse(spec), collapse = " ")
    s <- gsub("\\s+", " ", s)

    # If embedded "| group" exists, split it
    if (.has_bar(s)) {
      parts <- strsplit(s, "\\|", fixed = FALSE)[[1]]
      left  <- .trim(parts[1])
      right <- .trim(parts[2])
      group <- if (!is.null(group_name) && nzchar(group_name)) group_name else right

      # detect pdDiag => independent components => use '||'
      use_double_bar <- grepl("^pdDiag\\s*\\(", left)

      # pull inside (...) and then inside "~ ..."
      inside <- sub("^\\w+\\s*\\((.*)\\)\\s*$", "\\1", left)   # gets "~ 1 + Time" for pdDiag/pdSymm
      inside <- .trim(inside)
      inside <- sub("^~\\s*", "", inside)
      re <- if (inside == "") "1" else inside

      bar <- if (use_double_bar) "||" else "|"
      return(paste0("(", re, " ", bar, " ", group, ")"))
    }

    # No embedded group: must come from list name
    if (is.null(group_name) || !nzchar(group_name)) {
      stop("pdMat random specification provided without a grouping factor. Use pdDiag(~... ) | id or list(id = pdDiag(~...)).",
           call. = FALSE)
    }

    use_double_bar <- grepl("^pdDiag\\s*\\(", s)
    inside <- sub("^\\w+\\s*\\((.*)\\)\\s*$", "\\1", s)
    inside <- .trim(inside)
    inside <- sub("^~\\s*", "", inside)
    re <- if (inside == "") "1" else inside

    bar <- if (use_double_bar) "||" else "|"
    return(paste0("(", re, " ", bar, " ", group_name, ")"))
  }

  stop("Unsupported `random` specification. Supply a formula, a pdMat call, or a list of these.", call. = FALSE)
}

# Public: convert `random=` (as used in lme) to lme4 random terms
lme_random_to_lme4_terms <- function(random) {
  terms <- character(0)

  if (inherits(random, "formula") || is.call(random)) {
    terms <- c(terms, .lme_random_one_to_lme4(random, group_name = NULL))
    return(terms)
  }

  if (is.list(random)) {
    nms <- names(random)
    for (k in seq_along(random)) {
      gname <- NULL
      if (!is.null(nms) && nzchar(nms[k])) gname <- nms[k]
      terms <- c(terms, .lme_random_one_to_lme4(random[[k]], group_name = gname))
    }
    return(terms)
  }

  stop("`random` must be a formula/call or a list (as accepted by nlme::lme).", call. = FALSE)
}


#' Internal helper: convert lme format for random effect to lmer format
#' @keywords internal
# Public: build lmer() formula and a runnable call string
lme_to_lmer <- function(fixed, random) {
  if (!inherits(fixed, "formula")) stop("`fixed` must be a formula.", call. = FALSE)

  re_terms <- lme_random_to_lme4_terms(random)
  fixed_str <- paste(deparse(fixed), collapse = " ")
  rhs_add <- paste(re_terms, collapse = " + ")

  out_formula <- as.formula(paste0(fixed_str, " + ", rhs_add))

  

  list(formula = out_formula, re_terms = re_terms)
}



# ---- helpers ----

.trim <- function(x) gsub("^\\s+|\\s+$", "", x)

.parse_rand <- function(f) {
  if (!inherits(f, "formula")) stop("Random effect must be a formula.", call. = FALSE)
  s <- gsub("\\s+", " ", paste(deparse(f), collapse = " "))
  s <- sub("^\\s*~\\s*", "", s)
  parts <- strsplit(s, "\\|", fixed = FALSE)[[1]]
  if (length(parts) != 2) stop("Need random formula of the form `~ ... | group`.", call. = FALSE)
  list(lhs = .trim(parts[1]), group = .trim(parts[2]))
}

# Extract: (a) whether intercept is present, (b) whether intercept removed, (c) non-intercept terms
.canon_lhs <- function(lhs) {
  lhs <- gsub("\\s+", " ", .trim(lhs))

  removes_int <- grepl("(^|\\s)0(\\s|$)", lhs) || grepl("-\\s*1", lhs)
  has_explicit_1 <- grepl("(^|\\s)1(\\s|$)", lhs)

  # If intercept is not explicitly removed, it is present by default
  has_int <- if (removes_int) FALSE else TRUE
  if (has_explicit_1) has_int <- TRUE

  # Normalize "- 1" into "+ -1" so we can split on '+'
  lhs2 <- gsub("-\\s*1", "+ -1", lhs)
  pieces <- unlist(strsplit(lhs2, "\\+"))
  pieces <- .trim(pieces)
  pieces <- pieces[pieces != ""]
  pieces <- pieces[!pieces %in% c("1", "0", "-1", "- 1")]
  if (length(pieces) > 1) pieces <- pieces[!duplicated(pieces)]

  list(has_int = has_int, removes_int = removes_int, terms = pieces)
}

.combine_two_same_group <- function(r1, r2, keep_intercept = TRUE) {
  a <- .parse_rand(r1)
  b <- .parse_rand(r2)
  if (a$group != b$group) stop("Grouping factors differ; cannot combine.", call. = FALSE)

  ca <- .canon_lhs(a$lhs)
  cb <- .canon_lhs(b$lhs)

  if (keep_intercept) {
    int_keep <- ca$has_int || cb$has_int
  } else {
    int_keep <- !(ca$removes_int || cb$removes_int)
  }

  terms <- unique(c(ca$terms, cb$terms))

  lhs_new <- if (int_keep) {
    if (length(terms) == 0) "1" else paste(c("1", terms), collapse = " + ")
  } else {
    if (length(terms) == 0) "0" else paste(c("0", terms), collapse = " + ")
  }

  as.formula(paste0("~ ", lhs_new, " | ", a$group))
}

# Main: r_list can be a formula OR a list of formulas
merge_random_list_by_group <- function(r_list, r_new, keep_intercept = TRUE) {
  if (inherits(r_list, "formula")) r_list <- list(r_list)
  if (is.null(r_list)) r_list <- list()
  if (!is.list(r_list)) stop("`r_list` must be a list or a formula.", call. = FALSE)
  if (!inherits(r_new, "formula")) stop("`r_new` must be a formula.", call. = FALSE)

  target_group <- .parse_rand(r_new)$group
  same_idx <- which(vapply(r_list, function(x) .parse_rand(x)$group == target_group, logical(1)))

  if (length(same_idx) == 0) {
    r_list[[length(r_list) + 1]] <- r_new
    return(r_list)
  }

  combined <- r_new
  for (i in same_idx) {
    combined <- .combine_two_same_group(r_list[[i]], combined, keep_intercept = keep_intercept)
  }

  keep <- setdiff(seq_along(r_list), same_idx)
  out <- r_list[keep]

  insert_pos <- min(same_idx)
  out <- append(out, list(combined), after = insert_pos - 1)
  out
}


