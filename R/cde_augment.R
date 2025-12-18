\
#' Augment long-format data with CDE-style missingness features
#'
#' This function generalizes the core augmentation used in your simulation scripts:
#' it identifies timepoints with missingness, creates timepoint dummies `u_*`,
#' creates missingness indicators for baseline and time-dependent covariates, and
#' replaces missing covariate values with 0 (leaving the response unchanged).
#'
#' @param data A long-format `data.frame`.
#' @param id Column name for subject ID.
#' @param time Column name for time variable.
#' @param baseline_missing Baseline covariates that may be missing (character vector).
#' @param time_dependent_missing Time-dependent covariates that may be missing (character vector).
#' @param covariates Optional: the full covariate set to consider for `mis_any`. If `NULL`,
#'   it is taken as the union of `baseline_missing` and `time_dependent_missing`.
#' @param replace_na_with Value to replace `NA` for covariates (default 0).
#' @return A `data.frame` with additional columns:
#'   `mis_any`, `mis_base_any`, `mis_td_any`, per-variable flags, timepoint dummies `u_*`,
#'   and `ct` (row index).
#' @export
cde_augment_data <- function(
  data,
  id,
  time,
  baseline_missing = character(0),
  time_dependent_missing = character(0),
  covariates = NULL,
  replace_na_with = 0
) {
  id <- .as_chr1(id, "id")
  time <- .as_chr1(time, "time")
  baseline_missing <- .as_chrv(baseline_missing, "baseline_missing")
  time_dependent_missing <- .as_chrv(time_dependent_missing, "time_dependent_missing")

  if (!is.data.frame(data)) stop("`data` must be a data.frame.", call. = FALSE)
  if (!all(c(id, time) %in% names(data))) {
    stop("`id` and/or `time` not found in `data` column names.", call. = FALSE)
  }

  dat <- data
  dat[[id]] <- as.factor(dat[[id]])

  covs <- covariates %||% unique(c(baseline_missing, time_dependent_missing))
  covs <- covs[covs %in% names(dat)]

  # mis_any across selected covariates
  if (length(covs) > 0) {
    mis_mat <- is.na(dat[, covs, drop = FALSE])
    dat$mis_any <- as.integer(rowSums(mis_mat) > 0)
  } else {
    dat$mis_any <- 0L
  }

  # baseline flags
  dat$mis_base_any <- 0L
  for (v in baseline_missing) {
    if (!v %in% names(dat)) next
    nm <- paste0("mis_base_", v)
    dat[[nm]] <- as.integer(is.na(dat[[v]]))
    dat$mis_base_any <- pmax(dat$mis_base_any, dat[[nm]])
  }

  # time-dependent flags + timepoint dummies
  dat$mis_td_any <- 0L
  for (v in time_dependent_missing) {
    if (!v %in% names(dat)) next
    nm <- paste0("mis_td_", v)
    mis_v <- is.na(dat[[v]])
    dat[[nm]] <- as.integer(mis_v)
    dat$mis_td_any <- pmax(dat$mis_td_any, dat[[nm]])

    # Determine which timepoints have missingness for this variable
    t_miss <- sort(unique(dat[[time]][mis_v]))
    if (length(t_miss) > 0) {
      for (tt in t_miss) {
        u_nm <- paste0("u_", v, "_t", tt)
        dat[[u_nm]] <- as.integer(mis_v & (dat[[time]] == tt))
      }
    }
  }

  # Replace NA with 0 (or replace_na_with) ONLY for covariates listed in covs
  for (v in covs) {
    if (!v %in% names(dat)) next
    idx <- is.na(dat[[v]])
    if (any(idx)) dat[[v]][idx] <- replace_na_with
  }

  dat$ct <- seq_len(nrow(dat))
  dat
}
