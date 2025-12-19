#' Fit CDE-style LMM/GLMM with missing covariates
#'
#' This is the main user-facing function for **cdeMissing**.
#' It (1) augments the input long-format dataset with missingness indicators and
#' timepoint dummies, then (2) fits an LMM (via `nlme::lme`) for Gaussian outcomes or
#' a GLMM (via `lme4::glmer`) for non-Gaussian outcomes.
#'
#' By default, it mirrors the structure used in your simulation scripts:
#' - LMM: `random = ~ 1 + mis_base_any + mis_td_any | id` and
#'   `weights = varIdent(form = ~ 1 | mis_base_any * mis_td_any)`.
#' - GLMM: random terms `(mis_any - 1 | id) + (mis_td_any - 1 | ct)` and `nAGQ = 0`.
#'
#' @param formula A model formula for the *substantive* fixed effects (e.g., `y ~ x1 + Treatment + Time + x2`).
#'   Missingness terms will be appended automatically.
#' @param data Long-format `data.frame`.
#' @param id Column name for subject ID.
#' @param time Column name for time variable.
#' @param missing A list with elements `baseline` and `time_dependent`, each a character vector.
#' @param family A GLM family. Use `gaussian()` for LMM.
#' @param engine One of `"auto"`, `"nlme"`, `"lme4"`.
#' @param control A list of control arguments. For LMM, passed to `nlme::lmeControl()`.
#'   For GLMM, currently supports `nAGQ` (default 0).
#' @param return_augmented_data Whether to store the augmented data in the returned object.
#'
#' @return An object of class `cde_fit` with fields:
#'   - `fit`: the underlying `lme` or `glmerMod` object
#'   - `engine`, `family`
#'   - `data_aug` (optional)
#'   - `call`
#' @export
cde_fit <- function(
  formula,
  data,
  id,
  time,
  missing = list(baseline = character(0), time_dependent = character(0)),
  family = stats::gaussian(),
  engine = c("auto", "nlme", "lme4"),
  control = list(),
  return_augmented_data = TRUE
) {
  engine <- match.arg(engine)
  id <- .as_chr1(id, "id")
  time <- .as_chr1(time, "time")

  baseline_missing <- .as_chrv(missing$baseline %||% character(0), "missing$baseline")
  td_missing <- .as_chrv(missing$time_dependent %||% character(0), "missing$time_dependent")

  if (engine == "auto") {
    engine <- if (identical(family$family, "gaussian")) "nlme" else "lme4"
  }

  dat_aug <- cde_augment_data(
    data = data,
    id = id,
    time = time,
    baseline_missing = baseline_missing,
    time_dependent_missing = td_missing
  )

  y_name <- .cde_response(formula)
  rhs <- .cde_fixed_rhs(formula, dat_aug)

  if (engine == "nlme") {
    form_fix <- stats::as.formula(paste(y_name, "~", rhs))

    # Default mirrors your LMM CDE: random slopes for missingness + varIdent by missingness pattern
    fit <- nlme::lme(
      fixed = form_fix,
      data  = dat_aug,
      random = stats::as.formula(paste("~ 1 + mis_any |", id)),
      weights = nlme::varIdent(form = ~ 1 | mis_td_any),
      control = nlme::lmeControl(opt = "nlminb")),silent = TRUE)
    )

    out <- list(
      fit = fit,
      engine = "nlme",
      family = family,
      data_aug = if (isTRUE(return_augmented_data)) dat_aug else NULL,
      call = match.call()
    )
    class(out) <- "cde_fit"
    return(out)
  }

  if (engine == "lme4") {
    # Default mirrors your GLMM CDE: (mis_id-1|id) + (mis_t_id-1|ct)
    # Here we generalize to mis_any and mis_td_any.
    glmm_rhs <- paste0(
      rhs,
      " + (mis_any - 1 | ", id, ")",
      " + (mis_td_any - 1 | ct)"
    )
    form_glmm <- stats::as.formula(paste(y_name, "~", glmm_rhs))

    nAGQ <- control$nAGQ %||% 0
    fit <- lme4::glmer(
      form_glmm,
      data = dat_aug,
      family = family,
      nAGQ = nAGQ
    )

    out <- list(
      fit = fit,
      engine = "lme4",
      family = family,
      data_aug = if (isTRUE(return_augmented_data)) dat_aug else NULL,
      call = match.call()
    )
    class(out) <- "cde_fit"
    return(out)
  }

  stop("Unknown `engine`.", call. = FALSE)
}

#' Extract fixed effects from a `cde_fit`
#' @param object A `cde_fit`.
#' @param ... Unused.
#' @return A named numeric vector of fixed effect estimates.
#' @export
cde_fixef <- function(object, ...) {
  stopifnot(inherits(object, "cde_fit"))
  if (object$engine == "nlme") return(nlme::fixef(object$fit))
  if (object$engine == "lme4") return(lme4::fixef(object$fit))
  stop("Unknown engine.", call. = FALSE)
}

#' Retrieve augmented data (if stored)
#' @param object A `cde_fit`.
#' @export
cde_augmented_data <- function(object) {
  stopifnot(inherits(object, "cde_fit"))
  object$data_aug
}
