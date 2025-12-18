test_that("cde_augment_data creates expected columns", {
  skip_if_not_installed("nlme")
  dat <- data.frame(
    id = rep(1:3, each = 3),
    Time = rep(1:3, times = 3),
    y = rnorm(9),
    x1 = c(1, NA, 3, 4, 5, 6, 7, 8, 9),
    x2 = c(NA, 2, 3, 4, 5, NA, 7, 8, 9)
  )

  aug <- cde_augment_data(
    data = dat,
    id = "id",
    time = "Time",
    baseline_missing = c("x1"),
    time_dependent_missing = c("x2")
  )

  expect_true(all(c("mis_any", "mis_base_any", "mis_td_any", "ct") %in% names(aug)))
  expect_true(any(grepl("^u_x2_t", names(aug))))
  expect_true(all(!is.na(aug$x1))) # replaced NA with 0
  expect_true(all(!is.na(aug$x2)))
})

test_that("cde_fit runs for gaussian with nlme engine", {
  skip_if_not_installed("nlme")

  set.seed(1)
  dat <- data.frame(
    id = rep(1:10, each = 3),
    Time = rep(1:3, times = 10)
  )
  dat$x1 <- rnorm(nrow(dat))
  dat$x2 <- rnorm(nrow(dat))
  dat$y <- 1 + 0.5 * dat$x1 + 0.2 * dat$Time + rnorm(nrow(dat))

  # introduce missing
  dat$x1[sample.int(nrow(dat), 5)] <- NA
  dat$x2[sample.int(nrow(dat), 5)] <- NA

  fit <- cde_fit(
    y ~ x1 + Time + x2,
    data = dat,
    id = "id",
    time = "Time",
    missing = list(baseline = "x1", time_dependent = "x2"),
    family = gaussian(),
    engine = "nlme",
    control = list(msMaxIter = 50)
  )

  expect_s3_class(fit, "cde_fit")
  expect_true(!is.null(cde_fixef(fit)))
})
