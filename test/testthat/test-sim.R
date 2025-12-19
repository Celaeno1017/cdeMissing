test_that("sim_cde_data returns format compatible with cde_fit", {
  skip_if_not_installed("nlme")

  sim <- sim_cde_data(model = "lmm", m = 20, t = 4, q = 0.2, seed = 123)
  dat <- sim$data_mis

  expect_true(is.data.frame(dat))
  expect_true(all(c("y", "id", "x1", "Treatment", "Time", "x2") %in% names(dat)))

  fit <- cde_fit(
    y ~ x1 + Treatment + Time + x2,
    data = dat,
    id = "id",
    time = "Time",
    missing = list(baseline = "x1", time_dependent = "x2"),
    family = gaussian(),
    engine = "nlme",
    control = list(msMaxIter = 50)
  )
  expect_s3_class(fit, "cde_fit")
})
