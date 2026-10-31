test_that("fit_ar2 rejects series it cannot fit", {
  expect_error(fit_ar2(1:5), "Need >= 6 timepoints")
  expect_error(fit_ar2(c(1, 2, NA, 4, 5, 6, 7)), "missing values")
})

test_that("fit_ar2 recovers known coefficients", {
  set.seed(4)
  x <- as.numeric(stats::arima.sim(list(ar = c(0.6, -0.3)), n = 4000))
  fit <- fit_ar2(x)
  expect_equal(fit$phi1, 0.6, tolerance = 0.05)
  expect_equal(fit$phi2, -0.3, tolerance = 0.05)
  expect_equal(fit$eigenvalue, sqrt(0.3), tolerance = 0.05)
  expect_identical(fit$root_type, "Complex")
})

test_that("a pure oscillation gives complex roots and the right eigenperiod", {
  t <- seq_len(48)
  x <- sin(2 * pi * t / 24)
  fit <- fit_ar2(x)
  expect_identical(fit$root_type, "Complex")
  expect_equal(fit$eigenperiod, 24, tolerance = 1e-3)
  expect_gt(fit$eigenvalue, 0.99)
  expect_equal(fit$r2, 1, tolerance = 1e-9)
})

test_that("an overdamped series gives real roots and no eigenperiod", {
  x <- 2^(-seq_len(20) / 3)
  fit <- fit_ar2(x)
  expect_identical(fit$root_type, "Real")
  expect_true(is.na(fit$eigenperiod))
})

test_that("centring makes the fit invariant to a shift in level", {
  set.seed(9)
  x <- as.numeric(stats::arima.sim(list(ar = c(0.4, -0.2)), n = 60))
  expect_equal(fit_ar2(x)$eigenvalue, fit_ar2(x + 1000)$eigenvalue)
})

test_that("half_life and eigenperiod agree with the values carried in the fit", {
  set.seed(11)
  x <- as.numeric(stats::arima.sim(list(ar = c(0.7, -0.25)), n = 96))
  fit <- fit_ar2(x)
  expect_equal(fit$half_life, round(half_life(fit$eigenvalue, 1), 4))
  expect_equal(fit$eigenperiod, round(eigenperiod(fit$phi1, fit$phi2, 1), 4))
})

test_that("half_life and eigenperiod handle their boundaries", {
  expect_identical(half_life(0), 0)
  expect_identical(half_life(1), Inf)
  expect_true(is.na(eigenperiod(1.2, -0.2)))
  expect_equal(half_life(0.8, 2), 2 * half_life(0.8, 1))
})

test_that("classify_dynamics separates persistence from oscillation", {
  expect_identical(classify_dynamics(1.02, "Complex"), "Unstable")
  expect_identical(classify_dynamics(0.9, "Complex"), "Sustained oscillator")
  expect_identical(classify_dynamics(0.5, "Complex"), "Damped oscillator")
  expect_identical(classify_dynamics(0.9, "Real"), "Overdamped decay")
  expect_identical(classify_dynamics(0.1, "Real"), "Rapid decay")
})

test_that("bootstrap intervals are reproducible and bracket the estimate", {
  set.seed(3)
  x <- as.numeric(stats::arima.sim(list(ar = c(0.6, -0.3)), n = 24))
  a <- fit_ar2(x, n_bootstrap = 300, seed = 42)
  b <- fit_ar2(x, n_bootstrap = 300, seed = 42)
  expect_identical(a$eigenvalue_ci, b$eigenvalue_ci)
  expect_lt(a$eigenvalue_ci[1L], a$eigenvalue_ci[2L])
  expect_identical(a$n_timepoints, 24L)
  expect_false(identical(
    fit_ar2(x, n_bootstrap = 300, seed = 7)$eigenvalue_ci,
    a$eigenvalue_ci
  ))
})

test_that("without a seed, bootstrapping draws from the caller's stream", {
  set.seed(3)
  x <- as.numeric(stats::arima.sim(list(ar = c(0.5, -0.2)), n = 24))
  set.seed(123)
  a <- fit_ar2(x, n_bootstrap = 50)
  set.seed(123)
  b <- fit_ar2(x, n_bootstrap = 50)
  expect_identical(a$eigenvalue_ci, b$eigenvalue_ci)
  expect_identical(fit_ar2(x, n_bootstrap = 50, seed = 123)$eigenvalue_ci,
                   a$eigenvalue_ci)
})

test_that("print shows the fit without altering it", {
  set.seed(5)
  x <- as.numeric(stats::arima.sim(list(ar = c(0.6, -0.3)), n = 48))
  fit <- fit_ar2(x)
  expect_output(print(fit), "\\|lambda\\|")
  expect_identical(print(fit), fit)
})

test_that("an unstable fit reports its spectral radius without clipping", {
  fit <- fit_ar2(3^(0:11))
  expect_equal(fit$phi1, 4, tolerance = 1e-6)
  expect_equal(fit$phi2, -3, tolerance = 1e-6)
  expect_equal(fit$eigenvalue, 3, tolerance = 1e-6)
  expect_true(is.na(fit$half_life))
})

test_that("timescale converts |lambda| to time units", {
  expect_equal(timescale(exp(-1), sampling_interval = 2), 2)
  expect_equal(timescale(0.5, 1) * log(2), half_life(0.5, 1))
  expect_identical(timescale(0), 0)
  expect_identical(timescale(1), Inf)
})

test_that("fit_ar2 adds timescale only when sampling_interval is given", {
  set.seed(6)
  x <- rnorm(48)
  expect_null(fit_ar2(x)$timescale)
  one <- fit_ar2(x, sampling_interval = 1)
  two <- fit_ar2(x, sampling_interval = 2)
  expect_equal(one$timescale, timescale(one$eigenvalue, 1), tolerance = 1e-4)
  expect_equal(two$timescale, 2 * one$timescale, tolerance = 1e-3)
  expect_error(fit_ar2(x, sampling_interval = 0), "sampling_interval")
})

test_that("unstable fits report no timescale", {
  expect_true(is.na(fit_ar2(3^(0:11), sampling_interval = 2)$timescale))
})

test_that("a precise fit gives an ordered interval without a warning", {
  set.seed(8)
  x <- as.numeric(stats::arima.sim(list(ar = c(0.5, -0.3)), n = 400))
  expect_warning(
    fit <- fit_ar2(x, n_bootstrap = 300, seed = 1, sampling_interval = 1),
    NA
  )
  expect_lt(fit$timescale_ci[1L], fit$timescale)
  expect_gt(fit$timescale_ci[2L], fit$timescale)
  expect_false(fit$timescale_imprecise)
})

test_that("a wide timescale interval is flagged with a warning", {
  set.seed(9)
  x <- cos(2 * pi * (0:11) / 12) + 0.5 * rnorm(12)
  expect_warning(
    fit <- fit_ar2(x, n_bootstrap = 300, seed = 1, sampling_interval = 2),
    "too short or too noisy"
  )
  expect_true(fit$timescale_imprecise)
})

test_that("fit_ar2_batch carries a timescale column on request", {
  set.seed(10)
  m <- matrix(rnorm(72), nrow = 3)
  expect_false("timescale" %in% names(fit_ar2_batch(m)))
  expect_true("timescale" %in% names(fit_ar2_batch(m, sampling_interval = 2)))
})
