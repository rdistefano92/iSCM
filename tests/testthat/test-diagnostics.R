test_that("omega_diagnostics matches the paper's Germany numbers", {
  Omega <- matrix(c(1, -0.42, -0.33, 1), 2, 2, byrow = TRUE)
  d <- omega_diagnostics(Omega)
  expect_equal(d$det, 1 - 0.42 * 0.33, tolerance = 1e-12)      # 0.8614
  expect_equal(d$alpha, 1 / (1 - 0.42 * 0.33), tolerance = 1e-12)  # ~1.16
  expect_true(d$diag_dominant)
  expect_equal(unname(d$sigma_min), 0.63, tolerance = 0.02)
})

test_that("amplification equals 1 / sigma_min", {
  Omega <- matrix(c(1, -0.42, -0.33, 1), 2, 2, byrow = TRUE)
  d <- omega_diagnostics(Omega)
  expect_equal(d$amplification, 1 / d$sigma_min, tolerance = 1e-12)
})

test_that("a diagonally dominant system can still be poorly conditioned", {
  # off-diagonal 0.9 < 1, so Omega IS strictly diagonally dominant, yet the
  # amplification alpha = 1/(1 - 0.9*0.9) ~ 5.26 is large: dominance and
  # conditioning are distinct.
  Omega <- matrix(c(1, -0.9, -0.9, 1), 2, 2, byrow = TRUE)
  d <- omega_diagnostics(Omega)
  expect_true(d$diag_dominant)
  expect_gt(d$alpha, 5)
})

test_that("diagonal dominance flags a genuinely non-dominant system", {
  # each row's off-diagonal absolute sum is 1.2 > 1, so Omega is NOT strictly
  # diagonally dominant; sigma_min = 0.2 gives amplification 1/sigma_min = 5.
  Omega <- matrix(c(1, -0.6, -0.6,
                    -0.6, 1, -0.6,
                    -0.6, -0.6, 1), 3, 3, byrow = TRUE)
  d <- omega_diagnostics(Omega)
  expect_false(d$diag_dominant)
  expect_true(is.na(d$alpha))          # alpha is defined only for m = 2
  expect_gt(d$amplification, 3)
})

test_that("rmspe computes the root mean squared prediction error", {
  expect_equal(rmspe(c(3, 4)), sqrt((9 + 16) / 2))
})
