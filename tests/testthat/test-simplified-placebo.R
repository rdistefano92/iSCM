test_that("simplified iSCM gives a triangular Omega and no inversion instability", {
  sim <- simulate_iscm(seed = 7, spill = 3)
  fit <- iscm(sim$data, "unit", "time", "Y", sim$treated, sim$affected,
              sim$t0, simplified = TRUE)
  # affected-unit rows of Omega are unit vectors
  expect_equal(fit$Omega[2, 1], 0)
  expect_equal(mean(fit$effects[1, ]), 2, tolerance = 0.3)
})

test_that("full and simplified iSCM agree when l1 is small", {
  sim <- simulate_iscm(seed = 7, spill = 3)
  f_full <- iscm(sim$data, "unit", "time", "Y", sim$treated, sim$affected, sim$t0)
  f_simp <- iscm(sim$data, "unit", "time", "Y", sim$treated, sim$affected,
                 sim$t0, simplified = TRUE)
  expect_lt(abs(mean(f_full$effects[1, ]) - mean(f_simp$effects[1, ])), 0.3)
})

test_that("placebo ranks the treated unit as the most extreme", {
  sim <- simulate_iscm(seed = 7, spill = 3, n_pure = 6)
  fit <- iscm(sim$data, "unit", "time", "Y", sim$treated, sim$affected, sim$t0)
  pb <- iscm_placebo(fit, sim$data, "unit", "time", "Y", target = sim$treated)
  expect_s3_class(pb, "iscm_placebo")
  expect_equal(pb$rank, 1L)          # largest post/pre RMSPE ratio
  expect_true(pb$p_value <= 1 / (length(sim$pure)))
})

test_that("placebo p-value counts ties conservatively", {
  ratios <- c(target = 2, p1 = 2, p2 = 1)
  r_target <- ratios[["target"]]
  placebo_ratios <- ratios[names(ratios) != "target"]
  p_value <- (1 + sum(placebo_ratios >= r_target)) /
    (length(placebo_ratios) + 1)
  expect_equal(p_value, 2 / 3)
})

test_that("placebo reuses the fitted backend specification", {
  sim <- simulate_iscm(seed = 7, spill = 3, n_pure = 3)
  seen <- new.env(parent = emptyenv())
  backend <- function(X1, X0, v = NULL, marker = NULL, ...) {
    seen$marker <- marker
    seen$v <- v
    rep(1 / ncol(X0), ncol(X0))
  }
  fit <- iscm(sim$data, "unit", "time", "Y", sim$treated, sim$affected,
              sim$t0, estimator = backend, v = 2, marker = "original")
  seen$marker <- NULL
  seen$v <- NULL
  iscm_placebo(fit, sim$data, "unit", "time", "Y")
  expect_identical(seen$marker, "original")
  expect_equal(seen$v, 2)
})

test_that("iscm object has a working print and summary", {
  sim <- simulate_iscm(seed = 7)
  fit <- iscm(sim$data, "unit", "time", "Y", sim$treated, sim$affected, sim$t0)
  expect_output(print(fit), "Inclusive synthetic control")
  s <- summary(fit)
  expect_s3_class(s, "summary.iscm")
  expect_equal(nrow(s$effects), length(sim$tau))
})
