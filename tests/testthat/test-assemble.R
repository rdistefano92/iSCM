test_that("iscm_assemble matches the data-level iscm() on a simulated panel", {
  sim <- simulate_iscm(seed = 7, spill = 3)
  fit <- iscm(sim$data, "unit", "time", "Y", sim$treated, sim$affected, sim$t0)

  gaps <- list(fit$gaps_unrestricted[, sim$treated],
               fit$gaps_unrestricted[, sim$affected[1]])
  names(gaps) <- c(sim$treated, sim$affected)
  weights <- list(fit$fits[[1]]$w, fit$fits[[2]]$w)
  names(weights) <- c(sim$treated, sim$affected)

  fit2 <- iscm_assemble(gaps, weights, t0 = sim$t0, times = fit$times)
  expect_equal(unname(fit2$Omega), unname(fit$Omega), tolerance = 1e-10)
  expect_equal(rowMeans(fit2$effects), rowMeans(fit$effects), tolerance = 1e-10)
})

test_that("iscm_assemble reproduces the m = 2 closed form", {
  set.seed(2)
  gap_tr <- c(rnorm(10, 0, 0.1), rnorm(4, 1.8, 0.1))
  gap_a  <- c(rnorm(10, 0, 0.1), rnorm(4, 1.0, 0.1))
  w2 <- 0.42; l1 <- 0.33
  fit <- iscm_assemble(
    gaps    = list(TR = gap_tr, A = gap_a),
    weights = list(TR = c(A = w2, C1 = 1 - w2),
                   A  = c(TR = l1, C1 = 1 - l1)),
    t0 = 10, times = 1:14)
  post <- 11:14
  den <- 1 - w2 * l1
  theta_cf <- (gap_tr[post] + w2 * gap_a[post]) / den
  expect_equal(as.numeric(fit$effects[1, ]), theta_cf, tolerance = 1e-10)
})

test_that("germany dataset is available and well formed", {
  data(germany, package = "iscm", envir = environment())
  expect_true(is.data.frame(germany))
  expect_equal(nrow(germany), 748L)
  expect_true(all(c("index","country","year","gdp") %in% names(germany)))
  expect_true("West Germany" %in% germany$country)
})
