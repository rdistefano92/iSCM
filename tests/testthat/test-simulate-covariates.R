test_that("simulate_iscm_covariates returns a balanced panel with covariates", {
  sim <- simulate_iscm_covariates(seed = 1)
  expect_equal(nrow(sim$data), 10 * 28)
  expect_true(all(c("unit", "time", "Y", "x1", "x2") %in% names(sim$data)))
  expect_equal(sim$treated, "TR")
  expect_equal(sim$affected, "A1")
  expect_length(sim$pure, 8)
  expect_equal(unname(sim$u[c("TR", "A1", "C2", "C3")]), c(1, 1, 1, 0))
})

test_that("in the covariate design iSCM is accurate and the rule selects it", {
  sim <- simulate_iscm_covariates(seed = 4)
  args <- list(sim$data, "unit", "time", "Y", "TR", "A1", sim$t0,
               predictors = c("x1", "x2"), estimator = scm_penalized)
  fit <- do.call(iscm, args)
  expect_equal(mean(fit$effects["TR", ]), 2, tolerance = 0.1)
  expect_equal(mean(fit$effects["A1", ]), 3, tolerance = 0.1)
  expect_lt(fit$diagnostics$alpha, 1.25)               # similar to the German case
  cmp <- do.call(iscm_compare, args)
  expect_equal(cmp$recommendation, "iscm_full")
})

test_that("over replications iSCM beats the restricted SCM", {
  err <- t(vapply(1:10, function(s) {
    sim <- simulate_iscm_covariates(seed = s)
    args <- list(sim$data, "unit", "time", "Y", "TR", "A1", sim$t0,
                 predictors = c("x1", "x2"), estimator = scm_penalized)
    fit  <- do.call(iscm, args)
    rest <- do.call(restricted_scm, args)
    post <- fit$times > sim$t0
    abs(c(iscm = mean(fit$effects["TR", ]), restricted = mean(rest$gap[post])) - 2)
  }, numeric(2)))
  expect_lt(mean(err[, "iscm"]), 0.15)
  expect_gt(mean(err[, "restricted"]), 4 * mean(err[, "iscm"]))
})
