test_that("simulate_iscm_covariates returns a balanced panel with covariates", {
  sim <- simulate_iscm_covariates(seed = 1)
  expect_equal(nrow(sim$data), 10 * 28)
  expect_true(all(c("unit", "time", "Y", "x1", "x2") %in% names(sim$data)))
  expect_equal(sim$treated, "TR")
  expect_equal(sim$affected, "A1")
  expect_length(sim$pure, 8)
  expect_equal(unname(sim$u[c("TR", "A1", "C3")]), c(1, 1, 0))
})

test_that("in the covariate design iSCM beats the restricted SCM", {
  sim <- simulate_iscm_covariates(seed = 3)
  args <- list(sim$data, "unit", "time", "Y", "TR", "A1", sim$t0,
               predictors = c("x1", "x2"), estimator = scm_penalized)
  fit  <- do.call(iscm, args)
  rest <- do.call(restricted_scm, args)
  post <- fit$times > sim$t0
  err_iscm <- abs(mean(fit$effects["TR", ]) - 2)
  err_rest <- abs(mean(rest$gap[post]) - 2)
  expect_lt(err_iscm, 0.3)
  expect_gt(err_rest, 1)
  expect_equal(mean(fit$effects["A1", ]), 3, tolerance = 0.15)

  cmp <- do.call(iscm_compare, args)
  expect_equal(cmp$recommendation, "iscm_full")
})
