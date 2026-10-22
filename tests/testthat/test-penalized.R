test_that("scm_penalized returns simplex weights", {
  set.seed(1)
  X0 <- matrix(runif(15), nrow = 3, dimnames = list(NULL, paste0("d", 1:5)))
  X1 <- X0 %*% c(0.3, 0.3, 0.4, 0, 0)
  w <- scm_penalized(X1, X0, lambda = 0.1)
  expect_named(w, paste0("d", 1:5))
  expect_true(all(w >= 0))
  expect_equal(sum(w), 1, tolerance = 1e-8)
})

test_that("scm_penalized with lambda = 0 matches scm_base", {
  set.seed(2)
  X0 <- matrix(runif(20), nrow = 4)
  X1 <- X0 %*% c(0.5, 0.2, 0.3, 0, 0)
  expect_equal(as.numeric(scm_penalized(X1, X0, lambda = 0,
                                        control = list(maxit = 2000, reltol = 1e-12))),
               as.numeric(scm_base(X1, X0)), tolerance = 1e-6)
})

test_that("the penalty selects the nearby donors among exact solutions", {
  X0 <- cbind(d1 = c(-1, 0), d2 = c(1, 0), d3 = c(-3, 0), d4 = c(3, 0))
  w <- scm_penalized(c(0, 0), X0, lambda = 0.1, standardize = FALSE)
  expect_gt(w[["d1"]] + w[["d2"]], 0.99)
})

test_that("scm_penalized rejects an invalid lambda", {
  X0 <- diag(2)
  expect_error(scm_penalized(c(0.5, 0.5), X0, lambda = -1), "non-negative")
})

test_that("scm_penalized works as an iscm backend and passes lambda", {
  sim <- simulate_iscm(n_affected = 1, n_pure = 6, seed = 7)
  f0 <- iscm(sim$data, "unit", "time", "Y", sim$treated, sim$affected, sim$t0,
             estimator = scm_penalized, lambda = 0)
  f1 <- iscm(sim$data, "unit", "time", "Y", sim$treated, sim$affected, sim$t0,
             estimator = scm_penalized, lambda = 5)
  expect_s3_class(f1, "iscm")
  expect_false(isTRUE(all.equal(f0$fits$TR$w, f1$fits$TR$w)))
})
