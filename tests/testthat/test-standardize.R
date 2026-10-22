test_that("iscm_compare standardises Step 2 imbalance by predictor SD", {
  set.seed(1)
  units <- c("TR", "A", "C1", "C2", "C3", "C4"); times <- 1:12; t0 <- 8
  rows <- expand.grid(unit = units, time = times, stringsAsFactors = FALSE)
  base <- setNames(c(100, 98, 101, 99, 102, 97), units)
  rows$Y <- base[rows$unit] + rnorm(nrow(rows), 0, 1) +
    ifelse(rows$unit == "TR" & rows$time > t0, 5, 0)
  # TR is extreme on BOTH predictors, so its synthetic control cannot reproduce
  # it exactly and the Step 2 imbalance is strictly positive.
  pbig   <- setNames(c(16000, 15000, 14000, 13500, 14200, 13800), units)  # ~1e4
  psmall <- setNames(c(7.0, 5.0, 3.0, 2.5, 3.2, 2.8), units)              # ~1
  rows$pbig   <- pbig[rows$unit]   + rnorm(nrow(rows), 0, 20)
  rows$psmall <- psmall[rows$unit] + rnorm(nrow(rows), 0, 0.05)
  rows$pbig_k <- rows$pbig / 1000                                          # rescaled

  imb <- function(pred, std) {
    cmp <- iscm_compare(rows, "unit", "time", "Y", "TR", "A", t0,
                        predictors = pred, standardize = std)
    m <- cmp$comparison
    m[m$unit == "TR", c("imbalance_unres", "imbalance_restr")]
  }
  std   <- imb(c("pbig", "psmall"),   TRUE)
  raw   <- imb(c("pbig", "psmall"),   FALSE)
  std_k <- imb(c("pbig_k", "psmall"), TRUE)
  raw_k <- imb(c("pbig_k", "psmall"), FALSE)

  # standardised imbalance is positive and differs from the raw sum
  expect_gt(std$imbalance_unres, 0.05)
  expect_false(isTRUE(all.equal(std$imbalance_unres, raw$imbalance_unres)))

  # standardising makes the imbalance invariant to the units of a predictor ...
  expect_equal(std$imbalance_unres, std_k$imbalance_unres, tolerance = 1e-4)
  expect_equal(std$imbalance_restr, std_k$imbalance_restr, tolerance = 1e-4)

  # ... whereas the raw sum is not (it is dominated by the large-scale predictor)
  expect_false(isTRUE(all.equal(raw$imbalance_unres, raw_k$imbalance_unres,
                                tolerance = 1e-4)))
})

test_that("standardize = TRUE is the default in iscm_compare", {
  expect_identical(formals(iscm_compare)$standardize, TRUE)
})
