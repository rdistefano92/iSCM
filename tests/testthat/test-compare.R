test_that("iscm_decision reproduces the paper's Germany recommendation", {
  d <- iscm_decision(rmspe_ratio = 260.74 / 119.07, alpha = 1.16,
                     diag_dominant = TRUE)
  expect_equal(d$recommendation, "iscm_full")
})

test_that("iscm_decision covers all four branches", {
  expect_equal(iscm_decision(2.3, alpha = 3.6, diag_dominant = FALSE)$recommendation,
               "iscm_simplified")
  expect_equal(iscm_decision(1.05, alpha = 1.1, diag_dominant = TRUE)$recommendation,
               "restricted")
  expect_equal(iscm_decision(2.3, alpha = 1.1, diag_dominant = TRUE,
                             substantial = FALSE)$recommendation, "restricted")
  expect_equal(iscm_decision(2.3, alpha = 1.1, diag_dominant = TRUE)$recommendation,
               "iscm_full")
})

test_that("iscm_decision uses sigma_min when alpha is NA (m > 2)", {
  d <- iscm_decision(rmspe_ratio = 2, sigma_min = 0.1, diag_dominant = TRUE)
  # 1/sigma_min = 10 > alpha_tol -> simplified
  expect_equal(d$recommendation, "iscm_simplified")
})

test_that("iscm_compare returns a full comparison and a recommendation", {
  sim <- simulate_iscm(seed = 7, spill = 3)
  cmp <- iscm_compare(sim$data, "unit", "time", "Y", sim$treated, sim$affected,
                      sim$t0, validation = 8)
  expect_s3_class(cmp, "iscm_compare")
  expect_true(cmp$recommendation %in% c("restricted", "iscm_full", "iscm_simplified"))
  expect_true(all(c("imbalance_unres", "rmspe_unres", "rmspe_val_unres")
                  %in% names(cmp$comparison)))
  expect_equal(nrow(cmp$comparison), 1L + length(sim$affected))
})

test_that("iscm_compare_fits applies the rule to pre-computed gaps", {
  set.seed(3)
  # restricted fits much worse than unrestricted for the treated unit
  gU_tr <- c(rnorm(20, 0, 0.1),  rnorm(6, 1.9, 0.1))
  gR_tr <- c(rnorm(20, 0, 0.5),  rnorm(6, 2.3, 0.1))   # worse pre-fit
  gU_a  <- c(rnorm(20, 0, 0.1),  rnorm(6, 1.0, 0.1))
  gR_a  <- c(rnorm(20, 0, 0.12), rnorm(6, 1.0, 0.1))
  cmp <- iscm_compare_fits(
    gaps_unres = list(TR = gU_tr, A = gU_a),
    gaps_restr = list(TR = gR_tr, A = gR_a),
    weights    = list(TR = c(A = 0.42, C1 = 0.58), A = c(TR = 0.33, C1 = 0.67)),
    t0 = 20, times = 1:26, rule = "alpha")
  expect_s3_class(cmp, "iscm_compare")
  expect_equal(cmp$recommendation, "iscm_full")   # well conditioned, restricted worse
  expect_true(cmp$rmspe_ratio > 1.25)
})

test_that("iscm_compare_fits reproduces the paper's recommendation shape", {
  # feed RMSPE-equivalent gaps consistent with the Germany figures
  cmp <- iscm_compare_fits(
    gaps_unres = list("West Germany" = c(rep(0.119, 30), rep(-2, 14)),
                      "Austria"      = c(rep(0.195, 30), rep(0.3, 14))),
    gaps_restr = list("West Germany" = c(rep(0.271, 30), rep(-3, 14)),
                      "Austria"      = c(rep(0.181, 30), rep(1.3, 14))),
    weights    = list("West Germany" = c("Austria" = 0.42, C = 0.58),
                      "Austria"      = c("West Germany" = 0.33, C = 0.67)),
    t0 = 30, times = 1:44, rule = "alpha")
  expect_equal(cmp$recommendation, "iscm_full")
  expect_equal(unname(round(cmp$alpha, 2)), 1.16)
})

test_that("the implied-RMSPE rule matches the closed forms for m = 2", {
  w2 <- 0.42; l1 <- 0.33; alpha <- 1 / (1 - w2 * l1)
  e1 <- rep(0.119, 30); e2 <- rep(0.195, 30); et2 <- rep(0.181, 30)
  cmp <- iscm_compare_fits(
    gaps_unres = list(TR = c(e1, rep(-2, 14)), A = c(e2, rep(0.3, 14))),
    gaps_restr = list(TR = c(rep(0.271, 30), rep(-3, 14)), A = c(et2, rep(1.3, 14))),
    weights    = list(TR = c(A = w2, C = 1 - w2), A = c(TR = l1, C = 1 - l1)),
    t0 = 30, times = 1:44)
  expect_equal(cmp$implied_rmspe[["full"]], alpha * (0.119 + w2 * 0.195))
  expect_equal(cmp$implied_rmspe[["simplified"]], 0.119 + w2 * 0.181)
  expect_equal(cmp$implied_rmspe[["restricted"]], 0.271)
  # same-sign residuals add up, so the amplification costs more than it gains
  expect_equal(cmp$recommendation, "iscm_simplified")
})

test_that("the implied-RMSPE rule prefers the full iSCM when residuals offset", {
  e1 <- rep(c(0.1, -0.1), 15); e2 <- -e1
  cmp <- iscm_compare_fits(
    gaps_unres = list(TR = c(e1, rep(-2, 14)), A = c(e2, rep(0.3, 14))),
    gaps_restr = list(TR = c(rep(0.5, 30), rep(-3, 14)),
                      A = c(rep(0.3, 30), rep(1.3, 14))),
    weights    = list(TR = c(A = 0.42, C = 0.58), A = c(TR = 0.33, C = 0.67)),
    t0 = 30, times = 1:44)
  expect_lt(cmp$implied_rmspe[["full"]], cmp$implied_rmspe[["simplified"]])
  expect_equal(cmp$recommendation, "iscm_full")
})

test_that("the implied full-iSCM gap equals the gap of the implied weights", {
  sim <- simulate_iscm(n_affected = 2, n_pure = 6, spill = c(3, 1), seed = 5)
  fit <- iscm(sim$data, "unit", "time", "Y", sim$treated, sim$affected, sim$t0)
  pre <- fit$times <= sim$t0
  direct <- fit$Y[pre, sim$treated] -
    fit$Y[pre, sim$pure, drop = FALSE] %*% fit$implied_weights[sim$pure]
  cmp <- iscm_compare(sim$data, "unit", "time", "Y", sim$treated, sim$affected,
                      sim$t0)
  expect_equal(cmp$implied_rmspe[["full"]], rmspe(as.numeric(direct)),
               tolerance = 1e-8)
})

test_that("iscm_decision uses the implied RMSPEs when supplied", {
  d <- iscm_decision(rmspe_ratio = 260.74 / 119.07,
                     rmspe_full = 138.54, rmspe_simplified = 144.73)
  expect_equal(d$recommendation, "iscm_full")
  d2 <- iscm_decision(rmspe_ratio = 2, rmspe_full = 150, rmspe_simplified = 120)
  expect_equal(d2$recommendation, "iscm_simplified")
})

test_that("iscm_sensitivity returns a table across specifications of M", {
  sim <- simulate_iscm(n_affected = 2, n_pure = 6, tau = 2, spill = c(3, 1),
                       seed = 5)
  sens <- iscm_sensitivity(sim$data, "unit", "time", "Y", sim$treated,
    affected_sets = list(minimal = "A1", extended = c("A1", "A2")), t0 = sim$t0)
  expect_s3_class(sens, "iscm_sensitivity")
  expect_equal(nrow(sens$table), 2L)
  expect_true(all(c("theta", "alpha", "sigma_min") %in% names(sens$table)))
})
