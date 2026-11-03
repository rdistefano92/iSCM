test_that("iscm_decision reproduces the paper's Germany recommendation", {
  # Section 7: validation RMSPE 358.09 vs 93.08, imbalance 1.50 vs 0.26, alpha 1.16
  d <- iscm_decision(rmspe_ratio = 358.09 / 93.08, imbalance_ratio = 1.50 / 0.26,
                     alpha = 1.16)
  expect_equal(d$recommendation, "iscm_full")
})

test_that("iscm_decision covers every branch of Section 5.1", {
  expect_equal(iscm_decision(2, alpha = 1.1, substantial = FALSE)$recommendation,
               "restricted")
  # comparable fits
  expect_equal(iscm_decision(1.10, imbalance_ratio = 1.05, alpha = 1.1)$recommendation,
               "restricted")
  # intermediate (between 25% and 50%) is not substantially worse
  d <- iscm_decision(1.35, alpha = 1.1)
  expect_equal(d$recommendation, "restricted")
  expect_match(d$rationale, "25% and 50%")
  # substantially worse on either indicator warrants iSCM
  expect_equal(iscm_decision(1.1, imbalance_ratio = 2, alpha = 1.1)$recommendation,
               "iscm_full")
  expect_equal(iscm_decision(1.6, alpha = 1.2)$recommendation, "iscm_full")
  # large amplification: simplified if the affected units' restricted fits are adequate
  expect_equal(iscm_decision(2, alpha = 1.4, affected_rmspe_ratio = 1.1)$recommendation,
               "iscm_simplified")
  expect_equal(iscm_decision(2, alpha = 1.4, affected_rmspe_ratio = 1.8)$recommendation,
               "none_reliable")
  # the alpha threshold is strict: alpha must be below 1.25
  expect_equal(iscm_decision(2, alpha = 1.25, affected_rmspe_ratio = 1)$recommendation,
               "iscm_simplified")
})

test_that("iscm_decision uses 1/sigma_min when alpha is NA (m > 2)", {
  d <- iscm_decision(rmspe_ratio = 2, sigma_min = 0.5, affected_rmspe_ratio = c(1, 1.1))
  expect_equal(d$recommendation, "iscm_simplified")       # 1/0.5 = 2 >= 1.25
  expect_equal(iscm_decision(2, sigma_min = 0.9)$recommendation, "iscm_full")
  expect_error(iscm_decision(2), "alpha")
})

test_that("iscm_compare returns a full comparison with a validation window", {
  sim <- simulate_iscm(seed = 7, spill = 3)
  cmp <- iscm_compare(sim$data, "unit", "time", "Y", sim$treated, sim$affected, sim$t0)
  expect_s3_class(cmp, "iscm_compare")
  expect_true(cmp$recommendation %in%
                c("restricted", "iscm_full", "iscm_simplified", "none_reliable"))
  expect_true(all(c("imbalance_unres", "rmspe_unres", "rmspe_val_unres")
                  %in% names(cmp$comparison)))
  expect_equal(nrow(cmp$comparison), 1L + length(sim$affected))
  # default window: last third of the 24 pre-intervention periods
  expect_equal(cmp$validation$valid, c(17, 24))
  expect_equal(cmp$validation$train, c(1, 16))
  cmp0 <- iscm_compare(sim$data, "unit", "time", "Y", sim$treated, sim$affected,
                       sim$t0, validation = NULL)
  expect_null(cmp0$comparison$rmspe_val_unres)
  expect_equal(cmp0$rmspe_ratio, cmp0$rmspe_ratio_in_sample)
})

test_that("the validation RMSPE takes precedence over the in-sample one", {
  # in-sample: restricted better; validation: restricted much worse
  gU <- list(TR = c(rep(0.2, 30), rep(-2, 14)), A = c(rep(0.2, 30), rep(1, 14)))
  gR <- list(TR = c(rep(0.1, 30), rep(-3, 14)), A = c(rep(0.2, 30), rep(1, 14)))
  w  <- list(TR = c(A = 0.42, C = 0.58), A = c(TR = 0.33, C = 0.67))
  a <- iscm_compare_fits(gU, gR, w, t0 = 30, times = 1:44)
  expect_equal(a$recommendation, "restricted")
  b <- iscm_compare_fits(gU, gR, w, t0 = 30, times = 1:44,
                         rmspe_val_unres = c(TR = 90, A = 200),
                         rmspe_val_restr = c(TR = 360, A = 210))
  expect_equal(b$recommendation, "iscm_full")
  expect_equal(unname(round(b$alpha, 2)), 1.16)
})

test_that("Germany-like figures lead to the full iSCM", {
  cmp <- iscm_compare_fits(
    gaps_unres = list("West Germany" = c(rep(0.119, 30), rep(-2, 14)),
                      "Austria"      = c(rep(0.195, 30), rep(0.3, 14))),
    gaps_restr = list("West Germany" = c(rep(0.261, 30), rep(-3, 14)),
                      "Austria"      = c(rep(0.181, 30), rep(1.3, 14))),
    weights    = list("West Germany" = c("Austria" = 0.42, C = 0.58),
                      "Austria"      = c("West Germany" = 0.33, C = 0.67)),
    imbalance_unres = c("West Germany" = 0.26, "Austria" = 0.53),
    imbalance_restr = c("West Germany" = 1.50, "Austria" = 2.06),
    t0 = 30, times = 1:44)
  expect_equal(cmp$recommendation, "iscm_full")
})

test_that("the implied counterfactuals match the closed forms for m = 2", {
  w2 <- 0.42; l1 <- 0.33; alpha <- 1 / (1 - w2 * l1)
  e1 <- rep(0.119, 30); e2 <- rep(0.195, 30); et2 <- rep(0.181, 30)
  args <- list(
    gaps_unres = list(TR = c(e1, rep(-2, 14)), A = c(e2, rep(0.3, 14))),
    gaps_restr = list(TR = c(rep(0.271, 30), rep(-3, 14)), A = c(et2, rep(1.3, 14))),
    weights    = list(TR = c(A = w2, C = 1 - w2), A = c(TR = l1, C = 1 - l1)),
    t0 = 30, times = 1:44)
  cmp <- do.call(iscm_compare_fits, args)
  expect_equal(cmp$implied_rmspe[["full"]], alpha * (0.119 + w2 * 0.195))
  expect_equal(cmp$implied_rmspe[["simplified"]], 0.119 + w2 * 0.181)
  expect_equal(cmp$implied_rmspe[["restricted"]], 0.271)
  expect_equal(cmp$recommendation, "iscm_full")               # paper rule: alpha < 1.25
  alt <- do.call(iscm_compare_fits, c(args, rule = "implied_rmspe"))
  expect_equal(alt$recommendation, "iscm_simplified")         # same-sign residuals add up
})

test_that("the implied full-iSCM gap equals the gap of the implied weights", {
  sim <- simulate_iscm(n_affected = 2, n_pure = 6, spill = c(3, 1), seed = 5)
  fit <- iscm(sim$data, "unit", "time", "Y", sim$treated, sim$affected, sim$t0)
  pre <- fit$times <= sim$t0
  direct <- fit$Y[pre, sim$treated] -
    fit$Y[pre, sim$pure, drop = FALSE] %*% fit$implied_weights[sim$pure]
  cmp <- iscm_compare(sim$data, "unit", "time", "Y", sim$treated, sim$affected,
                      sim$t0, validation = NULL)
  expect_equal(cmp$implied_rmspe[["full"]], rmspe(as.numeric(direct)),
               tolerance = 1e-8)
})

test_that("iscm_decision supports the implied-RMSPE alternative", {
  d <- iscm_decision(rmspe_ratio = 2, rule = "implied_rmspe",
                     rmspe_full = 138.54, rmspe_simplified = 144.73)
  expect_equal(d$recommendation, "iscm_full")
  d2 <- iscm_decision(rmspe_ratio = 2, rule = "implied_rmspe",
                      rmspe_full = 150, rmspe_simplified = 120)
  expect_equal(d2$recommendation, "iscm_simplified")
  expect_error(iscm_decision(2, rule = "implied_rmspe"), "rmspe_full")
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
