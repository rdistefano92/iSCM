sim <- simulate_iscm(n_affected = 1, n_pure = 6, tau = 2, spill = 3, seed = 7)
fit <- iscm(sim$data, "unit", "time", "Y", sim$treated, sim$affected, sim$t0)

test_that("iscm_purge subtracts the iSCM effects from the focal units only", {
  purged <- iscm_purge(fit, sim$data, "unit", "time", "Y")
  post_tr <- sim$data$unit == "TR" & sim$data$time > sim$t0
  pre_tr  <- sim$data$unit == "TR" & sim$data$time <= sim$t0
  expect_equal(purged$Y[post_tr], sim$data$Y[post_tr] - fit$effects["TR", ])
  expect_equal(purged$Y[pre_tr], sim$data$Y[pre_tr])
  pure <- sim$data$unit %in% sim$pure
  expect_equal(purged$Y[pure], sim$data$Y[pure])
  kept <- iscm_purge(fit, sim$data, "unit", "time", "Y", exclude = "A1")
  a1 <- sim$data$unit == "A1"
  expect_equal(kept$Y[a1], sim$data$Y[a1])
})

test_that("iscm_placebo_fits reproduces iscm_placebo with the same weights", {
  ids <- unique(sim$data$unit)
  Y <- with(sim$data, tapply(Y, list(time, unit), sum))[, ids]
  pre <- seq_len(sim$t0)
  pw <- lapply(setNames(sim$pure, sim$pure), function(j) {
    donors <- setdiff(ids, j)
    scm_base(Y[pre, j], Y[pre, donors, drop = FALSE], standardize = FALSE)
  })
  for (tg in c("TR", "A1")) {
    a <- iscm_placebo(fit, sim$data, "unit", "time", "Y", target = tg)
    b <- iscm_placebo_fits(fit, sim$data, "unit", "time", "Y",
                           placebo_weights = pw, target = tg)
    expect_equal(b$table, a$table, tolerance = 1e-6)
    expect_equal(b$p_value, a$p_value)
    expect_equal(b$rank, a$rank)
  }
})

test_that("the target's post-intervention gap equals the iSCM estimate", {
  pb <- iscm_placebo(fit, sim$data, "unit", "time", "Y")
  row <- pb$table[pb$table$unit == "TR", ]
  expect_equal(row$post_rmspe, rmspe(fit$effects["TR", ]), tolerance = 1e-10)
})

test_that("iscm_placebo stops on assembled objects and points to the new tools", {
  asm <- iscm_assemble(
    gaps = list(TR = fit$fits$TR$gap, A1 = fit$fits$A1$gap),
    weights = list(TR = fit$fits$TR$w, A1 = fit$fits$A1$w),
    t0 = sim$t0, pure = sim$pure)
  expect_error(iscm_placebo(asm, sim$data, "unit", "time", "Y"),
               "iscm_placebo_fits")
})

test_that("iscm_placebo_fits validates its inputs", {
  pw <- list(C1 = c(C2 = 1))
  expect_error(iscm_placebo_fits(fit, sim$data, "unit", "time", "Y",
                                 placebo_weights = pw), "no entry for")
  bad <- setNames(lapply(sim$pure, function(j) c(ZZ = 1)), sim$pure)
  expect_error(iscm_placebo_fits(fit, sim$data, "unit", "time", "Y",
                                 placebo_weights = bad), "not found: ZZ")
})

test_that("the printed placebo result reports the caveats of Online Appendix B", {
  pb <- iscm_placebo(fit, sim$data, "unit", "time", "Y")
  expect_equal(pb$min_p_value, 1 / 7)
  out <- capture.output(print(pb))
  expect_true(any(grepl("Smallest attainable p-value with 6 pure controls: 0.143",
                        out, fixed = TRUE)))
  expect_true(any(grepl("alpha", out)))
})
