sim  <- simulate_iscm(n_affected = 1, n_pure = 6, tau = 2, spill = 3, seed = 7)
fit  <- iscm(sim$data, "unit", "time", "Y", sim$treated, sim$affected, sim$t0)
rest <- restricted_scm(sim$data, "unit", "time", "Y", sim$treated, sim$affected, sim$t0)

test_that("plot.iscm returns ggplot objects for every type", {
  for (tp in c("gaps", "trends", "weights")) {
    p <- plot(fit, type = tp)
    expect_s3_class(p, "ggplot")
    expect_silent(ggplot2::ggplot_build(p))
  }
  p <- plot(fit, restricted = rest)
  expect_true("Restricted SCM" %in% unique(p$data$series))
  expect_s3_class(plot(fit, unit = c("TR", "A1")), "ggplot")
})

test_that("gaps plot uses the stored gaps", {
  p <- plot(fit)
  d <- p$data[p$data$series == "iSCM", ]
  expect_equal(d$value, as.numeric(fit$gaps_iscm[, "TR"]))
})

test_that("plot.iscm validates its inputs", {
  expect_error(plot(fit, unit = "C1"), "focal units")
  expect_error(plot(fit, unit = c("TR", "A1"), restricted = rest), "named list")
  asm <- iscm_assemble(gaps = list(TR = fit$fits$TR$gap, A1 = fit$fits$A1$gap),
                       weights = list(TR = fit$fits$TR$w, A1 = fit$fits$A1$w),
                       t0 = sim$t0, pure = sim$pure)
  expect_error(plot(asm, type = "trends"), "iscm_assemble")
  expect_s3_class(plot(asm), "ggplot")
})

test_that("plot.iscm_placebo draws ratios and gaps", {
  pb <- iscm_placebo(fit, sim$data, "unit", "time", "Y")
  expect_equal(dim(pb$gaps), c(length(fit$times), 7L))
  expect_equal(as.numeric(pb$gaps[fit$times > sim$t0, "TR"]),
               as.numeric(fit$effects["TR", ]), tolerance = 1e-10)
  for (tp in c("ratios", "gaps")) {
    p <- plot(pb, type = tp)
    expect_s3_class(p, "ggplot")
    expect_silent(ggplot2::ggplot_build(p))
  }
  p <- plot(pb, type = "gaps", prune = 0)
  expect_equal(unique(p$data$unit), "TR")
})
