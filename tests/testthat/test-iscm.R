# Build a perfect-fit panel with KNOWN cross-weights, then use an oracle backend
# that returns those weights, so iSCM must recover the true effects exactly.
make_perfect <- function(w2 = 0.42, l1 = 0.33, tau = 2, gamma = -0.7,
                         Tt = 30, T0 = 24, k = 4, seed = 3) {
  set.seed(seed)
  P <- matrix(NA_real_, Tt, k)
  for (j in seq_len(k))
    P[, j] <- 1 + 0.1 * (1:Tt) + sin((1:Tt) / (1 + j)) + rnorm(Tt, 0, 0.3)
  wc <- c(0.16, 0.09, 0.11, 0.22); wc <- wc / sum(wc) * (1 - w2)
  lc <- c(0.21, 0.31, 0.10, 0.06); lc <- lc / sum(lc) * (1 - l1)
  den <- 1 - w2 * l1
  YN_TR <- P %*% ((wc + w2 * lc) / den)
  YN_A1 <- P %*% ((lc + l1 * wc) / den)
  Y_TR <- YN_TR; Y_TR[(T0 + 1):Tt] <- Y_TR[(T0 + 1):Tt] + tau
  Y_A1 <- YN_A1; Y_A1[(T0 + 1):Tt] <- Y_A1[(T0 + 1):Tt] + gamma
  ids <- c("TR", "A1", paste0("C", seq_len(k)))
  Ym <- cbind(Y_TR, Y_A1, P); colnames(Ym) <- ids
  data <- data.frame(unit = rep(ids, each = Tt),
                     time = rep(seq_len(Tt), times = k + 2),
                     Y = as.vector(Ym), stringsAsFactors = FALSE)
  truw <- list(TR = c(A1 = w2, setNames(wc, paste0("C", seq_len(k)))),
               A1 = c(TR = l1, setNames(lc, paste0("C", seq_len(k)))))
  list(data = data, truw = truw, T0 = T0, tau = tau, gamma = gamma,
       w2 = w2, l1 = l1)
}

oracle_backend <- function(truw) {
  function(X1, X0, v = NULL, Z1 = NULL, Z0 = NULL, ...) {
    donors <- colnames(X0)
    who <- setdiff(c("TR", "A1"), donors)
    who <- who[who %in% names(truw)][1]
    w <- truw[[who]][donors]; w[is.na(w)] <- 0
    as.numeric(w)
  }
}

test_that("iSCM recovers the true effects exactly under a perfect fit", {
  s <- make_perfect()
  fit <- iscm(s$data, "unit", "time", "Y", "TR", "A1", t0 = s$T0,
              estimator = oracle_backend(s$truw))
  expect_equal(-fit$Omega[1, 2], s$w2, tolerance = 1e-8)
  expect_equal(-fit$Omega[2, 1], s$l1, tolerance = 1e-8)
  expect_equal(mean(fit$effects[1, ]), s$tau, tolerance = 1e-8)
  expect_equal(mean(fit$effects[2, ]), s$gamma, tolerance = 1e-8)
})

test_that("naive estimate equals tau - w2 * gamma (contamination)", {
  s <- make_perfect()
  fit <- iscm(s$data, "unit", "time", "Y", "TR", "A1", t0 = s$T0,
              estimator = oracle_backend(s$truw))
  naive <- mean(fit$gaps_unrestricted[fit$times > s$T0, "TR"])
  expect_equal(naive, s$tau - s$w2 * s$gamma, tolerance = 1e-8)
})

test_that("implied pure-control weights sum to one (Lemma 3.1)", {
  s <- make_perfect()
  fit <- iscm(s$data, "unit", "time", "Y", "TR", "A1", t0 = s$T0,
              estimator = oracle_backend(s$truw))
  expect_equal(sum(fit$implied_weights), 1, tolerance = 1e-8)
})

test_that("scm_base recovers effects on a well-conditioned simulated panel", {
  sim <- simulate_iscm(seed = 7, spill = 3)
  fit <- iscm(sim$data, "unit", "time", "Y", sim$treated, sim$affected, sim$t0)
  # iSCM is much closer to the truth than the contaminated naive estimate
  naive <- mean(fit$gaps_unrestricted[fit$times > sim$t0, sim$treated])
  expect_lt(abs(mean(fit$effects[1, ]) - 2), abs(naive - 2))
  expect_equal(mean(fit$effects[2, ]), 3, tolerance = 0.15)
  expect_true(fit$diagnostics$diag_dominant)
})
