## Monte Carlo: a design in which iSCM outperforms the restricted SCM.
## Units are matched on observed covariates (x1, x2) with the penalized SCM;
## TR and A1 share an unobserved characteristic that the pure controls closest
## in covariates lack. See ?simulate_iscm_covariates. True effects: 2 on TR,
## 3 on A1. Runs in about a minute.

library(iscm)

one_rep <- function(seed) {
  sim  <- simulate_iscm_covariates(seed = seed)
  args <- list(sim$data, "unit", "time", "Y", "TR", "A1", sim$t0,
               predictors = c("x1", "x2"), estimator = scm_penalized)
  fit   <- do.call(iscm, args)
  fit_s <- do.call(iscm, c(args, simplified = TRUE))
  rest  <- do.call(restricted_scm, args)
  cmp   <- do.call(iscm_compare, args)
  post  <- fit$times > sim$t0
  c(unrestricted = mean(fit$gaps_unrestricted[post, "TR"]),
    restricted   = mean(rest$gap[post]),
    full_iSCM    = mean(fit$effects["TR", ]),
    simpl_iSCM   = mean(fit_s$effects["TR", ]),
    gamma_iSCM   = mean(fit$effects["A1", ]),
    alpha        = fit$diagnostics$alpha,
    recommends_full = cmp$recommendation == "iscm_full")
}

res <- t(vapply(1:100, one_rep, numeric(7)))

est <- res[, c("unrestricted", "restricted", "full_iSCM", "simpl_iSCM")]
print(round(rbind(mean = colMeans(est),
                  sd = apply(est, 2, sd),
                  mean_abs_error = colMeans(abs(est - 2))), 2))
cat(sprintf("\nSpillover on A1 (true 3): mean %.2f, sd %.2f\n",
            mean(res[, "gamma_iSCM"]), sd(res[, "gamma_iSCM"])))
cat(sprintf("Full iSCM closer to the truth than restricted SCM: %.0f%% of replications\n",
            100 * mean(abs(est[, "full_iSCM"] - 2) < abs(est[, "restricted"] - 2))))
cat(sprintf("iscm_compare() recommends the full iSCM: %.0f%% of replications\n",
            100 * mean(res[, "recommends_full"])))
