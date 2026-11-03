## Monte Carlo: a design in which iSCM outperforms the restricted SCM.
## Units are matched on observed covariates (x1, x2) with the penalized SCM;
## TR and A1 share an unobserved characteristic that some of the pure controls
## closest in covariates lack. See ?simulate_iscm_covariates. True effects: 2 on
## TR, 3 on A1. For each replication the decision rule of Section 5.1
## (iscm_compare()) selects a specification. Runs in about two minutes.

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
  est <- c(unrestricted = mean(fit$gaps_unrestricted[post, "TR"]),
           restricted   = mean(rest$gap[post]),
           full_iSCM    = mean(fit$effects["TR", ]),
           simpl_iSCM   = mean(fit_s$effects["TR", ]))
  chosen <- switch(cmp$recommendation, restricted = est[["restricted"]],
                   iscm_full = est[["full_iSCM"]], iscm_simplified = est[["simpl_iSCM"]],
                   NA_real_)
  data.frame(t(est), chosen = chosen, gamma = mean(fit$effects["A1", ]),
             alpha = cmp$amplification, recommendation = cmp$recommendation)
}

res <- do.call(rbind, lapply(1:100, one_rep))

est <- res[, c("unrestricted", "restricted", "full_iSCM", "simpl_iSCM", "chosen")]
print(round(rbind(mean = colMeans(est, na.rm = TRUE),
                  sd = apply(est, 2, sd, na.rm = TRUE),
                  mean_abs_error = colMeans(abs(est - 2), na.rm = TRUE)), 2))
cat(sprintf("\nSpillover on A1 (true 3): mean %.2f, sd %.2f\n",
            mean(res$gamma), sd(res$gamma)))
cat(sprintf("Amplification factor alpha: mean %.2f, range %.2f-%.2f\n",
            mean(res$alpha), min(res$alpha), max(res$alpha)))
cat(sprintf("Full iSCM closer to the truth than restricted SCM: %.0f%% of replications\n",
            100 * mean(abs(res$full_iSCM - 2) < abs(res$restricted - 2))))
cat("\nSpecification recommended by iscm_compare() (Section 5.1):\n")
print(table(res$recommendation))
