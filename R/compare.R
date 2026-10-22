# Internal: pre-intervention RMSPE of the implied pure-control counterfactuals of
# the full and the simplified iSCM (Step 4 of the decision rule).
#
# By the identity in the proof of Lemma 3.1, beta_t = Omega (Y_1t..Y_mt)' - P Y_Ct
# holds in every period, so the first component of Omega^{-1} e_t, where e_t
# stacks the pre-intervention gaps of the m unrestricted fits, is the gap of the
# implied weights, Y_1t - Y_Ct' a. For m = 2 it equals alpha (e_1t + w_2 e_2t)
# (Section 3.1). The simplified iSCM replaces the fits of the affected units with
# their restricted fits and has a triangular Omega, so its implied gap is
# e_1t + sum_j w_j e~_jt (Section 4.2).
.implied_pre_rmspe <- function(Omega, E_unres, E_restr, w_aff) {
  full <- tryCatch(
    rmspe(as.numeric(E_unres %*% solve(Omega)[1, ])),
    error = function(e) Inf)
  simpl <- E_unres[, 1]
  if (length(w_aff)) simpl <- simpl + as.numeric(E_restr[, -1, drop = FALSE] %*% w_aff)
  c(full = full, simplified = rmspe(simpl))
}

# Internal: the Section 5.1 recommendation logic, shared by iscm_compare(),
# iscm_compare_fits() and iscm_decision().
.iscm_decide <- function(substantial, worse_restricted, rule = "implied_rmspe",
                         implied = c(full = NA_real_, simplified = NA_real_),
                         diag_dominant = NA, amp = NA, alpha_tol = 2) {
  if (!substantial) {
    return(list(recommendation = "restricted",
         rationale = "Potentially affected units receive negligible weight, so all specifications are expected to agree (Remark 2.1); the restricted SCM avoids amplification."))
  }
  if (!worse_restricted) {
    return(list(recommendation = "restricted",
         rationale = "The restricted specification fits about as well as the unrestricted one, so it is preferable: it avoids the amplification of approximation errors."))
  }
  if (rule == "implied_rmspe") {
    if (isTRUE(implied[["full"]] <= implied[["simplified"]])) {
      list(recommendation = "iscm_full",
           rationale = sprintf("The affected units are essential for the fit. The implied pure-control counterfactual of the full iSCM fits the pre-intervention period better than that of the simplified iSCM (RMSPE %s vs %s), so the gain in fit outweighs the amplification.",
                               format(signif(implied[["full"]], 4)),
                               format(signif(implied[["simplified"]], 4))))
    } else {
      list(recommendation = "iscm_simplified",
           rationale = sprintf("The affected units are essential for the fit. The implied pure-control counterfactual of the simplified iSCM fits the pre-intervention period better than that of the full iSCM (RMSPE %s vs %s), so avoiding the inversion step is preferable.",
                               format(signif(implied[["simplified"]], 4)),
                               format(signif(implied[["full"]], 4))))
    }
  } else if (isTRUE(diag_dominant) && isTRUE(amp <= alpha_tol)) {
    list(recommendation = "iscm_full",
         rationale = "The affected units are essential for the fit and the system is well conditioned, so the full iSCM can be applied safely.")
  } else {
    list(recommendation = "iscm_simplified",
         rationale = "The affected units are essential but the system is poorly conditioned; the simplified iSCM removes the contamination without the inversion step.")
  }
}

#' Apply the iSCM decision rule to pre-computed comparison quantities
#'
#' The recommendation logic of Section 5.1 of Di Stefano and Mellace,
#' separated from any estimation so that it can be applied to fits produced by any
#' estimator, including a fully customised \code{Synth} specification. Given the
#' ratio of the restricted to the unrestricted pre-intervention fit and the
#' pre-intervention RMSPE of the implied counterfactuals of the full and the
#' simplified iSCM, it returns the recommended specification; see
#' \code{\link{iscm_compare}} for the rule.
#'
#' If \code{rmspe_full} and \code{rmspe_simplified} are not supplied, the choice
#' between full and simplified iSCM falls back to the amplification threshold
#' (\code{rule = "alpha"} in \code{\link{iscm_compare}}).
#'
#' @param rmspe_ratio Ratio of the restricted to the unrestricted pre-intervention
#'   (or validation) RMSPE of the main treated unit.
#' @param imbalance_ratio Optional ratio of the restricted to the unrestricted
#'   predictor imbalance. Defaults to \code{NA} (ignored).
#' @param alpha Amplification factor \eqn{\alpha} for \eqn{m = 2}; use
#'   \code{sigma_min} for the general case.
#' @param sigma_min Smallest singular value of \eqn{\hat\Omega}; its inverse is the
#'   general amplification bound. Used when \code{alpha} is \code{NA}.
#' @param rmspe_full,rmspe_simplified Pre-intervention RMSPE of the implied
#'   pure-control counterfactual of the main treated unit under the full and the
#'   simplified iSCM (reported by \code{\link{iscm_compare}} and
#'   \code{\link{iscm_compare_fits}} as \code{implied_rmspe}).
#' @param diag_dominant Logical: is \eqn{\hat\Omega} strictly diagonally dominant?
#' @param substantial Logical: do the affected units receive substantial weight in
#'   the unrestricted fit (Step 1)? Defaults to \code{TRUE}.
#' @param ratio_tol,alpha_tol Thresholds; see \code{\link{iscm_compare}}.
#'
#' @return A list of class \code{"iscm_decision"} with \code{recommendation} in
#'   \code{c("restricted", "iscm_full", "iscm_simplified")} and a \code{rationale}.
#'
#' @examples
#' # German-reunification figures: the restricted SCM fits much worse
#' # (RMSPE 260.74 vs 119.07); the implied counterfactual of the full iSCM has
#' # pre-intervention RMSPE 138.54, that of the simplified iSCM 144.73.
#' iscm_decision(rmspe_ratio = 260.74 / 119.07,
#'               rmspe_full = 138.54, rmspe_simplified = 144.73)
#'
#' # Previous rule, based on the amplification factor only
#' iscm_decision(rmspe_ratio = 260.74 / 119.07, alpha = 1.16,
#'               diag_dominant = TRUE)
#' @export
iscm_decision <- function(rmspe_ratio, imbalance_ratio = NA, alpha = NA,
                          sigma_min = NA, diag_dominant = NA, substantial = TRUE,
                          ratio_tol = 1.25, alpha_tol = 2,
                          rmspe_full = NA, rmspe_simplified = NA) {
  amp <- if (!is.na(alpha)) alpha else if (!is.na(sigma_min)) 1 / sigma_min else NA
  worse <- (!is.na(imbalance_ratio) && imbalance_ratio >= ratio_tol) ||
           (!is.na(rmspe_ratio) && rmspe_ratio >= ratio_tol)
  rule <- if (!is.na(rmspe_full) && !is.na(rmspe_simplified)) "implied_rmspe" else "alpha"
  out <- .iscm_decide(substantial, worse, rule,
                      implied = c(full = rmspe_full, simplified = rmspe_simplified),
                      diag_dominant = diag_dominant, amp = amp, alpha_tol = alpha_tol)
  structure(out, class = "iscm_decision")
}

#' @export
print.iscm_decision <- function(x, ...) {
  cat("Recommendation:", switch(x$recommendation,
      restricted = "restricted SCM", iscm_full = "full iSCM",
      iscm_simplified = "simplified iSCM"), "\n")
  cat(strwrap(x$rationale, exdent = 0), sep = "\n"); cat("\n")
  invisible(x)
}

# Internal: refit a unit's SC on a training sub-window and return the gap over a
# validation window (Step 3 variant of Section 5.1).
.refit_window <- function(focal_id, donor_ids, data, id, time, outcome,
                          predictors, ids, times, train, valid, Y, v,
                          estimator, zvar, ...) {
  Xtr <- .predictor_matrix(data, id, time, outcome, predictors, ids, times, train)
  zcol <- if (is.null(zvar)) outcome else zvar
  Ztr <- .panel_matrix(data, id, time, zcol)$Y[which(times %in% train), , drop = FALSE]
  fit <- .fit_one(focal_id, donor_ids, Xtr, Y, Ztr, v, estimator, ...)
  vi <- which(times %in% valid)
  rmspe(fit$gap[vi])
}

#' Decision rule from pre-computed fits
#'
#' Version of \code{\link{iscm_compare}} that applies the Section 5.1 decision rule
#' to synthetic controls estimated externally (e.g. the exact \code{Synth}
#' specification of the empirical application), rather than re-estimating them with
#' a package backend. The user supplies, for the main treated unit and each
#' potentially affected unit, the full gap series of the unrestricted and
#' restricted fits and the unrestricted donor weights; predictor imbalances and
#' out-of-sample (validation) RMSPEs may be supplied as well.
#'
#' @param gaps_unres,gaps_restr Named lists of full gap series (unrestricted and
#'   restricted fits), one element per focal unit, main treated first.
#' @param weights Named list of unrestricted donor-weight vectors, parallel to
#'   \code{gaps_unres}, used to build \eqn{\hat\Omega} and to read the weights on
#'   the affected units (Step 1).
#' @param t0 Last pre-intervention period.
#' @param times Optional period labels aligned with the gap series.
#' @param imbalance_unres,imbalance_restr Optional named numeric vectors of
#'   predictor imbalances per focal unit (Step 2). When the predictors are on
#'   very different scales, standardise each discrepancy by the predictor's
#'   cross-sectional standard deviation before aggregating (Section 5.1, Step 2),
#'   e.g. \code{s <- apply(cbind(X1, X0), 1, sd); sum(abs(X1 - X0 \%*\% w) / s)}
#'   from \code{Synth::dataprep}, using the same \code{s} for the unrestricted and
#'   restricted fits. Pass \code{sum(abs(X1 - X0 \%*\% w))} for the raw sum.
#' @param rmspe_val_unres,rmspe_val_restr Optional named numeric vectors of
#'   out-of-sample (validation) RMSPEs per focal unit (Step 3 variant).
#' @param weight_tol,ratio_tol,alpha_tol,rule Thresholds and rule; see
#'   \code{\link{iscm_compare}}.
#'
#' @return An object of class \code{"iscm_compare"}.
#' @seealso \code{\link{iscm_compare}}, \code{\link{iscm_decision}},
#'   \code{\link{iscm_assemble}}.
#'
#' @examples
#' # Unrestricted and restricted fits estimated outside the package
#' sim <- simulate_iscm(n_affected = 1, n_pure = 6, tau = 2, spill = 3, seed = 7)
#' Y <- with(sim$data, tapply(Y, list(time, unit), sum))
#' pre <- seq_len(sim$t0)
#' sc <- function(unit, donors) {
#'   w <- scm_base(Y[pre, unit], Y[pre, donors])
#'   list(w = w, gap = as.numeric(Y[, unit] - Y[, donors] %*% w))
#' }
#' all_but <- function(u) setdiff(colnames(Y), u)
#' un <- list(TR = sc("TR", all_but("TR")), A1 = sc("A1", all_but("A1")))
#' re <- list(TR = sc("TR", sim$pure),      A1 = sc("A1", sim$pure))
#'
#' iscm_compare_fits(
#'   gaps_unres = lapply(un, `[[`, "gap"),
#'   gaps_restr = lapply(re, `[[`, "gap"),
#'   weights    = lapply(un, `[[`, "w"),
#'   t0 = sim$t0)
#' @export
iscm_compare_fits <- function(gaps_unres, gaps_restr, weights, t0, times = NULL,
                              imbalance_unres = NULL, imbalance_restr = NULL,
                              rmspe_val_unres = NULL, rmspe_val_restr = NULL,
                              weight_tol = 0.05, ratio_tol = 1.25, alpha_tol = 2,
                              rule = c("implied_rmspe", "alpha")) {
  rule <- match.arg(rule)
  focal <- names(gaps_unres)
  treated <- focal[1]; affected <- focal[-1]
  len <- length(gaps_unres[[1]])
  if (is.null(times)) times <- seq_len(len)
  pre <- if (t0 %in% times) times <= t0 else seq_along(times) <= t0

  cross <- lapply(focal, function(fi) weights[[fi]][names(weights[[fi]]) %in% focal])
  Omega <- build_omega(cross, focal)
  diag <- omega_diagnostics(Omega)
  amp <- if (!is.na(diag$alpha)) diag$alpha else diag$amplification

  rows <- lapply(focal, function(fi) {
    out <- data.frame(unit = fi,
      role = if (fi == treated) "treated" else "affected",
      rmspe_unres = rmspe(gaps_unres[[fi]][pre]),
      rmspe_restr = rmspe(gaps_restr[[fi]][pre]),
      stringsAsFactors = FALSE)
    if (!is.null(imbalance_unres)) out$imbalance_unres <- imbalance_unres[[fi]]
    if (!is.null(imbalance_restr)) out$imbalance_restr <- imbalance_restr[[fi]]
    if (!is.null(rmspe_val_unres)) out$rmspe_val_unres <- rmspe_val_unres[[fi]]
    if (!is.null(rmspe_val_restr)) out$rmspe_val_restr <- rmspe_val_restr[[fi]]
    out
  })
  tab <- do.call(rbind, rows)

  w_main <- weights[[treated]]
  aff_weights <- vapply(affected, function(a)
    if (a %in% names(w_main)) w_main[[a]] else 0, numeric(1))
  substantial <- any(aff_weights >= weight_tol)
  main_in_affected <- vapply(affected, function(a) {
    wa <- weights[[a]]; if (treated %in% names(wa)) wa[[treated]] else 0
  }, numeric(1))

  main <- tab[tab$unit == treated, ]
  imb_ratio <- if (!is.null(imbalance_unres))
    main$imbalance_restr / main$imbalance_unres else NA_real_
  rmspe_ratio <- if (!is.null(rmspe_val_unres))
    main$rmspe_val_restr / main$rmspe_val_unres else main$rmspe_restr / main$rmspe_unres
  worse <- (!is.na(imb_ratio) && imb_ratio >= ratio_tol) || (rmspe_ratio >= ratio_tol)

  E_unres <- vapply(focal, function(fi) as.numeric(gaps_unres[[fi]][pre]),
                    numeric(sum(pre)))
  E_restr <- vapply(focal, function(fi) as.numeric(gaps_restr[[fi]][pre]),
                    numeric(sum(pre)))
  implied <- .implied_pre_rmspe(Omega, matrix(E_unres, ncol = length(focal)),
                                matrix(E_restr, ncol = length(focal)), aff_weights)
  dec <- .iscm_decide(substantial, worse, rule, implied = implied,
                      diag_dominant = diag$diag_dominant, amp = amp,
                      alpha_tol = alpha_tol)

  structure(list(
    treated = treated, affected = affected,
    step1_affected_weights = aff_weights, step1_substantial = substantial,
    comparison = tab, amplification = amp, alpha = diag$alpha,
    sigma_min = diag$sigma_min, diag_dominant = diag$diag_dominant,
    main_weight_in_affected = main_in_affected,
    imbalance_ratio = imb_ratio, rmspe_ratio = rmspe_ratio,
    implied_rmspe = c(implied, restricted = main$rmspe_restr), rule = rule,
    validation = if (is.null(rmspe_val_unres)) NULL else TRUE,
    recommendation = dec$recommendation, rationale = dec$rationale
  ), class = "iscm_compare")
}

#' Data-driven choice among restricted SCM and full/simplified iSCM
#'
#' Runs the decision rule of Section 5.1 of Di Stefano and Mellace and
#' returns a recommendation. For the main treated unit it (Step 1) reports the
#' weights assigned to the potentially affected units, (Step 2) compares the
#' predictor imbalance of the unrestricted and restricted specifications,
#' (Step 3) compares their pre-intervention outcome RMSPE, optionally with a
#' training/validation split, and (Step 4) chooses between the full and the
#' simplified iSCM. Steps 1-3 are then repeated for each potentially affected unit,
#' treating it as if it were the main treated unit, to assess Assumption 2.3 and
#' whether the main treated unit is needed in its donor pool (relevant for the
#' simplified iSCM).
#'
#' @inheritParams iscm
#' @param validation Optional integer: number of final pre-intervention periods to
#'   hold out as a validation window (the earlier pre-periods are used for
#'   training). If \code{NULL}, only the in-sample RMSPE comparison is reported.
#' @param weight_tol Weight below which a potentially affected unit is treated as
#'   receiving negligible weight (Step 1 / Remark 2.1). Default \code{0.05}.
#' @param ratio_tol Factor above which the restricted specification is judged
#'   substantially worse than the unrestricted one on predictor imbalance or
#'   RMSPE, warranting the iSCM correction. Default \code{1.25}.
#' @param rule How Step 4 chooses between the full and the simplified iSCM.
#'   \code{"implied_rmspe"} (default) compares the pre-intervention RMSPE of
#'   their implied pure-control counterfactuals; \code{"alpha"} uses the
#'   amplification threshold \code{alpha_tol} (the rule of earlier versions).
#' @param alpha_tol Used only with \code{rule = "alpha"}: amplification above
#'   which the simplified iSCM is preferred. Default \code{2}. For \eqn{m > 2}
#'   the general bound \code{1/sigma_min} is used in place of \code{alpha}.
#' @param standardize Logical; if \code{TRUE} (default), each predictor's absolute
#'   discrepancy is divided by the predictor's cross-sectional standard deviation
#'   (computed across all units and held fixed across the unrestricted and
#'   restricted fits) before the Step 2 imbalances are aggregated, as recommended
#'   in Section 5.1 when the predictors are on very different scales. Set to
#'   \code{FALSE} for the unstandardised sum of absolute discrepancies.
#'
#' @details \strong{Step 4.} Both versions of iSCM estimate the counterfactual of
#'   the main treated unit with a weighted average of pure controls (Lemma 3.1,
#'   Section 4.2). Their errors are \eqn{-\alpha(B^1_{sc} + \hat w_2 B^2_{sc})}
#'   and \eqn{-(B^1_{sc} + \hat w_2 \tilde B^2_{sc})} for \eqn{m = 2}
#'   (equations 4.6 and 4.7): the full iSCM amplifies the errors of the
#'   unrestricted fits by \eqn{\alpha}, the simplified iSCM avoids the
#'   amplification but uses the restricted fits of the affected units. The default
#'   rule evaluates both expressions where they are observable, before the
#'   intervention: it computes the pre-intervention RMSPE of the two implied
#'   counterfactuals and recommends the one with the smaller value. For the
#'   German application these are 138.54 (full) and 144.73 (simplified). This is
#'   a diagnostic, not a guarantee: a good pre-intervention fit is neither
#'   necessary nor sufficient for a good post-intervention counterfactual.
#'   The values are returned in \code{implied_rmspe}, together with the
#'   pre-intervention RMSPE of the restricted SCM, which can be compared with
#'   them as in Section 3.1.
#'
#'   \strong{Thresholds.} \code{weight_tol} and \code{ratio_tol} are package
#'   defaults, not values prescribed by the paper. Section 5.1 regards fits as
#'   comparable when their RMSPEs differ by less than 25\% and the restricted
#'   fit as substantially worse above 50\%; with the default
#'   \code{ratio_tol = 1.25} the range in between counts as "worse". The
#'   recommendation is meant to select the specification to emphasise: Section
#'   5.1 recommends reporting all four (unrestricted, restricted, full and
#'   simplified iSCM).
#'
#' @return An object of class \code{"iscm_compare"} with the per-step comparison
#'   tables, the diagnostics, and a \code{recommendation} in
#'   \code{c("restricted", "iscm_full", "iscm_simplified")} with a rationale.
#'   Printing it shows the four steps in order. The \code{comparison} table has
#'   one row per focal unit, so the rows of the affected units show whether
#'   their own fits are adequate (Assumption 2.3) and how much they lose when
#'   the main treated unit is excluded (relevant for the simplified iSCM).
#'
#' @seealso \code{\link{iscm_compare_fits}} and \code{\link{iscm_decision}} for
#'   fits estimated outside the package; \code{\link{iscm_sensitivity}} for the
#'   choice of the affected units.
#'
#' @examples
#' sim <- simulate_iscm(n_affected = 1, n_pure = 6, tau = 2, spill = 3, seed = 7)
#' cmp <- iscm_compare(sim$data, "unit", "time", "Y",
#'                     treated = sim$treated, affected = sim$affected,
#'                     t0 = sim$t0, validation = 8)
#' cmp                       # Steps 1-4 and the recommendation
#' cmp$recommendation
#' cmp$comparison            # imbalance and RMSPE for every focal unit
#' @export
iscm_compare <- function(data, id, time, outcome, treated, affected, t0,
                         predictors = NULL, pure = NULL, estimator = scm_base,
                         v = NULL, zvar = NULL, validation = NULL,
                         weight_tol = 0.05, ratio_tol = 1.25, alpha_tol = 2,
                         standardize = TRUE, rule = c("implied_rmspe", "alpha"),
                         ...) {
  rule <- match.arg(rule)
  pm <- .panel_matrix(data, id, time, outcome)
  ids <- pm$ids; times <- pm$times; Y <- pm$Y
  treated <- as.character(treated); affected <- as.character(affected)
  focal <- c(treated, affected)
  if (is.null(pure)) pure <- setdiff(ids, focal) else pure <- as.character(pure)
  pre <- times[times <= t0]
  Xmat <- .predictor_matrix(data, id, time, outcome, predictors, ids, times, pre)
  zcol <- if (is.null(zvar)) outcome else zvar
  Zmat <- .panel_matrix(data, id, time, zcol)$Y[which(times %in% pre), , drop = FALSE]

  # Step 2 predictor imbalance. When the predictors are on very different scales
  # (as in the German-reunification application), Section 5.1 Step 2 prescribes
  # standardising each predictor by its cross-sectional standard deviation before
  # aggregating the absolute discrepancies. The scale is fixed across all units
  # and reused for the unrestricted and restricted fits, so the two imbalances
  # remain comparable. Set standardize = FALSE for the raw sum.
  if (standardize) {
    sdv <- apply(Xmat, 1L, stats::sd)
    sdv[sdv == 0 | !is.finite(sdv)] <- 1
  } else {
    sdv <- rep(1, nrow(Xmat))
  }
  imbalance <- function(focal_id, donor_ids) {
    w <- .fit_one(focal_id, donor_ids, Xmat, Y, Zmat, v, estimator, ...)$w
    sum(abs(Xmat[, focal_id] - Xmat[, donor_ids, drop = FALSE] %*% w) / sdv)
  }
  fit_gap <- function(focal_id, donor_ids) {
    .fit_one(focal_id, donor_ids, Xmat, Y, Zmat, v, estimator, ...)$gap
  }

  # ---- the iSCM fit for Omega / amplification (Step 4) ------------------------
  fit <- iscm(data, id, time, outcome, treated, affected, t0,
              predictors = predictors, pure = pure, estimator = estimator,
              v = v, zvar = zvar, ...)
  diag <- fit$diagnostics
  amp <- if (!is.na(diag$alpha)) diag$alpha else diag$amplification

  # ---- per-unit comparison (main treated first, then each affected) ----------
  if (!is.null(validation)) {
    valid <- utils::tail(pre, validation); train <- setdiff(pre, valid)
  }
  gaps_re <- list()
  rows <- lapply(focal, function(fi) {
    donors_un <- setdiff(ids, fi)
    donors_re <- pure
    g_un <- fit_gap(fi, donors_un); g_re <- fit_gap(fi, donors_re)
    gaps_re[[fi]] <<- g_re
    out <- data.frame(
      unit = fi,
      role = if (fi == treated) "treated" else "affected",
      imbalance_unres = imbalance(fi, donors_un),
      imbalance_restr = imbalance(fi, donors_re),
      rmspe_unres = rmspe(g_un[times <= t0]),
      rmspe_restr = rmspe(g_re[times <= t0]),
      stringsAsFactors = FALSE)
    if (!is.null(validation)) {
      out$rmspe_val_unres <- .refit_window(fi, donors_un, data, id, time, outcome,
        predictors, ids, times, train, valid, Y, v, estimator, zvar, ...)
      out$rmspe_val_restr <- .refit_window(fi, donors_re, data, id, time, outcome,
        predictors, ids, times, train, valid, Y, v, estimator, zvar, ...)
    }
    out
  })
  tab <- do.call(rbind, rows)

  # ---- Step 1: weights on affected units in the main treated fit -------------
  w_main <- fit$fits[[1]]$w
  aff_weights <- vapply(affected, function(a)
    if (a %in% names(w_main)) w_main[[a]] else 0, numeric(1))
  substantial <- any(aff_weights >= weight_tol)

  # For each affected unit: weight it gives to the main treated unit (whether the
  # simplified iSCM is costless for it, Section 5.1 / Corollary 4.2).
  main_in_affected <- vapply(affected, function(a) {
    wa <- fit$fits[[match(a, focal)]]$w
    if (treated %in% names(wa)) wa[[treated]] else 0
  }, numeric(1))

  # ---- recommendation --------------------------------------------------------
  main <- tab[tab$unit == treated, ]
  imb_ratio <- main$imbalance_restr / main$imbalance_unres
  rmspe_ratio <- if (!is.null(validation))
    main$rmspe_val_restr / main$rmspe_val_unres
  else main$rmspe_restr / main$rmspe_unres
  worse_restricted <- (imb_ratio >= ratio_tol) || (rmspe_ratio >= ratio_tol)

  # ---- Step 4: implied counterfactuals of full vs simplified iSCM -----------
  pre_idx <- which(times <= t0)
  E_unres <- fit$gaps_unrestricted[pre_idx, focal, drop = FALSE]
  E_restr <- vapply(focal, function(fi) gaps_re[[fi]][pre_idx],
                    numeric(length(pre_idx)))
  implied <- .implied_pre_rmspe(fit$Omega, E_unres,
                                matrix(E_restr, ncol = length(focal)), aff_weights)
  dec <- .iscm_decide(substantial, worse_restricted, rule, implied = implied,
                      diag_dominant = diag$diag_dominant, amp = amp,
                      alpha_tol = alpha_tol)
  rec <- dec$recommendation; why <- dec$rationale

  structure(list(
    treated = treated, affected = affected,
    step1_affected_weights = aff_weights,
    step1_substantial = substantial,
    comparison = tab,
    amplification = amp, alpha = diag$alpha, sigma_min = diag$sigma_min,
    diag_dominant = diag$diag_dominant,
    main_weight_in_affected = main_in_affected,
    imbalance_ratio = imb_ratio, rmspe_ratio = rmspe_ratio,
    implied_rmspe = c(implied, restricted = main$rmspe_restr), rule = rule,
    validation = validation,
    recommendation = rec, rationale = why
  ), class = "iscm_compare")
}

#' @export
print.iscm_compare <- function(x, digits = 3, ...) {
  cat("iSCM specification comparison (Section 5.1)\n")
  cat("  Main treated unit:", x$treated, "\n\n")
  cat("Step 1 - weights on affected units in the unrestricted main fit:\n")
  print(round(x$step1_affected_weights, digits))
  if (!x$step1_substantial)
    cat("  -> all below weight_tol: affected units carry little contamination.\n")
  cat("\nSteps 2-3 - predictor imbalance and pre-intervention RMSPE:\n")
  print(format(x$comparison, digits = digits))
  cat(sprintf("\nRatios (restricted / unrestricted, main unit): imbalance = %s, RMSPE = %.2f\n",
              if (is.na(x$imbalance_ratio)) "n/a" else sprintf("%.2f", x$imbalance_ratio),
              x$rmspe_ratio))
  cat(sprintf("\nStep 4 - conditioning: %s = %.3f | sigma_min = %.3f | diag. dominant = %s\n",
              if (!is.na(x$alpha)) "alpha" else "1/sigma_min",
              x$amplification, x$sigma_min, x$diag_dominant))
  if (!is.null(x$implied_rmspe)) {
    ir <- x$implied_rmspe
    cat(sprintf("  Pre-intervention RMSPE of the implied counterfactual: full iSCM = %s | simplified iSCM = %s | restricted SCM = %s\n",
                format(signif(ir[["full"]], 4)), format(signif(ir[["simplified"]], 4)),
                format(signif(ir[["restricted"]], 4))))
    if (identical(x$rule, "alpha"))
      cat("  (rule = \"alpha\": the choice uses the amplification threshold)\n")
  }
  cat("\nRecommendation:", switch(x$recommendation,
      restricted = "restricted SCM",
      iscm_full = "full iSCM",
      iscm_simplified = "simplified iSCM"), "\n")
  cat(strwrap(x$rationale, indent = 2, exdent = 2), sep = "\n")
  cat("\n")
  invisible(x)
}

#' Sensitivity to the choice of potentially affected units
#'
#' Implements the sensitivity analysis of Section 5.2 of Di Stefano and Mellace:
#' re-estimates iSCM under several candidate definitions of the set of
#' potentially affected units \eqn{M} (for example a minimal specification with
#' only the most clearly exposed units and an extended one with all plausibly
#' exposed units) and reports, for each, the estimated effect on the main treated
#' unit together with the conditioning diagnostics. It also reports the weights
#' the units receive in the fully unrestricted fit, to inform which units are
#' empirically influential (Remark 2.1).
#'
#' @inheritParams iscm
#' @param affected_sets A named list of candidate \code{affected} vectors. Each
#'   element defines one specification of \eqn{M}.
#' @param ... Passed to \code{\link{iscm}}.
#'
#' @return An object of class \code{"iscm_sensitivity"} with a table of the treated
#'   effect and diagnostics across specifications.
#'
#' @seealso \code{\link{iscm}}, \code{\link{iscm_compare}}.
#'
#' @examples
#' sim <- simulate_iscm(n_affected = 2, n_pure = 6, tau = 2, spill = c(3, 1),
#'                      seed = 5)
#' iscm_sensitivity(sim$data, "unit", "time", "Y", treated = sim$treated,
#'                  affected_sets = list(minimal  = "A1",
#'                                       extended = c("A1", "A2")),
#'                  t0 = sim$t0)
#' @export
iscm_sensitivity <- function(data, id, time, outcome, treated, affected_sets, t0,
                             predictors = NULL, estimator = scm_base, v = NULL,
                             zvar = NULL, ...) {
  treated <- as.character(treated)
  if (is.null(names(affected_sets)))
    names(affected_sets) <- paste0("spec", seq_along(affected_sets))

  rows <- lapply(names(affected_sets), function(nm) {
    aff <- as.character(affected_sets[[nm]])
    fit <- iscm(data, id, time, outcome, treated, aff, t0,
                predictors = predictors, estimator = estimator, v = v,
                zvar = zvar, ...)
    data.frame(
      spec = nm,
      M = paste(aff, collapse = ", "),
      m = length(aff) + 1L,
      theta = mean(fit$effects[1, ]),
      alpha = fit$diagnostics$alpha,
      sigma_min = fit$diagnostics$sigma_min,
      diag_dominant = fit$diagnostics$diag_dominant,
      stringsAsFactors = FALSE)
  })
  tab <- do.call(rbind, rows)

  # weights every candidate unit receives in the fully unrestricted main fit
  all_aff <- unique(unlist(affected_sets))
  full <- iscm(data, id, time, outcome, treated, all_aff, t0,
               predictors = predictors, estimator = estimator, v = v,
               zvar = zvar, ...)
  wmain <- full$fits[[1]]$w
  unit_weights <- vapply(all_aff, function(a)
    if (a %in% names(wmain)) wmain[[a]] else 0, numeric(1))

  structure(list(treated = treated, table = tab, unit_weights = unit_weights),
            class = "iscm_sensitivity")
}

#' @export
print.iscm_sensitivity <- function(x, digits = 3, ...) {
  cat("iSCM sensitivity to the set of potentially affected units (Section 5.2)\n")
  cat("  Main treated unit:", x$treated, "\n\n")
  print(format(x$table, digits = digits))
  cat("\nWeights in the fully unrestricted main fit:\n")
  print(round(x$unit_weights, digits))
  invisible(x)
}
