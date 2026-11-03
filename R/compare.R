# Internal: pre-intervention RMSPE of the implied pure-control counterfactuals of
# the full and the simplified iSCM (an additional Step 4 diagnostic).
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

# Internal: classify a restricted / unrestricted ratio (Section 5.1, Step 3):
# comparable below `comparable_tol`, substantially worse from `worse_tol` on.
.classify_ratio <- function(r, comparable_tol, worse_tol) {
  if (is.null(r) || length(r) == 0L || is.na(r)) return(NA_character_)
  if (r < comparable_tol) "comparable" else if (r >= worse_tol) "worse" else "intermediate"
}

# Internal: the decision rule of Section 5.1, shared by iscm_compare(),
# iscm_compare_fits() and iscm_decision().
#   substantial   Step 1: do the affected units receive substantial weight?
#   imb_ratio     Step 2: restricted / unrestricted predictor imbalance (or NA)
#   rmspe_ratio   Step 3: restricted / unrestricted RMSPE of the main unit; the
#                 validation-split ratio when available, which takes precedence
#                 over the in-sample one (Section 5.1)
#   aff_ratios    restricted / unrestricted RMSPE of each affected unit
#   amp           alpha (m = 2) or 1 / sigma_min
.iscm_decide <- function(substantial, imb_ratio, rmspe_ratio, aff_ratios, amp,
                         rule = "paper", implied = c(full = NA, simplified = NA),
                         comparable_tol = 1.25, worse_tol = 1.5, alpha_tol = 1.25) {
  out <- function(rec, why) list(recommendation = rec, rationale = why)
  if (!substantial)
    return(out("restricted", "Step 1: the potentially affected units receive negligible weight, so all specifications are expected to perform similarly (Remark 2.1), and the restricted SCM can be adopted directly."))

  cls <- c(imbalance = .classify_ratio(imb_ratio, comparable_tol, worse_tol),
           rmspe = .classify_ratio(rmspe_ratio, comparable_tol, worse_tol))
  cls <- cls[!is.na(cls)]
  if (!any(cls == "worse")) {
    why <- "Steps 2-3: the restricted specification is not substantially worse than the unrestricted one, so the restricted SCM is preferable: it avoids the amplification of approximation errors (Section 4.4)."
    if (any(cls == "intermediate"))
      why <- paste(why, sprintf("Some ratios fall between %d%% and %d%%, a range Section 5.1 treats as neither comparable nor substantially worse: report the iSCM estimates alongside.",
                                round(100 * (comparable_tol - 1)), round(100 * (worse_tol - 1))))
    return(out("restricted", why))
  }

  if (rule == "implied_rmspe") {
    if (isTRUE(implied[["full"]] <= implied[["simplified"]]))
      return(out("iscm_full", sprintf("The restricted specification is substantially worse, so the iSCM correction is warranted. The implied pure-control counterfactual of the full iSCM fits the pre-intervention period better than that of the simplified iSCM (RMSPE %s vs %s).",
                 format(signif(implied[["full"]], 4)), format(signif(implied[["simplified"]], 4)))))
    return(out("iscm_simplified", sprintf("The restricted specification is substantially worse, so the iSCM correction is warranted. The implied pure-control counterfactual of the simplified iSCM fits the pre-intervention period better than that of the full iSCM (RMSPE %s vs %s).",
               format(signif(implied[["simplified"]], 4)), format(signif(implied[["full"]], 4)))))
  }

  if (isTRUE(amp < alpha_tol))
    return(out("iscm_full", sprintf("The restricted specification is substantially worse, so the iSCM correction is warranted; the amplification factor (%.3f) is below %s, so the full iSCM can be applied safely (Section 5.1, Step 4).",
               amp, format(alpha_tol))))
  aff_cls <- vapply(aff_ratios, .classify_ratio, character(1),
                    comparable_tol = comparable_tol, worse_tol = worse_tol)
  if (!any(aff_cls == "worse", na.rm = TRUE))
    return(out("iscm_simplified", sprintf("The restricted specification is substantially worse, so the iSCM correction is warranted; the amplification factor (%.3f) is not below %s, and the restricted fits of the potentially affected units remain adequate, so the simplified iSCM is the more reliable choice (Section 5.1, Step 4).",
               amp, format(alpha_tol))))
  out("none_reliable", sprintf("The restricted specification is substantially worse, but the amplification factor (%.3f) is not below %s and the restricted fits of the potentially affected units are also poor: no variant can be considered reliable. Reconsider the composition of M (Section 5.2).",
      amp, format(alpha_tol)))
}

.rec_label <- function(r) switch(r, restricted = "restricted SCM", iscm_full = "full iSCM",
                                 iscm_simplified = "simplified iSCM",
                                 none_reliable = "no variant can be considered reliable")

#' Apply the iSCM decision rule to summary figures
#'
#' The decision rule of Section 5.1 of Di Stefano and Mellace, separated from
#' any estimation, so that it can be applied to fits produced by any estimator,
#' including a fully customised \code{Synth} specification. See
#' \code{\link{iscm_compare}} for the rule.
#'
#' @param rmspe_ratio Ratio of the restricted to the unrestricted RMSPE of the main
#'   treated unit: on the validation window if a training/validation split was
#'   used (it then takes precedence, as in Section 5.1), otherwise in-sample.
#' @param imbalance_ratio Optional ratio of the restricted to the unrestricted
#'   predictor imbalance (Step 2).
#' @param alpha Amplification factor \eqn{\alpha} for \eqn{m = 2}.
#' @param sigma_min Smallest singular value of \eqn{\hat\Omega}; its inverse is the
#'   general amplification bound, used when \code{alpha} is \code{NA}.
#' @param affected_rmspe_ratio Optional ratios of the restricted to the
#'   unrestricted RMSPE of the potentially affected units, used when the
#'   amplification is not below \code{alpha_tol} to choose between the
#'   simplified iSCM and "no reliable variant".
#' @param substantial Logical: do the affected units receive substantial weight in
#'   the unrestricted fit (Step 1)? Defaults to \code{TRUE}.
#' @param rule,comparable_tol,worse_tol,alpha_tol See \code{\link{iscm_compare}}.
#' @param rmspe_full,rmspe_simplified Used with \code{rule = "implied_rmspe"}:
#'   pre-intervention RMSPE of the implied counterfactuals of the full and the
#'   simplified iSCM.
#'
#' @return A list of class \code{"iscm_decision"} with a \code{recommendation} in
#'   \code{c("restricted", "iscm_full", "iscm_simplified", "none_reliable")} and a
#'   \code{rationale}.
#'
#' @examples
#' # German reunification (Section 7): restricted RMSPE 358.09 against 93.08 on
#' # the validation window, imbalance 1.50 against 0.26, alpha = 1.16.
#' iscm_decision(rmspe_ratio = 358.09 / 93.08, imbalance_ratio = 1.50 / 0.26,
#'               alpha = 1.16)
#'
#' # Large amplification and poor restricted fits of the affected unit
#' iscm_decision(rmspe_ratio = 2, alpha = 1.8, affected_rmspe_ratio = 1.9)
#' @export
iscm_decision <- function(rmspe_ratio, imbalance_ratio = NA, alpha = NA,
                          sigma_min = NA, affected_rmspe_ratio = NA,
                          substantial = TRUE, rule = c("paper", "implied_rmspe"),
                          comparable_tol = 1.25, worse_tol = 1.5, alpha_tol = 1.25,
                          rmspe_full = NA, rmspe_simplified = NA) {
  rule <- match.arg(rule)
  amp <- if (!is.na(alpha)) alpha else if (!is.na(sigma_min)) 1 / sigma_min else NA
  if (rule == "paper" && is.na(amp))
    stop("Supply `alpha` (m = 2) or `sigma_min`.", call. = FALSE)
  if (rule == "implied_rmspe" && (is.na(rmspe_full) || is.na(rmspe_simplified)))
    stop("rule = \"implied_rmspe\" needs `rmspe_full` and `rmspe_simplified`.", call. = FALSE)
  out <- .iscm_decide(substantial, imbalance_ratio, rmspe_ratio, affected_rmspe_ratio,
                      amp, rule, c(full = rmspe_full, simplified = rmspe_simplified),
                      comparable_tol, worse_tol, alpha_tol)
  structure(out, class = "iscm_decision")
}

#' @export
print.iscm_decision <- function(x, ...) {
  cat("Recommendation:", .rec_label(x$recommendation), "\n")
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

# Internal: number of validation periods from the `validation` argument.
.n_validation <- function(validation, n_pre) {
  if (is.null(validation) || isFALSE(validation)) return(NULL)
  if (identical(validation, "auto")) validation <- 1 / 3
  if (!is.numeric(validation) || length(validation) != 1L || validation <= 0)
    stop("`validation` must be \"auto\", a share in (0, 1), a number of periods, or NULL.",
         call. = FALSE)
  k <- if (validation < 1) round(n_pre * validation) else as.integer(validation)
  if (k < 1L || k > n_pre - 2L)
    stop("The validation window must leave at least two training periods.", call. = FALSE)
  k
}

# Internal: assemble the "iscm_compare" object (shared by both front ends).
.compare_object <- function(treated, affected, aff_weights, substantial, tab, diag,
                            main_in_affected, implied, rule, validation, thresholds) {
  main <- tab[tab$unit == treated, ]
  use_val <- !is.null(main$rmspe_val_unres)
  imb_ratio <- if (!is.null(main$imbalance_unres))
    main$imbalance_restr / main$imbalance_unres else NA_real_
  rmspe_ratio <- if (use_val) main$rmspe_val_restr / main$rmspe_val_unres else
    main$rmspe_restr / main$rmspe_unres
  aff <- tab[tab$unit != treated, , drop = FALSE]
  aff_ratios <- if (use_val) aff$rmspe_val_restr / aff$rmspe_val_unres else
    aff$rmspe_restr / aff$rmspe_unres
  names(aff_ratios) <- aff$unit
  amp <- if (!is.na(diag$alpha)) diag$alpha else diag$amplification
  dec <- .iscm_decide(substantial, imb_ratio, rmspe_ratio, aff_ratios, amp, rule,
                      implied, thresholds$comparable_tol, thresholds$worse_tol,
                      thresholds$alpha_tol)
  structure(list(
    treated = treated, affected = affected,
    step1_affected_weights = aff_weights, step1_substantial = substantial,
    comparison = tab,
    amplification = amp, alpha = diag$alpha, sigma_min = diag$sigma_min,
    diag_dominant = diag$diag_dominant,
    main_weight_in_affected = main_in_affected,
    imbalance_ratio = imb_ratio, rmspe_ratio = rmspe_ratio,
    rmspe_ratio_in_sample = main$rmspe_restr / main$rmspe_unres,
    affected_rmspe_ratio = aff_ratios,
    implied_rmspe = c(implied, restricted = main$rmspe_restr),
    rule = rule, validation = validation, thresholds = thresholds,
    recommendation = dec$recommendation, rationale = dec$rationale
  ), class = "iscm_compare")
}

#' Decision rule from pre-computed fits
#'
#' Version of \code{\link{iscm_compare}} that applies the decision rule of
#' Section 5.1 to synthetic controls estimated externally (e.g. the exact
#' \code{Synth} specification of the empirical application), rather than
#' re-estimating them with a package backend. The user supplies, for the main
#' treated unit and each potentially affected unit, the full gap series of the
#' unrestricted and restricted fits and the unrestricted donor weights; predictor
#' imbalances and out-of-sample (validation) RMSPEs may be supplied as well, and
#' the validation RMSPEs then take precedence over the in-sample ones.
#'
#' @param gaps_unres,gaps_restr Named lists of full gap series (unrestricted and
#'   restricted fits), one element per focal unit, main treated first.
#' @param weights Named list of unrestricted donor-weight vectors, parallel to
#'   \code{gaps_unres}, used to build \eqn{\hat\Omega} and to read the weights on
#'   the affected units (Step 1).
#' @param t0 Last pre-intervention period.
#' @param times Optional period labels aligned with the gap series.
#' @param imbalance_unres,imbalance_restr Optional named numeric vectors of
#'   predictor imbalances per focal unit (Step 2): by default the sum of the
#'   absolute discrepancies, \code{sum(abs(X1 - X0 \%*\% w))}. When the
#'   predictors are on very different scales, standardise each discrepancy by
#'   the predictor's cross-sectional standard deviation before aggregating, e.g.
#'   \code{s <- apply(cbind(X1, X0), 1, sd); sum(abs(X1 - X0 \%*\% w) / s)} from
#'   \code{Synth::dataprep}, using the same \code{s} for both fits.
#' @param rmspe_val_unres,rmspe_val_restr Optional named numeric vectors of
#'   out-of-sample (validation) RMSPEs per focal unit (Step 3 variant).
#' @param weight_tol,comparable_tol,worse_tol,alpha_tol,rule See
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
#'   w <- scm_base(Y[pre, unit], Y[pre, donors], standardize = FALSE)
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
                              weight_tol = 0.05, comparable_tol = 1.25,
                              worse_tol = 1.5, alpha_tol = 1.25,
                              rule = c("paper", "implied_rmspe")) {
  rule <- match.arg(rule)
  focal <- names(gaps_unres)
  treated <- focal[1]; affected <- focal[-1]
  len <- length(gaps_unres[[1]])
  if (is.null(times)) times <- seq_len(len)
  pre <- if (t0 %in% times) times <= t0 else seq_along(times) <= t0

  cross <- lapply(focal, function(fi) weights[[fi]][names(weights[[fi]]) %in% focal])
  Omega <- build_omega(cross, focal)
  diag <- omega_diagnostics(Omega)

  rows <- lapply(focal, function(fi) {
    out <- data.frame(unit = fi,
      role = if (fi == treated) "treated" else "affected",
      stringsAsFactors = FALSE)
    if (!is.null(imbalance_unres)) out$imbalance_unres <- imbalance_unres[[fi]]
    if (!is.null(imbalance_restr)) out$imbalance_restr <- imbalance_restr[[fi]]
    out$rmspe_unres <- rmspe(gaps_unres[[fi]][pre])
    out$rmspe_restr <- rmspe(gaps_restr[[fi]][pre])
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

  E_unres <- vapply(focal, function(fi) as.numeric(gaps_unres[[fi]][pre]),
                    numeric(sum(pre)))
  E_restr <- vapply(focal, function(fi) as.numeric(gaps_restr[[fi]][pre]),
                    numeric(sum(pre)))
  implied <- .implied_pre_rmspe(Omega, matrix(E_unres, ncol = length(focal)),
                                matrix(E_restr, ncol = length(focal)), aff_weights)

  .compare_object(treated, affected, aff_weights, substantial, tab, diag,
                  main_in_affected, implied, rule,
                  validation = if (is.null(rmspe_val_unres)) NULL else "supplied",
                  thresholds = list(weight_tol = weight_tol, comparable_tol = comparable_tol,
                                    worse_tol = worse_tol, alpha_tol = alpha_tol))
}

#' Data-driven choice among restricted SCM and full/simplified iSCM
#'
#' Runs the decision rule of Section 5.1 of Di Stefano and Mellace and returns
#' a recommendation, together with every quantity the paper recommends
#' reporting.
#'
#' \strong{The rule.}
#' \describe{
#'   \item{Step 1}{Weights of the potentially affected units in the unrestricted
#'     fit of the main treated unit. If all are below \code{weight_tol}, the
#'     specifications are expected to agree (Remark 2.1): restricted SCM.}
#'   \item{Step 2}{Predictor imbalance of the unrestricted and restricted fits,
#'     aggregated as the sum of the absolute discrepancies (standardised if
#'     \code{standardize = TRUE}).}
#'   \item{Step 3}{Pre-intervention RMSPE of the two fits, in-sample and, by
#'     default, on a validation window: the weights are estimated on the first
#'     two thirds of the pre-intervention period and evaluated on the last third.
#'     A ratio below \code{comparable_tol} (25\%) means the fits are comparable; a
#'     ratio of at least \code{worse_tol} (50\%) means the restricted fit is
#'     substantially worse. When Steps 2 and 3 disagree, precedence goes to the
#'     imbalance and to the validation comparison, so the in-sample RMSPE is used
#'     only when no validation window is available.}
#'   \item{Decision}{If neither the imbalance nor the RMSPE of the restricted fit
#'     is substantially worse: restricted SCM. Otherwise the iSCM correction is
#'     warranted, and Step 4 chooses the variant.}
#'   \item{Step 4}{Amplification \eqn{\alpha} (\eqn{m = 2}) or
#'     \eqn{1/\sigma_{\min}(\hat\Omega)}. Below \code{alpha_tol} (1.25): full iSCM.
#'     Otherwise: simplified iSCM if the restricted fits of the potentially
#'     affected units are not substantially worse than their unrestricted fits;
#'     if they are, no variant can be considered reliable and the composition of
#'     \eqn{M} should be reconsidered (Section 5.2).}
#' }
#' Steps 2-3 are computed for each potentially affected unit as well, treating it
#' as if it were the main treated unit (\code{comparison} table).
#'
#' Ratios between 25\% and 50\% are neither comparable nor substantially worse in
#' the terms of Section 5.1; the rule then keeps the restricted SCM and says so in
#' the rationale. The recommendation selects the specification to emphasise:
#' Section 5.1 recommends reporting all four (unrestricted, restricted, full and
#' simplified iSCM).
#'
#' \strong{Alternative Step 4.} With \code{rule = "implied_rmspe"}, the choice
#' between the full and the simplified iSCM compares instead the
#' pre-intervention RMSPE of their implied pure-control counterfactuals, the
#' observable counterparts of equations (4.6) and (4.7). These values are always
#' reported in \code{implied_rmspe}, as an additional diagnostic.
#'
#' @inheritParams iscm
#' @param validation Validation window for Step 3: \code{"auto"} (default, the last
#'   third of the pre-intervention period, as in Section 5.1), a share in (0, 1),
#'   a number of final pre-intervention periods, or \code{NULL} for the in-sample
#'   comparison only.
#' @param weight_tol Weight below which a potentially affected unit is treated as
#'   receiving negligible weight (Step 1). Default \code{0.05}.
#' @param comparable_tol,worse_tol Restricted / unrestricted ratios below which
#'   the fits are comparable (default \code{1.25}) and from which the restricted
#'   fit is substantially worse (default \code{1.5}), as in Section 5.1.
#' @param alpha_tol Amplification below which the full iSCM can be applied safely.
#'   Default \code{1.25}, as in Section 5.1. For \eqn{m > 2} the general bound
#'   \code{1/sigma_min} is used in place of \code{alpha}.
#' @param standardize Logical; if \code{TRUE}, each predictor's absolute
#'   discrepancy is divided by the predictor's cross-sectional standard deviation
#'   (computed across all units and held fixed across the unrestricted and
#'   restricted fits) before the Step 2 imbalances are aggregated. Default
#'   \code{FALSE}, the plain sum of Section 5.1; set it to \code{TRUE} when the
#'   predictors are measured on very different scales, as in the German
#'   application.
#' @param rule \code{"paper"} (default) for the rule of Section 5.1, or
#'   \code{"implied_rmspe"} for the alternative Step 4 described above.
#'
#' @return An object of class \code{"iscm_compare"} with the per-step comparison
#'   table, the diagnostics, and a \code{recommendation} in
#'   \code{c("restricted", "iscm_full", "iscm_simplified", "none_reliable")} with a
#'   rationale.
#'
#' @seealso \code{\link{iscm_compare_fits}} and \code{\link{iscm_decision}} for
#'   fits estimated outside the package; \code{\link{iscm_sensitivity}} for the
#'   choice of the affected units.
#'
#' @examples
#' sim <- simulate_iscm(n_affected = 1, n_pure = 6, tau = 2, spill = 3, seed = 7)
#' cmp <- iscm_compare(sim$data, "unit", "time", "Y",
#'                     treated = sim$treated, affected = sim$affected,
#'                     t0 = sim$t0)
#' cmp                       # Steps 1-4 and the recommendation
#' cmp$recommendation
#' cmp$comparison            # imbalance and RMSPE for every focal unit
#' @export
iscm_compare <- function(data, id, time, outcome, treated, affected, t0,
                         predictors = NULL, pure = NULL, estimator = scm_base,
                         v = NULL, zvar = NULL, validation = "auto",
                         weight_tol = 0.05, comparable_tol = 1.25, worse_tol = 1.5,
                         alpha_tol = 1.25, standardize = FALSE,
                         rule = c("paper", "implied_rmspe"), ...) {
  rule <- match.arg(rule)
  pm <- .panel_matrix(data, id, time, outcome)
  ids <- pm$ids; times <- pm$times; Y <- pm$Y
  treated <- as.character(treated); affected <- as.character(affected)
  focal <- c(treated, affected)
  if (is.null(pure)) pure <- setdiff(ids, focal) else pure <- as.character(pure)
  pre <- times[times <= t0]
  n_val <- .n_validation(validation, length(pre))
  Xmat <- .predictor_matrix(data, id, time, outcome, predictors, ids, times, pre)
  zcol <- if (is.null(zvar)) outcome else zvar
  Zmat <- .panel_matrix(data, id, time, zcol)$Y[which(times %in% pre), , drop = FALSE]

  # Step 2 predictor imbalance: sum of absolute discrepancies, optionally
  # standardised by each predictor's cross-sectional standard deviation, held
  # fixed across the unrestricted and restricted fits.
  if (standardize) {
    sdv <- apply(Xmat, 1L, stats::sd)
    sdv[sdv == 0 | !is.finite(sdv)] <- 1
  } else {
    sdv <- rep(1, nrow(Xmat))
  }
  fit_unit <- function(focal_id, donor_ids)
    .fit_one(focal_id, donor_ids, Xmat, Y, Zmat, v, estimator, ...)
  imbalance <- function(focal_id, donor_ids, w)
    sum(abs(Xmat[, focal_id] - Xmat[, donor_ids, drop = FALSE] %*% w) / sdv)

  # ---- the iSCM fit for Omega / amplification (Step 4) ------------------------
  fit <- iscm(data, id, time, outcome, treated, affected, t0,
              predictors = predictors, pure = pure, estimator = estimator,
              v = v, zvar = zvar, ...)

  # ---- per-unit comparison (main treated first, then each affected) ----------
  if (!is.null(n_val)) {
    valid <- utils::tail(pre, n_val); train <- setdiff(pre, valid)
  }
  gaps_re <- list()
  rows <- lapply(focal, function(fi) {
    donors_un <- setdiff(ids, fi)
    f_un <- fit_unit(fi, donors_un); f_re <- fit_unit(fi, pure)
    gaps_re[[fi]] <<- f_re$gap
    out <- data.frame(
      unit = fi,
      role = if (fi == treated) "treated" else "affected",
      imbalance_unres = imbalance(fi, donors_un, f_un$w),
      imbalance_restr = imbalance(fi, pure, f_re$w),
      rmspe_unres = rmspe(f_un$gap[times <= t0]),
      rmspe_restr = rmspe(f_re$gap[times <= t0]),
      stringsAsFactors = FALSE)
    if (!is.null(n_val)) {
      out$rmspe_val_unres <- .refit_window(fi, donors_un, data, id, time, outcome,
        predictors, ids, times, train, valid, Y, v, estimator, zvar, ...)
      out$rmspe_val_restr <- .refit_window(fi, pure, data, id, time, outcome,
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
  # simplified iSCM is costless for it, Corollary 4.2).
  main_in_affected <- vapply(affected, function(a) {
    wa <- fit$fits[[match(a, focal)]]$w
    if (treated %in% names(wa)) wa[[treated]] else 0
  }, numeric(1))

  # ---- additional Step 4 diagnostic: implied counterfactuals -----------------
  pre_idx <- which(times <= t0)
  E_unres <- fit$gaps_unrestricted[pre_idx, focal, drop = FALSE]
  E_restr <- vapply(focal, function(fi) gaps_re[[fi]][pre_idx],
                    numeric(length(pre_idx)))
  implied <- .implied_pre_rmspe(fit$Omega, E_unres,
                                matrix(E_restr, ncol = length(focal)), aff_weights)

  .compare_object(treated, affected, aff_weights, substantial, tab, fit$diagnostics,
                  main_in_affected, implied, rule,
                  validation = if (is.null(n_val)) NULL else
                    list(train = range(train), valid = range(valid)),
                  thresholds = list(weight_tol = weight_tol, comparable_tol = comparable_tol,
                                    worse_tol = worse_tol, alpha_tol = alpha_tol))
}

#' @export
print.iscm_compare <- function(x, digits = 3, ...) {
  th <- x$thresholds
  cat("iSCM specification comparison (Section 5.1)\n")
  cat("  Main treated unit:", x$treated, "\n\n")
  cat("Step 1 - weights on affected units in the unrestricted main fit:\n")
  print(round(x$step1_affected_weights, digits))
  if (!x$step1_substantial)
    cat("  -> all below weight_tol: affected units carry little contamination.\n")
  cat("\nSteps 2-3 - predictor imbalance and pre-intervention RMSPE:\n")
  print(format(x$comparison, digits = digits))
  if (is.list(x$validation))
    cat(sprintf("  (validation: weights estimated on %s-%s, RMSPE evaluated on %s-%s)\n",
                x$validation$train[1], x$validation$train[2],
                x$validation$valid[1], x$validation$valid[2]))
  cat(sprintf("\nRatios restricted / unrestricted, main unit: imbalance = %s, RMSPE%s = %.2f\n",
              if (is.na(x$imbalance_ratio)) "n/a" else sprintf("%.2f", x$imbalance_ratio),
              if (is.null(x$validation)) "" else " (validation)", x$rmspe_ratio))
  cat(sprintf("  comparable below %.2f, substantially worse from %.2f\n",
              th$comparable_tol, th$worse_tol))
  cat(sprintf("\nStep 4 - amplification: %s = %.3f (full iSCM below %s) | sigma_min = %.3f\n",
              if (!is.na(x$alpha)) "alpha" else "1/sigma_min",
              x$amplification, format(th$alpha_tol), x$sigma_min))
  cat("  Restricted / unrestricted RMSPE of the affected units:",
      paste(sprintf("%s = %.2f", names(x$affected_rmspe_ratio), x$affected_rmspe_ratio),
            collapse = ", "), "\n")
  ir <- x$implied_rmspe
  cat(sprintf("  Pre-intervention RMSPE of the implied counterfactual: full iSCM = %s | simplified iSCM = %s | restricted SCM = %s\n",
              format(signif(ir[["full"]], 4)), format(signif(ir[["simplified"]], 4)),
              format(signif(ir[["restricted"]], 4))))
  if (identical(x$rule, "implied_rmspe"))
    cat("  (rule = \"implied_rmspe\": Step 4 compares the implied counterfactuals)\n")
  cat("\nRecommendation:", .rec_label(x$recommendation), "\n")
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
