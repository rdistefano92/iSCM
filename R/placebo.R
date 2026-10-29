#' In-space placebo test for iSCM
#'
#' Implements the placebo (permutation) inference of Online Appendix B of Di
#' Stefano and Mellace. Because the observed post-intervention outcomes of the
#' focal units embed the effects of the intervention, the test statistics are
#' computed from outcomes purged of the estimated iSCM effects. The statistic is
#' the ratio of the post- to the pre-intervention RMSPE; the intervention is then
#' reassigned to each pure control in turn and the same statistic is recomputed,
#' yielding a placebo distribution and a rank-based p-value.
#'
#' Run it once for the main treated unit and once for each affected unit
#' (\code{target}). As discussed in Online Appendix B, treat the result as a
#' transparent diagnostic rather than a test with formal size guarantees: with
#' \eqn{|C|} pure controls the smallest attainable p-value is \eqn{1/(|C|+1)},
#' and for large \eqn{\alpha} the statistic of the target is tilted upwards.
#' Printing the result reports both, together with the conditioning of
#' \eqn{\hat\Omega}, as the paper recommends.
#'
#' \code{iscm_placebo()} re-estimates the placebo synthetic controls with the
#' backend stored in \code{object}, so \code{object} should come from
#' \code{\link{iscm}}. For synthetic controls estimated outside the package
#' (objects built with \code{\link{iscm_assemble}}), use
#' \code{\link{iscm_purge}} and \code{\link{iscm_placebo_fits}}.
#'
#' @param object An object of class \code{"iscm"}.
#' @param data,id,time,outcome As in \code{\link{iscm}} (the same data used to fit
#'   \code{object}).
#' @param target The focal unit whose effect is being tested. Defaults to the main
#'   treated unit.
#' @param predictors As in \code{\link{iscm}}. Defaults to the value stored in
#'   \code{object}.
#' @param estimator,v,zvar,... As in \code{\link{iscm}}. By default the backend,
#'   predictor weights, outcome-fit variable and backend arguments stored in
#'   \code{object} are reused. Explicit values override the stored specification.
#'
#' @return An object of class \code{"iscm_placebo"}: a data frame
#'   (\code{table}) of per-unit \code{pre_rmspe}, \code{post_rmspe} and
#'   \code{ratio}, plus the target's \code{ratio}, \code{rank}, rank-based
#'   \code{p_value}, the smallest attainable p-value \code{min_p_value} and the
#'   conditioning diagnostics of \code{object} (\code{alpha}, \code{sigma_min},
#'   \code{amplification}), and the gaps of the target and of every placebo
#'   (\code{gaps}). \code{\link[=plot.iscm_placebo]{plot()}} draws the ratios
#'   as in Figure C.1 of the paper, or the gap paths.
#'
#' @seealso \code{\link{iscm_placebo_fits}} for synthetic controls estimated
#'   outside the package; \code{\link{iscm}}.
#'
#' @examples
#' sim <- simulate_iscm(n_affected = 1, n_pure = 6, tau = 2, spill = 3, seed = 7)
#' fit <- iscm(sim$data, "unit", "time", "Y", sim$treated, sim$affected, sim$t0)
#'
#' pb_tr <- iscm_placebo(fit, sim$data, "unit", "time", "Y")     # main unit
#' pb_tr
#' plot(pb_tr)
#'
#' pb_a1 <- iscm_placebo(fit, sim$data, "unit", "time", "Y", target = "A1")
#' pb_a1$table
#' @export
iscm_placebo <- function(object, data, id, time, outcome, target = NULL,
                         predictors = NULL, estimator = NULL, v = NULL,
                         zvar = NULL, ...) {
  stopifnot(inherits(object, "iscm"))
  target <- .placebo_target(object, target)
  if (is.null(predictors)) predictors <- object$predictors

  supplied_args <- list(...)
  if (!is.null(object$backend)) {
    if (is.null(estimator)) estimator <- object$backend$estimator
    if (is.null(v)) v <- object$backend$v
    if (is.null(zvar)) zvar <- object$backend$zvar
    estimator_args <- utils::modifyList(object$backend$args, supplied_args)
  } else {
    if (is.null(estimator)) {
      stop("`object` was not produced by iscm(), so its estimator is unknown. ",
           "Estimate the placebo synthetic controls with your own estimator on ",
           "the data returned by iscm_purge() and pass their weights to ",
           "iscm_placebo_fits().", call. = FALSE)
    }
    estimator_args <- supplied_args
  }

  pm <- .panel_matrix(data, id, time, outcome)
  ids <- pm$ids; times <- pm$times; Y <- pm$Y
  pre <- times <= object$t0
  Xmat <- .predictor_matrix(data, id, time, outcome, predictors, ids, times, times[pre])
  zcol <- if (is.null(zvar)) outcome else zvar
  Zmat <- .panel_matrix(data, id, time, zcol)$Y[which(pre), , drop = FALSE]

  # Placebo synthetic controls: each pure control on all other units. iscm()
  # builds predictors from pre-intervention periods only, so these weights do
  # not depend on the purged post-intervention outcomes.
  placebo_weights <- lapply(object$pure, function(j) {
    do.call(.fit_one, c(list(
      focal_id = j, donor_ids = setdiff(ids, j), Xmat = Xmat, Ymat = Y,
      Zmat = Zmat, v = v, estimator = estimator
    ), estimator_args))$w
  })
  names(placebo_weights) <- object$pure

  .placebo_core(object, Y, times, target, placebo_weights)
}

#' Purge the estimated effects from the outcomes of the focal units
#'
#' Returns a copy of \code{data} in which the post-intervention outcomes of the
#' main treated unit and of the potentially affected units are replaced by the
#' observed outcomes minus their iSCM estimates, \eqn{Y_{1t} - \hat\theta^{iSCM}_{1t}}
#' and \eqn{Y_{it} - \hat\gamma^{iSCM}_{it}}. These are the corrected outcomes on
#' which the placebo analysis of Online Appendix B is computed.
#'
#' Use it when the synthetic controls were estimated outside the package: estimate
#' the synthetic control of each pure control on the purged data, with your own
#' estimator and specification, and pass the weights to
#' \code{\link{iscm_placebo_fits}}. Using the purged data matters whenever the
#' specification uses periods after \code{t0}, for example predictors averaged
#' over a window that includes the intervention year.
#'
#' @param object An object of class \code{"iscm"}, from \code{\link{iscm}} or
#'   \code{\link{iscm_assemble}}.
#' @param data,id,time,outcome The long panel and its column names. The unit
#'   identifiers in \code{id} must match the focal-unit names used in
#'   \code{object}.
#' @param exclude Optional focal units whose outcomes are left untouched.
#'
#' @return \code{data} with the purged outcome in column \code{outcome}.
#'
#' @seealso \code{\link{iscm_placebo_fits}}.
#'
#' @examples
#' sim <- simulate_iscm(seed = 7)
#' fit <- iscm(sim$data, "unit", "time", "Y", sim$treated, sim$affected, sim$t0)
#' purged <- iscm_purge(fit, sim$data, "unit", "time", "Y")
#' # post-intervention outcomes of TR: observed and purged of the iSCM effect
#' post_tr <- sim$data$unit == "TR" & sim$data$time > sim$t0
#' cbind(observed = sim$data$Y[post_tr], purged = purged$Y[post_tr],
#'       effect = fit$effects["TR", ])
#' @export
iscm_purge <- function(object, data, id, time, outcome, exclude = character(0)) {
  stopifnot(inherits(object, "iscm"))
  units <- as.character(data[[id]])
  if (!all(object$focal %in% units))
    stop("The unit identifiers in `data[[id]]` must include the focal units of ",
         "`object`: ", paste(object$focal, collapse = ", "), ".", call. = FALSE)
  post_times <- as.character(object$post)
  for (u in setdiff(object$focal, as.character(exclude))) {
    eff <- object$effects[match(u, object$focal), ]
    rows <- which(units == u & as.character(data[[time]]) %in% post_times)
    k <- match(as.character(data[[time]][rows]), post_times)
    data[[outcome]][rows] <- data[[outcome]][rows] - eff[k]
  }
  data
}

#' Placebo test for iSCM from synthetic controls estimated externally
#'
#' Version of \code{\link{iscm_placebo}} for synthetic controls estimated outside
#' the package, for example with the exact \code{Synth} specification of the
#' paper. The user supplies, for each pure control, the donor weights of its own
#' synthetic control (the intervention reassigned to it, all other units as
#' donors), estimated on the purged data returned by \code{\link{iscm_purge}}.
#' The function computes the statistics on the purged outcomes, exactly as
#' \code{iscm_placebo()} does.
#'
#' @param object An object of class \code{"iscm"}, typically from
#'   \code{\link{iscm_assemble}}. Its \code{fits} must contain the donor weights
#'   of the focal units.
#' @param data,id,time,outcome The original (not purged) long panel and its column
#'   names. The unit identifiers must match the names used in \code{object} and in
#'   \code{placebo_weights}.
#' @param placebo_weights A named list with one element per pure control: the
#'   named vector of donor weights of that unit's synthetic control.
#' @param target The focal unit whose effect is being tested. Defaults to the main
#'   treated unit.
#'
#' @return An object of class \code{"iscm_placebo"}; see \code{\link{iscm_placebo}}.
#'
#' @seealso \code{\link{iscm_purge}}, \code{\link{iscm_placebo}},
#'   \code{\link{iscm_assemble}}.
#'
#' @examples
#' # Synthetic controls estimated "outside" the package, here with scm_base()
#' sim <- simulate_iscm(n_affected = 1, n_pure = 6, tau = 2, spill = 3, seed = 7)
#' Y <- with(sim$data, tapply(Y, list(time, unit), sum))
#' pre <- seq_len(sim$t0)
#' sc <- function(unit, Ymat) {
#'   donors <- setdiff(colnames(Ymat), unit)
#'   w <- scm_base(Ymat[pre, unit], Ymat[pre, donors], standardize = FALSE)
#'   list(w = w, gap = as.numeric(Ymat[, unit] - Ymat[, donors] %*% w))
#' }
#' f_tr <- sc("TR", Y); f_a1 <- sc("A1", Y)
#' fit <- iscm_assemble(gaps = list(TR = f_tr$gap, A1 = f_a1$gap),
#'                      weights = list(TR = f_tr$w, A1 = f_a1$w),
#'                      t0 = sim$t0, pure = sim$pure)
#'
#' # 1. purge the estimated effects; 2. fit each pure control on the purged data
#' purged <- iscm_purge(fit, sim$data, "unit", "time", "Y")
#' Yp <- with(purged, tapply(Y, list(time, unit), sum))
#' pw <- lapply(setNames(sim$pure, sim$pure), function(j) sc(j, Yp)$w)
#'
#' # 3. compute the placebo statistics
#' iscm_placebo_fits(fit, sim$data, "unit", "time", "Y", placebo_weights = pw)
#' @export
iscm_placebo_fits <- function(object, data, id, time, outcome, placebo_weights,
                              target = NULL) {
  stopifnot(inherits(object, "iscm"))
  target <- .placebo_target(object, target)
  if (!length(object$pure))
    stop("`object` has no pure controls; supply `pure` to iscm_assemble().",
         call. = FALSE)
  if (!is.list(placebo_weights) || is.null(names(placebo_weights)))
    stop("`placebo_weights` must be a named list, one element per pure control.",
         call. = FALSE)
  missing_pure <- setdiff(object$pure, names(placebo_weights))
  if (length(missing_pure))
    stop("`placebo_weights` has no entry for: ",
         paste(missing_pure, collapse = ", "), ".", call. = FALSE)

  pm <- .panel_matrix(data, id, time, outcome)
  if (length(pm$times) != length(object$times) ||
      any(as.character(pm$times) != as.character(object$times)))
    stop("The periods in `data` must coincide with those of `object`.",
         call. = FALSE)
  .placebo_core(object, pm$Y, pm$times, target, placebo_weights[object$pure])
}

# Internal: validate the target unit.
.placebo_target <- function(object, target) {
  if (is.null(target)) target <- object$treated
  target <- as.character(target)
  if (!target %in% object$focal) stop("`target` must be a focal unit.", call. = FALSE)
  target
}

# Internal: placebo statistics of Online Appendix B, shared by iscm_placebo()
# and iscm_placebo_fits(). Y is the T x J matrix of observed outcomes (columns
# named by unit); placebo_weights holds the donor weights of each pure control.
.placebo_core <- function(object, Y, times, target, placebo_weights) {
  pre <- times <= object$t0; post <- !pre
  check_donors <- function(w, who) {
    miss <- setdiff(names(w), colnames(Y))
    if (is.null(names(w)) || length(miss))
      stop(sprintf("The donor weights of %s must be named by units present in the data%s.",
                   who, if (length(miss)) paste0("; not found: ",
                                                 paste(miss, collapse = ", ")) else ""),
           call. = FALSE)
  }

  # outcomes purged of the estimated iSCM effects (focal units, post period)
  purify <- function(exclude = character(0)) {
    Yp <- Y
    for (u in setdiff(object$focal, exclude)) {
      Yp[which(post), u] <- Y[which(post), u] - object$effects[match(u, object$focal), ]
    }
    Yp
  }

  # --- statistic for the target, with its own estimated weights ----------------
  # Purging the other focal units makes the post-intervention gap equal to the
  # iSCM estimate (Online Appendix B); the denominator is the pre-intervention
  # RMSPE of the target's unrestricted fit.
  tf <- object$fits[[match(target, object$focal)]]
  check_donors(tf$w, target)
  Yp_t <- purify(exclude = target)
  gap_t <- as.numeric(Y[, target]) -
    as.numeric(Yp_t[, names(tf$w), drop = FALSE] %*% tf$w)
  pre_t <- rmspe(as.numeric(tf$gap)[pre])
  r_target <- rmspe(gap_t[post]) / pre_t

  # --- placebo statistics: each pure control as treated -------------------------
  Yp_all <- purify()
  gaps <- matrix(NA_real_, length(times), length(placebo_weights) + 1L,
                 dimnames = list(as.character(times), c(target, names(placebo_weights))))
  gaps[, target] <- gap_t
  rows <- lapply(names(placebo_weights), function(j) {
    w <- placebo_weights[[j]]
    check_donors(w, j)
    if (j %in% names(w)) stop(sprintf("The weights of %s include %s itself.", j, j),
                              call. = FALSE)
    gap <- as.numeric(Yp_all[, j]) - as.numeric(Yp_all[, names(w), drop = FALSE] %*% w)
    gaps[, j] <<- gap
    data.frame(unit = j, pre_rmspe = rmspe(gap[pre]), post_rmspe = rmspe(gap[post]),
               ratio = rmspe(gap[post]) / rmspe(gap[pre]), stringsAsFactors = FALSE)
  })
  tab <- rbind(
    data.frame(unit = target, pre_rmspe = pre_t, post_rmspe = rmspe(gap_t[post]),
               ratio = r_target, stringsAsFactors = FALSE),
    do.call(rbind, rows)
  )
  placebo_ratios <- tab$ratio[-1]
  rank_target <- 1L + sum(placebo_ratios > r_target)
  p_value <- (1 + sum(placebo_ratios >= r_target)) / (length(placebo_ratios) + 1)
  tab <- tab[order(-tab$ratio), ]
  rownames(tab) <- NULL

  d <- object$diagnostics
  structure(list(table = tab, target = target, ratio = r_target,
                 rank = rank_target, n = nrow(tab), p_value = p_value,
                 min_p_value = 1 / (length(placebo_ratios) + 1),
                 alpha = d$alpha, sigma_min = d$sigma_min,
                 amplification = d$amplification,
                 gaps = gaps, times = times, t0 = object$t0),
            class = "iscm_placebo")
}

#' @export
print.iscm_placebo <- function(x, digits = 3, ...) {
  cat("iSCM in-space placebo test (Online Appendix B)\n")
  cat("  Target unit:", x$target, "\n")
  cat(sprintf("  RMSPE ratio: %.3f (rank %d of %d)\n", x$ratio, x$rank, x$n))
  cat(sprintf("  Placebo p-value: %.3f\n", x$p_value))
  if (!is.null(x$min_p_value)) {
    cat("\nTo report alongside the p-value:\n")
    cat(sprintf("  Smallest attainable p-value with %d pure controls: %.3f\n",
                x$n - 1L, x$min_p_value))
    if (!is.null(x$sigma_min)) {
      cat(sprintf("  Conditioning of Omega: %ssigma_min = %.3f (amplification 1/sigma_min = %.3f)\n",
                  if (!is.null(x$alpha) && !is.na(x$alpha))
                    sprintf("alpha = %.3f, ", x$alpha) else "",
                  x$sigma_min, x$amplification))
      cat(strwrap(paste(
        "The target's statistic amplifies the fitting errors by the conditioning",
        "of Omega, unlike the placebos, so it is tilted upwards when alpha is far",
        "from one. Read the result as a diagnostic, not as a test with exact size."),
        indent = 2, exdent = 2), sep = "\n")
    }
  }
  invisible(x)
}
