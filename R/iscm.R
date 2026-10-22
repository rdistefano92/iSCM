# Internal: reshape long panel into a T x J outcome matrix (times x units) --------
.panel_matrix <- function(data, id, time, outcome) {
  ids <- unique(data[[id]])
  times <- sort(unique(data[[time]]))
  Y <- matrix(NA_real_, nrow = length(times), ncol = length(ids),
              dimnames = list(as.character(times), as.character(ids)))
  idx_t <- match(as.character(data[[time]]), rownames(Y))
  idx_u <- match(as.character(data[[id]]), colnames(Y))
  Y[cbind(idx_t, idx_u)] <- data[[outcome]]
  if (anyNA(Y)) stop("The panel is unbalanced: some unit-time cells are missing.")
  list(Y = Y, ids = as.character(ids), times = times)
}

# Internal: predictor matrix (k x J), one column per unit ------------------------
# If predictors is NULL, uses the pre-intervention outcome path (Y-only matching).
.predictor_matrix <- function(data, id, time, outcome, predictors, ids, times, pre) {
  if (is.null(predictors)) {
    Ywide <- .panel_matrix(data, id, time, outcome)$Y
    X <- Ywide[which(times %in% pre), , drop = FALSE]
    rownames(X) <- paste0("Y", times[times %in% pre])
    X <- X[, ids, drop = FALSE]
    attr(X, "outcome_path") <- TRUE
    return(X)
  }
  # aggregate each predictor over the pre-intervention window (mean) per unit
  d <- data[data[[time]] %in% pre, , drop = FALSE]
  X <- vapply(predictors, function(p) {
    tapply(d[[p]], factor(as.character(d[[id]]), levels = ids), mean, na.rm = TRUE)
  }, numeric(length(ids)))
  X <- t(X)                      # k x J
  colnames(X) <- ids
  X
}

# Internal: is `estimator` one of the package's own backends?
.is_builtin_backend <- function(estimator) {
  identical(estimator, scm_base) || identical(estimator, scm_penalized)
}

# Internal: fit one SC-type synthetic control -----------------------------------
# Returns weights (named over donors) and the full-series counterfactual & gap.
# When matching on the pre-intervention outcome path, every period is a
# predictor; standardising each one by its cross-sectional SD would give the
# most weight to the periods in which units are most alike (where the noise
# dominates), so the built-in backends do not standardise in that case unless
# the user asks for it explicitly.
.fit_one <- function(focal_id, donor_ids, Xmat, Ymat, Zmat, v, estimator, ...) {
  X1 <- Xmat[, focal_id]
  X0 <- Xmat[, donor_ids, drop = FALSE]
  Z1 <- if (is.null(Zmat)) NULL else Zmat[, focal_id]
  Z0 <- if (is.null(Zmat)) NULL else Zmat[, donor_ids, drop = FALSE]
  args <- list(...)
  if (isTRUE(attr(Xmat, "outcome_path")) && is.null(args$standardize) &&
      .is_builtin_backend(estimator)) {
    args$standardize <- FALSE
  }
  w <- do.call(estimator, c(list(X1 = X1, X0 = X0, v = v, Z1 = Z1, Z0 = Z0), args))
  w <- as.numeric(w)
  names(w) <- donor_ids
  yhat <- as.numeric(Ymat[, donor_ids, drop = FALSE] %*% w)
  gap  <- as.numeric(Ymat[, focal_id]) - yhat
  list(w = w, yhat = yhat, gap = gap)
}

#' Inclusive synthetic control method
#'
#' Estimates the effect of an intervention on a main treated unit while keeping
#' units potentially affected by the intervention (other treated units or units
#' exposed to spillovers) in the donor pool, following Di Stefano and Mellace. A
#' synthetic control is constructed for the main treated unit and for every
#' potentially affected unit; the resulting contaminated estimates satisfy a
#' linear system whose solution removes the contamination exactly (Theorem 2.1).
#'
#' This is the main entry point of the package. It carries out the four
#' implementation steps of Section 2.3: (1) fit the main treated unit on all other
#' units, (2) fit every potentially affected unit on all other units, (3) build
#' \eqn{\hat\Omega}, (4) solve the system. With \code{simplified = TRUE} the
#' affected units are fit on the pure controls only (Section 4.2).
#'
#' The function covers specifications in which every unit is matched on the same
#' predictors averaged over the whole pre-intervention period (or on the
#' pre-intervention outcome path). For richer specifications, such as
#' \code{Synth::dataprep} with special predictors and different windows, estimate
#' the fits yourself and use \code{\link{iscm_assemble}}.
#'
#' @param data A long-format \code{data.frame} with one row per unit-time.
#' @param id,time,outcome Column names (character) identifying the unit, the time
#'   period and the outcome.
#' @param treated Identifier of the main treated unit (unit 1).
#' @param affected Vector of identifiers of the potentially affected units (the set
#'   \eqn{M}). May be empty, in which case iSCM reduces to standard SCM.
#' @param t0 Last pre-intervention period; the intervention is active from
#'   \code{t0 + 1} onward.
#' @param predictors Optional character vector of predictor columns; each is
#'   averaged over the pre-intervention periods for every unit. If \code{NULL}
#'   (default), the pre-intervention outcome path is used for matching; in that
#'   case the built-in backends do not standardise the periods (pass
#'   \code{standardize = TRUE} to override).
#' @param pure Optional vector of pure-control identifiers. Defaults to every unit
#'   that is neither treated nor affected.
#' @param estimator A backend function; see \code{\link{backends}}. Defaults to
#'   \code{\link{scm_base}}.
#' @param v Optional predictor-weighting vector passed to the backend.
#' @param simplified Logical; if \code{TRUE}, computes the simplified iSCM, which
#'   fits each affected unit on the pure controls only and avoids the matrix
#'   inversion (Section 4.2). More robust when the cross-weights are large.
#' @param zvar Optional column name for a separate outcome-fit variable used only by
#'   \code{\link{scm_synth}} as \code{Z}; defaults to \code{outcome}.
#' @param ... Further arguments passed to the backend.
#'
#' @return An object of class \code{"iscm"}. Use \code{print()},
#'   \code{\link[=summary.iscm]{summary()}} and \code{\link[=plot.iscm]{plot()}}
#'   on it. Its most useful components are:
#'   \describe{
#'     \item{\code{effects}}{\code{m x P} matrix of iSCM estimates; row 1 is the
#'       effect on the main treated unit (\eqn{\hat\theta_{1t}}), the other rows
#'       the effects on the affected units (\eqn{\hat\gamma_{it}}); columns are
#'       post-intervention periods.}
#'     \item{\code{gaps_iscm}, \code{gaps_unrestricted}}{\code{T x m} matrices of
#'       gaps over all periods; the two coincide before \code{t0}.
#'       \code{gaps_unrestricted} is the uncorrected (contaminated) SCM.}
#'     \item{\code{Omega}, \code{diagnostics}}{The weight matrix and its
#'       conditioning diagnostics (see \code{\link{omega_diagnostics}}).}
#'     \item{\code{fits}}{For each focal unit, the donor weights \code{w}, the
#'       synthetic outcome \code{yhat} and the \code{gap}.}
#'     \item{\code{implied_weights}}{Weights of the implied pure-control synthetic
#'       control of the main treated unit (Lemma 3.1).}
#'     \item{\code{beta_post}}{The contaminated post-intervention gaps
#'       \eqn{\hat\beta_t}.}
#'   }
#'
#' @references Di Stefano, R. and Mellace, G. The inclusive synthetic control
#'   method.
#'
#' @seealso \code{\link{restricted_scm}} for the benchmark without the affected
#'   units, \code{\link{iscm_compare}} to choose among the specifications,
#'   \code{\link{iscm_placebo}} for inference, \code{\link{iscm_assemble}} for
#'   synthetic controls estimated outside the package.
#'
#' @examples
#' # Simulated panel: true effect on TR = 2, spillover on A1 = 3
#' sim <- simulate_iscm(n_affected = 1, n_pure = 6, tau = 2, spill = 3, seed = 7)
#'
#' # Full iSCM
#' fit <- iscm(sim$data, id = "unit", time = "time", outcome = "Y",
#'             treated = sim$treated, affected = sim$affected, t0 = sim$t0)
#' fit                         # average effects and Omega diagnostics
#' summary(fit)                # per-period effects
#' fit$effects                 # the estimates as a matrix
#' plot(fit)                   # inclusive vs unrestricted gaps
#'
#' # Simplified iSCM: the affected unit is fit on the pure controls only
#' fit_s <- iscm(sim$data, "unit", "time", "Y", sim$treated, sim$affected,
#'               sim$t0, simplified = TRUE)
#' rowMeans(fit_s$effects)
#'
#' \donttest{
#' # Classic SCM of Abadie et al. through the Synth package
#' if (requireNamespace("Synth", quietly = TRUE)) {
#'   iscm(sim$data, "unit", "time", "Y", sim$treated, sim$affected, sim$t0,
#'        estimator = scm_synth)
#' }
#' }
#' @export
iscm <- function(data, id, time, outcome, treated, affected, t0,
                 predictors = NULL, pure = NULL, estimator = scm_base,
                 v = NULL, simplified = FALSE, zvar = NULL, ...) {
  stopifnot(is.data.frame(data))
  for (nm in c(id, time, outcome)) if (!nm %in% names(data))
    stop(sprintf("Column '%s' not found in data.", nm))

  pm <- .panel_matrix(data, id, time, outcome)
  Y <- pm$Y; ids <- pm$ids; times <- pm$times
  treated <- as.character(treated)
  affected <- as.character(affected)
  focal <- c(treated, affected)                 # unit 1, then M
  if (!all(focal %in% ids)) stop("Some treated/affected ids are not in the data.")
  if (is.null(pure)) pure <- setdiff(ids, focal) else pure <- as.character(pure)
  m <- length(focal)
  estimator_args <- list(...)

  pre  <- times[times <= t0]
  post <- times[times >  t0]
  if (length(post) == 0) stop("No post-intervention periods (check t0).")

  Xmat <- .predictor_matrix(data, id, time, outcome, predictors, ids, times, pre)
  # Z (outcome path over pre) for the Synth backend's SSR optimisation
  zcol <- if (is.null(zvar)) outcome else zvar
  Zmat <- .panel_matrix(data, id, time, zcol)$Y[which(times %in% pre), , drop = FALSE]

  # ---- Step 1 & 2: fit a synthetic control for each focal unit -----------------
  fits <- vector("list", m); names(fits) <- focal
  cross <- vector("list", m)                      # cross-weights among focal units
  Pmat <- matrix(0, nrow = m, ncol = length(pure),
                 dimnames = list(focal, pure))     # pure-control weights (Lemma 3.1)

  for (i in seq_len(m)) {
    fi <- focal[i]
    if (simplified && fi != treated) {
      donors <- pure                              # affected units on pure controls only
    } else {
      donors <- setdiff(ids, fi)                  # unrestricted: everyone else
    }
    fit <- do.call(.fit_one, c(list(
      focal_id = fi, donor_ids = donors, Xmat = Xmat, Ymat = Y,
      Zmat = Zmat, v = v, estimator = estimator
    ), estimator_args))
    fits[[i]] <- fit
    w <- fit$w
    cross[[i]] <- w[names(w) %in% focal]
    Pmat[i, ] <- w[pure][match(pure, names(w[pure]))]
    Pmat[i, is.na(Pmat[i, ])] <- 0
  }

  # ---- Step 3: assemble Omega and beta ----------------------------------------
  if (simplified) {
    # affected-unit rows collapse to identity: triangular system, no inversion
    Omega <- diag(1, m); dimnames(Omega) <- list(focal, focal)
    w1 <- cross[[1]]
    for (j in names(w1)) if (j %in% focal && j != treated) Omega[1, j] <- -w1[[j]]
  } else {
    Omega <- build_omega(cross, focal)
  }

  post_idx <- which(times %in% post)
  beta_post <- vapply(fits, function(f) f$gap[post_idx], numeric(length(post_idx)))
  beta_post <- t(beta_post)                       # m x P
  rownames(beta_post) <- focal

  # ---- Step 4: correct ---------------------------------------------------------
  theta_post <- solve_iscm(beta_post, Omega)      # m x P
  rownames(theta_post) <- focal

  # inclusive gaps: identical to unrestricted pre-intervention, corrected post
  gaps_iscm <- vapply(seq_len(m), function(i) {
    g <- fits[[i]]$gap
    g[post_idx] <- theta_post[i, ]
    g
  }, numeric(length(times)))
  colnames(gaps_iscm) <- focal
  rownames(gaps_iscm) <- as.character(times)

  gaps_unrestricted <- vapply(fits, function(f) f$gap, numeric(length(times)))
  colnames(gaps_unrestricted) <- focal
  rownames(gaps_unrestricted) <- as.character(times)

  # implied pure-control weights a = (e1' Omega^{-1} P)'  (Lemma 3.1 / Prop 4.1)
  a_impl <- as.numeric(solve(t(Omega), diag(1, m)[, 1])) %*% Pmat
  a_impl <- as.numeric(a_impl); names(a_impl) <- pure

  structure(list(
    call = match.call(),
    focal = focal, treated = treated, affected = affected, pure = pure,
    times = times, t0 = t0, pre = pre, post = post,
    Omega = Omega,
    fits = fits,
    beta_post = beta_post,
    effects = theta_post,                 # m x P : row1 theta, rows>1 gamma
    gaps_iscm = gaps_iscm,
    gaps_unrestricted = gaps_unrestricted,
    P = Pmat,
    implied_weights = a_impl,
    Y = Y,
    simplified = simplified,
    predictors = predictors,
    diagnostics = omega_diagnostics(Omega),
    backend = list(estimator = estimator, v = v, zvar = zvar,
                   args = estimator_args)
  ), class = "iscm")
}
