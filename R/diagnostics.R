#' Diagnostics for the iSCM weight matrix
#'
#' Computes the quantities that govern the reliability of the correction
#' (Assumption 2.4 and Section 5.1 of Di Stefano and Mellace): the
#' determinant, the smallest singular value \eqn{\sigma_{\min}(\hat\Omega)} (whose
#' inverse bounds the amplification of the approximation errors), whether
#' \eqn{\hat\Omega} is strictly diagonally dominant by rows, and, when
#' \eqn{m = 2}, the closed-form amplification factor
#' \eqn{\alpha = 1/(1 - \hat w_2 \hat l_1)}.
#'
#' @param Omega The \code{m x m} weight matrix produced by \code{\link{build_omega}}.
#'
#' @return A list with elements \code{det}, \code{sigma_min},
#'   \code{amplification} (equal to \code{1 / sigma_min}, the general bound on
#'   error propagation), \code{alpha} (the \eqn{m = 2} factor, or \code{NA}),
#'   \code{diag_dominant} (logical) and \code{row_offdiag_sums}.
#'
#' @details How to read the output (Section 5.1, Step 4): \code{amplification}
#'   bounds how much the approximation errors of the underlying fits can be
#'   magnified by the correction (Theorem 2.1(a)); values close to one mean the
#'   correction is stable. \code{diag_dominant = TRUE} guarantees that
#'   \eqn{\hat\Omega} is invertible (Lemma 2.3). Every \code{"iscm"} object
#'   stores this output in \code{$diagnostics}.
#'
#' @seealso \code{\link{iscm_compare}}, which uses these quantities to choose
#'   between the full and the simplified iSCM.
#'
#' @examples
#' # German reunification: alpha = 1.16, sigma_min = 0.62 (Section 6)
#' Omega <- matrix(c(1, -0.333, -0.421, 1), 2, 2)
#' omega_diagnostics(Omega)
#'
#' # Large mutual weights: still invertible, but errors are amplified
#' omega_diagnostics(matrix(c(1, -0.9, -0.9, 1), 2, 2))
#' @export
omega_diagnostics <- function(Omega) {
  Omega <- as.matrix(Omega)
  m <- nrow(Omega)
  d <- svd(Omega, nu = 0, nv = 0)$d
  smin <- min(d)
  off <- rowSums(abs(Omega)) - abs(diag(Omega))
  # Omega[1,2] = -w2, Omega[2,1] = -l1  =>  w2*l1 = Omega[1,2]*Omega[2,1]
  alpha <- if (m == 2L) 1 / (1 - Omega[1, 2] * Omega[2, 1]) else NA_real_
  list(
    det = det(Omega),
    sigma_min = smin,
    amplification = 1 / smin,
    alpha = alpha,
    diag_dominant = all(off < abs(diag(Omega))),
    row_offdiag_sums = off
  )
}

#' Restricted SCM
#'
#' Standard synthetic control for the main treated unit that excludes all
#' potentially affected units from the donor pool (uses the pure controls only).
#' Provided as the natural benchmark for the decision rule of Section 5.1.
#'
#' To get the restricted fit of a potentially affected unit, call the function
#' with that unit as \code{treated} and the main treated unit in
#' \code{affected}. The restricted fit of each affected unit is also what the
#' simplified iSCM uses (\code{iscm(..., simplified = TRUE)}).
#'
#' @inheritParams iscm
#' @return A list with the donor weights \code{w}, the counterfactual \code{yhat},
#'   the \code{gap} over all periods, and the pre-intervention \code{rmspe}.
#'   It can be passed to \code{\link[=plot.iscm]{plot}(fit, restricted = ...)}.
#'
#' @seealso \code{\link{iscm}}, \code{\link{iscm_compare}}.
#'
#' @examples
#' sim <- simulate_iscm(n_affected = 1, n_pure = 6, tau = 2, spill = 3, seed = 7)
#' rest <- restricted_scm(sim$data, "unit", "time", "Y",
#'                        treated = sim$treated, affected = sim$affected,
#'                        t0 = sim$t0)
#' rest$rmspe                               # pre-intervention fit
#' mean(rest$gap[-seq_len(sim$t0)])         # average post-intervention effect
#'
#' fit <- iscm(sim$data, "unit", "time", "Y", sim$treated, sim$affected, sim$t0)
#' plot(fit, restricted = rest)
#' @export
restricted_scm <- function(data, id, time, outcome, treated, affected, t0,
                           predictors = NULL, pure = NULL, estimator = scm_base,
                           v = NULL, zvar = NULL, ...) {
  pm <- .panel_matrix(data, id, time, outcome)
  ids <- pm$ids; times <- pm$times; Y <- pm$Y
  treated <- as.character(treated); affected <- as.character(affected)
  if (is.null(pure)) pure <- setdiff(ids, c(treated, affected)) else pure <- as.character(pure)
  pre <- times[times <= t0]
  Xmat <- .predictor_matrix(data, id, time, outcome, predictors, ids, times, pre)
  zcol <- if (is.null(zvar)) outcome else zvar
  Zmat <- .panel_matrix(data, id, time, zcol)$Y[which(times %in% pre), , drop = FALSE]
  fit <- .fit_one(treated, pure, Xmat, Y, Zmat, v, estimator, ...)
  fit$rmspe <- rmspe(fit$gap[times <= t0])
  fit$w <- fit$w
  fit
}

#' Pre-intervention root mean squared prediction error
#'
#' @param gap Numeric vector of pre-intervention gaps.
#' @return The RMSPE, \eqn{(T_0^{-1}\sum_t \text{gap}_t^2)^{1/2}} (equation 5.1).
#' @examples
#' sim <- simulate_iscm(seed = 7)
#' fit <- iscm(sim$data, "unit", "time", "Y", sim$treated, sim$affected, sim$t0)
#' pre <- fit$times <= fit$t0
#' rmspe(fit$gaps_unrestricted[pre, "TR"])
#' @export
rmspe <- function(gap) sqrt(mean(gap^2))
