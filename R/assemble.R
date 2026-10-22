#' Assemble an iSCM estimate from pre-computed synthetic-control fits
#'
#' Low-level entry point for applying the iSCM correction to synthetic controls
#' that were estimated externally, for example with a fully customised
#' \code{Synth::dataprep}/\code{Synth::synth} specification (custom predictors,
#' special predictors, a cross-validated \code{V} matrix, or a penalized
#' estimator). The user supplies, for the main treated unit and for each
#' potentially affected unit, the full-length gap series and the estimated donor
#' weights; this function builds \eqn{\hat\Omega}, solves the system and returns an
#' object of class \code{"iscm"} with the same methods as \code{\link{iscm}}.
#'
#' This generalises to \eqn{m} units the two-unit, determinant-ratio computation
#' used in the replication scripts of the paper.
#'
#' What you need for each focal unit (main treated first, then each affected
#' unit):
#' \enumerate{
#'   \item its gap \eqn{Y_{it} - \hat Y^N_{it}} over all periods, from a
#'     synthetic control whose donor pool contains the other focal units;
#'   \item the named vector of donor weights of that synthetic control.
#' }
#' The returned object works with \code{print()}, \code{summary()} and
#' \code{plot()}. For the placebo analysis, estimate the synthetic control of
#' each pure control with your own estimator on the data returned by
#' \code{\link{iscm_purge}} and pass the weights to
#' \code{\link{iscm_placebo_fits}}; \code{system.file("examples",
#' "germany_replication.R", package = "iscm")} does this for the paper's
#' \pkg{Synth} specification.
#'
#' @param gaps A named list of numeric vectors, one per focal unit, each the full
#'   pre- and post-intervention gap \eqn{Y_{it} - \hat Y^N_{it}}. Names are the
#'   focal-unit identifiers; the first element must be the main treated unit.
#' @param weights A named list, parallel to \code{gaps}: element \code{i} is the
#'   named donor-weight vector of focal unit \code{i}'s synthetic control (names
#'   are donor identifiers and must include the other focal units they weight).
#' @param t0 Last pre-intervention period (index or value in \code{times}).
#' @param times Optional vector of period labels aligned with the gap vectors.
#'   Defaults to \code{seq_along} the first gap.
#' @param pure Optional vector of pure-control identifiers. If supplied together
#'   with pure-control weights present in \code{weights}, the implied pure-control
#'   weights (Lemma 3.1) are computed.
#' @param simplified Logical; set \code{TRUE} if the affected units were fit on the
#'   pure controls only (their cross-weights on other focal units are then ignored,
#'   giving a triangular system).
#'
#' @return An object of class \code{"iscm"}.
#'
#' @seealso \code{\link{iscm}} for the data-level interface,
#'   \code{\link{iscm_compare_fits}} for the decision rule on external fits,
#'   \code{\link{solve_iscm}} and \code{\link{build_omega}} for the bare linear
#'   algebra.
#'
#' @examples
#' # Estimate the two synthetic controls "outside" iscm(), here with scm_base()
#' # on the pre-intervention outcome path, then apply the correction.
#' sim <- simulate_iscm(n_affected = 1, n_pure = 6, tau = 2, spill = 3, seed = 7)
#' Y <- with(sim$data, tapply(Y, list(time, unit), sum))   # T x J matrix
#' pre <- seq_len(sim$t0)
#' sc <- function(unit) {
#'   donors <- setdiff(colnames(Y), unit)
#'   w <- scm_base(Y[pre, unit], Y[pre, donors])
#'   list(w = w, gap = Y[, unit] - Y[, donors] %*% w)
#' }
#' f_tr <- sc("TR"); f_a1 <- sc("A1")
#'
#' fit <- iscm_assemble(
#'   gaps    = list(TR = as.numeric(f_tr$gap), A1 = as.numeric(f_a1$gap)),
#'   weights = list(TR = f_tr$w, A1 = f_a1$w),
#'   t0 = sim$t0, pure = sim$pure)
#' fit
#'
#' \dontrun{
#' # The paper's application with the exact Synth specification:
#' file.edit(system.file("examples", "germany_replication.R", package = "iscm"))
#' }
#' @export
iscm_assemble <- function(gaps, weights, t0, times = NULL, pure = NULL,
                          simplified = FALSE) {
  if (!is.list(gaps) || is.null(names(gaps))) stop("`gaps` must be a named list.")
  if (!is.list(weights) || is.null(names(weights)))
    stop("`weights` must be a named list.")
  focal <- names(gaps)
  if (!all(focal %in% names(weights)))
    stop("`weights` must have an entry for every focal unit in `gaps`.")
  m <- length(focal)
  len <- length(gaps[[1]])
  if (is.null(times)) times <- seq_len(len)
  if (length(times) != len) stop("`times` must match the length of the gaps.")

  # split pre/post
  if (t0 %in% times) {
    post_idx <- which(times > t0); pre_idx <- which(times <= t0)
  } else {
    post_idx <- which(seq_along(times) > t0); pre_idx <- which(seq_along(times) <= t0)
  }
  post <- times[post_idx]; pre <- times[pre_idx]

  treated <- focal[1]; affected <- focal[-1]

  # cross-weights among focal units -> Omega
  cross <- lapply(focal, function(fi) {
    w <- weights[[fi]]
    w[names(w) %in% focal]
  })
  if (simplified) {
    Omega <- diag(1, m); dimnames(Omega) <- list(focal, focal)
    w1 <- cross[[1]]
    for (j in names(w1)) if (j %in% focal && j != treated) Omega[1, j] <- -w1[[j]]
  } else {
    Omega <- build_omega(cross, focal)
  }

  # contaminated post gaps -> correct
  beta_post <- vapply(focal, function(fi) gaps[[fi]][post_idx],
                      numeric(length(post_idx)))
  beta_post <- t(beta_post); rownames(beta_post) <- focal
  effects <- solve_iscm(beta_post, Omega); rownames(effects) <- focal

  gaps_unrestricted <- vapply(focal, function(fi) gaps[[fi]], numeric(len))
  colnames(gaps_unrestricted) <- focal; rownames(gaps_unrestricted) <- as.character(times)
  gaps_iscm <- gaps_unrestricted
  for (i in seq_len(m)) gaps_iscm[post_idx, i] <- effects[i, ]

  # optional pure-control representation
  implied <- NULL; Pmat <- NULL
  if (!is.null(pure)) {
    pure <- as.character(pure)
    Pmat <- matrix(0, m, length(pure), dimnames = list(focal, pure))
    for (i in seq_len(m)) {
      w <- weights[[focal[i]]]
      wp <- w[pure][match(pure, names(w[pure]))]; wp[is.na(wp)] <- 0
      Pmat[i, ] <- wp
    }
    implied <- as.numeric(solve(t(Omega), diag(1, m)[, 1]) %*% Pmat)
    names(implied) <- pure
  }

  fits <- lapply(focal, function(fi) list(w = weights[[fi]], gap = gaps[[fi]]))
  names(fits) <- focal

  structure(list(
    call = match.call(),
    focal = focal, treated = treated, affected = affected,
    pure = if (is.null(pure)) character(0) else pure,
    times = times, t0 = t0, pre = pre, post = post,
    Omega = Omega, fits = fits, beta_post = beta_post,
    effects = effects, gaps_iscm = gaps_iscm, gaps_unrestricted = gaps_unrestricted,
    P = Pmat, implied_weights = implied, Y = NULL,
    simplified = simplified, predictors = NULL,
    diagnostics = omega_diagnostics(Omega)
  ), class = "iscm")
}
