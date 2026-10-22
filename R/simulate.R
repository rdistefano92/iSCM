#' Simulate a panel with treatment and spillover effects
#'
#' Generates a long-format panel from a linear factor model with two common
#' factors and one "special" factor, in the spirit of the Monte Carlo designs used
#' to study iSCM. By default the main treated unit and the affected units share the
#' special factor, while the affected units differ from the treated unit in the
#' common factors and the pure controls bracket the treated unit's common-factor
#' loadings. This yields a well-conditioned system with moderate cross-weights (as
#' in the German-reunification application) in which keeping the affected units in
#' the donor pool improves the pre-intervention fit but contaminates the naive
#' estimate.
#'
#' @param n_affected Number of potentially affected units.
#' @param n_pure Number of pure-control units (\eqn{\ge 3} recommended).
#' @param t0,t1 Numbers of pre- and post-intervention periods.
#' @param tau Effect on the main treated unit (scalar or length-\code{t1} vector).
#' @param spill Effects on the affected units (scalar, length-\code{t1} vector, or
#'   \code{n_affected x t1} matrix).
#' @param special_loading Loadings on the shared special factor for the treated and
#'   affected units (scalar or length \code{1 + n_affected}). Pure controls load
#'   near zero on it, so they cannot reproduce it alone.
#' @param noise_sd Standard deviation of the idiosyncratic noise.
#' @param seed Optional random seed.
#'
#' @return A list with \code{data} (long \code{data.frame}: \code{unit}, \code{time},
#'   \code{Y}), the identifiers \code{treated}, \code{affected}, \code{pure}, the
#'   value \code{t0}, and the true \code{tau}/\code{spill} for benchmarking.
#'
#' @examples
#' sim <- simulate_iscm(n_affected = 1, n_pure = 6, tau = 2, spill = 3, seed = 7)
#' fit <- iscm(sim$data, "unit", "time", "Y", sim$treated, sim$affected, sim$t0)
#' summary(fit)
#' @export
simulate_iscm <- function(n_affected = 1, n_pure = 6, t0 = 24, t1 = 6,
                          tau = 2, spill = 3, special_loading = NULL,
                          noise_sd = 0.015, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  Tt <- t0 + t1
  tt <- seq_len(Tt)
  m <- 1L + n_affected
  n <- m + n_pure

  if (is.null(special_loading)) {
    special_loading <- rep(0.6, m)
  } else if (length(special_loading) == 1L) {
    special_loading <- rep(special_loading, m)
  }

  # Common-factor loadings (deterministic, for a reproducible well-conditioned
  # design): treated anchored high on factor 1; affected units shifted toward
  # factor 2 (so pure controls beat them on the common factors); pure controls
  # bracket the treated unit's common loadings and load near zero on the special
  # factor (so they cannot reproduce it alone, making the affected units useful).
  aff1 <- if (n_affected > 0) seq(0.66, 0.62, length.out = n_affected) else numeric(0)
  aff2 <- if (n_affected > 0) seq(0.64, 0.58, length.out = n_affected) else numeric(0)
  pur1 <- seq(0.60, 0.95, length.out = n_pure)
  pur2 <- seq(0.66, 0.40, length.out = n_pure)
  lam1 <- c(0.90, aff1, pur1)
  lam2 <- c(0.45, aff2, pur2)
  lam3 <- c(special_loading, rep(0.02, n_pure))

  f1 <- 0.10 * tt + 0.02 * tt^1.2
  f2 <- sin(tt / 2.4)
  g  <- 0.22 * sin(tt / 1.3) + 0.14 * cos(tt / 2.0)

  Y <- matrix(NA_real_, nrow = Tt, ncol = n)
  for (j in seq_len(n)) {
    upre <- c(stats::rnorm(t0, 0, 0.015), rep(0, t1))
    Y[, j] <- 1 + lam1[j] * f1 + lam2[j] * f2 + lam3[j] * g + upre +
      stats::rnorm(Tt, 0, noise_sd)
  }

  tau_v <- if (length(tau) == 1L) rep(tau, t1) else tau
  if (is.matrix(spill)) {
    spill_m <- spill
  } else {
    sv <- if (length(spill) == 1L) rep(spill, t1) else spill
    spill_m <- matrix(sv, nrow = max(n_affected, 0), ncol = t1, byrow = TRUE)
  }
  Y[(t0 + 1):Tt, 1] <- Y[(t0 + 1):Tt, 1] + tau_v
  if (n_affected > 0) {
    for (a in seq_len(n_affected))
      Y[(t0 + 1):Tt, 1 + a] <- Y[(t0 + 1):Tt, 1 + a] + spill_m[a, ]
  }

  unit_names <- c("TR", if (n_affected > 0) paste0("A", seq_len(n_affected)),
                  paste0("C", seq_len(n_pure)))
  data <- data.frame(
    unit = rep(unit_names, each = Tt),
    time = rep(tt, times = n),
    Y = as.vector(Y),
    stringsAsFactors = FALSE
  )
  list(
    data = data,
    treated = "TR",
    affected = if (n_affected > 0) paste0("A", seq_len(n_affected)) else character(0),
    pure = paste0("C", seq_len(n_pure)),
    t0 = t0,
    tau = tau_v,
    spill = if (n_affected > 0) spill_m else NULL
  )
}

#' Simulate a panel in which the affected unit is essential for the fit
#'
#' Generates a panel in which iSCM outperforms the restricted SCM, illustrating
#' the mechanism of Section 3.1 of Di Stefano and Mellace: when the synthetic
#' controls are chosen by matching on observed predictors, excluding the
#' potentially affected unit can lead the restricted optimizer to donors that
#' match the observed predictors but differ in unobserved characteristics.
#'
#' The outcome follows a linear factor model with two observed covariates
#' (\code{x1}, \code{x2}) and one unobserved loading \eqn{u}:
#' \deqn{Y_{jt} = 10 + 3x_{1j}f_{1t} + 3x_{2j}f_{2t} + u_j f_{3t} + \varepsilon_{jt}.}
#' The units are placed in covariate space as follows:
#' \describe{
#'   \item{\code{TR}, \code{A1}}{The main treated unit and the potentially
#'     affected unit are close to each other and share \eqn{u = 1} (think of
#'     West Germany and Austria).}
#'   \item{\code{C1}, \code{C2}}{Pure controls on either side of the pair, also
#'     with \eqn{u = 1}, but farther from \code{TR}.}
#'   \item{\code{C3}, \code{C4}}{Pure controls that are the closest to \code{TR}
#'     in covariates but have \eqn{u = 0}.}
#'   \item{\code{C5}-\code{C8}}{Remote pure controls.}
#' }
#' A synthetic control matched on \code{x1} and \code{x2} reproduces \code{TR}
#' from \code{A1} and \code{C1}, and \code{A1} from \code{TR} and \code{C2}, so the
#' unrestricted fits match \eqn{u} and the iSCM correction is accurate. Without
#' \code{A1}, the closest feasible combination relies on \code{C3}-\code{C4}, whose
#' unobserved loading is wrong, so the restricted SCM is biased; the simplified
#' iSCM inherits the same problem through the restricted fit of \code{A1}.
#'
#' Use it with a predictor-based estimator, typically
#' \code{iscm(..., predictors = c("x1", "x2"), estimator = scm_penalized)}. The
#' penalty makes the solution unique: with exact matching on two predictors,
#' many weight vectors attain a perfect fit.
#'
#' @inheritParams simulate_iscm
#' @param x_sd Standard deviation of the noise added to each unit's covariates.
#' @param u_scale Scale of the unobserved component \eqn{u_j f_{3t}}; larger
#'   values increase the bias of the restricted SCM.
#'
#' @return A list as in \code{\link{simulate_iscm}}; \code{data} additionally
#'   contains the time-invariant covariates \code{x1} and \code{x2}, and the list
#'   contains the unobserved loadings \code{u}.
#'
#' @seealso \code{\link{simulate_iscm}} for a design in which the restricted SCM
#'   performs well; \code{\link{scm_penalized}}.
#'
#' @examples
#' sim <- simulate_iscm_covariates(seed = 1)     # true effects: 2 on TR, 3 on A1
#' args <- list(sim$data, "unit", "time", "Y", "TR", "A1", sim$t0,
#'              predictors = c("x1", "x2"), estimator = scm_penalized)
#'
#' fit  <- do.call(iscm, args)
#' fit_s <- do.call(iscm, c(args, simplified = TRUE))
#' rest <- do.call(restricted_scm, args)
#' post <- fit$times > sim$t0
#' c(unrestricted = mean(fit$gaps_unrestricted[post, "TR"]),
#'   restricted   = mean(rest$gap[post]),
#'   full_iSCM    = mean(fit$effects["TR", ]),
#'   simpl_iSCM   = mean(fit_s$effects["TR", ]))
#' @export
simulate_iscm_covariates <- function(t0 = 20, t1 = 8, tau = 2, spill = 3,
                                     noise_sd = 0.1, x_sd = 0.005, u_scale = 1,
                                     seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  Tt <- t0 + t1
  tt <- seq_len(Tt)
  # columns: x1, x2, u
  geo <- rbind(TR = c( 0.0,  0.0, 1), A1 = c( 0.3,  0.0, 1),
               C1 = c(-0.3,  0.0, 1), C2 = c( 0.6,  0.0, 1),
               C3 = c( 0.2,  0.3, 0), C4 = c( 0.2, -0.3, 0),
               C5 = c(-0.6,  0.5, 1), C6 = c(-0.6, -0.5, 1),
               C7 = c( 0.9,  0.6, 0.5), C8 = c( 0.9, -0.6, 0.5))
  n <- nrow(geo)
  X <- geo[, 1:2] + matrix(stats::rnorm(2 * n, 0, x_sd), n)
  f1 <- 1 + 0.05 * tt
  f2 <- 0.5 * cos(tt / 3)
  f3 <- 1 + 0.1 * tt
  Y <- vapply(seq_len(n), function(j)
    10 + 3 * X[j, 1] * f1 + 3 * X[j, 2] * f2 + u_scale * geo[j, 3] * f3 +
      stats::rnorm(Tt, 0, noise_sd), numeric(Tt))

  tau_v <- if (length(tau) == 1L) rep(tau, t1) else tau
  spill_v <- if (length(spill) == 1L) rep(spill, t1) else spill
  Y[(t0 + 1):Tt, 1] <- Y[(t0 + 1):Tt, 1] + tau_v
  Y[(t0 + 1):Tt, 2] <- Y[(t0 + 1):Tt, 2] + spill_v

  ids <- rownames(geo)
  data <- data.frame(
    unit = rep(ids, each = Tt),
    time = rep(tt, times = n),
    Y = as.vector(Y),
    x1 = rep(X[, 1], each = Tt),
    x2 = rep(X[, 2], each = Tt),
    stringsAsFactors = FALSE
  )
  list(
    data = data,
    treated = "TR",
    affected = "A1",
    pure = ids[-(1:2)],
    t0 = t0,
    tau = tau_v,
    spill = matrix(spill_v, nrow = 1),
    u = stats::setNames(geo[, 3], ids)
  )
}
