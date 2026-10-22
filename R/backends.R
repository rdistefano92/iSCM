#' SC-type backends for iSCM
#'
#' An iSCM backend is a function that, given the pre-intervention predictors of a
#' treated unit (\code{X1}) and of the donor pool (\code{X0}), returns a vector of
#' synthetic-control weights over the donors. Any synthetic-control-type estimator
#' whose counterfactual is a weighted average of donor outcomes can be used
#' (Di Stefano and Mellace); the iSCM correction is applied on top of it and
#' does not depend on which backend produced the weights.
#'
#' A backend must have the signature \code{function(X1, X0, v = NULL, ...)} and
#' return a numeric vector \code{w} of length \code{ncol(X0)} such that the columns
#' of \code{X0} (the donors) are in the same order as supplied. When the estimator
#' imposes the usual constraints, \code{w >= 0} and \code{sum(w) == 1}.
#'
#' The package calls the backend with \code{X1}, \code{X0}, \code{v} and, in
#' addition, \code{Z1} and \code{Z0} (the pre-intervention outcome paths of the
#' focal unit and of the donors); backends that do not need them can absorb them
#' in \code{...}. Extra arguments given to \code{\link{iscm}} through \code{...}
#' are passed on to the backend.
#'
#' Available backends:
#' \describe{
#'   \item{\code{\link{scm_base}}}{Classic simplex-constrained SCM with fixed
#'     predictor weights. No dependencies; the default.}
#'   \item{\code{\link{scm_penalized}}}{Penalized SCM of Abadie and L'Hour
#'     (2021), which favours donors close to the focal unit. No dependencies.}
#'   \item{\code{\link{scm_synth}}}{Wrapper around \code{Synth::synth}, which also
#'     optimises the predictor weights \code{V} unless \code{v} is supplied.}
#'   \item{Your own function}{Any SC-type estimator, see the example.}
#' }
#' If an estimator cannot be written in this form (e.g. it needs the whole panel),
#' estimate the synthetic controls yourself and use \code{\link{iscm_assemble}}.
#'
#' @examples
#' sim <- simulate_iscm(n_affected = 1, n_pure = 6, seed = 7)
#'
#' # A custom backend: equal weights on the 3 donors closest in predictors
#' nearest3 <- function(X1, X0, v = NULL, ...) {
#'   d <- colSums((X0 - X1)^2)
#'   w <- as.numeric(rank(d, ties.method = "first") <= 3)
#'   w / sum(w)
#' }
#' fit <- iscm(sim$data, "unit", "time", "Y", sim$treated, sim$affected,
#'             sim$t0, estimator = nearest3)
#' fit$Omega
#'
#' # Passing options to the built-in backend through iscm()
#' iscm(sim$data, "unit", "time", "Y", sim$treated, sim$affected, sim$t0,
#'      estimator = scm_base, standardize = TRUE)   # behaviour of version 0.1.1
#' @name backends
NULL

#' Dependency-free simplex SCM backend
#'
#' Solves the classic synthetic-control problem
#' \eqn{\min_w \|X_1 - X_0 w\|_V^2} subject to \eqn{w \ge 0} and
#' \eqn{\sum_j w_j = 1}, using a softmax reparametrisation so that no external
#' quadratic-programming dependency is required. This makes the package usable and
#' testable out of the box. It is a genuine SC-type estimator, but unlike
#' \code{\link{scm_synth}} it does not re-optimise the predictor-weighting matrix
#' \code{V}; supply \code{v} to weight predictors, otherwise predictors are
#' standardised and weighted equally.
#'
#' @param X1 Numeric vector of length \code{k}: predictors of the treated unit.
#' @param X0 Numeric \code{k x J0} matrix: predictors of the \code{J0} donors.
#' @param v Optional numeric vector of length \code{k} of predictor weights
#'   (the diagonal of \code{V}). Defaults to equal weights on standardised
#'   predictors.
#' @param standardize Logical; if \code{TRUE} (default) predictors are scaled to
#'   unit variance before matching. When \code{\link{iscm}} and the other
#'   high-level functions match on the pre-intervention outcome path
#'   (\code{predictors = NULL}), they call the built-in backends with
#'   \code{standardize = FALSE} unless \code{standardize} is given explicitly:
#'   standardising each period would give the most weight to the periods in which
#'   units are most alike.
#' @param control List of optimiser controls passed to \code{stats::optim}.
#' @param ... Ignored.
#'
#' @return Numeric vector of donor weights, named with the column names of
#'   \code{X0} when available.
#'
#' @seealso \code{\link{backends}}, \code{\link{scm_synth}}.
#'
#' @examples
#' set.seed(1)
#' X0 <- matrix(runif(12), nrow = 3, dimnames = list(NULL, paste0("d", 1:4)))
#' X1 <- X0 %*% c(0.5, 0.5, 0, 0)          # inside the convex hull
#' round(scm_base(X1, X0), 3)
#' @export
scm_base <- function(X1, X0, v = NULL, standardize = TRUE,
                     control = list(maxit = 2000, reltol = 1e-12), ...) {
  p <- .scale_predictors(X1, X0, v, standardize)
  w <- .simplex_fit(p$X1, p$X0, lambda = 0, control = control)
  names(w) <- colnames(p$X0)
  w
}

# Internal: standardise and weight predictors (shared by the built-in backends).
.scale_predictors <- function(X1, X0, v, standardize) {
  X0 <- as.matrix(X0)
  cn <- colnames(X0)
  X1 <- as.numeric(X1)
  k <- length(X1)
  if (nrow(X0) != k) stop("X0 must have as many rows as length(X1).")
  if (standardize) {
    s <- apply(cbind(X1, X0), 1L, stats::sd)
    s[s == 0 | !is.finite(s)] <- 1
    X1 <- X1 / s
    X0 <- X0 / s
  }
  if (is.null(v)) v <- rep(1, k)
  vs <- sqrt(as.numeric(v) / sum(v))
  X0 <- X0 * vs
  colnames(X0) <- cn
  list(X1 = X1 * vs, X0 = X0)
}

# Internal: minimise ||X1 - X0 w||^2 + lambda * sum_j w_j ||X1 - X0_j||^2 over the
# simplex, via a softmax reparametrisation (no QP dependency).
.simplex_fit <- function(X1, X0, lambda, control) {
  J0 <- ncol(X0)
  if (J0 == 1L) return(1)
  d <- colSums((X0 - X1)^2)
  # softmax(z) -> simplex; optimise over z (J0 - 1 free, first fixed at 0)
  obj <- function(z) {
    zz <- c(0, z)
    w <- exp(zz - max(zz))
    w <- w / sum(w)
    r <- X1 - X0 %*% w
    sum(r * r) + lambda * sum(w * d)
  }
  opt <- stats::optim(rep(0, J0 - 1L), obj, method = "BFGS", control = control)
  zz <- c(0, opt$par)
  w <- exp(zz - max(zz))
  w / sum(w)
}

#' Dependency-free penalized SCM backend (Abadie and L'Hour, 2021)
#'
#' Solves the penalized synthetic-control problem of Abadie and L'Hour (2021),
#' \deqn{\min_w \|X_1 - X_0 w\|_V^2 + \lambda \sum_j w_j \|X_1 - X_j\|_V^2}
#' subject to \eqn{w \ge 0} and \eqn{\sum_j w_j = 1}. The penalty favours donors
#' that are individually close to the focal unit, which reduces interpolation
#' bias and makes the solution unique when the focal unit lies inside the convex
#' hull of the donors. With \code{lambda = 0} it coincides with
#' \code{\link{scm_base}}.
#'
#' Pass \code{lambda} through \code{\link{iscm}} and the other high-level
#' functions, e.g. \code{iscm(..., estimator = scm_penalized, lambda = 0.1)}. The
#' penalty is not cross-validated here; for the cross-validated estimator used in
#' Online Appendix C of the paper, see the \pkg{pensynth} example in
#' \code{system.file("examples", "germany_penalized.R", package = "iscm")}.
#'
#' @inheritParams scm_base
#' @param lambda Non-negative penalty parameter. Default \code{0.1}.
#'
#' @return Numeric vector of donor weights, named with the column names of
#'   \code{X0} when available.
#'
#' @references Abadie, A. and L'Hour, J. (2021). A penalized synthetic control
#'   estimator for disaggregated data. \emph{Journal of the American Statistical
#'   Association} 116(536), 1817-1834.
#'
#' @seealso \code{\link{backends}}, \code{\link{scm_base}},
#'   \code{\link{simulate_iscm_covariates}}.
#'
#' @examples
#' # Two exact solutions: d1 + d2 (close to X1) or d3 + d4 (far from X1).
#' X0 <- cbind(d1 = c(-1, 0), d2 = c(1, 0), d3 = c(-3, 0), d4 = c(3, 0))
#' X1 <- c(0, 0)
#' round(scm_penalized(X1, X0, lambda = 0.1, standardize = FALSE), 3)
#' @export
scm_penalized <- function(X1, X0, v = NULL, lambda = 0.1, standardize = TRUE,
                          control = list(maxit = 5000, reltol = 1e-14), ...) {
  if (!is.numeric(lambda) || length(lambda) != 1L || lambda < 0)
    stop("`lambda` must be a single non-negative number.")
  p <- .scale_predictors(X1, X0, v, standardize)
  w <- .simplex_fit(p$X1, p$X0, lambda = lambda, control = control)
  names(w) <- colnames(p$X0)
  w
}

#' Synth backend (classic SCM of Abadie et al.)
#'
#' Wraps \code{Synth::synth} to compute synthetic-control weights, using the
#' predictor-weighting matrix \code{V} supplied in \code{v} (e.g. one obtained by
#' cross-validation as in Abadie et al., 2015). Requires the \pkg{Synth} package.
#' This backend reproduces the estimator used in the empirical application of
#' Di Stefano and Mellace.
#'
#' @inheritParams scm_base
#' @param Z1 Optional numeric vector: treated unit's outcomes over the SSR-fit
#'   window (passed as \code{Z1} to \code{Synth::synth}). Defaults to \code{X1}.
#' @param Z0 Optional matrix: donors' outcomes over the SSR-fit window (\code{Z0}).
#'   Defaults to \code{X0}.
#' @param ... Additional arguments passed to \code{Synth::synth} (e.g.
#'   \code{Margin.ipop}, \code{Sigf.ipop}, \code{Bound.ipop}).
#'
#' @return Numeric vector of donor weights.
#'
#' @seealso \code{\link{backends}}; \code{\link{iscm_assemble}} when the paper's
#'   full \code{Synth::dataprep} specification is needed.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("Synth", quietly = TRUE)) {
#'   sim <- simulate_iscm(n_affected = 1, n_pure = 6, seed = 7)
#'   fit <- iscm(sim$data, "unit", "time", "Y", sim$treated, sim$affected,
#'               sim$t0, estimator = scm_synth)
#'   fit
#' }
#' }
#' @export
scm_synth <- function(X1, X0, v = NULL, Z1 = NULL, Z0 = NULL, ...) {
  if (!requireNamespace("Synth", quietly = TRUE)) {
    stop("scm_synth() requires the 'Synth' package. Install it, or use ",
         "estimator = scm_base.", call. = FALSE)
  }
  X0 <- as.matrix(X0)
  X1 <- matrix(as.numeric(X1), ncol = 1L)
  if (is.null(Z1)) Z1 <- X1
  if (is.null(Z0)) Z0 <- X0
  Z1 <- matrix(as.numeric(Z1), ncol = 1L)
  Z0 <- as.matrix(Z0)
  args <- list(X1 = X1, X0 = X0, Z1 = Z1, Z0 = Z0, ...)
  if (!is.null(v)) args$custom.v <- as.numeric(v)
  out <- do.call(Synth::synth, args)
  w <- as.numeric(out$solution.w)
  names(w) <- colnames(X0)
  w
}
