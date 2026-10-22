#' iscm: the inclusive synthetic control method
#'
#' Implements the inclusive synthetic control method (iSCM) of Di Stefano and
#' Mellace. iSCM keeps units that are \emph{potentially affected} by an
#' intervention (other treated units, or controls exposed to spillovers) in the
#' donor pool and removes the resulting contamination exactly, on top of any
#' synthetic-control-type estimator.
#'
#' @section Before you start:
#' You need a balanced long-format panel (one row per unit and period) and three
#' choices:
#' \itemize{
#'   \item \code{treated}: the main treated unit;
#'   \item \code{affected}: the potentially affected units (the set \eqn{M});
#'     every other unit is treated as a pure control;
#'   \item \code{t0}: the last pre-intervention period.
#' }
#' See Section 5.2 of the paper for guidance on choosing \code{affected}.
#'
#' @section Two ways to use the package:
#' \describe{
#'   \item{Let the package estimate the synthetic controls}{Call
#'     \code{\link{iscm}} on your data. All other high-level functions
#'     (\code{\link{restricted_scm}}, \code{\link{iscm_compare}},
#'     \code{\link{iscm_placebo}}, \code{\link{iscm_sensitivity}}) take the same
#'     arguments. The estimator is chosen with \code{estimator =}: the built-in
#'     \code{\link{scm_base}} (default, no dependencies),
#'     \code{\link{scm_penalized}} (penalized SCM, no dependencies),
#'     \code{\link{scm_synth}}
#'     (the \pkg{Synth} package), or any function you write (see
#'     \code{\link{backends}}).}
#'   \item{Bring your own synthetic controls}{If you need a specification the
#'     data-level interface cannot express (e.g. \code{Synth::dataprep} with
#'     special predictors and different time windows, as in the paper's
#'     application, or a penalized estimator), estimate the synthetic controls
#'     yourself and pass their gaps and weights to \code{\link{iscm_assemble}}.
#'     Use \code{\link{iscm_compare_fits}} for the decision rule and
#'     \code{\link{iscm_purge}} with \code{\link{iscm_placebo_fits}} for the
#'     placebo analysis on those fits.}
#' }
#'
#' @section What can I do, and with which command?:
#' \tabular{ll}{
#'   \strong{Task} \tab \strong{Command} \cr
#'   Estimate the full iSCM \tab \code{iscm(...)} \cr
#'   Estimate the simplified iSCM \tab \code{iscm(..., simplified = TRUE)} \cr
#'   Estimate the restricted SCM (affected units dropped) \tab \code{restricted_scm(...)} \cr
#'   Get the unrestricted (uncorrected) SCM gaps \tab \code{fit$gaps_unrestricted} \cr
#'   Average and per-period effects \tab \code{print(fit)}, \code{summary(fit)} \cr
#'   Plot the gaps (Figure 1 style) \tab \code{plot(fit, unit = ..., restricted = ...)} \cr
#'   Conditioning of the correction (det, \eqn{\sigma_{\min}}, \eqn{\alpha}) \tab \code{fit$diagnostics}, \code{omega_diagnostics()} \cr
#'   Implied weights on pure controls (Lemma 3.1) \tab \code{fit$implied_weights} \cr
#'   Choose restricted / full / simplified (Section 5.1) \tab \code{iscm_compare(...)} \cr
#'   Same, for fits estimated outside the package \tab \code{iscm_compare_fits()}, \code{iscm_decision()} \cr
#'   Placebo inference (Online Appendix B) \tab \code{iscm_placebo(fit, ...)}, \code{plot()} \cr
#'   Same, for fits estimated outside the package \tab \code{iscm_purge()}, \code{iscm_placebo_fits()} \cr
#'   Sensitivity to the choice of \eqn{M} (Section 5.2) \tab \code{iscm_sensitivity(...)} \cr
#'   Apply the correction to external fits \tab \code{iscm_assemble()} \cr
#'   Only the linear algebra (Theorem 2.1) \tab \code{build_omega()}, \code{solve_iscm()} \cr
#'   Simulate data with known effects \tab \code{simulate_iscm()} \cr
#'   Simulate a case where iSCM beats the restricted SCM \tab \code{simulate_iscm_covariates()} \cr
#'   Use the German-reunification data \tab \code{data(germany)} \cr
#' }
#'
#' @section Replicating the paper:
#' The complete scripts for the empirical application ship with the package:
#' \preformatted{
#' system.file("examples", "germany_replication.R", package = "iscm")
#' system.file("examples", "germany_penalized.R",   package = "iscm")
#' }
#' The first reproduces Tables 1 and C.1-C.2, Figure 1 and the placebo analysis
#' (Figure C.1); the second the inclusive penalized SCM (Figure C.2). Both require
#' \pkg{Synth}; the second also \pkg{pensynth} and \pkg{clarabel}. The vignette
#' (\code{vignette("iscm")}) walks through the whole workflow.
#'
#' @references Di Stefano, R. and Mellace, G. The inclusive synthetic control
#'   method.
#'
#' @examples
#' sim <- simulate_iscm(n_affected = 1, n_pure = 6, tau = 2, spill = 3, seed = 7)
#' fit <- iscm(sim$data, id = "unit", time = "time", outcome = "Y",
#'             treated = sim$treated, affected = sim$affected, t0 = sim$t0)
#' fit           # average effects (true values: 2 on TR, 3 on A1) and diagnostics
#' @importFrom stats optim sd rnorm runif
#' @importFrom utils tail
#' @importFrom graphics abline legend lines dotchart
"_PACKAGE"
