#' @export
print.iscm <- function(x, ...) {
  cat("Inclusive synthetic control method",
      if (x$simplified) "(simplified)" else "", "\n", sep = " ")
  cat("  Main treated unit :", x$treated, "\n")
  cat("  Affected units    :", if (length(x$affected)) paste(x$affected, collapse = ", ") else "(none)", "\n")
  cat("  Pure controls     :", length(x$pure), "\n")
  cat("  Pre / post periods :", length(x$pre), "/", length(x$post), "\n")
  eff <- rowMeans(x$effects)
  cat("\nAverage post-intervention effects:\n")
  lab <- paste0(ifelse(x$focal == x$treated, "theta[", "gamma["), x$focal, "]")
  for (i in seq_along(eff)) cat(sprintf("  %-14s % .4f\n", lab[i], eff[i]))
  d <- x$diagnostics
  cat(sprintf("\nOmega: det = %.3f, sigma_min = %.3f, amplification = %.3f%s\n",
              d$det, d$sigma_min, d$amplification,
              if (!is.na(d$alpha)) sprintf(", alpha = %.3f", d$alpha) else ""))
  invisible(x)
}

#' Summarise an iSCM fit
#'
#' @param object An object of class \code{"iscm"}.
#' @param ... Ignored.
#' @return An object of class \code{"summary.iscm"} with a per-period effects table
#'   (\code{$effects}), the average effects (\code{$average}), \code{Omega}, the
#'   diagnostics and the implied pure-control weights.
#' @examples
#' sim <- simulate_iscm(seed = 7)
#' fit <- iscm(sim$data, "unit", "time", "Y", sim$treated, sim$affected, sim$t0)
#' s <- summary(fit)
#' s
#' s$effects            # per-period effects as a matrix, ready to export
#' @export
summary.iscm <- function(object, ...) {
  eff <- object$effects
  tab <- t(eff)
  colnames(tab) <- paste0(ifelse(object$focal == object$treated, "theta_", "gamma_"),
                          object$focal)
  rownames(tab) <- as.character(object$post)
  out <- list(
    treated = object$treated, affected = object$affected,
    simplified = object$simplified,
    effects = tab,
    average = colMeans(tab),
    Omega = object$Omega,
    diagnostics = object$diagnostics,
    implied_weights = object$implied_weights
  )
  class(out) <- "summary.iscm"
  out
}

#' @export
print.summary.iscm <- function(x, digits = 4, ...) {
  cat("Inclusive synthetic control method",
      if (x$simplified) "(simplified)\n" else "\n")
  cat("\nPer-period effects:\n")
  print(round(x$effects, digits))
  cat("\nAverage effect(s):\n")
  print(round(x$average, digits))
  cat("\nWeight matrix Omega:\n")
  print(round(x$Omega, digits))
  d <- x$diagnostics
  cat(sprintf("\ndet = %.4f | sigma_min = %.4f | amplification (1/sigma_min) = %.4f%s\n",
              d$det, d$sigma_min, d$amplification,
              if (!is.na(d$alpha)) sprintf(" | alpha = %.4f", d$alpha) else ""))
  cat("Strictly diagonally dominant:", d$diag_dominant, "\n")
  invisible(x)
}

#' Plot iSCM gaps
#'
#' Plots the estimated gap for a focal unit under the inclusive, unrestricted and
#' (optionally) restricted synthetic control, in the style of Figure 1 of the
#' paper.
#'
#' @param x An object of class \code{"iscm"}.
#' @param unit Which focal unit to plot (identifier). Defaults to the main treated
#'   unit.
#' @param restricted Optional object returned by \code{\link{restricted_scm}} to
#'   overlay the restricted gap. Ignored for affected units unless supplied.
#' @param main,xlab,ylab Standard plotting labels.
#' @param ... Passed to \code{plot}.
#' @return Invisibly, a data frame of the plotted series.
#' @examples
#' sim <- simulate_iscm(seed = 7)
#' fit  <- iscm(sim$data, "unit", "time", "Y", sim$treated, sim$affected, sim$t0)
#' rest <- restricted_scm(sim$data, "unit", "time", "Y", sim$treated,
#'                        sim$affected, sim$t0)
#' op <- par(mfrow = c(1, 2))
#' plot(fit, restricted = rest)      # main treated unit
#' plot(fit, unit = "A1")            # affected unit
#' par(op)
#' @export
plot.iscm <- function(x, unit = NULL, restricted = NULL,
                      main = NULL, xlab = "Time",
                      ylab = "Gap", ...) {
  if (is.null(unit)) unit <- x$treated
  unit <- as.character(unit)
  if (!unit %in% x$focal) stop("`unit` must be one of the focal units.")
  ti <- x$times
  g_iscm <- x$gaps_iscm[, unit]
  g_unr  <- x$gaps_unrestricted[, unit]
  ylim <- range(c(g_iscm, g_unr,
                  if (!is.null(restricted)) restricted$gap), na.rm = TRUE)
  if (is.null(main)) main <- paste0(unit, if (x$simplified) " - simplified iSCM" else " - iSCM")
  plot(ti, g_iscm, type = "l", lwd = 2, col = "black",
       ylim = ylim, xlab = xlab, ylab = ylab, main = main, ...)
  graphics::lines(ti, g_unr, col = "red", lwd = 2)
  leg <- c("inclusive", "unrestricted"); col <- c("black", "red")
  if (!is.null(restricted)) {
    graphics::lines(ti, restricted$gap, col = "purple", lwd = 2)
    leg <- c(leg, "restricted"); col <- c(col, "purple")
  }
  graphics::abline(v = x$t0 + 0.5, lty = "dotted")
  graphics::abline(h = 0, lty = "dotted")
  graphics::legend("topleft", legend = leg, col = col, lwd = 2, bty = "n", cex = 0.8)
  invisible(data.frame(time = ti, iscm = g_iscm, unrestricted = g_unr))
}
