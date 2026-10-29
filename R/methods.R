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
