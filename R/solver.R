#' Solve the iSCM linear system
#'
#' Core linear-algebra step of the inclusive synthetic control method. Given the
#' matrix of contaminated estimates (gaps) for the \code{m} focal units over the
#' post-intervention periods and the weight matrix \eqn{\hat\Omega}, it returns the
#' disentangled effects \eqn{\hat\vartheta_t = \hat\Omega^{-1}\hat\beta_t}
#' (Theorem 2.1 of Di Stefano and Mellace).
#'
#' This function contains no estimation; it only inverts the contamination. It is
#' exposed so that the correction can be applied to gaps produced by any SC-type
#' estimator.
#'
#' @param beta Numeric \code{m x P} matrix of contaminated estimates: row \code{i}
#'   is the gap of focal unit \code{i}, columns are post-intervention periods.
#'   May also be a length-\code{m} vector (single period).
#' @param Omega Numeric \code{m x m} weight matrix with ones on the diagonal and
#'   \eqn{-\hat w} cross-weights off the diagonal (see \code{\link{iscm}}).
#' @param tol Numeric; if \code{sigma_min(Omega)} is below \code{tol} the system is
#'   treated as numerically singular and an error is thrown.
#'
#' @return Numeric matrix (or vector) of the same shape as \code{beta} containing
#'   the corrected effects; row 1 is \eqn{\hat\theta_{1t}}, rows 2..m are
#'   \eqn{\hat\gamma_{it}}.
#'
#' @seealso \code{\link{build_omega}} to construct \code{Omega},
#'   \code{\link{iscm_assemble}} for the same correction with diagnostics and
#'   methods, \code{\link{iscm}} for the full estimator.
#'
#' @examples
#' # German reunification (Section 6): West Germany gives Austria weight 0.42,
#' # Austria gives West Germany weight 0.33.
#' Omega <- build_omega(list(c(Austria = 0.42), c("West Germany" = 0.33)),
#'                      focal = c("West Germany", "Austria"))
#' # Contaminated gaps in two post-intervention years (unrestricted SCM)
#' beta <- rbind("West Germany" = c(-1217, -2645), "Austria" = c(249, 894))
#' solve_iscm(beta, Omega)
#' # Same as the closed form (4.1) for m = 2
#' (beta[1, ] + 0.42 * beta[2, ]) / (1 - 0.42 * 0.33)
#' @export
solve_iscm <- function(beta, Omega, tol = 1e-10) {
  Omega <- as.matrix(Omega)
  m <- nrow(Omega)
  if (ncol(Omega) != m) stop("Omega must be square.")
  vec_in <- is.null(dim(beta))
  B <- if (vec_in) matrix(beta, ncol = 1L) else as.matrix(beta)
  if (nrow(B) != m) stop("beta must have as many rows as Omega has dimensions.")

  smin <- min(svd(Omega, nu = 0, nv = 0)$d)
  if (!is.finite(smin) || smin < tol) {
    stop(sprintf(paste0("Omega is numerically singular (sigma_min = %.3g). ",
                        "The focal units assign near-degenerate weights to one ",
                        "another; consider the simplified iSCM (simplified = TRUE)."),
                 smin), call. = FALSE)
  }
  Theta <- solve(Omega, B)
  if (vec_in) as.numeric(Theta) else Theta
}

#' Build the iSCM weight matrix Omega
#'
#' Assembles \eqn{\hat\Omega} from the per-unit synthetic-control weights. Entry
#' \eqn{(i,i)} equals one and entry \eqn{(i,j)}, \eqn{i \ne j}, equals minus the
#' weight that focal unit \code{i}'s synthetic control assigns to focal unit
#' \code{j}.
#'
#' @param cross A named list of length \code{m}; element \code{i} is a named
#'   numeric vector giving the weights that focal unit \code{i} assigns to the
#'   other focal units (names are focal-unit identifiers). Missing focal units are
#'   treated as zero weight.
#' @param focal Character (or coercible) vector of the \code{m} focal-unit
#'   identifiers, in the order they should appear in \code{Omega} (main treated
#'   first).
#'
#' @return An \code{m x m} matrix with dimnames given by \code{focal}.
#'
#' @seealso \code{\link{solve_iscm}}, \code{\link{omega_diagnostics}}.
#'
#' @examples
#' # Weights that each focal unit's synthetic control assigns to the others
#' cross <- list(TR = c(A1 = 0.40, A2 = 0.10),
#'               A1 = c(TR = 0.30),
#'               A2 = c(TR = 0.05, A1 = 0.20))
#' Omega <- build_omega(cross, focal = c("TR", "A1", "A2"))
#' Omega
#' omega_diagnostics(Omega)
#' @export
build_omega <- function(cross, focal) {
  focal <- as.character(focal)
  m <- length(focal)
  Omega <- diag(1, m)
  dimnames(Omega) <- list(focal, focal)
  for (i in seq_len(m)) {
    wi <- cross[[i]]
    if (is.null(wi)) next
    nm <- intersect(names(wi), focal)
    for (j in nm) {
      if (j == focal[i]) next
      Omega[i, j] <- -wi[[j]]
    }
  }
  Omega
}
