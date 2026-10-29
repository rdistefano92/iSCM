# Colours: categorical slots 1-3 of a palette validated for colour-vision
# deficiencies; each series also has its own line type, so the plots remain
# readable in black and white.
.iscm_cols <- c("iSCM" = "#2a78d6", "Unrestricted SCM" = "#eb6834",
                "Restricted SCM" = "#1baf7a", "Observed" = "#0b0b0b")
.iscm_ltys <- c("iSCM" = "solid", "Unrestricted SCM" = "dashed",
                "Restricted SCM" = "dotted", "Observed" = "solid")
.ink_muted <- "#52514e"
.placebo_grey <- "#b8b7b2"

#' ggplot2 theme for iscm plots
#'
#' A light theme used by all plots of the package: recessive grid, legend at the
#' bottom, no minor grid lines. Add it to your own plots, or override it on the
#' package's plots with the usual \code{+ theme(...)}.
#'
#' @param base_size Base font size, in points.
#' @return A \code{ggplot2} theme.
#' @examples
#' sim <- simulate_iscm(seed = 7)
#' fit <- iscm(sim$data, "unit", "time", "Y", sim$treated, sim$affected, sim$t0)
#' plot(fit) + theme_iscm(base_size = 14)
#' @export
theme_iscm <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "#e8e7e3", linewidth = 0.3),
      axis.text = ggplot2::element_text(colour = .ink_muted),
      axis.title = ggplot2::element_text(colour = .ink_muted),
      plot.title = ggplot2::element_text(face = "bold", colour = "#0b0b0b"),
      plot.title.position = "plot",
      plot.subtitle = ggplot2::element_text(colour = .ink_muted),
      plot.caption = ggplot2::element_text(colour = .ink_muted, hjust = 0),
      strip.text = ggplot2::element_text(face = "bold", hjust = 0),
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      legend.key.width = ggplot2::unit(2, "lines")
    )
}

# Internal: position of the intervention line, halfway between t0 and the first
# post-intervention period.
.intervention_x <- function(times, t0) {
  times <- as.numeric(times)
  nxt <- times[times > t0]
  if (length(nxt)) (as.numeric(t0) + min(nxt)) / 2 else as.numeric(t0)
}

# Internal: restricted gaps for the requested units. `restricted` may be one
# object (a list with $gap, e.g. from restricted_scm()) for a single unit, or a
# named list of such objects, one per unit.
.restricted_gaps <- function(restricted, units, n_times) {
  if (is.null(restricted)) return(list())
  one <- function(r, u) {
    g <- if (is.list(r)) r$gap else r
    if (is.null(g) || length(g) != n_times)
      stop(sprintf("The restricted fit for '%s' must contain a gap over all %d periods.",
                   u, n_times), call. = FALSE)
    as.numeric(g)
  }
  if (is.list(restricted) && !is.null(restricted$gap)) {
    if (length(units) > 1L)
      stop("With several units, pass `restricted` as a named list, one restricted ",
           "fit per unit.", call. = FALSE)
    return(stats::setNames(list(one(restricted, units)), units))
  }
  out <- lapply(intersect(units, names(restricted)), function(u) one(restricted[[u]], u))
  stats::setNames(out, intersect(units, names(restricted)))
}

#' Plot an iSCM fit
#'
#' Draws the results of \code{\link{iscm}} or \code{\link{iscm_assemble}} with
#' \pkg{ggplot2}. Three types of plot are available:
#' \describe{
#'   \item{\code{"gaps"}}{(default) the estimated effects over time: the
#'     inclusive gap, the unrestricted SCM gap and, if supplied, the restricted
#'     SCM gap, as in Figure 1 of the paper. Before the intervention the
#'     inclusive and unrestricted gaps coincide.}
#'   \item{\code{"trends"}}{the observed outcome against the synthetic
#'     counterfactuals. Needs the observed outcomes, which are stored by
#'     \code{iscm()} but not by \code{iscm_assemble()}.}
#'   \item{\code{"weights"}}{the donor weights of the unit's synthetic control,
#'     with the potentially affected units highlighted; for the main treated unit
#'     also the implied weights on the pure controls (Lemma 3.1).}
#' }
#'
#' The result is a \code{ggplot} object, so it can be modified with the usual
#' \pkg{ggplot2} functions, e.g. \code{plot(fit) + labs(y = "GDP per capita")}.
#'
#' @param x An object of class \code{"iscm"}.
#' @param unit Focal unit(s) to plot. Defaults to the main treated unit. Several
#'   units are drawn in separate panels.
#' @param type One of \code{"gaps"}, \code{"trends"}, \code{"weights"}.
#' @param restricted Optional restricted fit to overlay, as returned by
#'   \code{\link{restricted_scm}} (or any list with a \code{gap} element over all
#'   periods). With several units, a named list with one fit per unit.
#' @param title,xlab,ylab Optional labels; sensible defaults are used.
#' @param min_weight For \code{type = "weights"}: donors with smaller weights are
#'   not shown. Default \code{0.001}.
#' @param ... Ignored.
#'
#' @return A \code{ggplot} object.
#'
#' @seealso \code{\link{plot.iscm_placebo}}, \code{\link{theme_iscm}}.
#'
#' @examples
#' sim  <- simulate_iscm(seed = 7)
#' fit  <- iscm(sim$data, "unit", "time", "Y", sim$treated, sim$affected, sim$t0)
#' rest <- restricted_scm(sim$data, "unit", "time", "Y", sim$treated,
#'                        sim$affected, sim$t0)
#'
#' plot(fit, restricted = rest)                   # gaps of the main treated unit
#' plot(fit, unit = c("TR", "A1"))                # treated and affected unit
#' plot(fit, type = "trends", restricted = rest)  # observed vs synthetic
#' plot(fit, type = "weights")                    # donor and implied weights
#'
#' # ggplot2 customisation
#' plot(fit) + ggplot2::labs(y = "Effect on Y", title = NULL)
#' @export
plot.iscm <- function(x, unit = NULL, type = c("gaps", "trends", "weights"),
                      restricted = NULL, title = NULL, xlab = "Time", ylab = NULL,
                      min_weight = 0.001, ...) {
  type <- match.arg(type)
  if (is.null(unit)) unit <- x$treated
  unit <- as.character(unit)
  bad <- setdiff(unit, x$focal)
  if (length(bad))
    stop("`unit` must be among the focal units: ", paste(x$focal, collapse = ", "),
         call. = FALSE)
  switch(type,
         gaps = .plot_gaps(x, unit, restricted, title, xlab, ylab),
         trends = .plot_trends(x, unit, restricted, title, xlab, ylab),
         weights = .plot_weights(x, unit, title, min_weight))
}

.plot_series <- function(df, x, unit, title, subtitle, xlab, ylab, zero_line) {
  df$series <- factor(df$series, levels = intersect(names(.iscm_cols), unique(df$series)))
  df$unit <- factor(df$unit, levels = unit)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$time, y = .data$value,
                                        colour = .data$series,
                                        linetype = .data$series))
  if (zero_line)
    p <- p + ggplot2::geom_hline(yintercept = 0, colour = .ink_muted, linewidth = 0.3)
  p <- p +
    ggplot2::geom_vline(xintercept = .intervention_x(x$times, x$t0),
                        colour = .ink_muted, linetype = "dashed", linewidth = 0.4) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::scale_colour_manual(values = .iscm_cols) +
    ggplot2::scale_linetype_manual(values = .iscm_ltys) +
    ggplot2::labs(title = title, subtitle = subtitle, x = xlab, y = ylab) +
    theme_iscm()
  if (length(unit) > 1L)
    p <- p + ggplot2::facet_wrap(~ unit, scales = "free_y")
  p
}

.plot_gaps <- function(x, unit, restricted, title, xlab, ylab) {
  times <- as.numeric(x$times)
  rg <- .restricted_gaps(restricted, unit, length(times))
  df <- do.call(rbind, lapply(unit, function(u) {
    out <- rbind(
      data.frame(time = times, value = as.numeric(x$gaps_iscm[, u]),
                 series = "iSCM", unit = u),
      data.frame(time = times, value = as.numeric(x$gaps_unrestricted[, u]),
                 series = "Unrestricted SCM", unit = u))
    if (!is.null(rg[[u]]))
      out <- rbind(out, data.frame(time = times, value = rg[[u]],
                                   series = "Restricted SCM", unit = u))
    out
  }))
  if (is.null(title))
    title <- if (length(unit) == 1L) paste("Estimated effect on", unit) else "Estimated effects"
  sub <- paste0(if (x$simplified) "Simplified iSCM. " else "",
                "Gap between observed and synthetic outcome; dashed line: intervention.")
  .plot_series(df, x, unit, title, sub, xlab, if (is.null(ylab)) "Gap" else ylab,
               zero_line = TRUE)
}

.plot_trends <- function(x, unit, restricted, title, xlab, ylab) {
  if (is.null(x$Y))
    stop("type = \"trends\" needs the observed outcomes, which are stored by ",
         "iscm() but not by iscm_assemble(). Use type = \"gaps\".", call. = FALSE)
  times <- as.numeric(x$times)
  rg <- .restricted_gaps(restricted, unit, length(times))
  df <- do.call(rbind, lapply(unit, function(u) {
    y <- as.numeric(x$Y[, u])
    out <- rbind(
      data.frame(time = times, value = y, series = "Observed", unit = u),
      data.frame(time = times, value = y - as.numeric(x$gaps_iscm[, u]),
                 series = "iSCM", unit = u),
      data.frame(time = times, value = y - as.numeric(x$gaps_unrestricted[, u]),
                 series = "Unrestricted SCM", unit = u))
    if (!is.null(rg[[u]]))
      out <- rbind(out, data.frame(time = times, value = y - rg[[u]],
                                   series = "Restricted SCM", unit = u))
    out
  }))
  if (is.null(title))
    title <- if (length(unit) == 1L) paste(unit, "and its synthetic controls") else
      "Observed and synthetic outcomes"
  .plot_series(df, x, unit, title,
               "Synthetic lines are the estimated outcomes without the intervention.",
               xlab, if (is.null(ylab)) "Outcome" else ylab, zero_line = FALSE)
}

.plot_weights <- function(x, unit, title, min_weight) {
  df <- do.call(rbind, lapply(unit, function(u) {
    w <- x$fits[[u]]$w
    out <- data.frame(donor = names(w), weight = as.numeric(w),
                      panel = paste("Synthetic", u),
                      role = ifelse(names(w) %in% x$focal, "Potentially affected or treated",
                                    "Pure control"))
    if (u == x$treated && !is.null(x$implied_weights)) {
      a <- x$implied_weights
      out <- rbind(out, data.frame(donor = names(a), weight = as.numeric(a),
                                   panel = paste("Implied weights of iSCM for", u),
                                   role = "Pure control"))
    }
    out
  }))
  df <- df[abs(df$weight) >= min_weight, , drop = FALSE]
  df$panel <- factor(df$panel, levels = unique(df$panel))
  df$role <- factor(df$role, levels = c("Pure control", "Potentially affected or treated"))
  df$donor <- stats::reorder(df$donor, df$weight, FUN = max)
  if (is.null(title)) title <- "Synthetic control weights"
  ggplot2::ggplot(df, ggplot2::aes(x = .data$weight, y = .data$donor, fill = .data$role)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", .data$weight)),
                       hjust = -0.15, size = 3, colour = .ink_muted) +
    ggplot2::scale_fill_manual(values = c("Pure control" = "#2a78d6",
                                          "Potentially affected or treated" = "#eb6834")) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) +
    ggplot2::facet_wrap(~ panel, scales = "free_y") +
    ggplot2::labs(title = title, x = "Weight", y = NULL,
                  subtitle = paste("Donors with weight below", min_weight, "are not shown.")) +
    theme_iscm() +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
}

#' Plot an iSCM placebo test
#'
#' Draws the result of \code{\link{iscm_placebo}} or
#' \code{\link{iscm_placebo_fits}} with \pkg{ggplot2}.
#' \describe{
#'   \item{\code{"ratios"}}{(default) the post- to pre-intervention RMSPE ratios
#'     of the target and of every placebo, as in Figure C.1 of the paper.}
#'   \item{\code{"gaps"}}{the gap paths of the target and of the placebos, on the
#'     purged outcomes.}
#' }
#'
#' @param x An object of class \code{"iscm_placebo"}.
#' @param type One of \code{"ratios"}, \code{"gaps"}.
#' @param prune For \code{type = "gaps"}: if a number \eqn{k}, placebos whose
#'   pre-intervention RMSPE exceeds \eqn{k} times that of the target are not
#'   drawn, as is common practice (Abadie et al., 2010). Default \code{NULL}
#'   (draw all). The p-value is not affected.
#' @param title,xlab,ylab Optional labels.
#' @param ... Ignored.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' sim <- simulate_iscm(seed = 7)
#' fit <- iscm(sim$data, "unit", "time", "Y", sim$treated, sim$affected, sim$t0)
#' pb  <- iscm_placebo(fit, sim$data, "unit", "time", "Y")
#' plot(pb)
#' plot(pb, type = "gaps")
#' @export
plot.iscm_placebo <- function(x, type = c("ratios", "gaps"), prune = NULL,
                              title = NULL, xlab = NULL, ylab = NULL, ...) {
  type <- match.arg(type)
  lev <- c("Target", "Placebo")
  cols <- stats::setNames(c("#2a78d6", .placebo_grey), lev)
  sub <- sprintf("Target: %s. Placebo p-value %.3f (smallest attainable %.3f).",
                 x$target, x$p_value, x$min_p_value)
  if (type == "ratios") {
    df <- x$table
    df$group <- factor(ifelse(df$unit == x$target, "Target", "Placebo"), levels = lev)
    df$unit <- stats::reorder(df$unit, df$ratio)
    return(
      ggplot2::ggplot(df, ggplot2::aes(x = .data$ratio, y = .data$unit, fill = .data$group)) +
        ggplot2::geom_col(width = 0.7) +
        ggplot2::scale_fill_manual(values = cols) +
        ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
        ggplot2::labs(title = if (is.null(title)) "Placebo test: post/pre RMSPE ratios" else title,
                      subtitle = sub,
                      x = if (is.null(xlab)) "Post- to pre-intervention RMSPE ratio" else xlab,
                      y = ylab) +
        theme_iscm() +
        ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
    )
  }
  if (is.null(x$gaps))
    stop("This placebo object does not store the gaps; recompute it with the ",
         "current version of the package.", call. = FALSE)
  keep <- colnames(x$gaps)
  if (!is.null(prune)) {
    tab <- x$table
    lim <- prune * tab$pre_rmspe[tab$unit == x$target]
    keep <- c(x$target, tab$unit[tab$unit != x$target & tab$pre_rmspe <= lim])
  }
  times <- as.numeric(x$times)
  df <- do.call(rbind, lapply(keep, function(u)
    data.frame(time = times, gap = as.numeric(x$gaps[, u]), unit = u,
               group = if (u == x$target) "Target" else "Placebo")))
  df$group <- factor(df$group, levels = rev(lev))       # draw the target last
  df <- df[order(df$group), ]
  ggplot2::ggplot(df, ggplot2::aes(x = .data$time, y = .data$gap, group = .data$unit,
                                   colour = .data$group, linewidth = .data$group)) +
    ggplot2::geom_hline(yintercept = 0, colour = .ink_muted, linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = .intervention_x(times, x$t0),
                        colour = .ink_muted, linetype = "dashed", linewidth = 0.4) +
    ggplot2::geom_line() +
    ggplot2::scale_colour_manual(values = cols, breaks = lev) +
    ggplot2::scale_linewidth_manual(values = c(Target = 1, Placebo = 0.5), breaks = lev) +
    ggplot2::labs(title = if (is.null(title)) "Placebo test: gaps of target and placebos" else title,
                  subtitle = if (is.null(prune)) sub else
                    paste(sub, sprintf("Placebos with pre-RMSPE > %g x target not drawn.", prune)),
                  x = if (is.null(xlab)) "Time" else xlab,
                  y = if (is.null(ylab)) "Gap (purged outcomes)" else ylab) +
    theme_iscm()
}
