## German-reunification application of Di Stefano and Mellace
## Reproduces the paper's iSCM estimates using the exact Synth specification.
## Requires: Synth, iscm.

library(Synth)
library(iscm)

data(germany)
x <- germany
## Tables and figures are written to the folder given by
## options(iscm.output_dir = "..."), or by the first command-line argument
## (Rscript germany_replication.R <folder>), or else to the subfolder
## "iscm_germany_replication" of the working directory. Runs in about 4 minutes.
args <- commandArgs(trailingOnly = TRUE)
figure_dir <- getOption("iscm.output_dir",
  if (length(args)) args[[1]] else file.path(getwd(), "iscm_germany_replication"))
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
key <- unique(x[, c("index", "country")])
named <- function(w) setNames(as.numeric(w),
  key$country[match(as.integer(rownames(w)), key$index)])
gap_of <- function(dp, so) as.numeric(dp$Y1plot - dp$Y0plot %*% so$solution.w)

## helpers for the two special-predictor specifications used in the paper
prep <- function(tr, controls, sp, tpp, ssr)
  dataprep(x, predictors = c("gdp","trade","infrate"), dependent = "gdp",
           unit.variable = 1, time.variable = 3, special.predictors = sp,
           treatment.identifier = tr, controls.identifier = controls,
           time.predictors.prior = tpp, time.optimize.ssr = ssr,
           unit.names.variable = 2, time.plot = 1960:2003)
sp_wg_train <- list(list("industry",1971:1980,"mean"),
                    list("schooling",c(1970,1975),"mean"), list("invest70",1980,"mean"))
sp_wg_main  <- list(list("industry",1981:1990,"mean"),
                    list("schooling",c(1980,1985),"mean"), list("invest80",1980,"mean"))
sp_au       <- list(list("industry",1971:1990,"mean"),
                    list("schooling",c(1970,1985),"mean"), list("invest80",1980,"mean"))
all_idx <- unique(x$index); WG <- 7; AU <- 3

## ---------------------------------------------------------------------------
## West Germany: V cross-validated on 1971-80, weights optimised on 1960-89
## ---------------------------------------------------------------------------
so_wg_tr <- synth(prep(WG, all_idx[-WG], sp_wg_train, 1971:1980, 1981:1990),
                  Margin.ipop = .005, Sigf.ipop = 7, Bound.ipop = 6)
dp_wg <- prep(WG, all_idx[-WG], sp_wg_main, 1981:1990, 1960:1989)
so_wg <- synth(dp_wg, custom.v = as.numeric(so_wg_tr$solution.v))
gap_wg <- gap_of(dp_wg, so_wg)

## West Germany restricted (Austria excluded)
so_wg_r_tr <- synth(prep(WG, all_idx[-c(WG,AU)], sp_wg_train, 1971:1980, 1981:1990),
                    Margin.ipop = .005, Sigf.ipop = 7, Bound.ipop = 6)
dp_wg_r <- prep(WG, all_idx[-c(WG,AU)], sp_wg_main, 1981:1990, 1960:1989)
so_wg_r <- synth(dp_wg_r, custom.v = as.numeric(so_wg_r_tr$solution.v))
gap_wg_r <- gap_of(dp_wg_r, so_wg_r)

## ---------------------------------------------------------------------------
## Austria: data-driven V
## ---------------------------------------------------------------------------
dp_au <- prep(AU, all_idx[-AU], sp_au, 1971:1990, 1960:1989)
so_au <- synth(dp_au, Margin.ipop = .005, Sigf.ipop = 7, Bound.ipop = 6)
gap_au <- gap_of(dp_au, so_au)

## Austria restricted (West Germany excluded)
dp_au_r <- prep(AU, all_idx[-c(AU,WG)], sp_au, 1971:1990, 1960:1989)
so_au_r <- synth(dp_au_r, Margin.ipop = .005, Sigf.ipop = 7, Bound.ipop = 6)
gap_au_r <- gap_of(dp_au_r, so_au_r)

## ---------------------------------------------------------------------------
## iSCM correction (generalises the paper's 2x2 determinant computation)
## ---------------------------------------------------------------------------
fit <- iscm_assemble(
  gaps    = list("West Germany" = gap_wg, "Austria" = gap_au),
  weights = list("West Germany" = named(so_wg$solution.w),
                 "Austria"      = named(so_au$solution.w)),
  t0 = 1989, times = 1960:2003,
  pure = setdiff(key$country, c("West Germany", "Austria")))

print(fit)                 # Omega, diagnostics, average effects
summary(fit)

## Gaps plotted against the unrestricted and restricted synthetic controls
rest_wg <- list(gap = gap_wg_r)      # minimal object accepted by plot.iscm
rest_au <- list(gap = gap_au_r)
print(plot(fit, unit = c("West Germany", "Austria"),
           restricted = list("West Germany" = rest_wg, "Austria" = rest_au),
           xlab = "Year", ylab = "Gap in GDP per capita (2002 USD)"))

## Decision rule (Section 5.1) on the actual Synth fits: predictor imbalance and
## pre-intervention RMSPE, unrestricted vs restricted, for both units.
## Step 2: the predictors are on very different scales, so we standardise each by
## its cross-sectional SD before aggregating (Section 5.1). The scale is fixed per
## treated unit from the unrestricted dataprep (all units) and reused for the
## restricted fit. imb_raw is kept only to report the old (unstandardised) figure.
imb_raw <- function(dp, so)    sum(abs(dp$X1 - dp$X0 %*% so$solution.w))
sd_ref  <- function(dp)        apply(cbind(dp$X1, dp$X0), 1, stats::sd)
imb_std <- function(dp, so, s) sum(abs(dp$X1 - dp$X0 %*% so$solution.w) / s)
s_wg <- sd_ref(dp_wg); s_au <- sd_ref(dp_au)
cmp <- iscm_compare_fits(
  gaps_unres = list("West Germany" = gap_wg,   "Austria" = gap_au),
  gaps_restr = list("West Germany" = gap_wg_r, "Austria" = gap_au_r),
  weights    = list("West Germany" = named(so_wg$solution.w),
                    "Austria"      = named(so_au$solution.w)),
  imbalance_unres = c("West Germany" = imb_std(dp_wg,   so_wg,   s_wg),
                      "Austria"      = imb_std(dp_au,   so_au,   s_au)),
  imbalance_restr = c("West Germany" = imb_std(dp_wg_r, so_wg_r, s_wg),
                      "Austria"      = imb_std(dp_au_r, so_au_r, s_au)),
  t0 = 1989, times = 1960:2003)
print(cmp)     # reproduces RMSPE 119 vs 261 for West Germany -> full iSCM

## Simplified iSCM: Austria is fit on the pure controls only. This is needed for
## Table C.1 and complements the full-iSCM object used above.
fit_simplified <- iscm_assemble(
  gaps = list("West Germany" = gap_wg, "Austria" = gap_au_r),
  weights = list("West Germany" = named(so_wg$solution.w),
                 "Austria" = named(so_au_r$solution.w)),
  t0 = 1989, times = 1960:2003,
  pure = setdiff(key$country, c("West Germany", "Austria")),
  simplified = TRUE)

## ---------------------------------------------------------------------------
## Export Figure 1: gaps and outcome trajectories
## ---------------------------------------------------------------------------
years <- 1960:2003
observed <- function(unit_name)
  x$gdp[x$country == unit_name][order(x$year[x$country == unit_name])]

draw_series <- function(y, series, cols, labels, ylab, panel, ylim = NULL) {
  if (is.null(ylim)) ylim <- range(c(y, unlist(series)), finite = TRUE)
  plot(years, y, type = "l", lwd = 2, col = cols[[1]], ylim = ylim,
       xlab = "Year", ylab = ylab)
  for (i in seq_along(series)) lines(years, series[[i]], lwd = 2, col = cols[[i + 1]])
  abline(v = 1990, lty = 3); if (grepl("Gap", ylab)) abline(h = 0, lty = 3)
  legend("topleft", legend = labels, col = cols, lwd = 2, bty = "n", cex = .72)
  mtext(panel, side = 3, line = .5, font = 2)
}

draw_figure1 <- function() {
  op <- par(mfrow = c(2, 2), mar = c(4.5, 5, 2.5, 1))
  on.exit(par(op), add = TRUE)
  for (u in c("West Germany", "Austria")) {
    gu <- if (u == "West Germany") gap_wg else gap_au
    gr <- if (u == "West Germany") gap_wg_r else gap_au_r
    gi <- as.numeric(fit$gaps_iscm[, u])
    draw_series(gi, list(gu, gr), c("black", "red", "purple"),
      c(paste("iSCM", u), paste("unrestricted SCM", u), paste("restricted SCM", u)),
      "Gap in per capita GDP (PPP, 2002 USD)",
      paste0(if (u == "West Germany") "(a) Gaps, " else "(b) Gaps, ", u),
      ylim = c(-4500, 4500))
  }
  for (u in c("West Germany", "Austria")) {
    yy <- observed(u)
    gu <- if (u == "West Germany") gap_wg else gap_au
    gr <- if (u == "West Germany") gap_wg_r else gap_au_r
    gi <- as.numeric(fit$gaps_iscm[, u])
    draw_series(yy, list(yy - gi, yy - gu, yy - gr),
      c("black", "blue", "red", "purple"),
      c(u, "inclusive synthetic", "unrestricted synthetic", "restricted synthetic"),
      "Per capita GDP (PPP, 2002 USD)",
      paste0(if (u == "West Germany") "(c) Trends, " else "(d) Trends, ", u),
      ylim = c(0, 32000))
  }
}

pdf(file.path(figure_dir, "figure_1_estimated_effects.pdf"),
    width = 11, height = 9, useDingbats = FALSE)
draw_figure1(); dev.off()
png(file.path(figure_dir, "figure_1_estimated_effects.png"),
    width = 2200, height = 1800, res = 200)
draw_figure1(); dev.off()

## ---------------------------------------------------------------------------
## Export Tables 1-3 and reusable numerical results
## ---------------------------------------------------------------------------
table1_block <- function(unit, dp_u, so_u, dp_r, so_r) {
  obs <- as.numeric(dp_u$X1)
  unres <- as.numeric(dp_u$X0 %*% so_u$solution.w)
  restr <- as.numeric(dp_r$X0 %*% so_r$solution.w)
  data.frame(unit = unit, predictor = rownames(dp_u$X1), observed = obs,
             unrestricted_synthetic = unres, restricted_synthetic = restr,
             unrestricted_discrepancy = abs(obs - unres),
             restricted_discrepancy = abs(obs - restr), check.names = FALSE)
}
table1 <- rbind(
  table1_block("West Germany", dp_wg, so_wg, dp_wg_r, so_wg_r),
  table1_block("Austria", dp_au, so_au, dp_au_r, so_au_r))
write.csv(table1, file.path(figure_dir, "table_1_predictors.csv"), row.names = FALSE)

country_names <- key$country
weight_col <- function(w) {
  z <- setNames(rep(0, length(country_names)), country_names)
  z[names(named(w))] <- named(w)
  z
}
table2 <- data.frame(
  country = country_names,
  west_germany_unrestricted = weight_col(so_wg$solution.w),
  west_germany_restricted = weight_col(so_wg_r$solution.w),
  austria_unrestricted = weight_col(so_au$solution.w),
  austria_restricted = weight_col(so_au_r$solution.w),
  row.names = NULL, check.names = FALSE)
table2[table2$country == "West Germany", c("west_germany_unrestricted", "west_germany_restricted")] <- NA
table2[table2$country == "Austria", c("austria_unrestricted", "austria_restricted")] <- NA
write.csv(table2, file.path(figure_dir, "table_C2_weights.csv"), row.names = FALSE)

effects_frame <- function(unit) {
  unres <- if (unit == "West Germany") gap_wg[31:44] else gap_au[31:44]
  restr <- if (unit == "West Germany") gap_wg_r[31:44] else gap_au_r[31:44]
  data.frame(unit = unit, year = 1990:2003, unrestricted_SCM = unres,
             restricted_SCM = restr, full_iSCM = as.numeric(fit$effects[unit, ]),
             simplified_iSCM = as.numeric(fit_simplified$effects[unit, ]),
             check.names = FALSE)
}
effects_annual <- rbind(effects_frame("West Germany"), effects_frame("Austria"))
write.csv(effects_annual, file.path(figure_dir, "effects_annual_1990_2003.csv"),
          row.names = FALSE)

selected <- effects_annual[effects_annual$year %in% c(1990, 1995, 2000, 2003), ]
averages <- aggregate(cbind(unrestricted_SCM, restricted_SCM, full_iSCM,
                            simplified_iSCM) ~ unit,
                      effects_annual, mean)
averages$year <- "Average 1990-2003"
selected$year <- as.character(selected$year)
table3 <- rbind(selected, averages[names(selected)])
write.csv(table3, file.path(figure_dir, "table_C1_effects.csv"),
          row.names = FALSE)

diagnostics <- data.frame(
  determinant = fit$diagnostics$det,
  sigma_min = fit$diagnostics$sigma_min,
  amplification = fit$diagnostics$amplification,
  alpha = fit$diagnostics$alpha,
  diagonally_dominant = fit$diagnostics$diag_dominant)
write.csv(diagnostics, file.path(figure_dir, "iscm_diagnostics.csv"), row.names = FALSE)
write.csv(cmp$comparison, file.path(figure_dir, "specification_comparison.csv"),
          row.names = FALSE)
saveRDS(list(fit = fit, fit_simplified = fit_simplified, comparison = cmp,
             table1 = table1, table_C2 = table2, table_C1 = table3,
             effects_annual = effects_annual),
        file.path(figure_dir, "germany_replication_results.rds"))

## In-space placebo (Online Appendix B). The intervention is reassigned to each
## pure control in turn, using the same Synth specification as the unit being
## tested, on outcomes purged of the estimated iSCM effects. The specification
## uses predictors averaged over windows that include 1990, so the placebo
## weights are estimated on the purged data returned by iscm_purge().
## iscm_placebo_fits() then computes the RMSPE ratios and the p-value. The
## reference distribution contains the 15 pure controls only.
pure_idx <- setdiff(all_idx, c(WG, AU))
x_purged <- iscm_purge(fit, x, id = "country", time = "year", outcome = "gdp")

## `prep_data` is the data-argument version of `prep`.
prep_data <- function(dat, tr, controls, sp, tpp, ssr)
  dataprep(dat, predictors = c("gdp", "trade", "infrate"), dependent = "gdp",
           unit.variable = 1, time.variable = 3, special.predictors = sp,
           treatment.identifier = tr, controls.identifier = controls,
           time.predictors.prior = tpp, time.optimize.ssr = ssr,
           unit.names.variable = 2, time.plot = 1960:2003)

## West-Germany specification for a reassigned treated unit: predictor weights
## selected with the same sample-split procedure, donor weights on 1960-89.
weights_placebo_wg <- function(tr) {
  controls <- setdiff(all_idx, tr)
  so_tr <- synth(
    prep_data(x_purged, tr, controls, sp_wg_train, 1971:1980, 1981:1990),
    Margin.ipop = .005, Sigf.ipop = 7, Bound.ipop = 6)
  dp <- prep_data(x_purged, tr, controls, sp_wg_main, 1981:1990, 1960:1989)
  named(synth(dp, custom.v = as.numeric(so_tr$solution.v))$solution.w)
}

## Austria specification, including its data-driven V selection.
weights_placebo_au <- function(tr) {
  controls <- setdiff(all_idx, tr)
  dp <- prep_data(x_purged, tr, controls, sp_au, 1971:1990, 1960:1989)
  named(synth(dp, Margin.ipop = .005, Sigf.ipop = 7, Bound.ipop = 6)$solution.w)
}

pure_names <- key$country[match(pure_idx, key$index)]

## Synth is verbose; retain its diagnostics in a log.
inference_log <- tempfile("germany-placebo-", fileext = ".log")
invisible(capture.output({
  pw_wg <- setNames(lapply(pure_idx, weights_placebo_wg), pure_names)
  pw_au <- setNames(lapply(pure_idx, weights_placebo_au), pure_names)
}, file = inference_log, type = "output"))

pb_wg <- iscm_placebo_fits(fit, x, "country", "year", "gdp",
                           placebo_weights = pw_wg, target = "West Germany")
pb_au <- iscm_placebo_fits(fit, x, "country", "year", "gdp",
                           placebo_weights = pw_au, target = "Austria")
print(pb_wg)
print(pb_au)

cat(sprintf("West Germany: RMSPE ratio %.4f, placebo p-value %d/16 = %.4f\n",
            pb_wg$ratio, round(pb_wg$p_value * 16), pb_wg$p_value))
cat(sprintf("Austria: RMSPE ratio %.4f, placebo p-value %d/16 = %.4f\n",
            pb_au$ratio, round(pb_au$p_value * 16), pb_au$p_value))

## named vectors of ratios for the figure and the exported table
pb_wg$ratios <- setNames(pb_wg$table$ratio, pb_wg$table$unit)
pb_au$ratios <- setNames(pb_au$table$ratio, pb_au$table$unit)

draw_placebo_figure <- function() {
  op <- par(mfrow = c(1, 2), mar = c(6, 7, 1, 1), oma = c(0, 0, 0, 0))
  on.exit(par(op), add = TRUE)
  for (z in list(pb_wg, pb_au)) {
    ord <- order(z$ratios)
    vals <- z$ratios[ord]
    labs <- names(vals)
    cols <- ifelse(labs == z$target, "black", "grey25")
    pch <- ifelse(labs == z$target, 19, 16)
    plot(vals, seq_along(vals), yaxt = "n", ylab = "", xlab = "",
         pch = pch, col = cols, xlim = c(0, max(vals) * 1.05))
    axis(2, at = seq_along(vals), labels = labs, las = 1, cex.axis = .72)
    abline(h = seq_along(vals), col = "grey88", lty = 3)
    points(vals, seq_along(vals), pch = pch, col = cols)
    mtext("RMSPE ratio", side = 1, line = 2.2)
    mtext(paste0("(", if (z$target == "West Germany") "a" else "b",
                 ") ", z$target), side = 1, line = 4.2)
  }
}

pdf(file.path(figure_dir, "figure_C1_placebo_RMSPE.pdf"),
    width = 10, height = 6.5, useDingbats = FALSE)
draw_placebo_figure()
dev.off()

png(file.path(figure_dir, "figure_C1_placebo_RMSPE.png"),
    width = 2000, height = 1300, res = 200)
draw_placebo_figure()
dev.off()

write.csv(data.frame(
  unit = c(names(pb_wg$ratios), names(pb_au$ratios)),
  target = rep(c("West Germany", "Austria"),
               c(length(pb_wg$ratios), length(pb_au$ratios))),
  rmspe_ratio = c(pb_wg$ratios, pb_au$ratios),
  included_in_p_value = c(names(pb_wg$ratios) %in%
                            key$country[match(pure_idx, key$index)],
                          names(pb_au$ratios) %in%
                            key$country[match(pure_idx, key$index)])),
  file.path(figure_dir, "figure_C1_placebo_RMSPE.csv"), row.names = FALSE)

message("Tables and figures saved in: ", normalizePath(figure_dir))
