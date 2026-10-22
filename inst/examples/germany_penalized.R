## ===========================================================================
##  German reunification with the INCLUSIVE PENALIZED SCM  (Online Appendix C)
##  Di Stefano and Mellace, "The inclusive synthetic control method"
##
##  iSCM is estimator-agnostic: it corrects the contamination on top of ANY
##  SC-type set of weights. Here the backend is the penalized SCM of Abadie and
##  L'Hour (2021), estimated with the 'pensynth' package (cross-validated penalty
##  via cv_pensynth(), clarabel solver), retaining the predictor-weighting matrix
##  V from the corresponding Synth training fit -- exactly the specification used
##  for Figure C.2 of the paper. We then assemble the inclusive penalized SCM with
##  iscm_assemble(), which builds Omega from the cross-weights and solves the
##  system (the m-unit generalisation of the 2x2 determinant-ratio computation in
##  the paper's replication script).
##
##  Requires: Synth, pensynth, clarabel, iscm.
## ===========================================================================

for (pkg in c("Synth", "pensynth", "clarabel")) if (!requireNamespace(pkg, quietly = TRUE))
  stop(sprintf("This example needs the '%s' package; install it first.", pkg), call. = FALSE)
library(Synth)
library(pensynth)
library(iscm)

data(germany)
x   <- germany
args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args)) args[[1]] else getwd()
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
key <- unique(x[, c("index", "country")])
WG <- 7; AU <- 3; t0 <- 1989; years <- 1960:2003
all_idx <- unique(x$index)

## same specification as germany_replication.R ------------------------------
prep <- function(tr, controls, sp, tpp, ssr)
  dataprep(x, predictors = c("gdp", "trade", "infrate"), dependent = "gdp",
           unit.variable = 1, time.variable = 3, special.predictors = sp,
           treatment.identifier = tr, controls.identifier = controls,
           time.predictors.prior = tpp, time.optimize.ssr = ssr,
           unit.names.variable = 2, time.plot = years)
sp_wg_train <- list(list("industry", 1971:1980, "mean"),
                    list("schooling", c(1970, 1975), "mean"), list("invest70", 1980, "mean"))
sp_wg_main  <- list(list("industry", 1981:1990, "mean"),
                    list("schooling", c(1980, 1985), "mean"), list("invest80", 1980, "mean"))
sp_au       <- list(list("industry", 1971:1990, "mean"),
                    list("schooling", c(1970, 1985), "mean"), list("invest80", 1980, "mean"))

## name a donor-weight vector by country, from the dataprep's X0 columns ------
named_w <- function(w, dp) {
  idx <- suppressWarnings(as.integer(colnames(dp$X0)))
  nm  <- if (all(is.na(idx))) colnames(dp$X0) else key$country[match(idx, key$index)]
  setNames(as.numeric(w), nm)
}
gap_of <- function(dp, w) as.numeric(dp$Y1plot - dp$Y0plot %*% as.numeric(w))

## penalized SCM with cross-validated lambda (pensynth), given a Synth V -------
pen_fit <- function(dp, v) {
  cv <- with(dp, pensynth::cv_pensynth(
    X1, X0, Z1, Z0, v = as.numeric(v), nlambda = 1000,
    opt_pars = clarabel::clarabel_control(), standardize = TRUE,
    return_solver_info = FALSE))
  list(w = named_w(cv$w_opt, dp), gap = gap_of(dp, cv$w_opt), lambda = cv$l_opt)
}

## ---------------------------------------------------------------------------
##  West Germany (Austria IN): sample-split V, penalized weights
## ---------------------------------------------------------------------------
v_wg <- synth(prep(WG, all_idx[-WG], sp_wg_train, 1971:1980, 1981:1990),
              Margin.ipop = .005, Sigf.ipop = 7, Bound.ipop = 6)$solution.v
dp_wg  <- prep(WG, all_idx[-WG], sp_wg_main, 1981:1990, 1960:1989)
f_wg   <- pen_fit(dp_wg, v_wg)

## West Germany (Austria OUT): its own sample-split V
v_wg_r <- synth(prep(WG, all_idx[-c(WG, AU)], sp_wg_train, 1971:1980, 1981:1990),
                Margin.ipop = .005, Sigf.ipop = 7, Bound.ipop = 6)$solution.v
dp_wg_r <- prep(WG, all_idx[-c(WG, AU)], sp_wg_main, 1981:1990, 1960:1989)
f_wg_r  <- pen_fit(dp_wg_r, v_wg_r)

## ---------------------------------------------------------------------------
##  Austria (West Germany IN): data-driven V, penalized weights
## ---------------------------------------------------------------------------
dp_au  <- prep(AU, all_idx[-AU], sp_au, 1971:1990, 1960:1989)
v_au   <- synth(dp_au, Margin.ipop = .005, Sigf.ipop = 7, Bound.ipop = 6)$solution.v
f_au   <- pen_fit(dp_au, v_au)

## Austria (West Germany OUT): data-driven V
dp_au_r <- prep(AU, all_idx[-c(AU, WG)], sp_au, 1971:1990, 1960:1989)
v_au_r  <- synth(dp_au_r, Margin.ipop = .005, Sigf.ipop = 7, Bound.ipop = 6)$solution.v
f_au_r  <- pen_fit(dp_au_r, v_au_r)

cat(sprintf("Selected penalties (cv_pensynth): lambda_WG = %.4g, lambda_AU = %.4g\n",
            f_wg$lambda, f_au$lambda))

## ---------------------------------------------------------------------------
##  Inclusive penalized SCM: assemble Omega from the cross-weights and solve
## ---------------------------------------------------------------------------
fit_pen <- iscm_assemble(
  gaps    = list("West Germany" = f_wg$gap, "Austria" = f_au$gap),
  weights = list("West Germany" = f_wg$w,   "Austria" = f_au$w),
  t0 = t0, times = years,
  pure = setdiff(key$country, c("West Germany", "Austria")))

cat("\n=================  INCLUSIVE PENALIZED SCM  =================\n")
print(fit_pen)

## Decision rule (Section 5.1) with the STANDARDISED Step 2 imbalance, as in the
## original-SCM replication: divide each predictor discrepancy by its
## cross-sectional SD before aggregating; fix the scale per treated unit.
sd_ref  <- function(dp)           apply(cbind(dp$X1, dp$X0), 1, stats::sd)
imb_std <- function(dp, w, s)     sum(abs(dp$X1 - dp$X0 %*% as.numeric(w)) / s)
s_wg <- sd_ref(dp_wg); s_au <- sd_ref(dp_au)
cmp_pen <- iscm_compare_fits(
  gaps_unres = list("West Germany" = f_wg$gap,   "Austria" = f_au$gap),
  gaps_restr = list("West Germany" = f_wg_r$gap, "Austria" = f_au_r$gap),
  weights    = list("West Germany" = f_wg$w,     "Austria" = f_au$w),
  imbalance_unres = c("West Germany" = imb_std(dp_wg,   f_wg$w,   s_wg),
                      "Austria"      = imb_std(dp_au,   f_au$w,   s_au)),
  imbalance_restr = c("West Germany" = imb_std(dp_wg_r, f_wg_r$w, s_wg),
                      "Austria"      = imb_std(dp_au_r, f_au_r$w, s_au)),
  t0 = t0, times = years)
print(cmp_pen)

## inclusive penalized gaps (pre-period = unrestricted, post = corrected)
gapWG_ipen <- f_wg$gap; gapWG_ipen[years > t0] <- fit_pen$effects["West Germany", ]
gapAU_ipen <- f_au$gap; gapAU_ipen[years > t0] <- fit_pen$effects["Austria", ]

## ---------------------------------------------------------------------------
##  FIGURE 2 (a)-(b): penalized gaps, inclusive vs unrestricted vs restricted
## ---------------------------------------------------------------------------
plot_gaps <- function(main, g_incl, g_unres, g_restr) {
  plot(years, g_incl, type = "l", col = "black", lwd = 2, ylim = c(-4500, 4500),
       xlab = "Year", ylab = "Gap in Per capita GDP (PPP, 2002 USD)", main = main)
  lines(years, g_unres, col = "red",    lwd = 2)
  lines(years, g_restr, col = "purple", lwd = 2)
  abline(v = 1990, lty = "dotted"); abline(h = 0, lty = "dotted")
  legend("topleft", bty = "n", cex = 0.8, lwd = 2,
         col = c("black", "red", "purple"),
         legend = c("inclusive penalized", "unrestricted penalized", "restricted penalized"))
}
op <- par(mfrow = c(1, 2))
plot_gaps("WEST GERMANY - Penalized SCM", gapWG_ipen, f_wg$gap, f_wg_r$gap)
plot_gaps("AUSTRIA - Penalized SCM",      gapAU_ipen, f_au$gap, f_au_r$gap)
par(op)

## ---------------------------------------------------------------------------
##  FIGURE 2 (c)-(d): inclusive PENALIZED vs inclusive iSCM (original SCM)
## ---------------------------------------------------------------------------
named_synth <- function(so, dp) named_w(so$solution.w, dp)
so_wg <- synth(dp_wg, custom.v = as.numeric(v_wg))
so_au <- synth(dp_au, Margin.ipop = .005, Sigf.ipop = 7, Bound.ipop = 6)
fit_iscm <- iscm_assemble(
  gaps    = list("West Germany" = gap_of(dp_wg, so_wg$solution.w),
                 "Austria"      = gap_of(dp_au, so_au$solution.w)),
  weights = list("West Germany" = named_synth(so_wg, dp_wg),
                 "Austria"      = named_synth(so_au, dp_au)),
  t0 = t0, times = years,
  pure = setdiff(key$country, c("West Germany", "Austria")))
gapWG_iscm <- fit_iscm$gaps_unrestricted[, "West Germany"]
gapWG_iscm[years > t0] <- fit_iscm$effects["West Germany", ]
gapAU_iscm <- fit_iscm$gaps_unrestricted[, "Austria"]
gapAU_iscm[years > t0] <- fit_iscm$effects["Austria", ]

op <- par(mfrow = c(1, 2))
plot(years, gapWG_ipen, type = "l", col = "black", lwd = 2, ylim = c(-4500, 4500),
     xlab = "Year", ylab = "Gap", main = "West Germany: penalized vs iSCM")
lines(years, gapWG_iscm, col = "red", lwd = 2)
abline(v = 1990, lty = "dotted"); abline(h = 0, lty = "dotted")
legend("topleft", bty = "n", cex = 0.8, lwd = 2, col = c("black", "red"),
       legend = c("inclusive penalized", "inclusive (original SCM)"))
plot(years, gapAU_ipen, type = "l", col = "black", lwd = 2, ylim = c(-4500, 4500),
     xlab = "Year", ylab = "Gap", main = "Austria: penalized vs iSCM")
lines(years, gapAU_iscm, col = "red", lwd = 2)
abline(v = 1990, lty = "dotted"); abline(h = 0, lty = "dotted")
legend("topleft", bty = "n", cex = 0.8, lwd = 2, col = c("black", "red"),
       legend = c("inclusive penalized", "inclusive (original SCM)"))
par(op)

## ---------------------------------------------------------------------------
## Export the complete four-panel Figure C.2 and numerical penalized results
## ---------------------------------------------------------------------------
draw_pen_panel <- function(gaps, cols, labels, panel) {
  plot(years, gaps[[1]], type = "l", col = cols[[1]], lwd = 2,
       ylim = c(-4500, 4500), xlab = "Year",
       ylab = "Gap in per capita GDP (PPP, 2002 USD)")
  if (length(gaps) > 1L)
    for (i in 2:length(gaps)) lines(years, gaps[[i]], col = cols[[i]], lwd = 2)
  abline(v = 1990, lty = 3); abline(h = 0, lty = 3)
  legend("topleft", bty = "n", cex = .7, lwd = 2, col = cols, legend = labels)
  mtext(panel, side = 3, line = .5, font = 2)
}

draw_figure2 <- function() {
  op <- par(mfrow = c(2, 2), mar = c(4.5, 5, 2.5, 1))
  on.exit(par(op), add = TRUE)
  draw_pen_panel(list(gapWG_ipen, f_wg$gap, f_wg_r$gap),
    c("black", "red", "purple"),
    c("inclusive penalized", "unrestricted penalized", "restricted penalized"),
    "(a) Gaps, West Germany, penalized SCM")
  draw_pen_panel(list(gapAU_ipen, f_au$gap, f_au_r$gap),
    c("black", "red", "purple"),
    c("inclusive penalized", "unrestricted penalized", "restricted penalized"),
    "(b) Gaps, Austria, penalized SCM")
  draw_pen_panel(list(gapWG_iscm, gapWG_ipen), c("red", "black"),
    c("iSCM", "inclusive penalized"),
    "(c) Penalized iSCM vs iSCM, West Germany")
  draw_pen_panel(list(gapAU_iscm, gapAU_ipen), c("red", "black"),
    c("iSCM", "inclusive penalized"),
    "(d) Penalized iSCM vs iSCM, Austria")
}

pdf(file.path(output_dir, "figure_2_penalized_SCM.pdf"),
    width = 11, height = 9, useDingbats = FALSE)
draw_figure2(); dev.off()
png(file.path(output_dir, "figure_2_penalized_SCM.png"),
    width = 2200, height = 1800, res = 200)
draw_figure2(); dev.off()

penalized_effects <- rbind(
  data.frame(unit = "West Germany", year = years,
             unrestricted_penalized = f_wg$gap,
             restricted_penalized = f_wg_r$gap,
             inclusive_penalized = gapWG_ipen,
             original_iSCM = gapWG_iscm),
  data.frame(unit = "Austria", year = years,
             unrestricted_penalized = f_au$gap,
             restricted_penalized = f_au_r$gap,
             inclusive_penalized = gapAU_ipen,
             original_iSCM = gapAU_iscm))
write.csv(penalized_effects,
          file.path(output_dir, "penalized_effects_1960_2003.csv"), row.names = FALSE)
write.csv(data.frame(unit = c("West Germany", "Austria"),
                     selected_lambda = c(f_wg$lambda, f_au$lambda)),
          file.path(output_dir, "penalized_selected_lambdas.csv"), row.names = FALSE)
write.csv(cmp_pen$comparison,
          file.path(output_dir, "penalized_specification_comparison.csv"),
          row.names = FALSE)
saveRDS(list(fit_penalized = fit_pen, comparison = cmp_pen,
             effects = penalized_effects,
             lambdas = c("West Germany" = f_wg$lambda, "Austria" = f_au$lambda)),
        file.path(output_dir, "germany_penalized_results.rds"))
