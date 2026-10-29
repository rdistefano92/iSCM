# iscm 0.1.3

## Plots

* All plots are now drawn with `ggplot2` (new required dependency) and return
  `ggplot` objects that can be customised with the usual `ggplot2` functions.
  Colours are colour-blind safe and every series also has its own line type.
* `plot.iscm()` gains `type`: `"gaps"` (default, Figure 1 style), `"trends"`
  (observed against synthetic outcomes) and `"weights"` (donor weights, with the
  implied pure-control weights of Lemma 3.1). Several focal units can be drawn
  at once in separate panels (`unit = c(...)`).
* `plot.iscm_placebo()` gains `type = "gaps"`, the gap paths of the target and
  of the placebos, with optional pruning of poorly fitted placebos (`prune`).
  Placebo objects now store these gaps.
* New `theme_iscm()`.

## Documentation

* README and vignette: the three groups of units (main treated, potentially
  affected, pure controls) are defined before they are used; the arguments
  `id`, `time`, `outcome` and `predictors` are listed with the others; full
  references (Abadie, Diamond and Hainmueller 2010, 2011, 2015; Abadie and
  L'Hour 2021) and a link to the paper on arXiv; the simulation design of
  `simulate_iscm_covariates()` is explained step by step.
* The `iscm_assemble()` example in the README is now self-contained and runs as
  it is; the "typical workflow" block no longer depends on objects created in
  earlier blocks.
* The replication scripts write their output to the folder set with
  `options(iscm.output_dir = ...)`, or to a new subfolder of the working
  directory, instead of the working directory itself, and say where. Output
  files are named after the paper (`table_C1_effects.csv`,
  `table_C2_weights.csv`, `figure_C1_placebo_RMSPE.*`, `figure_C2_penalized_SCM.*`).

# iscm 0.1.2

## New features

* `scm_penalized()`: dependency-free backend for the penalized SCM of Abadie and
  L'Hour (2021), with the penalty passed as `lambda`.
* `iscm_purge()` and `iscm_placebo_fits()`: placebo analysis (Online Appendix B)
  for synthetic controls estimated outside the package. `iscm_purge()` returns
  the data with the estimated effects removed from the focal units; the user
  estimates each pure control's synthetic control on it with any estimator, and
  `iscm_placebo_fits()` computes the RMSPE ratios, the p-value and the plot with
  the same code as `iscm_placebo()`. The German replication script now uses
  them and reproduces the published placebo results (1/16 and 15/16).
* Printed placebo results now report the smallest attainable p-value and the
  conditioning of Omega, as recommended in Online Appendix B; both are also
  returned (`min_p_value`, `alpha`, `sigma_min`, `amplification`).
* `simulate_iscm_covariates()`: a simulation design in which units are matched on
  observed covariates and the treated and affected units share an unobserved
  characteristic. Here iSCM clearly outperforms the restricted SCM (Section 3.1).
  `inst/examples/simulation_covariates.R` reproduces the Monte Carlo summary
  reported in the README and the vignette.

## Changes in behaviour

* `iscm_placebo()` on an object that does not store its estimator (from
  `iscm_assemble()`) now stops with a message pointing to `iscm_placebo_fits()`,
  instead of warning and silently re-estimating with `scm_base()`.
* **Decision rule, Step 4.** `iscm_compare()` and `iscm_compare_fits()` now
  choose between the full and the simplified iSCM by comparing the
  pre-intervention RMSPE of their implied pure-control counterfactuals, the
  observable counterparts of the error expressions (4.6) and (4.7), instead of
  comparing the amplification factor with the fixed threshold `alpha_tol = 2`.
  The values are returned in `implied_rmspe` and printed in Step 4. The previous
  rule is available with `rule = "alpha"`. `iscm_decision()` gains
  `rmspe_full` and `rmspe_simplified`; without them it uses the previous rule.
  For the German application the new rule gives 138.54 (full) against 144.73
  (simplified) and still recommends the full iSCM.
* **Matching on the outcome path.** When `predictors = NULL`, the built-in
  backends (`scm_base`, `scm_penalized`) are now called with
  `standardize = FALSE`. Standardising every pre-intervention period by its
  cross-sectional standard deviation gave the most weight to the periods in which
  units are most alike, where noise dominates. Pass `standardize = TRUE` to
  restore the previous behaviour. The results of the paper are unaffected: its
  application uses `Synth` through `iscm_assemble()`.

## Documentation

* New package help page (`package?iscm`) with a map from tasks to commands and the
  two ways of using the package (`iscm()` versus `iscm_assemble()`).
* Every exported function now has runnable examples and links to related
  functions. `?iscm` documents the components of the returned object.
* `?backends` explains how to plug in another SC-type estimator, with a worked
  custom backend.
* The vignette is reorganised as a step-by-step workflow (estimate, diagnostics,
  choice of specification, placebo, sensitivity, other estimators, German
  replication) and reports the four specifications side by side.
* The README is rewritten around what users can do and which command to use.
* `?iscm_compare` states that the thresholds are package defaults and how they
  relate to Section 5.1.
* `?iscm_placebo` and `?iscm_assemble` explain that the placebo analysis
  re-estimates fits with the stored backend, so objects from `iscm_assemble()`
  need the loop in `inst/examples/germany_replication.R`.
* References updated to the current manuscript (placebo: Online Appendix B;
  penalized SCM: Online Appendix C; Tables C.1-C.2, Figures C.1-C.2).

## Other changes

* `NAMESPACE` is now generated by roxygen2.
* `DESCRIPTION` lists Roberta Di Stefano once, as author and maintainer.
* Removed a redundant line in `omega_diagnostics()` (no change in results).

# iscm 0.1.1

## Step 2 imbalance standardisation (Section 5.1)

* `iscm_compare()` gains a `standardize` argument (default `TRUE`). The Step 2
  predictor imbalance is now computed by dividing each predictor's absolute
  discrepancy by the predictor's cross-sectional standard deviation before
  aggregation, as recommended in Section 5.1 when the predictors are on very
  different scales. The scale is computed across all units and held fixed across
  the unrestricted and restricted fits, so the two imbalances remain comparable.
  Set `standardize = FALSE` for the previous unstandardised sum of absolute
  discrepancies.
* The documentation of `imbalance_unres` / `imbalance_restr` in
  `iscm_compare_fits()` now shows the standardised formula
  (`s <- apply(cbind(X1, X0), 1, sd); sum(abs(X1 - X0 %*% w) / s)`), used with the
  same `s` for the unrestricted and restricted fits.
* The shipped example `inst/examples/germany_replication.R` now standardises the
  Step 2 imbalance for both West Germany and Austria.

## New example

* Added `inst/examples/germany_penalized.R`, the inclusive **penalized** SCM of
  Section 7.1. It reproduces the paper's Figure 2 specification: penalized SCM of
  Abadie and L'Hour (2021) via `pensynth::cv_pensynth()` (clarabel solver,
  cross-validated penalty, `standardize = TRUE`), retaining the `Synth` training
  `V`, with the inclusive correction assembled by `iscm_assemble()`. Adds
  `pensynth` and `clarabel` to Suggests.

## Fixes

* Fixed a parse error in `tests/testthat/test-compare.R` (the `iscm_sensitivity`
  block was missing its `test_that()` wrapper).
* Corrected `tests/testthat/test-diagnostics.R`: the matrix
  `[[1, -0.9], [-0.9, 1]]` is strictly diagonally dominant, so `diag_dominant`
  is `TRUE`; the test now checks that a diagonally dominant system can still be
  poorly conditioned (large `alpha`) and adds a genuinely non-dominant case.
* Added `tests/testthat/test-standardize.R` covering the new `standardize`
  argument (scale invariance and the raw-vs-standardised difference).
