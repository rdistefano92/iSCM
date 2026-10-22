# iscm — The Inclusive Synthetic Control Method

An R implementation of the **inclusive synthetic control method (iSCM)** of
Di Stefano and Mellace. iSCM lets units that are *potentially affected* by an
intervention (other treated units, or controls exposed to spillovers) remain in
the donor pool, and removes the resulting contamination exactly, on top of any
synthetic-control-type estimator.

## Installation

```r
# install.packages("remotes")
remotes::install_github("rdistefano92/iscm", build_vignettes = TRUE)
```

The package has no required dependencies beyond base R. To use the classic SCM of
Abadie et al. (as in the paper), also install `Synth`; the penalized-SCM example
additionally needs `pensynth` and `clarabel`.

## What you need

A balanced long-format panel (one row per unit and period) and three choices:

- `treated`: the main treated unit;
- `affected`: the units potentially affected by the intervention (the set *M*);
  every other unit is treated as a pure control;
- `t0`: the last pre-intervention period.

## Quick start

```r
library(iscm)

# simulated panel: true effect on TR = 2, spillover on A1 = 3
sim <- simulate_iscm(n_affected = 1, n_pure = 6, tau = 2, spill = 3, seed = 7)

fit <- iscm(sim$data, id = "unit", time = "time", outcome = "Y",
            treated = "TR", affected = "A1", t0 = sim$t0)

fit                    # average effects + conditioning of Omega
summary(fit)           # effects period by period
plot(fit)              # inclusive vs unrestricted gaps
fit$effects            # estimates as a matrix (row 1 = treated, then affected)
```

## What can I do, and with which command?

| I want to... | Command |
|---|---|
| estimate the **full iSCM** | `iscm(data, id, time, outcome, treated, affected, t0)` |
| estimate the **simplified iSCM** | `iscm(..., simplified = TRUE)` |
| estimate the **restricted SCM** (affected units dropped) | `restricted_scm(...)` (same arguments) |
| get the **unrestricted** (uncorrected) SCM | `fit$gaps_unrestricted` |
| see the effects | `print(fit)`, `summary(fit)`, `fit$effects` |
| plot the gaps (Figure 1 style) | `plot(fit, unit = "A1", restricted = rest)` |
| check the **conditioning** of the correction (det, σ<sub>min</sub>, α) | `fit$diagnostics`, `omega_diagnostics(Omega)` |
| see the implied weights on pure controls (Lemma 3.1) | `fit$implied_weights` |
| **choose** among restricted / full / simplified (Section 5.1) | `iscm_compare(...)` |
| same, for fits estimated outside the package | `iscm_compare_fits()`, `iscm_decision()` |
| run the **placebo** analysis (Online Appendix B) | `iscm_placebo(fit, data, id, time, outcome, target = ...)` |
| same, for fits estimated outside the package | `iscm_purge()`, then `iscm_placebo_fits(fit, ..., placebo_weights)` |
| check **sensitivity** to the set of affected units (Section 5.2) | `iscm_sensitivity(..., affected_sets = list(...))` |
| use **another SC-type estimator** | `iscm(..., estimator = scm_penalized)`, `scm_synth`, or your own function (`?backends`) |
| correct synthetic controls **estimated elsewhere** | `iscm_assemble(gaps, weights, t0)` |
| apply only the linear algebra (Theorem 2.1) | `build_omega()`, `solve_iscm()` |
| simulate data with known effects | `simulate_iscm()`, `simulate_iscm_covariates()` |

All data-level functions (`iscm`, `restricted_scm`, `iscm_compare`,
`iscm_sensitivity`) take the same first arguments, so a typical analysis is the
same call repeated with a different function name.

## When does iSCM beat the restricted SCM?

The iSCM counterfactual is itself a weighted average of pure controls
(Lemma 3.1), so iSCM gains when the restricted optimizer, guided by its own
criterion, picks a worse combination of pure controls than the one iSCM implies
(Section 3.1). `simulate_iscm_covariates()` builds such a case: units are matched
on observed covariates, the treated and the affected unit share an unobserved
characteristic, and the pure controls closest in covariates do not.

```r
sim  <- simulate_iscm_covariates(seed = 1)     # true effects: 2 on TR, 3 on A1
args <- list(sim$data, "unit", "time", "Y", "TR", "A1", sim$t0,
             predictors = c("x1", "x2"), estimator = scm_penalized)
fit  <- do.call(iscm, args)
rest <- do.call(restricted_scm, args)
```

Over 100 replications the average estimate of the effect on `TR` (true value 2) is
0.51 for the unrestricted SCM, 4.06 for the restricted SCM, 3.20 for the
simplified iSCM and 2.00 for the full iSCM, and `iscm_compare()` recommends the
full iSCM in every replication. The script
`system.file("examples", "simulation_covariates.R", package = "iscm")` reproduces
these figures.

## Two ways to use the package

**1. Let the package estimate the synthetic controls: `iscm()`.** Each unit is
matched on the pre-intervention outcome path (default) or on the columns given in
`predictors`, averaged over the pre-intervention period. Choose the estimator with
`estimator =`:

- `scm_base` (default): classic simplex-constrained SCM, no dependencies;
- `scm_penalized`: penalized SCM of Abadie and L'Hour (2021), no dependencies;
- `scm_synth`: `Synth::synth`, which also optimises the predictor weights;
- any function `function(X1, X0, v = NULL, ...)` returning donor weights.

**2. Bring your own synthetic controls: `iscm_assemble()`.** When you need a
specification `iscm()` cannot express, such as `Synth::dataprep` with special
predictors and different time windows, or a penalized estimator, estimate the
synthetic control of each focal unit yourself (with the other focal units in its
donor pool) and pass the gap series and donor weights:

```r
fit <- iscm_assemble(
  gaps    = list("West Germany" = gap_wg, "Austria" = gap_au),
  weights = list("West Germany" = w_wg,   "Austria" = w_au),
  t0 = 1989, times = 1960:2003)
```

The result has the same `print`, `summary` and `plot` methods. Use
`iscm_compare_fits()` for the decision rule on such fits. For the placebo
analysis, get the data purged of the estimated effects with `iscm_purge()`,
estimate the synthetic control of each pure control on it with your own
estimator, and pass the weights to `iscm_placebo_fits()`.

## A typical workflow

```r
args <- list(data = sim$data, id = "unit", time = "time", outcome = "Y",
             treated = "TR", affected = "A1", t0 = sim$t0)

# 1. Estimate the four specifications
fit   <- do.call(iscm, args)
fit_s <- do.call(iscm, c(args, simplified = TRUE))
rest  <- do.call(restricted_scm, args)
# unrestricted SCM: fit$gaps_unrestricted

# 2. Diagnostics and choice of specification (Section 5.1)
fit$diagnostics
do.call(iscm_compare, c(args, validation = 8))

# 3. Placebo inference for the treated and the affected unit
iscm_placebo(fit, sim$data, "unit", "time", "Y", target = "TR")
iscm_placebo(fit, sim$data, "unit", "time", "Y", target = "A1")

# 4. Sensitivity to the set of affected units (Section 5.2)
sim2 <- simulate_iscm(n_affected = 2, spill = c(3, 1), seed = 5)
iscm_sensitivity(sim2$data, "unit", "time", "Y", treated = "TR",
                 affected_sets = list(minimal = "A1", extended = c("A1", "A2")),
                 t0 = sim2$t0)
```

The vignette, `vignette("iscm")`, goes through each step with output and
explanations.

## Replicating the paper

The `germany` dataset is the panel of Abadie et al. (2015). Two scripts installed
with the package reproduce the empirical application with the exact `Synth`
specification of the paper:

```r
# Tables 1 and C.1-C.2, Figure 1, decision rule, placebo analysis (Figure C.1)
source(system.file("examples", "germany_replication.R", package = "iscm"))

# Inclusive penalized SCM (Online Appendix C, Figure C.2)
source(system.file("examples", "germany_penalized.R", package = "iscm"))
```

They recover the published figures: Austria's weight in synthetic West Germany is
`0.42` and West Germany's in synthetic Austria `0.33`, so Ω has determinant
`0.86`, smallest singular value `0.62` and α ≈ `1.16`; the average effect on West
Germany is `-1,677` dollars per year (unrestricted SCM: `-1,570`; restricted SCM:
`-2,387`), and the spillover on Austria is `-254`.

## Getting help

- `package?iscm` — overview of the package and the map of commands;
- `?iscm`, `?iscm_compare`, `?iscm_placebo`, ... — every function has runnable
  examples;
- `vignette("iscm")` — the full workflow.

## Notes for developers

The `.Rd` help files and `NAMESPACE` are generated from the roxygen comments in
`R/`. After editing, regenerate them and check the package with:

```r
devtools::document()
devtools::test()
devtools::check()
```

## Reference

Di Stefano, R. and Mellace, G. *The inclusive synthetic control method.*
