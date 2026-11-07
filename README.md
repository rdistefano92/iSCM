# iscm — The Inclusive Synthetic Control Method

An R implementation of the **inclusive synthetic control method (iSCM)** of
Di Stefano and Mellace. iSCM lets units that are *potentially affected* by an
intervention (other treated units, or controls exposed to spillovers) remain in
the donor pool, and removes the resulting contamination exactly, on top of any
synthetic-control-type estimator.

## Installation

```r
# install.packages("remotes")
remotes::install_github("rdistefano92/iSCM", build_vignettes = TRUE)
```

The only required dependency is `ggplot2`, used for the plots. To use the
classic SCM of Abadie, Diamond and Hainmueller (2010, 2015), as in the paper's
application, also install the `Synth` package (Abadie, Diamond and Hainmueller,
2011). The penalized-SCM replication script additionally needs `pensynth` and
`clarabel`.

## What you need

iSCM sorts the units of the panel into three groups:

- the **main treated unit**, whose effect you want to estimate;
- the **potentially affected units**, the set *M* in the paper: units that may
  also have been affected by the intervention, because they were treated as well
  or exposed to spillovers. Unlike in standard SCM, they stay in the donor pool;
- the **pure controls**: all remaining units, assumed unaffected.

You then need a balanced panel in long format (one row per unit and period) and
these arguments:

| Argument | What it is |
|---|---|
| `data` | the panel, a `data.frame` |
| `id`, `time`, `outcome` | names of the columns with the unit identifier, the period and the outcome |
| `treated` | identifier of the main treated unit |
| `affected` | identifiers of the potentially affected units; every unit not listed here or in `treated` is a pure control |
| `t0` | last pre-intervention period |
| `predictors` | optional: names of the columns used for matching, averaged over the pre-intervention period. If omitted, units are matched on the pre-intervention path of the outcome |

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
| plot observed vs synthetic outcomes, or the weights | `plot(fit, type = "trends")`, `plot(fit, type = "weights")` |
| add the true effect to the plot (simulated data) | `plot(fit, truth = sim$tau)` |
| check the **conditioning** of the correction (det, σ<sub>min</sub>, α) | `fit$diagnostics`, `omega_diagnostics(Omega)` |
| see the implied weights on pure controls (Lemma 3.1) | `fit$implied_weights` |
| **choose** among restricted / full / simplified (Section 5.1) | `iscm_compare(...)` |
| same, for fits estimated outside the package | `iscm_compare_fits()`, `iscm_decision()` |
| run the **placebo** analysis ([Section 6](https://arxiv.org/abs/2403.17624)) | `iscm_placebo(fit, data, id, time, outcome, target = ...)` |
| same, for fits estimated outside the package | `iscm_purge()`, then `iscm_placebo_fits(fit, ..., placebo_weights)` |
| check **sensitivity** to the set of affected units (Section 5.2) | `iscm_sensitivity(..., affected_sets = list(...))` |
| use **another SC-type estimator** | `iscm(..., estimator = scm_penalized)`, `scm_synth`, or your own function (`?backends`) |
| correct synthetic controls **estimated elsewhere** | `iscm_assemble(gaps, weights, t0)` |
| apply only the linear algebra (Theorem 2.1) | `build_omega()`, `solve_iscm()` |
| simulate data with known effects | `simulate_iscm()`, `simulate_iscm_covariates()` |

All data-level functions (`iscm`, `restricted_scm`, `iscm_compare`,
`iscm_sensitivity`) take the same first arguments, so a typical analysis is the
same call repeated with a different function name. All plots are `ggplot2`
objects, so they can be customised, e.g. `plot(fit) + ggplot2::labs(y = "GDP")`.

## When does iSCM beat the restricted SCM?

The iSCM counterfactual is itself a weighted average of pure controls
(Lemma 3.1), so iSCM gains when the restricted optimizer, guided by its own
criterion, picks a worse combination of pure controls than the one iSCM implies
(Section 3.1). `simulate_iscm_covariates()` builds such a case, in the spirit of
the German application:

- the outcome depends on two **observed covariates**, `x1` and `x2`, and on an
  **unobserved characteristic** *u* (for example, institutions or economic
  structure, which the covariates do not capture);
- the main treated unit `TR` and the affected unit `A1` are close to each other
  in the covariates and share the same *u*, like West Germany and Austria;
- two pure controls close to `TR` in the covariates, `C3` and `C4`, have a
  different *u*; the pure controls with the same *u*, `C1` and `C2`, lie on
  either side of the pair, farther from `TR` than `A1`.

All synthetic controls are matched on `x1` and `x2` only. With `A1` in the donor
pool, the synthetic `TR` is built mostly from `A1` and therefore has the right
*u*; the cross-weights (about 0.5 and 0.3) are similar to those of the German
application. Without `A1`, the restricted optimizer can match the covariates
equally well with `C1` and `C2` (right *u*) or with `C1`, `C3` and `C4` (wrong
*u*); which one it picks depends on the noise, so the restricted SCM is unstable
and biased on average. iSCM keeps `A1` and removes the contamination from its
spillover.

```r
sim  <- simulate_iscm_covariates(seed = 4)     # true effects: 2 on TR, 3 on A1
args <- list(sim$data, "unit", "time", "Y", "TR", "A1", sim$t0,
             predictors = c("x1", "x2"), estimator = scm_penalized)
fit  <- do.call(iscm, args)
rest <- do.call(restricted_scm, args)
```

Over 100 replications, with a true effect on `TR` of 2:

| | mean estimate | mean absolute error |
|---|---|---|
| unrestricted SCM | 0.51 | 1.49 |
| restricted SCM | 2.93 | 0.93 |
| full iSCM | 1.97 | 0.05 |
| simplified iSCM | 1.94 | 0.07 |
| specification chosen by `iscm_compare()` | 1.99 | 0.06 |

The decision rule of Section 5.1 recommends the full iSCM in 84 replications and
the restricted SCM in the other 16, in which the restricted fit is comparable;
the amplification factor α stays between 1.00 and 1.23. The script
`system.file("examples", "simulation_covariates.R", package = "iscm")` reproduces
these figures.

## Two ways to use the package

**1. Let the package estimate the synthetic controls: `iscm()`.** Each unit is
matched on the pre-intervention outcome path (default) or on the columns given in
`predictors`, averaged over the pre-intervention period. Choose the estimator with
`estimator =`:

- `scm_base` (default): classic simplex-constrained SCM, no extra dependencies;
- `scm_penalized`: penalized SCM of Abadie and L'Hour (2021), no extra dependencies;
- `scm_synth`: `Synth::synth`, which also optimises the predictor weights;
- any function `function(X1, X0, v = NULL, ...)` returning donor weights.

**2. Bring your own synthetic controls: `iscm_assemble()`.** When you need a
specification `iscm()` cannot express, such as `Synth::dataprep` with special
predictors and different time windows, or a penalized estimator, estimate the
synthetic control of the main treated unit and of each affected unit yourself,
with the other units (affected ones included) in the donor pool, and pass to
`iscm_assemble()`, for each of them, the gap series (observed minus synthetic
outcome, over all periods) and the named vector of donor weights.

The example below estimates the two synthetic controls by hand on simulated data,
so that it runs as it is:

```r
library(iscm)
sim <- simulate_iscm(n_affected = 1, n_pure = 6, tau = 2, spill = 3, seed = 7)
Y   <- with(sim$data, tapply(Y, list(time, unit), sum))   # periods x units
pre <- seq_len(sim$t0)

# your own estimator: here scm_base() on the pre-intervention outcome path
my_sc <- function(unit) {
  donors <- setdiff(colnames(Y), unit)
  w <- scm_base(Y[pre, unit], Y[pre, donors], standardize = FALSE)
  list(w = w, gap = as.numeric(Y[, unit] - Y[, donors] %*% w))
}
sc_tr <- my_sc("TR")
sc_a1 <- my_sc("A1")

fit <- iscm_assemble(
  gaps    = list(TR = sc_tr$gap, A1 = sc_a1$gap),
  weights = list(TR = sc_tr$w,   A1 = sc_a1$w),
  t0 = sim$t0, pure = sim$pure)
fit
```

The replication script of the paper (see [Replicating the paper](#replicating-the-paper))
does the same with the exact `Synth` specification for West Germany and Austria.

The result has the same `print`, `summary` and `plot` methods. Use
`iscm_compare_fits()` for the decision rule on such fits. For the placebo
analysis, get the data purged of the estimated effects with `iscm_purge()`,
estimate the synthetic control of each pure control on it with your own
estimator, and pass the weights to `iscm_placebo_fits()`.

## A typical workflow

```r
library(iscm)
sim  <- simulate_iscm(n_affected = 1, n_pure = 6, tau = 2, spill = 3, seed = 7)
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

The `germany` dataset is the panel of Abadie, Diamond and Hainmueller (2015).
Two scripts installed with the package reproduce the empirical application with
the exact `Synth` specification of the paper. Each one runs the whole analysis,
prints the main results and saves the tables (as `.csv`) and figures (as `.pdf`
and `.png`) in a folder of your choice:

```r
options(iscm.output_dir = "results")   # folder for tables and figures

# West Germany and Austria, about 4 minutes. Needs Synth.
# Saves Tables 1, 2 and 3, Figure 1 and the placebo ratios of Figure B.1,
# and prints Omega, the effects, the decision rule and the placebo p-values.
source(system.file("examples", "germany_replication.R", package = "iscm"))

# Inclusive penalized SCM (Figure 2). Needs Synth, pensynth and clarabel.
source(system.file("examples", "germany_penalized.R", package = "iscm"))
```

Without `options(iscm.output_dir = ...)`, the files go to a new subfolder of the
working directory (`iscm_germany_replication` and `iscm_germany_penalized`). To
read the code instead of running it, open the scripts with
`file.edit(system.file("examples", "germany_replication.R", package = "iscm"))`.

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

## References

Abadie, A., Diamond, A. and Hainmueller, J. (2010). Synthetic control methods for
comparative case studies: Estimating the effect of California's tobacco control
program. *Journal of the American Statistical Association* 105(490), 493–505.

Abadie, A., Diamond, A. and Hainmueller, J. (2011). Synth: An R package for
synthetic control methods in comparative case studies. *Journal of Statistical
Software* 42(13), 1–17.

Abadie, A., Diamond, A. and Hainmueller, J. (2015). Comparative politics and the
synthetic control method. *American Journal of Political Science* 59(2), 495–510.

Abadie, A. and L'Hour, J. (2021). A penalized synthetic control estimator for
disaggregated data. *Journal of the American Statistical Association* 116(536),
1817–1834.

Di Stefano, R. and Mellace, G. *The inclusive synthetic control method.*
[arXiv:2403.17624](https://arxiv.org/abs/2403.17624).
