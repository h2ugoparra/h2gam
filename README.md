# h2gam

Automated GAM covariate validation & selection for [mgcv](https://cran.r-project.org/package=mgcv), data-independent.

Given a response and a pool of candidate covariates, `select_gam_covariates()`:

1. **Profiles the response** and auto-suggests a family (count / proportion / positive / gaussian branches), ranked by AIC. Only families whose AIC is comparable are ranked against each other, so Tweedie is not auto-suggested for counts or for proportions with exact 0/1; pass `family = tw` to use it. Dispersion and a DHARMa residual-uniformity p-value are reported alongside for you to inspect, but do not affect the choice — and since both are measured on the structural-only model, a low `resid_ks_p` reflects not-yet-modelled covariates and structure as much as the family. Re-check it on the final fit before acting on it.
2. **Forward-selects** a decorrelated, prediction-relevant covariate subset, scored by random/spatial-block CV skill, with a concurvity gate and a 1-SE stopping rule. Fold construction is set by `cv_scheme`: `"spatial"` (block CV, the default) or `"stratified"` — random k-fold stratified on the response, for binomial / low-prevalence targets where spatial blocks would leave folds with too few positives.
3. Runs a **backward elimination re-check** (BIC), a final `select=TRUE` shrinkage fit, and adequacy diagnostics.

Every threshold is a logged, overridable argument; nothing is hardcoded to a particular dataset. Structural (spatial/temporal) terms are forced in and not screened — covariates must earn their place over and above them.

## Installation

```r
# install.packages("remotes")
remotes::install_github("h2ugoparra/h2gam")
```

Pin a release for reproducibility, e.g. `remotes::install_github("h2ugoparra/h2gam@v0.3.0")`.

## Usage

```r
library(h2gam)

df$year <- factor(df$year)   # random-effect terms (bs = 're') need a factor

res <- select_gam_covariates(
  data       = df,
  response   = "my_target",
  candidates = c("sst", "chl", "mld", ...),
  structural = c("s(year, bs='re')", "te(lon, lat)"),
  offset     = "nhooks",
  cv_k       = 5,
  outdir     = "output/my_target"
)

saveRDS(res$model, "output/my_target/gam_model.rds")
```

The variable in a random-effect structural term such as `s(year, bs='re')` must be a factor. On a numeric variable mgcv fits a penalised slope instead of random intercepts, so `select_gam_covariates()` stops with an error asking you to convert it. The default `structural` includes `s(year, bs='re')`, so this applies whenever you rely on the defaults.

Each project keeps its own short driver script (data loading, preprocessing, target list); the selection logic lives here. See `?select_gam_covariates` for the full argument reference and a runnable synthetic example, or `inst/examples/synthetic_demo.R` for a self-contained end-to-end demo.

## Dependencies

Required: `mgcv` (attached with the package).
Optional (feature-gated via `Suggests`): `DHARMa` (residual diagnostics in family suggestion), `sf` (only when reprojecting lon/lat for the `euclidean` spatial-CV metric, i.e. when you pass `cv_crs`). The default `haversine` metric, and `euclidean` on already-projected coordinates, need neither.

## Outputs

Written to `outdir`:

- `covariate_decision_table.csv` — per-candidate accept/reject decisions with scores
- `selection_diagnostics.png` — adequacy diagnostics of the final model
- `partial_effects.png` — partial effect plots of the kept smooths

The returned list carries the fitted `mgcv` model, the kept covariate set, the chosen family, the decision table, the selection path, the final formula, deviance explained, and the CV folds.
