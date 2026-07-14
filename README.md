# h2gam

Automated GAM covariate validation & selection for [mgcv](https://cran.r-project.org/package=mgcv), data-independent.

Given a response and a pool of candidate covariates, `select_gam_covariates()`:

1. **Profiles the response** and auto-suggests a family (count / proportion / positive / gaussian branches), ranked by AIC + residual diagnostics.
2. **Forward-selects** a decorrelated, prediction-relevant covariate subset, scored by random/spatial-block CV skill, with a concurvity gate and a BIC tie-break / 1-SE stopping rule. Fold construction is set by `cv_scheme`: `"spatial"` (block CV, the default) or `"stratified"` — random k-fold stratified on the response, for binomial / low-prevalence targets where spatial blocks would leave folds with too few positives.
3. Runs a **bidirectional backward re-check**, a final `select=TRUE` shrinkage fit, and adequacy diagnostics.

Every threshold is a logged, overridable argument; nothing is hardcoded to a particular dataset. Structural (spatial/temporal) terms are forced in and not screened — covariates must earn their place over and above them.

## Usage

This repo is meant to be cloned once and sourced by path from consuming projects. Sourcing loads only the functions (the example run at the bottom of the script executes only under `Rscript`, not when sourced):

```r
source("C:/Users/h2ugo/Documents/h2gam/covariate_selection.R")

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

Each project keeps its own short driver script (data loading, preprocessing, target list); the selection logic lives here.

## Dependencies

Required: `mgcv`, `cluster`.
Optional (feature-gated): `DHARMa` (residual diagnostics), `sf` (spatial CV blocking).

## Outputs

Written to `outdir`:

- `covariate_decision_table.csv` — per-candidate accept/reject decisions with scores
- `selection_diagnostics.png` — adequacy diagnostics of the final model
- `partial_effects.png` — partial effect plots of the kept smooths

The returned list carries the fitted `mgcv` model, the kept covariate set, the chosen family, the decision table, the selection path, the final formula, deviance explained, and the CV folds.
