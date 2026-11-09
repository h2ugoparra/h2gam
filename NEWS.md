# h2gam 0.3.0

## Breaking changes

* A random-effect structural term such as `s(year, bs='re')` now requires its
  variable to be a factor; a numeric variable stops with an error asking you to
  convert it. Previously every held-out CV prediction came back `NA`, so the
  search silently selected nothing. The default `structural` includes
  `s(year, bs='re')`.
* Spatial CV with `cv_metric = "euclidean"` and the default linkage now uses
  `"ward.D2"` (scikit-learn's Ward) instead of `hclust`'s legacy `"ward"`
  (`"ward.D"`), so fold assignments, and possibly selections, differ from
  0.2.0. Haversine runs are unchanged.

## Bug fixes

* Missing offset values no longer void every CV fold. Rows with an `NA`
  offset are dropped, and non-positive offsets are rejected.
* Forward selection now errors when the structural-only baseline cannot be
  scored on at least 2 folds, and warns about candidates scored on fewer than
  half the folds, instead of reporting "no CV improvement".
* `na_frac` in the decision table now reflects the data as supplied (it read
  0 for every screened-in candidate).
* A parametric structural term (e.g. `"lat"`) no longer crashes the final
  partial-effects plot.
* Fixed crashes with: no numeric candidates or a rare response under spatial
  CV; `Inf` covariate values; a final model with no smooth terms.
* `make_spatial_folds()` and `make_stratified_folds()` restore the caller's
  random-number state; folds remain reproducible from `seed`.
* A family object passed as `family` (e.g. `nb()`) has its theta reset before
  every fit, so fits no longer start from the previous fit's estimate.
* `verbose = FALSE` is now fully quiet, including `gam.check()` output and
  DHARMa's mgcViz notice.

## Documentation

* Removed claims of a "BIC tie-break" in forward selection and a
  "bidirectional" backward re-check; neither existed.
* Documented the O(n^2) memory use of `make_spatial_folds()`.
