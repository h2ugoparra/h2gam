# h2gam (development version)

## Bug fixes

* A family supplied through `family` with an estimated shape parameter
  (`nb`, `nb()`, `tw`, `betar`, `scat`) now has that parameter held fixed at
  its structural-only estimate during the search, as auto-suggested `nb` and
  `tw` already did; auto-suggested `betar` is now fixed too. Left free, the
  shape was re-estimated on every CV fit, which made held-out deviance
  incomparable between models: real drivers could be rejected (`nb`, `scat`)
  or noise accepted (`tw`). Selections with these families may change. The
  supplied link is kept. Families whose shape cannot be fixed (e.g. `ocat()`)
  give a warning.
* The final model is now fitted on every row complete for the kept
  covariates. It used to reuse the search data, which is restricted to rows
  complete for every screened-in candidate, so a candidate with missing
  values (up to `na_max`) silently removed those rows from the final model
  even when it was not selected. The search is unchanged; the number of rows
  it uses, and the final fit's when different, are now logged.
* `family = "auto"` no longer offers Tweedie for count responses or for
  proportions with exact 0/1. Its likelihood is a density, which cannot be
  ranked by AIC against Poisson / negative binomial probabilities or against
  a squeezed-beta fit; on low counts it beat the true Poisson or negative
  binomial family by hundreds of AIC units. Counts now rank Poisson against
  negative binomial only, and boundary proportions use squeezed beta. Auto
  runs that chose `tw` on counts will now choose `nb` or `poisson`; pass
  `family = tw` to keep Tweedie. Tweedie is still ranked for strictly
  positive continuous responses.

# h2gam 0.3.0

## Breaking changes

* A random-effect structural term such as `s(year, bs='re')` now requires its
  variable to be a factor; a numeric variable stops with an error asking you to
  convert it, and a missing column is reported as missing. Previously every
  held-out CV prediction came back `NA`, so the search silently selected
  nothing. The default `structural` includes `s(year, bs='re')`.
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
