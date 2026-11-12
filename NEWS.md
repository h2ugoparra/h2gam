# h2gam (development version)

## Breaking changes

* `make_spatial_folds()` now defaults to `linkage = "ward.D2"`, the default
  `select_gam_covariates()` already passed, instead of `"average"`. Called
  directly with `metric = "euclidean"` and no `linkage`, it returns different
  folds than before; pass `linkage = "average"` to keep the old ones. The
  default `"haversine"` metric is unaffected (Ward is downgraded to average
  there), and so is `select_gam_covariates()`.
* `outdir` now defaults to `tempdir()` instead of the working directory, so
  a call that does not set it no longer writes the decision table, PNGs and
  run log into wherever R happens to be running (CRAN policy). Those files
  are deleted when the session ends: pass `outdir` to keep them. Calls that
  already set `outdir` are unaffected.

## New features

* Verbose output, and so the run log, now opens with a header giving the
  h2gam, mgcv and R versions, the time, and every setting of the run
  (response, candidates, structural terms, family, thresholds, CV setup). A
  log now shows which code and thresholds produced it.
* The decision table gains `cv_gain` and `cv_se` columns (placed after
  `concurvity`): the cross-validated gain, and its standard error, behind
  each candidate's decision. Candidates that never entered the model are no
  longer all labelled `"no CV improvement"`: that reason now means a gain of
  zero or less, `"CV gain below 1 SE"` marks a positive gain that failed the
  one-standard-error rule, and `"CV scoring failed"` marks candidates whose
  fits or held-out predictions failed.
* The run log now also records warnings, messages and the error that stopped
  a failed run. It used to capture printed output only, so a log could look
  clean while the console showed warnings. Each distinct text is written
  once, where it first occurred.
* `structural = character()` is now accepted and gives an intercept-only
  baseline; it used to fail with a formula parse error.

## Input validation

* `select_gam_covariates()` now stops with a clear message when `data` is
  not a data.frame, when the response or (for spatial CV) the `coords`
  columns are missing, or when the response is listed among the candidates.
  Missing `coords` used to surface as "must have n >= 2 objects to cluster"
  after the family fits had already run. Duplicate candidates are dropped
  with a warning.
* `make_spatial_folds()` and `make_stratified_folds()` stop on missing
  columns and on missing or non-finite values. `make_stratified_folds()`
  used to return fold id 0 for rows with a missing response.
* A failure in the partial-effects plot no longer discards the finished
  selection; the message is written into `partial_effects.png` instead.

## Performance

* A candidate rejected by the concurvity gate is no longer re-scored at every
  later step of the forward search (one model fit per CV fold each time). Its
  worst-case concurvity against a growing set of accepted smooths cannot
  fall, so the outcome was always the same rejection. Selections are
  unchanged; the selection path lists each such candidate once.

## Bug fixes

* Forward selection now warns when the reference model (structural terms
  plus the covariates accepted so far) could not be scored on every CV fold.
  Such a fold is lost to every candidate compared against that model, so
  gains and standard errors silently rested on fewer folds than `cv_k`;
  only a baseline scorable on fewer than 2 folds was reported.
* The concurvity gate now fails closed. A candidate whose concurvity could
  not be computed (its full-data fit or `concurvity()` failed) used to be
  accepted, and logged with concurvity 0.00 in the second case. It is now
  rejected with a warning naming it, and its decision-table reason is
  `"concurvity check failed"`.
* Random-effect structural terms are now recognised in every form, not only
  `s(var, bs='re')`: random slopes (`s(vessel, sst, bs='re')`), `by =`
  variables, factor-smooth interactions (`bs='fs'`) and tensor terms with an
  `'re'` margin are all excluded from held-out CV prediction. Previously
  these stayed in the prediction (or were excluded under the wrong label),
  and a random slope written with the numeric variable first was rejected as
  "needs to be a factor". Selections with such terms may change. Only
  held-out levels unseen in training are remapped now, so a grouping factor
  that also appears in another term keeps its real levels there.
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
