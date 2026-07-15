# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

`h2gam` is an R package: automated GAM covariate validation and selection for mgcv. The single entry point is `select_gam_covariates()`; almost all logic lives in `R/covariate_selection.R`.

## Commands

Run from the package root with Rscript:

```sh
# All tests
Rscript -e "devtools::test()"

# A single test file (filter matches tests/testthat/test-<name>.R)
Rscript -e "devtools::test(filter = 'family')"

# Regenerate NAMESPACE and man/ after editing roxygen comments
Rscript -e "devtools::document()"

# Full check — CI runs R CMD check --as-cran and fails on WARNINGs, so keep this clean
Rscript -e "devtools::check(args = c('--no-manual', '--as-cran'))"
```

Documentation is roxygen2 (markdown enabled): edit the `#'` blocks in `R/`, never `man/*.Rd` or `NAMESPACE` directly, then run `devtools::document()`.

## CI gates (on PRs to dev/main)

- Conventional commit messages (`feat:`, `fix:`, `docs:`, `chore:`, `perf:`, `refactor:`, ...).
- Branch names must be prefixed `feature/`, `fix/`, `docs/`, `chore/`, `perf/`, `refactor/` (PRs to dev).
- `R CMD check --as-cran` with `error-on: warning`.

## Architecture

`R/covariate_selection.R` is a pipeline orchestrated by `select_gam_covariates()`:

1. **Family suggestion** (`profile_response` → `family_candidates` → `suggest_family`): profiles the response, fits each candidate family on the structural-only model, ranks by AIC alone. Dispersion / DHARMa columns are reported, never scored.
2. **CV folds**: `make_spatial_folds()` (SPCV — a faithful R port of the Python h2ml `SPCVSplitter`: AHC blocks then an HBGF cluster ensemble with spectral clustering) or `make_stratified_folds()` for rare-event responses. Both are exported.
3. **Forward selection** (`forward_select`): scored by paired per-fold held-out deviance with a 1-SE rule and a concurvity gate against already-accepted smooths (structural terms excluded from the gate).
4. **Backward re-check** (BIC) and a **final `select = TRUE` shrinkage fit**, plus a decision table CSV and diagnostic PNGs written to `outdir`.

Non-obvious invariants to preserve when editing:

- **Family generators are zero-arg closures, called fresh per fit** — `nb()`/`tw()` carry mutable theta/p state in their environment, so a shared family object leaks state across fits.
- **Search fits vs final fit are deliberately different**: the search uses `bam(discrete = TRUE, method = "fREML")` with the family shape parameter *fixed* (`make_select_gen`) for speed; the final fit is a free `gam(method = "REML", select = TRUE)`.
- **Random-effect structural terms (`bs = 're'`) are excluded from held-out CV prediction** (via `predict(..., exclude = ...)` plus a factor-level remap), matching deployment behavior.
- **The SV squeeze for betar is applied once, up front**, so the CV search and the final fit see the same (0,1) response.
- **Structural terms are forced in and never screened** — candidates must improve CV skill over and above them.
- mgcv is in `Depends` (not just Imports) on purpose, so `s()`/`te()` and family constructors are visible to users writing structural terms; the package's own code still binds via explicit `@importFrom` (see `R/h2gam-package.R`).
- `sf` and `DHARMa` are `Suggests` and feature-gated with `requireNamespace()` — keep new optional dependencies behind the same pattern.

Tests are testthat edition 3 in `tests/testthat/`; the end-to-end demo lives in `inst/examples/synthetic_demo.R`.
