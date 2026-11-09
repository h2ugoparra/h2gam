# Regression tests for inputs that used to void every CV fold silently, so the
# search reported "no CV improvement" instead of failing or selecting.

sim_counts <- function(n = 300, seed = 11) {
  set.seed(seed)
  df <- data.frame(lon = runif(n, -20, -5), lat = runif(n, 35, 50),
                   sst = rnorm(n), noise = rnorm(n),
                   year = sample(2000:2009, n, TRUE), eff = runif(n, 50, 150))
  df$catch <- rnbinom(n, mu = exp(1 + 0.8 * sin(2 * df$sst)) * df$eff / 100, size = 3)
  df
}

run <- function(df, ...) {
  select_gam_covariates(df, "catch", c("sst", "noise"), cv_k = 3,
                        cv_scheme = "stratified", outdir = tempdir(),
                        logfile = NULL, verbose = FALSE, ...)
}

test_that("NA offset rows are dropped instead of voiding the CV folds", {
  skip_on_cran()
  df <- sim_counts()
  df$eff[1:20] <- NA
  res <- run(df, structural = "te(lon, lat)", offset = "eff")
  expect_true("sst" %in% res$kept)
  expect_equal(nrow(res$model$model), nrow(df) - 20)
})

test_that("non-positive offset values are rejected", {
  df <- sim_counts()
  df$eff[1] <- 0
  expect_error(run(df, structural = "te(lon, lat)", offset = "eff"), "must be > 0")
})

test_that("a numeric random-effect variable is rejected up front", {
  df <- sim_counts()
  expect_error(run(df, structural = c("s(year, bs='re')", "te(lon, lat)")),
               "needs 'year' to be a factor")
})

test_that("a factor random-effect variable still selects the real driver", {
  skip_on_cran()
  df <- sim_counts()
  df$year <- factor(df$year)
  res <- run(df, structural = c("s(year, bs='re')", "te(lon, lat)"))
  expect_true("sst" %in% res$kept)
})

test_that("decision-table na_frac reflects the data as supplied", {
  skip_on_cran()
  df <- sim_counts()
  df$sst[1:30] <- NA
  res <- run(df, structural = "te(lon, lat)")
  expect_equal(res$decision$na_frac[res$decision$covariate == "sst"], 30 / nrow(df))
})

test_that("forward selection fails loudly when the baseline cannot be scored", {
  df <- sim_counts(n = 60)
  folds <- rep_len(1:3, nrow(df))
  expect_error(
    forward_select(df, "catch", c("sst"), "te(lon, lat)",
                   family_gen = function() stop("boom"), folds = folds,
                   knots = NULL, verbose = FALSE),
    "scored on only 0 of 3 CV folds")
})

test_that("a missing random-effect variable is reported as missing", {
  df <- sim_counts()
  df$year <- NULL
  expect_error(run(df, structural = c("s(year, bs='re')", "te(lon, lat)")),
               "column 'year' not found")
})
