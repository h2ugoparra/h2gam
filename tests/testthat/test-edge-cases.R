# Inputs at the edges that used to crash with cryptic errors.

test_that("make_spatial_folds works with no numeric candidates", {
  set.seed(21)
  df <- data.frame(lon = runif(150, -20, -5), lat = runif(150, 35, 50),
                   y = rnorm(150), f = sample(letters[1:3], 150, TRUE))
  for (cands in list(character(), "f")) {
    folds <- make_spatial_folds(df, "y", candidates = cands, k = 3, verbose = FALSE)
    expect_setequal(unique(folds), 1:3)
  }
})

test_that("make_spatial_folds works when few blocks have distinct responses", {
  set.seed(22)
  n  <- 200
  df <- data.frame(lon = runif(n, -20, -5), lat = runif(n, 35, 50), x1 = rnorm(n),
                   y = 0)
  df$y[1] <- 1   # rare event: almost every block has mean response exactly 0
  folds <- make_spatial_folds(df, "y", candidates = "x1", k = 4, verbose = FALSE)
  expect_length(folds, n)
  expect_setequal(unique(folds), 1:4)
})

test_that("non-finite covariate values are dropped like NAs", {
  skip_on_cran()
  set.seed(23)
  n  <- 200
  df <- data.frame(lon = runif(n, -20, -5), lat = runif(n, 35, 50),
                   sst = rnorm(n), noise = rnorm(n))
  df$catch <- rnbinom(n, mu = exp(1 + 0.8 * sin(2 * df$sst)), size = 3)
  df$noise[1:5] <- Inf
  res <- select_gam_covariates(df, "catch", c("sst", "noise"), structural = "te(lon, lat)",
                               cv_k = 3, cv_scheme = "stratified", outdir = tempdir(),
                               logfile = NULL, verbose = FALSE)
  expect_length(res$folds, n - 5)                  # dropped from the search
  # and from the final fit only if the covariate carrying them is kept
  expect_equal(nrow(res$model$model), if ("noise" %in% res$kept) n - 5 else n)
})

test_that("a sparse unselected candidate does not cost the final fit its rows", {
  skip_on_cran()
  set.seed(11)
  n  <- 300
  df <- data.frame(lon = runif(n, -20, -5), lat = runif(n, 35, 50),
                   sst = rnorm(n), noise = rnorm(n))
  df$catch <- rnbinom(n, mu = exp(1 + 0.8 * sin(2 * df$sst)), size = 3)
  df$noise[1:80] <- NA
  df$sst[81:90]  <- NA
  res <- select_gam_covariates(df, "catch", c("sst", "noise"), structural = "te(lon, lat)",
                               cv_k = 3, cv_scheme = "stratified", outdir = tempdir(),
                               logfile = NULL, verbose = FALSE)
  expect_identical(res$kept, "sst")
  expect_length(res$folds, n - 90)                 # search: complete for both candidates
  expect_equal(nrow(res$model$model), n - 10)      # final fit: complete for sst only
})

test_that("a final model with no smooth terms is handled", {
  skip_on_cran()
  set.seed(24)
  n  <- 150
  df <- data.frame(lon = runif(n, -20, -5), lat = runif(n, 35, 50), noise = rnorm(n))
  df$y <- 0.2 * df$lat + rnorm(n)
  res <- select_gam_covariates(df, "y", "noise", structural = "lat",
                               cv_k = 3, cv_scheme = "stratified", outdir = tempdir(),
                               logfile = NULL, verbose = FALSE)
  expect_length(res$kept, 0)
  expect_true(all(is.na(res$decision$final_edf)))
  expect_identical(res$cv$gain, 0)                 # nothing kept, nothing gained
  expect_identical(res$cv$structural, res$cv$selected)
})

test_that("a parametric structural term doesn't break the partial-effects plot", {
  skip_on_cran()
  set.seed(25)
  n  <- 200
  df <- data.frame(lon = runif(n, -20, -5), lat = runif(n, 35, 50), sst = rnorm(n))
  df$y <- 0.2 * df$lat + sin(2 * df$sst) + rnorm(n, sd = 0.3)
  dir <- file.path(tempdir(), "param_plot")
  res <- select_gam_covariates(df, "y", "sst", structural = "lat",
                               cv_k = 3, cv_scheme = "stratified", outdir = dir,
                               logfile = NULL, verbose = FALSE)
  expect_true("sst" %in% res$kept)
  expect_true(file.exists(file.path(dir, "partial_effects.png")))
})

# Input validation: clear errors up front instead of cryptic ones mid-run.

sim_valid <- function(n = 200, seed = 31) {
  set.seed(seed)
  df <- data.frame(lon = runif(n, -20, -5), lat = runif(n, 35, 50),
                   sst = rnorm(n), noise = rnorm(n))
  df$y <- sin(2 * df$sst) + rnorm(n, sd = 0.3)
  df
}
run_valid <- function(df, candidates = c("sst", "noise"), ...) {
  select_gam_covariates(df, "y", candidates, cv_k = 3, outdir = tempdir(),
                        logfile = NULL, verbose = FALSE, ...)
}

test_that("an empty structural gives an intercept-only baseline", {
  skip_on_cran()
  res <- run_valid(sim_valid(), structural = character(), cv_scheme = "stratified")
  expect_true("sst" %in% res$kept)

  # nothing kept: the final model is intercept-only and has no terms to plot
  dir <- file.path(tempdir(), "intercept_only")
  res <- select_gam_covariates(sim_valid(), "y", "noise", structural = character(),
                               cv_k = 3, cv_scheme = "stratified", outdir = dir,
                               logfile = NULL, verbose = FALSE)
  expect_length(res$kept, 0)
  expect_true(file.exists(file.path(dir, "partial_effects.png")))
})

test_that("bad inputs are rejected with a clear message", {
  df <- sim_valid()
  expect_error(run_valid(as.matrix(df), structural = "te(lon, lat)"), "must be a data.frame")
  expect_error(run_valid(df, structural = "te(lon, lat)", coords = c("x", "lat")),
               "coords column(s) not found in data: x", fixed = TRUE)
  expect_error(run_valid(df, structural = "te(lon, lat)", coords = "lon"),
               "exactly two columns")
  expect_error(run_valid(df, candidates = c("sst", "y"), structural = "te(lon, lat)"),
               "also listed in `candidates`")
  expect_error(select_gam_covariates(df, "nope", "sst", outdir = tempdir()),
               "response column(s) not found in data: nope", fixed = TRUE)
})

test_that("duplicate candidates are dropped with a warning", {
  skip_on_cran()
  expect_warning(
    res <- run_valid(sim_valid(), candidates = c("sst", "sst", "noise"),
                     structural = "te(lon, lat)", cv_scheme = "stratified"),
    "duplicate candidates dropped: sst")
  expect_identical(res$decision$covariate, c("sst", "noise"))
})

test_that("fold builders reject missing values instead of returning bad folds", {
  df <- sim_valid()
  df$y[1] <- NA
  expect_error(make_stratified_folds(df, "y", k = 3, verbose = FALSE), "missing values")
  expect_error(make_spatial_folds(df, "y", "sst", k = 3, verbose = FALSE),
               "coordinates or response")
  df <- sim_valid(); df$sst[2] <- NA
  expect_error(make_spatial_folds(df, "y", "sst", k = 3, verbose = FALSE),
               "numeric candidates")
  expect_error(make_spatial_folds(sim_valid(), "y", "sst", coords = c("x", "lat")),
               "coords column")
})
