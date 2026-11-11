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
