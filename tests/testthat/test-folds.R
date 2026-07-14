test_that("make_stratified_folds returns k balanced folds covering every row", {
  set.seed(1)
  df <- data.frame(y = rbinom(200, 1, 0.15))
  folds <- make_stratified_folds(df, "y", k = 5, verbose = FALSE)

  expect_length(folds, nrow(df))
  expect_setequal(unique(folds), 1:5)
  # cyclic assignment is balanced within each stratum, so across s strata the
  # fold sizes differ by at most s (here binary => 2 strata).
  expect_lte(diff(range(tabulate(folds, 5))), 2)
})

test_that("make_stratified_folds spreads rare positives across folds", {
  set.seed(2)
  df <- data.frame(y = rbinom(300, 1, 0.1))
  folds <- make_stratified_folds(df, "y", k = 5, verbose = FALSE)

  pos_per_fold <- tapply(df$y, folds, sum)
  # stratifying on the label => no fold is left without positives
  expect_true(all(pos_per_fold > 0))
})

test_that("make_spatial_folds (haversine) partitions rows into k spatial folds", {
  set.seed(3)
  n  <- 200
  df <- data.frame(
    lon = runif(n, -20, -5),
    lat = runif(n,  35,  50),
    x1  = rnorm(n),
    y   = rnorm(n)
  )
  folds <- make_spatial_folds(df, "y", candidates = "x1", k = 4,
                              metric = "haversine", verbose = FALSE)

  expect_length(folds, n)
  expect_setequal(unique(folds), 1:4)
})

test_that("make_spatial_folds (euclidean, crs=NULL) uses planar coords as-is", {
  set.seed(4)
  n  <- 200
  # already-projected coordinates (e.g. metres); no CRS / sf needed
  df <- data.frame(
    x  = runif(n, 0, 5e5),
    y  = runif(n, 0, 5e5),
    x1 = rnorm(n),
    z  = rnorm(n)
  )
  folds <- make_spatial_folds(df, "z", candidates = "x1", k = 4,
                              coords = c("x", "y"), metric = "euclidean",
                              crs = NULL, verbose = FALSE)

  expect_length(folds, n)
  expect_setequal(unique(folds), 1:4)
})
