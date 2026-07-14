test_that("select_gam_covariates recovers true drivers and a count family", {
  skip_on_cran()
  set.seed(42)
  n  <- 500
  df <- data.frame(
    lon    = runif(n, -20, -5),
    lat    = runif(n,  35,  50),
    sst    = rnorm(n),
    noise1 = rnorm(n)
  )
  eta      <- 1 + 0.6 * sin(df$lon / 4) + 0.8 * sin(df$sst)
  df$catch <- rnbinom(n, mu = exp(eta), size = 2)

  res <- select_gam_covariates(
    data       = df,
    response   = "catch",
    candidates = c("sst", "noise1"),
    structural = "te(lon, lat)",
    cv_k       = 3,
    cv_metric  = "haversine",
    outdir     = tempdir(),
    verbose    = FALSE
  )

  expect_true("sst" %in% res$kept)          # real driver kept
  expect_false("noise1" %in% res$kept)      # pure noise dropped
  expect_match(res$family, "nb|poisson|tw") # non-negative count family
  expect_s3_class(res$model, "gam")
})
