test_that("select_gam_covariates writes a run log, and logfile=NULL disables it", {
  skip_on_cran()
  set.seed(7)
  n  <- 300
  df <- data.frame(lon = runif(n, -20, -5), lat = runif(n, 35, 50),
                   sst = rnorm(n), z = rnorm(n))
  df$catch <- rnbinom(n, mu = exp(1 + 0.5 * sin(df$sst)), size = 2)

  fit <- function(dir, ...) {
    select_gam_covariates(df, "catch", candidates = "sst",
                          structural = "te(lon, lat)", cv_k = 3,
                          cv_metric = "haversine", outdir = dir,
                          verbose = TRUE, ...)
  }

  # default: log written into outdir and its path returned
  d1  <- file.path(tempdir(), "log_on")
  res <- fit(d1)
  logpath <- file.path(d1, "covariate_selection_log.txt")
  expect_true(file.exists(logpath))
  expect_gt(file.info(logpath)$size, 0)          # not empty
  expect_identical(res$logfile, logpath)

  # logfile = NULL: nothing written, NULL returned
  d2   <- file.path(tempdir(), "log_off")
  res2 <- fit(d2, logfile = NULL)
  expect_false(file.exists(file.path(d2, "covariate_selection_log.txt")))
  expect_null(res2$logfile)

  # sink stack is balanced afterwards (no dangling redirection)
  expect_identical(sink.number(), 0L)
})
