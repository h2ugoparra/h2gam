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

test_that("the run log records warnings and the error that stopped a run", {
  skip_on_cran()
  set.seed(7)
  n  <- 200
  df <- data.frame(lon = runif(n, -20, -5), lat = runif(n, 35, 50), sst = rnorm(n))
  df$y <- sin(2 * df$sst) + rnorm(n, sd = 0.3)
  fit <- function(dir, ...) {
    utils::capture.output(
      res <- select_gam_covariates(df, "y", structural = "te(lon, lat)", cv_k = 3,
                                   cv_scheme = "stratified", outdir = dir, ...))
    res
  }

  # a warning still reaches the caller, and is in the log once
  d1 <- file.path(tempdir(), "log_warn")
  expect_warning(res <- fit(d1, candidates = c("sst", "sst")), "duplicate candidates")
  log <- readLines(res$logfile)
  expect_length(grep("^Warning: duplicate candidates dropped: sst", log), 1)

  # an error still stops the run, is in the log, and the console is restored
  d2 <- file.path(tempdir(), "log_error")
  expect_error(fit(d2, candidates = "sst", offset = "nope"), "offset column 'nope'")
  log <- readLines(file.path(d2, "covariate_selection_log.txt"))
  expect_length(grep("^Error: offset column 'nope' not found", log), 1)
  expect_identical(sink.number(), 0L)
})

test_that("outdir defaults to the temporary directory, not the working directory", {
  expect_identical(eval(formals(select_gam_covariates)$outdir), tempdir())
})

test_that("verbose = FALSE prints nothing", {
  skip_on_cran()
  set.seed(8)
  n  <- 200
  df <- data.frame(lon = runif(n, -20, -5), lat = runif(n, 35, 50), sst = rnorm(n))
  df$catch <- rnbinom(n, mu = exp(1 + 0.5 * sin(df$sst)), size = 2)
  expect_silent(select_gam_covariates(df, "catch", candidates = "sst",
                                      structural = "te(lon, lat)", cv_k = 3,
                                      cv_scheme = "stratified", outdir = tempdir(),
                                      verbose = FALSE))
})
