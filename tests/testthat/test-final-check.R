# final_check$resid_moran_*: residual spatial autocorrelation of the final model.

sim_field <- function(n = 300, seed = 5) {
  set.seed(seed)
  df <- data.frame(lon = runif(n, -20, -5), lat = runif(n, 35, 50), x1 = rnorm(n))
  df$y <- 2 * sin(df$lon / 2) + 2 * cos(df$lat / 2) + rnorm(n, sd = 0.3)
  df
}

# R CMD check hides packages the DESCRIPTION does not declare. DHARMa then
# fails to simulate from a gam where mgcViz is installed but hidden, and
# every DHARMa result is NA: nothing to test in that environment.
skip_if_no_dharma_sim <- function() {
  skip_if_not_installed("DHARMa")
  skip_if(is.null(dharma_sim(gam(y ~ 1, data = data.frame(y = rnorm(20))))),
          "DHARMa cannot simulate from a gam here")
}

run_fc <- function(df, ...) {
  select_gam_covariates(df, "y", "x1", cv_k = 3, cv_scheme = "stratified",
                        outdir = tempdir(), logfile = NULL, verbose = FALSE, ...)$final_check
}

test_that("final_check flags spatial structure the model leaves in the residuals", {
  skip_on_cran()
  skip_if_no_dharma_sim()
  df <- sim_field()
  # no spatial term: the field stays in the residuals
  bare <- run_fc(df, structural = character())
  expect_gt(bare$resid_moran_i, 0.05)
  expect_lt(bare$resid_moran_p, 0.001)
  # with the spatial smooth it is absorbed
  smooth <- run_fc(df, structural = "te(lon, lat)")
  expect_lt(abs(smooth$resid_moran_i), bare$resid_moran_i / 5)
  # planar distances give the same verdict
  expect_lt(run_fc(df, structural = character(), cv_metric = "euclidean")$resid_moran_p, 0.001)
})

test_that("the Moran test aggregates rows that share a location", {
  skip_on_cran()
  skip_if_no_dharma_sim()
  df <- sim_field()
  df$lon <- round(df$lon); df$lat <- round(df$lat)       # many repeated locations
  df$y <- 2 * sin(df$lon / 2) + 2 * cos(df$lat / 2) + rnorm(nrow(df), sd = 0.3)
  fc <- run_fc(df, structural = character())
  expect_gt(fc$resid_moran_i, 0.05)
  expect_lt(fc$resid_moran_p, 0.001)
})

test_that("the Moran test is NA when it cannot be run", {
  skip_on_cran()
  df <- sim_field(n = 120)
  fc <- run_fc(df, structural = character(), coords = c("absent", "absent2"))
  expect_true(is.na(fc$resid_moran_i) && is.na(fc$resid_moran_p))
  expect_identical(dharma_moran(NULL, 1:5, 1:5), list(i = NA_real_, p = NA_real_))
})

test_that("the residual map and its data are produced, one row per location", {
  skip_on_cran()
  df <- sim_field(n = 150)
  df$lon <- round(df$lon); df$lat <- round(df$lat)       # repeated locations
  out <- file.path(tempdir(), "resid-map")
  res <- select_gam_covariates(df, "y", "x1", structural = character(), cv_k = 3,
                               cv_scheme = "stratified", outdir = out, logfile = NULL,
                               verbose = FALSE)
  rm <- res$resid_map
  expect_named(rm, c("lon", "lat", "n", "resid"))
  expect_identical(nrow(rm), nrow(unique(df[c("lon", "lat")])))
  expect_identical(sum(rm$n), nrow(df))
  # n counts the rows at each location
  one <- rm[1, ]
  expect_identical(one$n, sum(df$lon == one$lon & df$lat == one$lat))
  expect_true(all(is.finite(rm$resid)))
  expect_true(attr(rm, "type") %in% c("DHARMa scaled", "deviance"))
  if (attr(rm, "type") == "DHARMa scaled") expect_true(all(rm$resid >= 0 & rm$resid <= 1))
  expect_gt(file.size(file.path(out, "residual_map.png")), 0)
})

test_that("without coordinates there is no residual map", {
  skip_on_cran()
  out <- file.path(tempdir(), "resid-map-none")
  res <- select_gam_covariates(sim_field(n = 120), "y", "x1", structural = character(),
                               cv_k = 3, cv_scheme = "stratified", coords = c("absent", "absent2"),
                               outdir = out, logfile = NULL, verbose = FALSE)
  expect_null(res$resid_map)
  expect_false(file.exists(file.path(out, "residual_map.png")))
})
