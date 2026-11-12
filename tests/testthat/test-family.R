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
  expect_match(res$family, "nb|poisson")    # non-negative count family
  expect_s3_class(res$model, "gam")
})

test_that("auto candidates only include families with comparable AIC", {
  cands <- function(y) names(family_candidates(profile_response(y)))
  set.seed(9)
  # zero-heavy, overdispersed counts: no Tweedie (density vs probabilities)
  expect_setequal(cands(rnbinom(300, mu = 0.5, size = 0.5)), c("poisson", "nb"))
  expect_setequal(cands(rpois(300, 5)), c("poisson", "nb"))
  # proportions with exact 0/1: squeezed beta only
  expect_identical(cands(c(0, 1, runif(100))), "betar_sv")
  # strictly positive continuous: all densities on the same response
  expect_setequal(cands(rgamma(300, shape = 0.5)), c("Gamma_log", "tw"))
})

test_that("auto does not pick Tweedie on low Poisson counts", {
  # Tweedie's density at integer values used to beat the Poisson/nb
  # probabilities by hundreds of AIC units here.
  set.seed(1)
  x  <- rnorm(400)
  df <- data.frame(x = x, y = rpois(400, 0.5 * exp(0.5 * sin(2 * x))))
  fs <- suggest_family(df, "y", "s(x)", knots = NULL)
  expect_match(fs$chosen, "poisson|nb")
})

test_that("a user-supplied shape family still selects the real driver", {
  # With theta left free in the search, adding the driver raised theta and the
  # held-out deviance with it, so the driver was rejected.
  skip_on_cran()
  set.seed(11)
  n  <- 300
  df <- data.frame(lon = runif(n, -20, -5), lat = runif(n, 35, 50),
                   sst = rnorm(n), noise = rnorm(n))
  df$catch <- rnbinom(n, mu = exp(1 + 0.8 * sin(2 * df$sst)), size = 3)
  for (fam in list(nb, nb(), tw)) {
    res <- select_gam_covariates(df, "catch", c("sst", "noise"), structural = "te(lon, lat)",
                                 family = fam, cv_k = 3, cv_scheme = "stratified",
                                 outdir = tempdir(), logfile = NULL, verbose = FALSE)
    expect_identical(res$kept, "sst")
  }
})

test_that("make_select_gen fixes the shape and keeps the fitted link", {
  set.seed(3)
  x <- rnorm(200)
  y <- rnbinom(200, mu = exp(1 + sin(x)), size = 2)
  free <- function() nb(link = "sqrt")
  gen  <- make_select_gen(free, gam(y ~ s(x), family = free()))
  expect_identical(gen()$n.theta, 0)
  expect_identical(gen()$link, "sqrt")

  # nothing to fix: the generator passes through
  pois <- function() poisson()
  expect_identical(make_select_gen(pois, gam(y ~ s(x), family = pois())), pois)
  expect_identical(make_select_gen(free, NULL), free)

  # a shape that cannot be fixed warns and stays free
  z <- gam(pmin(y, 3) + 1 ~ s(x), family = ocat(R = 4))
  expect_warning(g <- make_select_gen(free, z), "cannot be held fixed")
  expect_identical(g, free)
})

test_that("a user family object starts every fit from its supplied theta", {
  set.seed(3)
  x <- rnorm(200)
  y <- rnbinom(200, mu = exp(1 + sin(x)), size = 2)
  f   <- nb()
  gen <- family_object_gen(f)
  th0 <- f$getTheta()

  gam(y ~ s(x), family = gen())          # fit mutates the shared object's theta
  expect_false(isTRUE(all.equal(f$getTheta(), th0)))
  expect_equal(gen()$getTheta(), th0)    # next fit starts fresh again

  # plain families pass through unchanged
  expect_identical(family_object_gen(poisson())(), poisson())
})
