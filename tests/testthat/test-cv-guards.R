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

test_that("every random-effect form is excluded from held-out prediction", {
  expect_true(is_re_smooth(s(v, bs = "re")))
  expect_true(is_re_smooth(s(x, v, bs = "fs")))
  expect_true(is_re_smooth(te(x, v, bs = c("tp", "re"))))
  expect_false(is_re_smooth(s(x)))
  expect_false(is_re_smooth(1))

  # A strong per-vessel effect with one vessel per fold: held-out vessels are
  # unseen, so a term left in the prediction costs the fold or its accuracy.
  df <- sim_counts(n = 240)
  df$vessel <- factor(rep(letters[1:6], each = 40))
  df$catch  <- rpois(240, exp(1 + 0.8 * sin(2 * df$sst) + rep(rnorm(6), each = 40)))
  folds <- rep(rep(1:3, each = 40), 2)
  for (term in c("s(vessel, bs='re')", "s(vessel, sst, bs='re')",
                 "s(vessel, bs='re', by=eff)", "s(lat, vessel, bs='fs', k=4)")) {
    fw <- forward_select(df, "catch", "sst", c("te(lon, lat, k=4)", term),
                         family_gen = function() poisson(), folds = folds,
                         knots = NULL, verbose = FALSE)
    expect_true(is.finite(fw$scores$cv_gain), info = term)
  }
})

test_that("a random-effect term with no factor is rejected whatever its form", {
  df <- sim_counts()
  expect_error(run(df, structural = "s(sst, year, bs='re')"),
               "needs 'sst' or 'year' to be a factor")
  expect_error(run(df, structural = "s(lat, nope, bs='fs')"), "column 'nope' not found")
})

test_that("decision-table na_frac reflects the data as supplied", {
  skip_on_cran()
  df <- sim_counts()
  df$sst[1:30] <- NA
  res <- run(df, structural = "te(lon, lat)")
  expect_equal(res$decision$na_frac[res$decision$covariate == "sst"], 30 / nrow(df))
})

test_that("the decision table carries the CV score and reason behind each decision", {
  skip_on_cran()
  df <- sim_counts()
  df$const <- 1
  res <- select_gam_covariates(df, "catch", c("sst", "noise", "const"), cv_k = 3,
                               structural = "te(lon, lat)", cv_scheme = "stratified",
                               outdir = tempdir(), logfile = NULL, verbose = FALSE)
  d <- res$decision
  row <- function(v) d[d$covariate == v, ]

  # kept: the score at its entry step
  acc <- res$path[res$path$var == "sst" & res$path$action == "accept", ]
  expect_equal(row("sst")$cv_gain, acc$cv_gain)
  expect_equal(row("sst")$cv_se, acc$cv_se)
  expect_gt(row("sst")$cv_gain, row("sst")$cv_se)

  # dropped in the search: scored, and the reason agrees with the score
  n <- row("noise")
  expect_true(is.finite(n$cv_gain) && is.finite(n$cv_se))
  expect_identical(n$reason, if (n$cv_gain <= 0) "no CV improvement" else "CV gain below 1 SE")
  if (n$cv_gain > 0) expect_lt(n$cv_gain, n$cv_se)

  # screened out before the search: no score
  expect_identical(row("const")$reason, "zero variance")
  expect_true(is.na(row("const")$cv_gain))
})

test_that("a candidate whose fits fail is reported as a scoring failure", {
  df <- sim_counts(n = 120)
  folds <- rep_len(1:3, nrow(df))
  # s(bad) cannot be fitted (non-numeric), so every fold is NA for it
  df$bad <- sample(letters, nrow(df), TRUE)
  fw <- suppressWarnings(
    forward_select(df, "catch", "bad", "te(lon, lat)", family_gen = function() poisson(),
                   folds = folds, knots = NULL, verbose = FALSE))
  expect_length(fw$accepted, 0)
  expect_identical(fw$scores$var, "bad")
  expect_false(is.finite(fw$scores$cv_gain))
})

test_that("the concurvity gate fails closed when concurvity cannot be computed", {
  set.seed(5)
  n  <- 240
  df <- data.frame(lon = runif(n, -20, -5), lat = runif(n, 35, 50),
                   a = rnorm(n), b = rnorm(n))
  df$y <- sin(2 * df$a) + 0.6 * sin(2 * df$b) + rnorm(n, sd = 0.3)
  fs <- function() forward_select(df, "y", c("a", "b"), "te(lon, lat, k=4)",
                                  family_gen = function() gaussian(),
                                  folds = rep_len(1:3, n), knots = NULL, verbose = FALSE)
  expect_setequal(fs()$accepted, c("a", "b"))        # both real drivers pass normally

  local_mocked_bindings(concurvity = function(...) stop("boom"))
  m <- gam(y ~ s(a) + s(b), data = df)
  expect_identical(concurvity_new(m, "b", character()), 0)   # nothing to compare against
  expect_true(is.na(concurvity_new(m, "b", "a")))

  # the second driver can no longer be gated: rejected and reported, not accepted
  expect_warning(fw <- fs(), "concurvity could not be computed for b")
  expect_identical(fw$accepted, "a")
  rej <- fw$path[fw$path$action == "reject-concurvity", ]
  expect_identical(unique(rej$var), "b")
  expect_true(all(is.na(rej$concurvity)))
})

test_that("a concurvity-rejected candidate is not scored again", {
  set.seed(6)
  n  <- 240
  df <- data.frame(lon = runif(n, -20, -5), lat = runif(n, 35, 50),
                   a = rnorm(n), b = rnorm(n), c = rnorm(n))
  df$y <- 1.5 * sin(2 * df$a) + sin(2 * df$b) + 0.6 * sin(2 * df$c) + rnorm(n, sd = 0.3)
  # report b as concurve with a, everything else as independent
  local_mocked_bindings(concurvity = function(model, ...) {
    labs <- vapply(model$smooth, `[[`, "", "label")
    W <- matrix(0, length(labs), length(labs), dimnames = list(labs, labs))
    if (all(c("s(a)", "s(b)") %in% labs)) W["s(b)", "s(a)"] <- 1
    list(worst = W)
  })
  fw <- forward_select(df, "y", c("a", "b", "c"), "te(lon, lat, k=4)",
                       family_gen = function() gaussian(),
                       folds = rep_len(1:3, n), knots = NULL, verbose = FALSE)
  expect_identical(fw$accepted, c("a", "c"))
  # rejected at step 2, before c is accepted; not revisited at step 3
  expect_identical(fw$path$action[fw$path$var == "b"], "reject-concurvity")
  expect_identical(max(fw$path$step), 2)
})

test_that("the result reports the CV gain of the kept set and final-model checks", {
  skip_on_cran()
  res <- run(sim_counts(), structural = "te(lon, lat)")
  expect_identical(res$kept, "sst")
  # one covariate kept: the summary is its accepted step
  acc <- res$path[res$path$action == "accept", ]
  expect_equal(res$cv$gain, acc$cv_gain)
  expect_equal(res$cv$se, acc$cv_se)
  expect_equal(res$cv$structural - res$cv$selected, res$cv$gain)
  expect_identical(res$cv$folds, 3L)
  # the per-fold table carries the scores behind those means
  bf <- res$cv$by_fold
  expect_identical(bf$fold, 1:3)
  expect_identical(bf$n, as.vector(table(res$folds)))
  expect_equal(mean(bf$structural), res$cv$structural)
  expect_equal(mean(bf$gain), res$cv$gain)

  expect_true(is.finite(res$final_check$dispersion))
  p <- res$final_check$resid_ks_p
  expect_true(is.na(p) || (p >= 0 && p <= 1))
})

test_that("a CV fold lost to the reference model is reported", {
  df <- sim_counts(n = 120)
  folds <- rep(1:3, c(30, 40, 50))          # training sizes 90, 80, 70
  real_fit <- fit_search
  # every fit that holds out fold 2 fails
  local_mocked_bindings(fit_search = function(form, data, ...) {
    if (nrow(data) == 80) NULL else real_fit(form, data, ...)
  })
  expect_warning(
    fw <- forward_select(df, "catch", "sst", "te(lon, lat, k=4)",
                         family_gen = function() poisson(), folds = folds,
                         knots = NULL, verbose = FALSE),
    "only 2 of 3 CV folds")
  expect_true(is.finite(fw$scores$cv_gain))
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

test_that("the backward re-check drops a covariate made redundant by a later one", {
  # per-fold held-out deviance by covariate set: b alone does as well as a + b
  dev <- list("a,b" = c(1.00, 1.10, 0.90), "b" = c(1.01, 1.09, 0.90),
              "a" = c(1.50, 1.60, 1.40), "a,b,c" = c(0.80, 0.90, 0.70),
              "b,c" = c(0.81, 0.89, 0.70), "a,c" = c(1.30, 1.45, 1.20),
              "c" = c(1.40, 1.55, 1.30))
  fake <- function(accepted) list(
    accepted = accepted, acc_dev = dev[[paste(sort(accepted), collapse = ",")]],
    fold_dev = function(terms) dev[[paste(sort(terms), collapse = ",")]])

  bw <- backward_recheck(fake(c("a", "b")), verbose = FALSE)
  expect_identical(bw$kept, "b")
  expect_identical(bw$dev, dev[["b"]])
  expect_identical(bw$dropped$var, "a")
  expect_lt(bw$dropped$cv_gain, bw$dropped$cv_se)

  # after a drop the survivors are re-tested, and both still earn their place
  bw <- backward_recheck(fake(c("a", "b", "c")), verbose = FALSE)
  expect_identical(bw$kept, c("b", "c"))

  # nothing accepted, or one covariate: nothing to re-test
  expect_identical(backward_recheck(fake(character()), verbose = FALSE)$kept, character())
  one <- backward_recheck(list(accepted = "a", acc_dev = dev[["a"]],
                               fold_dev = function(terms) stop("not called")), verbose = FALSE)
  expect_identical(one$kept, "a")
  expect_null(one$dropped)
})
