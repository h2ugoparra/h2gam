# =============================================================================
# Runnable demo of select_gam_covariates() on synthetic data.
#
# Self-contained: no external files, no project preprocessing. It builds a
# spatial dataset where only some candidates truly drive the response, then
# shows that the selection keeps the real drivers and drops the noise. Run with
#   Rscript inst/examples/synthetic_demo.R
# or step through it interactively.
# =============================================================================

library(h2gam)

set.seed(42)
n <- 800

# --- Build a synthetic survey -------------------------------------------------
# Coordinates over a rectangular footprint, plus candidate covariates. Two of
# them (sst, chl) genuinely drive the response; noise1/noise2 do not; chl_dup is
# a near-copy of chl to exercise the concurvity gate.
df <- data.frame(
  lon    = runif(n, -20, -5),
  lat    = runif(n,  35,  50),
  sst    = rnorm(n),
  chl    = rnorm(n),
  noise1 = rnorm(n),
  noise2 = rnorm(n)
)
df$chl_dup <- df$chl + rnorm(n, sd = 0.05)   # ~collinear with chl

# A smooth spatial field plus nonlinear effects of the true drivers. The
# response is an overdispersed non-negative count, so family="auto" should land
# on a count family (Poisson / negative binomial / Tweedie).
spatial <- 0.6 * sin(df$lon / 4) + 0.6 * cos(df$lat / 4)
eta     <- 1.0 + spatial + 0.8 * sin(df$sst) - 0.5 * df$chl^2
df$catch <- rnbinom(n, mu = exp(eta), size = 2)

candidates <- c("sst", "chl", "chl_dup", "noise1", "noise2")

# --- Run selection ------------------------------------------------------------
# Structural spatial term te(lon, lat) is forced in; candidates must improve
# cross-validated held-out deviance over and above it. Spatial-block folds
# (the default) need the sf package for the euclidean metric; use haversine to
# stay dependency-free in this demo.
res <- select_gam_covariates(
  data       = df,
  response   = "catch",
  candidates = candidates,
  structural = "te(lon, lat)",
  cv_k       = 4,
  cv_metric  = "haversine",
  outdir     = tempdir()
)

cat("\n--- Demo result ---\n")
cat("chosen family:", res$family, "\n")
cat("kept covariates:", paste(res$kept, collapse = ", "), "\n")
cat("deviance explained:", sprintf("%.1f%%", 100 * res$dev_expl), "\n")
cat("outputs written to:", tempdir(), "\n")
# Expected: sst and chl kept; noise1/noise2 dropped for no CV improvement; and
# chl_dup dropped (concurvity with chl, or redundant once chl is in).
