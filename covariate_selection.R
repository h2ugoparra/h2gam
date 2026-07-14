# =============================================================================
# covariate_selection.R
#
# Automated GAM covariate validation & selection for mgcv, data-independent.
#
# Given a response and a pool of candidate covariates it:
#   1. profiles the response and auto-suggests a family (count / proportion /
#      positive / gaussian branches), auto-ranked by AIC + residual diagnostics;
#   2. forward-selects a DECORRELATED, prediction-relevant covariate subset,
#      scored by spatial-block CV skill, with a concurvity gate and a BIC
#      tie-break / 1-SE stopping rule;
#   3. runs a bidirectional backward re-check, a final select=TRUE shrinkage
#      fit, and adequacy diagnostics.
#
# Every threshold is a logged, overridable argument; nothing is hardcoded to a
# particular dataset. Structural (spatial/temporal) terms are forced in and not
# screened -- covariates must earn their place OVER AND ABOVE them.
#
# Design decisions and rationale: see the approved plan.
# =============================================================================

suppressPackageStartupMessages({
  library(mgcv)
  library(cluster)
})

# ---------------------------------------------------------------------------
# 1. Response profiling & family cascade
# ---------------------------------------------------------------------------

profile_response <- function(y) {
  y <- y[is.finite(y)]
  list(
    n          = length(y),
    min        = min(y),
    max        = max(y),
    integer    = all(abs(y - round(y)) < 1e-8),
    binary     = all(y %in% c(0, 1)),
    nonneg     = all(y >= 0),
    positive   = all(y > 0),
    in01       = all(y >= 0 & y <= 1),
    frac_zero  = mean(y == 0),
    frac_one   = mean(y == 1),
    mean       = mean(y),
    var        = stats::var(y),
    dispersion = if (mean(y) > 0) stats::var(y) / mean(y) else NA_real_,
    skew       = mean((y - mean(y))^3) / stats::sd(y)^3
  )
}

# Ordered list of candidate family generators for the profiled response.
# Each element is a zero-arg closure returning a FRESH mgcv family (fresh
# matters: nb()/tw() carry a mutable theta/p estimate in their environment).
family_candidates <- function(p) {
  cand <- list()
  add <- function(name, gen) cand[[name]] <<- gen

  if (p$binary) {
    # ---- presence/absence ----
    # logit is the conventional default, but cloglog is the mechanistically
    # consistent link for effort-offset longline data: with offset(log(nhooks)),
    #   cloglog(p) = log(nhooks) + f  <=>  p = 1 - exp(-nhooks * exp(f)),
    # i.e. the probability of >=1 catch under a Poisson encounter rate ~ hooks
    # (the presence analogue of the log-link count model). probit is included as
    # a cheap symmetric alternative. AIC is comparable across links here (same
    # response and likelihood); the existing AIC + DHARMa ranking picks among them.
    add("binomial_logit",   function() binomial(link = "logit"))
    add("binomial_cloglog", function() binomial(link = "cloglog"))
    add("binomial_probit",  function() binomial(link = "probit"))

  } else if (p$integer && p$nonneg) {
    # ---- count branch ----
    if (!is.na(p$dispersion) && p$dispersion > 1.5) {
      add("nb",      function() nb())        # overdispersed -> negative binomial
      add("tw",      function() tw())        # Tweedie handles zero-mass + tail
      add("poisson", function() poisson())
    } else {
      add("poisson", function() poisson())
      add("nb",      function() nb())
    }
    if (p$frac_zero > 0.5) add("tw", function() tw())

  } else if (p$in01 && p$frac_zero == 0 && p$frac_one == 0) {
    # ---- proportion, open interval ----
    add("betar",    function() betar())
    add("gaussian", function() gaussian())

  } else if (p$in01) {
    # ---- proportion with 0/1 boundary (caller applies SV squeeze for betar) ----
    add("betar_sv", function() betar())
    add("tw",       function() tw())

  } else if (p$positive) {
    # ---- strictly positive continuous ----
    if (!is.na(p$skew) && p$skew > 1) {
      add("Gamma_log", function() Gamma(link = "log"))
      add("tw",        function() tw())
    } else {
      add("gaussian",  function() gaussian())
      add("Gamma_log", function() Gamma(link = "log"))
    }

  } else {
    # ---- real-valued fallback ----
    add("gaussian", function() gaussian())
  }
  cand
}

# Smithson-Verkuilen squeeze so betar can take exact 0/1 on a (0,1) response.
sv_squeeze <- function(y) {
  n <- sum(is.finite(y))
  (y * (n - 1) + 0.5) / n
}

# DHARMa uniformity p-value (NA if the family can't be simulated).
dharma_ks <- function(m) {
  if (!requireNamespace("DHARMa", quietly = TRUE)) return(NA_real_)
  out <- tryCatch({
    sim <- DHARMa::simulateResiduals(m, n = 250, plot = FALSE)
    DHARMa::testUniformity(sim, plot = FALSE)$p.value
  }, error = function(e) NA_real_)
  out
}

# Fit each candidate family on the structural-only model, rank by AIC (ok
# families first), and flag dispersion / residual uniformity. Returns the
# chosen family generator plus a diagnostic table.
suggest_family <- function(data, response, structural, knots, offset = NULL) {
  p     <- profile_response(data[[response]])
  cands <- family_candidates(p)
  rows  <- list()
  fits  <- list()

  # Rank families on the SAME offset model the pipeline deploys. This matters for
  # the binomial cloglog candidate (its log(nhooks) offset is what gives it the
  # "P(>=1 catch)" meaning) and keeps AIC comparable to the search/final fits.
  off_term <- if (!is.null(offset)) sprintf("offset(log(%s))", offset) else NULL

  for (nm in names(cands)) {
    fam  <- cands[[nm]]()
    resp <- response
    dat  <- data
    if (nm == "betar_sv") {                # squeeze the response for beta
      dat[[".y_sv"]] <- sv_squeeze(dat[[response]])
      resp <- ".y_sv"
    }
    form <- stats::as.formula(paste(resp, "~", paste(c(structural, off_term), collapse = " + ")))
    m <- tryCatch(gam(form, data = dat, family = fam, method = "REML", knots = knots),
                  error = function(e) e)
    if (inherits(m, "error")) {
      rows[[nm]] <- data.frame(family = nm, AIC = NA, dispersion = NA,
                               resid_ks_p = NA, ok = FALSE)
      next
    }
    disp <- sum(residuals(m, type = "pearson")^2) / m$df.residual
    rows[[nm]] <- data.frame(family = nm, AIC = AIC(m), dispersion = disp,
                             resid_ks_p = dharma_ks(m), ok = TRUE)
    fits[[nm]] <- m
  }

  tab <- do.call(rbind, rows)
  tab <- tab[order(!tab$ok, tab$AIC), ]
  chosen <- tab$family[which(tab$ok)[1]]
  if (is.na(chosen)) {
    print(tab)
    stop("suggest_family: no candidate family could be fitted (table above).")
  }
  rownames(tab) <- NULL

  list(profile = p, table = tab, chosen = chosen,
       family_gen = cands[[chosen]], squeeze = identical(chosen, "betar_sv"),
       fit = fits[[chosen]])
}

# Family generator with the shape parameter (Tweedie p / negbin theta) held
# FIXED at its structural-only estimate. Re-estimating it on every CV fit is the
# dominant cost of the search; fixing it gives a large speedup with negligible
# effect on ranking. The final fit still uses the free (re-estimating) family.
make_select_gen <- function(chosen, free_gen, fit) {
  if (is.null(fit)) return(free_gen)
  th <- tryCatch(fit$family$getTheta(TRUE), error = function(e) NA_real_)
  if (identical(chosen, "tw") && is.finite(th))
    return(function() Tweedie(p = th, link = power(0)))
  if (identical(chosen, "nb") && is.finite(th))
    return(function() negbin(theta = th, link = "log"))
  free_gen
}

# ---------------------------------------------------------------------------
# 2. Spatial-block folds -- SPCV (mirrors h2ml.features.spatial_cv.SPCVSplitter)
# ---------------------------------------------------------------------------
# Two-stage spatial CV (Wang et al. 2023), faithful R port of the h2ml method:
#   Stage 1  AHC blocks: agglomerative hierarchical clustering on coordinates
#            (euclidean, or haversine for lat/lon degrees) cut at a distance
#            threshold (default = 10th percentile of pairwise distances).
#   Stage 2  Cluster ensemble (HBGF): per-block means of location / covariates /
#            target -> KMeans on each -> one-hot memberships stacked into a
#            co-occurrence affinity -> spectral clustering into k folds.
# Folds are geographically separated AND representative of covariate+label space.

# Pairwise great-circle angular distance (radians), matching sklearn "haversine".
.haversine_dist <- function(lat, lon) {
  latr <- lat * pi / 180; lonr <- lon * pi / 180
  a <- sin(outer(latr, latr, "-") / 2)^2 +
       (cos(latr) %o% cos(latr)) * sin(outer(lonr, lonr, "-") / 2)^2
  root <- sqrt(a); root[root > 1] <- 1        # clamp in place (keeps matrix dims)
  2 * asin(root)
}

# Normalised spectral clustering (Ng-Jordan-Weiss) on a precomputed affinity,
# matching sklearn SpectralClustering(affinity="precomputed", assign_labels="kmeans").
.spectral_cluster <- function(S, k, seed) {
  dinv <- 1 / sqrt(pmax(rowSums(S), .Machine$double.eps))
  L    <- outer(dinv, dinv) * S                 # D^-1/2 S D^-1/2 (normalised affinity)
  U    <- eigen(L, symmetric = TRUE)$vectors[, seq_len(k), drop = FALSE]
  rn   <- sqrt(rowSums(U^2)); rn[rn == 0] <- 1
  set.seed(seed)
  stats::kmeans(U / rn, centers = k, nstart = 10, iter.max = 100)$cluster
}

make_spatial_folds <- function(data, response, candidates, coords = c("lon", "lat"),
                               k = 10, metric = "haversine", threshold = NULL,
                               linkage = "average", pca_var = 0.95, seed = 42,
                               verbose = TRUE) {
  lon <- data[[coords[1]]]; lat <- data[[coords[2]]]

  # -- Stage 1: AHC -> spatially coherent blocks --------------------------
  if (metric == "haversine") {
    dmat <- .haversine_dist(lat, lon); d <- stats::as.dist(dmat)
    link <- if (linkage == "ward") "average" else linkage    # ward invalid for haversine
    thr  <- if (is.null(threshold)) stats::quantile(d, 0.10) else threshold * pi / 180
    loc  <- cbind(lat, lon)                                  # stage-2 location view (degrees)
  } else {
    if (!requireNamespace("sf", quietly = TRUE))
      stop("metric='euclidean' needs the sf package to reproject lon/lat.")
    # reproject lon/lat -> Lambert Conformal Conic (metres) so euclidean distance is valid
    lcc <- "+proj=lcc +lat_1=25 +lat_2=45 +lat_0=35 +lon_0=-25 +datum=WGS84 +units=m +no_defs"
    loc  <- sf::sf_project(from = "+proj=longlat +datum=WGS84", to = lcc, pts = cbind(lon, lat))
    d    <- stats::dist(loc, method = "euclidean")           # stage-2 location view (metres)
    link <- linkage
    thr  <- if (is.null(threshold)) stats::quantile(d, 0.10) else threshold
  }
  blocks <- stats::cutree(stats::hclust(d, method = link), h = as.numeric(thr))
  blocks <- match(blocks, sort(unique(blocks)))              # compact 1..B
  B <- max(blocks)
  if (B < k) stop(sprintf("SPCV: AHC produced %d block(s) < k=%d; lower the threshold.", B, k))
  if (B < 2 * k)
    warning(sprintf("SPCV: only %d blocks for %d folds -- affinity may be weakly connected.", B, k))

  # -- Stage 2: cluster ensemble (HBGF) -> folds --------------------------
  num  <- candidates[vapply(candidates, function(v) is.numeric(data[[v]]), logical(1))]
  Xmat <- as.matrix(data[, num, drop = FALSE])
  bmean <- function(v) as.numeric(tapply(v, blocks, mean))
  block_coords <- cbind(bmean(loc[, 1]), bmean(loc[, 2]))    # projected coords when euclidean
  block_y      <- matrix(bmean(data[[response]]), ncol = 1)
  block_X      <- vapply(seq_len(ncol(Xmat)), function(j) bmean(Xmat[, j]), numeric(B))

  if (ncol(block_X) > 1) {                                    # PCA to pca_var variance
    pc   <- stats::prcomp(scale(block_X), center = FALSE, scale. = FALSE)
    npc  <- max(1, which(cumsum(pc$sdev^2) / sum(pc$sdev^2) >= pca_var)[1])
    block_Xp <- pc$x[, seq_len(npc), drop = FALSE]
  } else block_Xp <- block_X

  set.seed(seed)
  km  <- function(m) stats::kmeans(m, centers = k, nstart = 10, iter.max = 100)$cluster
  eye <- diag(k)
  H   <- cbind(eye[km(block_coords), ], eye[km(block_Xp), ], eye[km(block_y), ])  # (B, 3k)
  S   <- H %*% t(H)
  S   <- (S + 1 / (B + 1)); S <- S / max(S)                  # floor for connectivity, normalise

  folds <- .spectral_cluster(S, k, seed)[blocks]
  if (verbose) cat(sprintf("  SPCV (%s): %d AHC blocks -> %d folds (sizes %s)\n",
                           metric, B, k, paste(tabulate(folds, k), collapse = "/")))
  as.integer(folds)
}

# Stratified k-fold folds -- a non-spatial alternative to make_spatial_folds for
# rare-event responses. SPCV's spatially coherent blocks can leave a fold with
# zero positives when prevalence is low (~4%), making held-out binomial deviance
# degenerate. Stratifying on the response spreads the rare positives evenly, so
# every fold is evaluable. NB: this deliberately drops spatial separation --
# spatial autocorrelation can leak between train/test, so CV skill is optimistic
# relative to SPCV; the tradeoff buys stable fold estimates for a rare event.
# Returns the same integer fold vector (1..k) that forward_select consumes.
make_stratified_folds <- function(data, response, k = 5, seed = 1, verbose = TRUE) {
  y <- data[[response]]
  # Strata: exact label values for a few-valued response (binary / integer
  # label), else quantile bins so a continuous response stays balanced too.
  strata <- if (length(unique(y[is.finite(y)])) <= max(k, 10)) {
    as.integer(factor(y))
  } else {
    br <- unique(stats::quantile(y, seq(0, 1, length.out = k + 1), na.rm = TRUE))
    as.integer(cut(y, breaks = br, include.lowest = TRUE))
  }
  set.seed(seed)
  folds <- integer(length(y))
  for (s in unique(strata)) {                 # cyclic assignment within a stratum
    idx <- which(strata == s)                 # -> near-equal class counts / fold
    folds[sample(idx)] <- rep_len(seq_len(k), length(idx))
  }
  minc <- min(tabulate(strata))
  if (minc < k)
    warning(sprintf("stratified CV: rarest class has %d case(s) < k=%d; some folds will lack it.",
                    minc, k))
  if (verbose) {
    pos <- if (all(y %in% c(0, 1))) paste0("; pos ",
             paste(tapply(y, folds, function(v) sum(v == 1)), collapse = "/")) else ""
    cat(sprintf("  stratified %d-fold (sizes %s%s)\n", k,
                paste(tabulate(folds, k), collapse = "/"), pos))
  }
  as.integer(folds)
}

# ---------------------------------------------------------------------------
# 3. Univariate screen + Spearman clustering (annotation only)
# ---------------------------------------------------------------------------

univariate_screen <- function(data, candidates, na_max = 0.30) {
  keep <- character(); dropped <- character(); reason <- character()
  for (v in candidates) {
    x <- data[[v]]
    if (is.null(x)) { dropped <- c(dropped, v); reason <- c(reason, "absent"); next }
    if (!is.numeric(x)) { dropped <- c(dropped, v); reason <- c(reason, "non-numeric"); next }
    na_frac <- mean(is.na(x) | !is.finite(x))
    if (na_frac > na_max) {
      dropped <- c(dropped, v); reason <- c(reason, sprintf("NA %.0f%%", 100 * na_frac)); next
    }
    if (stats::sd(x, na.rm = TRUE) < sqrt(.Machine$double.eps)) {
      dropped <- c(dropped, v); reason <- c(reason, "zero variance"); next
    }
    keep <- c(keep, v)
  }
  list(keep = keep, dropped = dropped, reason = stats::setNames(reason, dropped))
}

# Single-linkage clusters on |Spearman rho| > thr. Used only to annotate the
# decision table; the forward loop's concurvity gate does the actual removal.
spearman_clusters <- function(data, vars, thr = 0.7) {
  num <- vars[sapply(vars, function(v) is.numeric(data[[v]]))]
  if (length(num) < 2) return(stats::setNames(seq_along(vars), vars))
  rho <- suppressWarnings(stats::cor(data[, num], method = "spearman", use = "pairwise"))
  adj <- abs(rho) > thr
  adj[is.na(adj)] <- FALSE               # NA rho (e.g. degenerate column) = no edge
  cl  <- stats::setNames(seq_along(num), num)          # start: each its own
  repeat {                                             # merge until stable
    changed <- FALSE
    for (i in seq_along(num)) for (j in seq_along(num)) {
      if (adj[i, j] && cl[i] != cl[j]) { cl[cl == cl[j]] <- cl[i]; changed <- TRUE }
    }
    if (!changed) break
  }
  out <- stats::setNames(rep(NA_integer_, length(vars)), vars)
  out[num] <- as.integer(factor(cl))
  out
}

# ---------------------------------------------------------------------------
# 4. Forward selection (CV-scored, concurvity-gated)
# ---------------------------------------------------------------------------

# Fast fitter for the search: bam(discrete=TRUE) is ~5x faster than gam on the
# repeated te()-smooth fits and predictions correlate ~0.99. Requires a fixed
# shape family (see make_select_gen). Falls back to gam if bam errors.
fit_search <- function(form, data, family_gen, knots) {
  m <- tryCatch(bam(form, data = data, family = family_gen(),
                    method = "fREML", discrete = TRUE, knots = knots),
                error = function(e) NULL)
  if (is.null(m))
    m <- tryCatch(gam(form, data = data, family = family_gen(),
                      method = "REML", knots = knots), error = function(e) NULL)
  m
}

# Concurvity of a newly added smooth against the already-accepted covariate
# smooths only (structural terms deliberately excluded). Returns the worst
# pairwise value, 0 when nothing to compare against.
concurvity_new <- function(model, newvar, accepted) {
  if (length(accepted) == 0) return(0)
  cc <- tryCatch(concurvity(model, full = FALSE), error = function(e) NULL)
  if (is.null(cc)) return(0)
  W    <- cc$worst
  rlab <- intersect(paste0("s(", newvar, ")"), rownames(W))
  clab <- intersect(paste0("s(", accepted, ")"), colnames(W))
  if (length(rlab) == 0 || length(clab) == 0) return(0)
  max(W[rlab, clab, drop = FALSE])
}

forward_select <- function(data, response, candidates, structural, family_gen,
                           folds, knots, concurvity_max = 0.8, offset = NULL,
                           re_vars = character(), re_labels = character(),
                           verbose = TRUE) {
  off_term <- if (!is.null(offset)) sprintf("offset(log(%s))", offset) else NULL
  fold_ids <- sort(unique(folds))

  build_form <- function(terms) stats::as.formula(
    paste(response, "~", paste(c(structural, off_term, sprintf("s(%s)", terms)), collapse = " + "))
  )
  fit_full <- function(terms) fit_search(build_form(terms), data, family_gen, knots)

  # Per-fold held-out mean deviance for a covariate set (lower is better).
  fold_dev <- function(terms) {
    form <- build_form(terms)
    d <- rep(NA_real_, length(fold_ids))
    for (i in seq_along(fold_ids)) {
      tr <- data[folds != fold_ids[i], ]; te <- data[folds == fold_ids[i], ]
      m  <- fit_search(form, tr, family_gen, knots)
      if (is.null(m)) next
      # Random-effect terms (e.g. s(year, bs='re')) are excluded from held-out
      # prediction -- their level-specific deviation is unknowable out of fold
      # (matches the deployment predict(..., exclude=...)). Remap the factor to a
      # training level so the design matrix builds; the excluded term zeroes it.
      for (v in re_vars) te[[v]] <- factor(as.character(tr[[v]][1]), levels = levels(tr[[v]]))
      mu <- tryCatch(predict(m, newdata = te, type = "response",
                             exclude = re_labels), error = function(e) NULL)
      if (is.null(mu) || length(mu) != nrow(te)) next
      dr <- tryCatch(sum(m$family$dev.resids(te[[response]], mu, rep(1, nrow(te)))) / nrow(te),
                     error = function(e) NA_real_)
      if (is.finite(dr)) d[i] <- dr
    }
    d
  }

  accepted  <- character()
  remaining <- candidates
  cur_dev   <- fold_dev(accepted)          # structural-only baseline
  path      <- list()
  step      <- 0

  repeat {
    step <- step + 1
    scored <- lapply(remaining, function(v) {
      cd   <- fold_dev(c(accepted, v))
      diff <- cur_dev - cd                 # paired per-fold improvement
      list(v = v, cd = cd,
           mean = mean(diff, na.rm = TRUE),
           se   = stats::sd(diff, na.rm = TRUE) / sqrt(sum(is.finite(diff))))
    })
    ord <- order(vapply(scored, function(s) -s$mean, numeric(1)))

    picked <- NULL
    for (j in ord) {
      s <- scored[[j]]
      if (!is.finite(s$mean) || s$mean <= 0) break     # ordered desc: nothing better remains
      if (!is.finite(s$se) || s$mean < s$se) next      # NA se (<2 scored folds) or fails 1-SE bar
      m  <- fit_full(c(accepted, s$v))
      cg <- if (is.null(m)) NA_real_ else concurvity_new(m, s$v, accepted)
      if (is.finite(cg) && cg > concurvity_max) {
        path[[length(path) + 1]] <- data.frame(step = step, var = s$v, action = "reject-concurvity",
                                                concurvity = cg, cv_gain = s$mean, cv_se = s$se)
        next
      }
      picked <- s; picked$concurvity <- cg; break
    }

    if (is.null(picked)) break
    path[[length(path) + 1]] <- data.frame(step = step, var = picked$v, action = "accept",
                                            concurvity = picked$concurvity,
                                            cv_gain = picked$mean, cv_se = picked$se)
    if (verbose) cat(sprintf("  [step %d] + %-18s  CV gain %.4f (se %.4f)  concurvity %.2f\n",
                             step, picked$v, picked$mean, picked$se, picked$concurvity))
    accepted  <- c(accepted, picked$v)
    remaining <- setdiff(remaining, picked$v)
    cur_dev   <- picked$cd
    if (length(remaining) == 0) break
  }

  list(accepted = accepted,
       path = if (length(path)) do.call(rbind, path) else NULL)
}

# ---------------------------------------------------------------------------
# 5. Backward re-check (BIC-based) & final shrinkage fit
# ---------------------------------------------------------------------------

backward_recheck <- function(data, response, accepted, structural, family_gen,
                             knots, offset = NULL, verbose = TRUE) {
  off_term <- if (!is.null(offset)) sprintf("offset(log(%s))", offset) else NULL
  build <- function(terms) stats::as.formula(
    paste(response, "~", paste(c(structural, off_term, sprintf("s(%s)", terms)), collapse = " + "))
  )
  fit <- function(terms) fit_search(build(terms), data, family_gen, knots)
  repeat {
    if (length(accepted) == 0) break
    full <- fit(accepted); if (is.null(full)) break
    bic_full <- BIC(full)
    bic_drop <- vapply(accepted, function(v) {
      r <- fit(setdiff(accepted, v)); if (is.null(r)) Inf else BIC(r)
    }, numeric(1))
    best <- which.min(bic_drop)
    if (is.finite(bic_drop[best]) && bic_drop[best] < bic_full) {
      if (verbose) cat(sprintf("  [backward] - %-18s  BIC %.1f -> %.1f\n",
                               accepted[best], bic_full, bic_drop[best]))
      accepted <- setdiff(accepted, accepted[best])
    } else break
  }
  accepted
}

final_fit <- function(data, response, accepted, structural, family_gen, knots,
                      offset = NULL, squeeze = FALSE) {
  off_term <- if (!is.null(offset)) sprintf("offset(log(%s))", offset) else NULL
  resp <- response
  if (squeeze) { data[[".y_sv"]] <- sv_squeeze(data[[response]]); resp <- ".y_sv" }
  cov_terms <- if (length(accepted)) sprintf("s(%s)", accepted) else NULL
  form <- stats::as.formula(paste(resp, "~", paste(c(structural, off_term, cov_terms), collapse = " + ")))
  # select=TRUE: double penalty shrinks any weak survivor toward edf ~ 0.
  gam(form, data = data, family = family_gen(), method = "REML", select = TRUE, knots = knots)
}

# ---------------------------------------------------------------------------
# 6. Orchestrator
# ---------------------------------------------------------------------------

select_gam_covariates <- function(
    data, response, candidates,
    structural     = c("s(month, bs='cc')", "s(year, bs='re')", "te(lon, lat)"),
    coords         = c("lon", "lat"),
    family         = "auto",
    concurvity_max = 0.8,
    spearman_pre   = 0.7,
    na_max         = 0.30,
    cv_k           = 5,
    cv_scheme      = "spatial",      # fold construction: "spatial" (SPCV) | "stratified"
                                    #   "stratified" = non-spatial k-fold stratified on the
                                    #   response; use for low-prevalence targets where SPCV
                                    #   leaves folds with too few positives (see make_stratified_folds).
    cv_metric      = "euclidean",   # SPCV distance: "haversine" (lat/lon deg) | "euclidean"
    cv_threshold   = NULL,          # AHC block threshold (NULL = 10th pct of pairwise dist)
    cv_linkage     = "ward",     # AHC linkage ("ward" downgraded to "average" for haversine)
    knots          = NULL,          # mgcv per-smooth knot positions, passed to every gam()/bam() fit.
                                    #   Mainly the cyclic boundary for s(month, bs='cc'): pass
                                    #   list(month = c(0.5, 12.5)) so the 12-month cycle wraps (Dec->Jan)
                                    #   instead of collapsing month 1 and 12 onto the same knot. NULL lets
                                    #   mgcv place knots automatically; entries for absent vars are ignored.
    offset         = NULL,
    outdir         = ".",
    seed           = 1,
    verbose        = TRUE) {

  stopifnot(response %in% names(data))
  cv_scheme <- match.arg(cv_scheme, c("spatial", "stratified"))
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  # -- family --------------------------------------------------------------
  if (identical(family, "auto")) {
    fs <- suggest_family(data, response, structural, knots, offset = offset)
    family_gen <- fs$family_gen; squeeze <- fs$squeeze
    if (verbose) { cat("\n== Family suggestion ==\n"); print(fs$table)
      cat(sprintf("-> chosen: %s%s\n", fs$chosen, if (squeeze) " (SV squeeze)" else "")) }
  } else {
    family_gen <- if (is.function(family)) function() family() else function() family
    squeeze <- FALSE; fs <- list(chosen = "user", table = NULL, fit = NULL,
                                 profile = profile_response(data[[response]]))
  }
  # fast fixed-shape family for the search; free family for the final fit
  select_gen <- make_select_gen(fs$chosen, family_gen, fs$fit)

  # -- univariate screen ---------------------------------------------------
  scr <- univariate_screen(data, candidates, na_max = na_max)
  if (verbose) cat(sprintf("\n== Univariate screen ==\n  kept %d, dropped %d\n",
                           length(scr$keep), length(scr$dropped)))

  # complete-case on the comparison columns so every CV model sees the same rows
  need <- c(response, coords, scr$keep,
            all.vars(stats::as.formula(paste("~", paste(structural, collapse = "+")))))
  data <- data[stats::complete.cases(data[, intersect(need, names(data))]), ]

  # SV squeeze applied HERE so the search sees the same (0,1) response as the
  # final fit -- raw 0/1 values give non-finite beta deviance in the CV scoring.
  search_resp <- response
  if (squeeze) { data[[".y_sv"]] <- sv_squeeze(data[[response]]); search_resp <- ".y_sv" }

  clusters <- spearman_clusters(data, scr$keep, thr = spearman_pre)
  if (verbose) cat(sprintf("\n== %s folds ==\n",
                   if (cv_scheme == "stratified") "Stratified k-fold" else "Spatial (SPCV)"))
  folds    <- if (cv_scheme == "stratified") {
    make_stratified_folds(data, response, k = cv_k, seed = seed, verbose = verbose)
  } else {
    make_spatial_folds(data, response, scr$keep, coords = coords, k = cv_k,
                       metric = cv_metric, threshold = cv_threshold,
                       linkage = cv_linkage, seed = seed, verbose = verbose)
  }

  # random-effect structural terms -> excluded from held-out CV prediction
  re_terms  <- structural[grepl("bs\\s*=\\s*['\"]re['\"]", structural)]
  re_vars   <- sub("^\\s*s\\(\\s*([A-Za-z0-9_.]+).*", "\\1", re_terms)
  re_labels <- if (length(re_vars)) paste0("s(", re_vars, ")") else character()

  # -- forward selection ---------------------------------------------------
  if (verbose) cat("\n== Forward selection (CV-scored, concurvity-gated) ==\n")
  fw <- forward_select(data, search_resp, scr$keep, structural, select_gen,
                       folds, knots, concurvity_max = concurvity_max,
                       offset = offset, re_vars = re_vars, re_labels = re_labels,
                       verbose = verbose)

  # -- backward re-check ---------------------------------------------------
  if (verbose) cat("\n== Backward re-check (BIC) ==\n")
  kept <- backward_recheck(data, search_resp, fw$accepted, structural, select_gen,
                           knots, offset = offset, verbose = verbose)

  # -- final fit + adequacy ------------------------------------------------
  model <- final_fit(data, response, kept, structural, family_gen, knots,
                     offset = offset, squeeze = squeeze)

  # -- decision table ------------------------------------------------------
  edf <- {
    s <- summary(model)$s.table
    stats::setNames(s[, "edf"], gsub("^s\\(|\\)$", "", rownames(s)))
  }
  decision <- data.frame(
    covariate       = candidates,
    na_frac         = sapply(candidates, function(v)
                        if (v %in% names(data)) mean(is.na(data[[v]])) else NA),
    spearman_cluster= clusters[candidates],
    entry_step      = NA_integer_,
    concurvity      = NA_real_,
    status          = NA_character_,
    reason          = NA_character_,
    final_edf       = NA_real_,
    row.names = NULL, stringsAsFactors = FALSE
  )
  # screened-out
  decision$status[match(scr$dropped, decision$covariate)] <- "dropped"
  decision$reason[match(scr$dropped, decision$covariate)] <- scr$reason[scr$dropped]
  # forward path
  if (!is.null(fw$path)) {
    acc <- fw$path[fw$path$action == "accept", ]
    rej <- fw$path[fw$path$action == "reject-concurvity", ]
    mi <- match(acc$var, decision$covariate)
    decision$entry_step[mi] <- acc$step; decision$concurvity[mi] <- acc$concurvity
    mr <- match(rej$var, decision$covariate)
    decision$status[mr] <- "dropped"; decision$reason[mr] <- "concurvity"
    decision$concurvity[mr] <- rej$concurvity
  }
  # backward drops vs final kept
  fwd_acc <- if (!is.null(fw$path)) fw$path$var[fw$path$action == "accept"] else character()
  back_dropped <- setdiff(fwd_acc, kept)
  decision$status[match(back_dropped, decision$covariate)] <- "dropped"
  decision$reason[match(back_dropped, decision$covariate)] <- "backward-BIC"
  ki <- match(kept, decision$covariate)
  decision$status[ki] <- "kept"; decision$reason[ki] <- "selected"
  decision$final_edf[ki] <- edf[kept]
  # anything still unlabelled = screened-in but never entered
  decision$status[is.na(decision$status)] <- "dropped"
  decision$reason[is.na(decision$reason)] <- "no CV improvement"

  utils::write.csv(decision, file.path(outdir, "covariate_decision_table.csv"), row.names = FALSE)

  # -- diagnostics plots ---------------------------------------------------
  grDevices::png(file.path(outdir, "selection_diagnostics.png"), width = 1400, height = 1100, res = 130)
  tryCatch({ par(mfrow = c(2, 2)); gam.check(model) },
           finally = grDevices::dev.off())

  grDevices::png(file.path(outdir, "partial_effects.png"), width = 1400, height = 1100, res = 130)
  tryCatch(plot(model, pages = 1, scale = 0, all.terms = TRUE, shade = TRUE),
           finally = grDevices::dev.off())

  final_formula <- paste(deparse(formula(model)), collapse = " ")
  dev_expl <- summary(model)$dev.expl

  if (verbose) {
    cat("\n== Result ==\n")
    cat("kept:", if (length(kept)) paste(kept, collapse = ", ") else "(none)", "\n")
    cat(sprintf("deviance explained: %.1f%%\n", 100 * dev_expl))
    cat("formula:", final_formula, "\n")
    cat("outputs ->", normalizePath(outdir), "\n")

    cat("\n== Final model summary ==\n")
    print(summary(model))
    cat("\n== Concurvity (worst) of final model ==\n")
    cc <- tryCatch(round(concurvity(model, full = FALSE)$worst, 3), error = function(e) NULL)
    if (!is.null(cc)) print(cc) else cat("(not available)\n")
  }

  invisible(list(model = model, kept = kept, family = fs$chosen,
                 decision = decision, path = fw$path, formula = final_formula,
                 dev_expl = dev_expl, folds = folds))
}

# =============================================================================
# Example run (executes only under Rscript, not when sourced)
# =============================================================================

if (sys.nframe() == 0) {
  
  TARGETS <- c("ASUP", "LNAS")
  
  FEATURES_DEFAULT <- c(
    "npp", "mnkc_epi", "zeu", "mnkc_mumeso", "zooc", "mnkc_hmlmeso",
    "sst", "sst_std", "sst_fdist", "chl_fdist", "chl",
    "ekman_anom_lag7", "ekman_anom", "ekman_anom_lag3", "ekman_7d",
    "n_upwell_events_14d", "n_upwell_events_3d", "n_upwell_events_7d",
    "tp", "ekman_pumping", "ekman_anom_lag14", "tisr", "slhf", "ssrd",
    "adt", "sla", "adt_std", "sla_std", "gke", "mld", "fsle_max",
    "ac_normdist", "c_normdist", "moon_phase", "bathy", "bathy_std",
    "o2_0", "o2_100", "o2_500"
  )
  
  data_path <- "C:/Users/h2ugo/Documents/COSTA/longline/data/processed"
  
  for (TARGET in TARGETS) {

    out_path <- sprintf("C:/Users/h2ugo/Documents/COSTA/longline/data_analysis/GAM/output/%s", TARGET)

    df <- read.csv(file.path(data_path, "LL_extracted.csv"))
    
    # FOR BINOMIAL: presence/absence
    df[[TARGET]] <- ifelse(df[[TARGET]] > 1, 1, df[[TARGET]])

    # Mirror the pipeline preprocessing (see project memory / run_pipeline_spatialCV).
    df <- df[df$embarcacao != "Arquipelago" & df$Year >= 2015, ]
    df <- df[is.finite(df$nhooks) & df$nhooks > 0, ]     # offset needs positive effort
    names(df)[names(df) == "x_centroid"] <- "lon"
    names(df)[names(df) == "y_centroid"] <- "lat"
    df$month <- df$Month
    df$year  <- factor(df$Year)

    dir.create(out_path, showWarnings = FALSE, recursive = TRUE)
    sink(file.path(out_path, "covariate_selection_log.txt"), split = TRUE)
    tryCatch({   # finally: never leave the console sinked if the run errors

      res <- select_gam_covariates(
        data       = df,
        response   = TARGET,
        candidates = FEATURES_DEFAULT,
        structural = c("s(year, bs='re')", "te(lon, lat)"),
        offset     = "nhooks",          # model catch per effort (log-hooks offset)
      #  knots      = list(month = c(0.5, 12.5)),
        cv_k       = 5,
        cv_scheme  = "stratified",      # rare-event target (~4% prevalence): stratify folds
                                        #   on the 0/1 label instead of SPCV so no fold lacks
                                        #   positives. Drop this arg to use the SPCV default.
        outdir     = out_path
      )

      # Persist for later inspection: readRDS(...) then summary()/gam.check()/plot().
      saveRDS(res$model, file.path(out_path, "gam_model.rds"))  # the mgcv gam object
      saveRDS(res,       file.path(out_path, "selection_result.rds"))  # + kept set, decision table, path

    }, finally = sink())
  }
}