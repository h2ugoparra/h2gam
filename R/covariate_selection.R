# =============================================================================
# covariate_selection.R
#
# Automated GAM covariate validation & selection for mgcv, data-independent.
#
# Given a response and a pool of candidate covariates it:
#   1. profiles the response and auto-suggests a family (count / proportion /
#      positive / gaussian branches), ranked by AIC alone; dispersion and
#      residual uniformity are reported for inspection, not scored;
#   2. forward-selects a DECORRELATED, prediction-relevant covariate subset,
#      scored by cross-validated skill (spatial-block or stratified folds),
#      with a concurvity gate and a 1-SE stopping rule;
#   3. runs a backward elimination re-check (BIC), a final select=TRUE
#      shrinkage fit, and adequacy diagnostics.
#
# Every threshold is a logged, overridable argument; nothing is hardcoded to a
# particular dataset. Structural (spatial/temporal) terms are forced in and not
# screened -- covariates must earn their place OVER AND ABOVE them.
# =============================================================================

# ---------------------------------------------------------------------------
# 1. Response profiling & family cascade
# ---------------------------------------------------------------------------

profile_response <- function(y) {
  y <- y[is.finite(y)]
  list(
    integer    = all(abs(y - round(y)) < 1e-8),
    binary     = all(y %in% c(0, 1)),
    nonneg     = all(y >= 0),
    positive   = all(y > 0),
    in01       = all(y >= 0 & y <= 1),
    frac_zero  = mean(y == 0),
    frac_one   = mean(y == 1),
    skew       = mean((y - mean(y))^3) / stats::sd(y)^3
  )
}

# Candidate family generators for the profiled response. Each element is a
# zero-arg closure returning a FRESH mgcv family (fresh matters: nb()/tw()
# carry a mutable theta/p estimate in their environment). Candidates within a
# branch are ranked by AIC, so they must share a likelihood type on the same
# response: all probabilities, or all densities.
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
    # response and likelihood), and AIC alone is what picks among them.
    add("binomial_logit",   function() binomial(link = "logit"))
    add("binomial_cloglog", function() binomial(link = "cloglog"))
    add("binomial_probit",  function() binomial(link = "probit"))

  } else if (p$integer && p$nonneg) {
    # ---- count branch ----
    # Tweedie is deliberately not a candidate: its likelihood is a density,
    # Poisson/nb give probabilities, and on low counts the density at integer
    # values is arbitrarily large -- it "wins" on AIC by hundreds of units on
    # truly Poisson/nb data. Pass family = tw explicitly to use it.
    add("poisson", function() poisson())
    add("nb",      function() nb())

  } else if (p$in01 && p$frac_zero == 0 && p$frac_one == 0) {
    # ---- proportion, open interval ----
    add("betar",    function() betar())
    add("gaussian", function() gaussian())

  } else if (p$in01) {
    # ---- proportion with 0/1 boundary (caller applies SV squeeze for betar) ----
    # The only candidate: the squeezed-beta AIC is set by the squeeze constant
    # (the density is unbounded at the boundary), so it cannot be ranked
    # against a model of the raw response such as Tweedie.
    add("betar_sv", function() betar())

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
  # suppressMessages: DHARMa message()s an "install mgcViz" notice on every gam
  # simulation when mgcViz is absent, which would break verbose = FALSE.
  out <- tryCatch(suppressMessages({
    sim <- DHARMa::simulateResiduals(m, n = 250, plot = FALSE)
    DHARMa::testUniformity(sim, plot = FALSE)$p.value
  }), error = function(e) NA_real_)
  out
}

# Fit each candidate family on the structural-only model, rank by AIC (ok
# families first), and flag dispersion / residual uniformity. `ok` means the fit
# converged without error -- it is not a diagnostic pass, and the dispersion /
# resid_ks_p columns do not enter the ranking. Returns the chosen family
# generator plus a diagnostic table.
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

  list(table = tab, chosen = chosen,
       family_gen = cands[[chosen]], squeeze = identical(chosen, "betar_sv"),
       fit = fits[[chosen]])
}

# Generator for a user-supplied family OBJECT. Unlike a constructor, the object
# is shared by every fit, and extended families (nb(), tw(), betar(), ...)
# store their estimated theta in its environment, so each fit would start from
# the previous fit's estimate. Reset theta to the supplied value before each use.
family_object_gen <- function(family) {
  if (!inherits(family, "extended.family") || !is.function(family$putTheta))
    return(function() family)
  th0 <- family$getTheta()
  function() { family$putTheta(th0); family }
}

# Family generator with the shape parameter (negbin theta / Tweedie p / beta
# precision / scaled-t df+scale) held FIXED at its structural-only estimate
# (`fit`). Two reasons: held-out deviance is only comparable between models
# that share the shape value -- left free, a real driver raises nb theta and
# with it the held-out deviance, so the search rejects it -- and re-estimating
# it on every CV fit is the dominant cost of the search. The final fit still
# uses the free (re-estimating) family. mgcv constructors take a positive
# `theta` as known; the fitted link is carried over.
make_select_gen <- function(free_gen, fit) {
  fam <- fit$family
  if (!inherits(fam, "extended.family") || !isTRUE(fam$n.theta > 0)) return(free_gen)
  ctor <- list("Negative Binomial" = nb, "Tweedie" = tw, "Beta regression" = betar,
               "Scaled t" = scat)[[sub("\\(.*", "", fam$family)]]
  th <- tryCatch(fam$getTheta(TRUE), error = function(e) NA_real_)
  if (is.null(ctor) || !all(is.finite(th)) || any(th <= 0)) {
    warning(sprintf(paste0("family '%s': its shape parameter cannot be held fixed during the ",
                           "search, so it is re-estimated on every CV fit and CV gains are ",
                           "not strictly comparable between models."), fam$family),
            call. = FALSE)
    return(free_gen)
  }
  # do.call passes the link by value: the constructors read it with
  # substitute(), and newer mgcv rejects a variable name there.
  args <- list(theta = th, link = fam$link)
  function() do.call(ctor, args)
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

# Stop with a clear message when named columns are missing from `data`.
.need_cols <- function(data, cols, what) {
  miss <- setdiff(cols, names(data))
  if (length(miss))
    stop(sprintf("%s column(s) not found in data: %s.", what, paste(miss, collapse = ", ")),
         call. = FALSE)
}

# Snapshot the caller's RNG state and return a function that restores it, so
# the fold builders can set.seed() for reproducibility without overwriting the
# user's random-number stream. Call it on entry, run the result on exit:
#   restore_rng <- .rng_restorer(); on.exit(restore_rng(), add = TRUE)
.rng_restorer <- function() {
  env <- globalenv()
  if (exists(".Random.seed", envir = env, inherits = FALSE)) {
    old <- get(".Random.seed", envir = env, inherits = FALSE)
    function() assign(".Random.seed", old, envir = env)
  } else {
    function() if (exists(".Random.seed", envir = env, inherits = FALSE))
      rm(".Random.seed", envir = env)
  }
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

#' Spatial-block cross-validation folds (SPCV)
#'
#' Two-stage spatial CV (Wang et al. 2023), a faithful R port of the h2ml
#' `SPCVSplitter`. Stage 1 builds spatially coherent blocks by agglomerative
#' hierarchical clustering on the coordinates, cut at a distance threshold.
#' Stage 2 groups the blocks into `k` folds by a cluster ensemble (HBGF):
#' per-block means of location, covariates, and response are each k-means
#' clustered, the one-hot memberships are stacked into a co-occurrence
#' affinity, and spectral clustering assigns folds. Folds are geographically
#' separated and representative of covariate and label space.
#'
#' Stage 1 builds the full pairwise distance matrix, so memory and time grow
#' with n^2 (roughly 8 * n^2 bytes for haversine: about 20 GB at n = 50,000).
#' For large data, thin or aggregate the rows before building folds.
#'
#' @param data A data.frame with the coordinate, response, and covariate columns.
#' @param response Name of the response column (block-mean view for stage 2).
#' @param candidates Character vector of covariate names (numeric ones feed the
#'   stage-2 covariate view; others are ignored).
#' @param coords Length-2 character vector naming the two coordinate columns
#'   (x/longitude then y/latitude). Default `c("lon", "lat")`. For
#'   `metric = "haversine"` these must be lon/lat degrees; for
#'   `metric = "euclidean"` with `crs = NULL` they are taken as already-planar
#'   coordinates.
#' @param k Number of folds.
#' @param metric `"haversine"` (great-circle distance on lon/lat degrees; no
#'   projection, valid anywhere on the globe) or `"euclidean"` (planar
#'   distance; see `crs`). Euclidean uses the C-optimised [stats::dist()] and
#'   scales better to large n than the haversine matrix.
#' @param threshold AHC cut height: degrees for haversine, projected units
#'   (e.g. metres) for euclidean. `NULL` (default) uses the 10th percentile of
#'   pairwise distances.
#' @param crs Only used when `metric = "euclidean"`. `NULL` (default) treats
#'   the `coords` columns as already projected. Otherwise an sf/PROJ CRS (a
#'   proj string or EPSG code) used to reproject lon/lat to a planar system
#'   before computing distances (needs the sf package). No projection is
#'   assumed by default, so this function makes no study-area assumption.
#' @param linkage AHC linkage method, passed to [stats::hclust()]. `"ward"`
#'   is taken as `"ward.D2"` (the Ward criterion on unsquared distances,
#'   matching scikit-learn). Ward variants are downgraded to `"average"` for
#'   haversine, where they are invalid.
#' @param pca_var Fraction of variance retained when the stage-2 covariate
#'   view is compressed by PCA.
#' @param seed Random seed for the k-means / spectral steps.
#' @param verbose Print a one-line fold summary.
#'
#' @return An integer vector of fold ids (1..k), one per row of `data`.
#'
#' @seealso [make_stratified_folds()] for the non-spatial alternative;
#'   [select_gam_covariates()], which builds these folds internally.
#' @export
make_spatial_folds <- function(data, response, candidates, coords = c("lon", "lat"),
                               k = 10, metric = "haversine", threshold = NULL,
                               crs = NULL, linkage = "average", pca_var = 0.95,
                               seed = 42, verbose = TRUE) {
  restore_rng <- .rng_restorer(); on.exit(restore_rng(), add = TRUE)
  if (length(coords) != 2) stop("`coords` must name exactly two columns.")
  .need_cols(data, coords, "coords"); .need_cols(data, response, "response")
  c1 <- data[[coords[1]]]; c2 <- data[[coords[2]]]           # coord 1 (x/lon), coord 2 (y/lat)
  if (!all(is.finite(c1) & is.finite(c2) & is.finite(data[[response]])))
    stop("make_spatial_folds: missing or non-finite values in the coordinates or response; drop those rows first.")
  # hclust's "ward" is the legacy alias for "ward.D"; scikit-learn's Ward (which
  # SPCV mirrors) is "ward.D2".
  if (identical(linkage, "ward")) linkage <- "ward.D2"

  # -- Stage 1: AHC -> spatially coherent blocks --------------------------
  if (metric == "haversine") {
    dmat <- .haversine_dist(c2, c1); d <- stats::as.dist(dmat)  # (lat, lon) degrees
    link <- if (startsWith(linkage, "ward")) "average" else linkage  # ward invalid for haversine
    thr  <- if (is.null(threshold)) stats::quantile(d, 0.10) else threshold * pi / 180
    loc  <- cbind(c2, c1)                                    # stage-2 location view (degrees)
  } else {
    # euclidean distance on PLANAR coordinates. This function assumes no
    # particular study area: with crs = NULL the given coords are taken as
    # already projected (project upstream, in whatever CRS fits your region);
    # pass a crs (an sf/PROJ string or EPSG code) to reproject lon/lat -> that
    # CRS here. Euclidean distance on raw lon/lat degrees is only valid near
    # the equator, so no default projection is baked in.
    if (is.null(crs)) {
      loc <- cbind(c1, c2)                                   # already planar (metres/units)
    } else {
      if (!requireNamespace("sf", quietly = TRUE))
        stop("make_spatial_folds(crs=...) needs the sf package to reproject lon/lat.")
      loc <- sf::sf_project(from = "+proj=longlat +datum=WGS84", to = crs, pts = cbind(c1, c2))
    }
    d    <- stats::dist(loc, method = "euclidean")           # stage-2 location view (planar)
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
  if (!all(is.finite(Xmat)))
    stop("make_spatial_folds: missing or non-finite values in the numeric candidates; drop those rows first.")
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
  # k-means labels for one view. A view with <= k distinct rows (e.g. a rare
  # response where most blocks have mean 0) can't take k centres: each distinct
  # row is its own cluster. With > k distinct but heavily duplicated rows, a
  # random start can still pick duplicate centres; fall back to clustering the
  # distinct rows and mapping back.
  km <- function(m) {
    key <- do.call(paste, as.data.frame(m))
    if (length(unique(key)) <= k) return(match(key, unique(key)))
    tryCatch(stats::kmeans(m, centers = k, nstart = 10, iter.max = 100)$cluster,
             error = function(e) {
               u <- !duplicated(key)
               stats::kmeans(m[u, , drop = FALSE], centers = k, nstart = 10,
                             iter.max = 100)$cluster[match(key, key[u])]
             })
  }
  eye   <- diag(k)
  views <- list(block_coords, block_Xp, block_y)
  views <- views[vapply(views, ncol, integer(1)) > 0]         # no numeric candidates -> no covariate view
  H     <- do.call(cbind, lapply(views, function(m) eye[km(m), , drop = FALSE]))  # (B, views*k)
  S   <- H %*% t(H)
  S   <- (S + 1 / (B + 1)); S <- S / max(S)                  # floor for connectivity, normalise

  folds <- .spectral_cluster(S, k, seed)[blocks]
  if (verbose) cat(sprintf("  SPCV (%s): %d AHC blocks -> %d folds (sizes %s)\n",
                           metric, B, k, paste(tabulate(folds, k), collapse = "/")))
  as.integer(folds)
}

#' Stratified k-fold cross-validation folds
#'
#' A non-spatial alternative to [make_spatial_folds()] for rare-event
#' responses. SPCV's spatially coherent blocks can leave a fold with zero
#' positives when prevalence is low, making held-out binomial deviance
#' degenerate. Stratifying on the response spreads the rare positives evenly,
#' so every fold is evaluable. Strata are the exact label values for a
#' few-valued response (binary / integer label), else quantile bins so a
#' continuous response stays balanced too.
#'
#' Note this deliberately drops spatial separation -- spatial autocorrelation
#' can leak between train and test, so CV skill is optimistic relative to
#' SPCV; the tradeoff buys stable fold estimates for a rare event.
#'
#' @param data A data.frame containing the response column.
#' @param response Name of the response column to stratify on.
#' @param k Number of folds.
#' @param seed Random seed for the within-stratum shuffle.
#' @param verbose Print a one-line fold summary (with per-fold positive counts
#'   for a binary response).
#'
#' @return An integer vector of fold ids (1..k), one per row of `data`.
#'
#' @export
make_stratified_folds <- function(data, response, k = 5, seed = 1, verbose = TRUE) {
  restore_rng <- .rng_restorer(); on.exit(restore_rng(), add = TRUE)
  .need_cols(data, response, "response")
  y <- data[[response]]
  if (anyNA(y)) stop("make_stratified_folds: the response has missing values; drop those rows first.")
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
    if (!isTRUE(stats::sd(x[is.finite(x)]) >= sqrt(.Machine$double.eps))) {
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
  m <- tryCatch(withCallingHandlers(
                  bam(form, data = data, family = family_gen(),
                      method = "fREML", discrete = TRUE, knots = knots),
                  warning = function(w)   # e.g. structural-only model with no smooths
                    if (grepl("no smooths, ignoring", conditionMessage(w)))
                      invokeRestart("muffleWarning")),
                error = function(e) NULL)
  if (is.null(m))
    m <- tryCatch(gam(form, data = data, family = family_gen(),
                      method = "REML", knots = knots), error = function(e) NULL)
  m
}

# Concurvity of a newly added smooth against the already-accepted covariate
# smooths only (structural terms deliberately excluded). Returns the worst
# pairwise value, 0 when nothing to compare against, and NA when it cannot be
# computed -- the caller must not read that as "no concurvity".
concurvity_new <- function(model, newvar, accepted) {
  if (length(accepted) == 0) return(0)
  cc <- tryCatch(concurvity(model, full = FALSE), error = function(e) NULL)
  if (is.null(cc)) return(NA_real_)
  W    <- cc$worst
  rlab <- paste0("s(", newvar, ")")
  clab <- paste0("s(", accepted, ")")
  if (!rlab %in% rownames(W) || !all(clab %in% colnames(W))) return(NA_real_)
  max(W[rlab, clab, drop = FALSE])
}

# A random-effect smooth or smooth spec: bs = 're', bs = 'fs', or a tensor
# with such a margin.
is_re_smooth <- function(s) {
  is.list(s) &&
    (inherits(s, c("random.effect", "fs.interaction", "re.smooth.spec", "fs.smooth.spec")) ||
       any(vapply(s$margin, is_re_smooth, logical(1))))
}

forward_select <- function(data, response, candidates, structural, family_gen,
                           folds, knots, concurvity_max = 0.8, offset = NULL,
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
      # (matches the deployment predict(..., exclude=...)). They are read off
      # the fit, so every form is covered (random slopes, by=, bs='fs', tensor
      # margins). Levels unseen in training are remapped to a training level so
      # the design matrix builds; the excluded term zeroes it.
      re <- Filter(is_re_smooth, m$smooth)
      for (v in unique(unlist(lapply(re, `[[`, "term")))) if (is.factor(tr[[v]])) {
        unseen <- !te[[v]] %in% tr[[v]]
        te[[v]][unseen] <- tr[[v]][1]
      }
      mu <- tryCatch(predict(m, newdata = te, type = "response",
                             exclude = vapply(re, `[[`, "", "label")),
                     error = function(e) NULL)
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
  # A fold is skipped (NA) when its fit, prediction, or deviance fails. Without
  # a scorable baseline every candidate would just look like "no CV improvement",
  # so fail loudly instead.
  if (sum(is.finite(cur_dev)) < 2)
    stop(sprintf(paste0("forward selection: the structural-only model could be scored on ",
                        "only %d of %d CV folds (fit/prediction failed or non-finite ",
                        "deviance); check for NA/invalid values in the model variables."),
                 sum(is.finite(cur_dev)), length(cur_dev)))
  path      <- list()
  step      <- 0
  low_folds <- character()                 # candidates scored on < half the folds
  last      <- list()                      # each candidate's most recent CV score
  ungated   <- character()                 # candidates whose concurvity check failed

  repeat {
    step <- step + 1
    scored <- lapply(remaining, function(v) {
      cd   <- fold_dev(c(accepted, v))
      if (sum(is.finite(cd)) < length(cd) / 2) low_folds <<- union(low_folds, v)
      diff <- cur_dev - cd                 # paired per-fold improvement
      list(v = v, cd = cd,
           mean = mean(diff, na.rm = TRUE),
           se   = stats::sd(diff, na.rm = TRUE) / sqrt(sum(is.finite(diff))))
    })
    for (s in scored) last[[s$v]] <- data.frame(var = s$v, cv_gain = s$mean, cv_se = s$se)
    ord <- order(vapply(scored, function(s) -s$mean, numeric(1)))

    picked <- NULL
    for (j in ord) {
      s <- scored[[j]]
      if (!is.finite(s$mean) || s$mean <= 0) break     # ordered desc: nothing better remains
      if (!is.finite(s$se) || s$mean < s$se) next      # NA se (<2 scored folds) or fails 1-SE bar
      m  <- fit_full(c(accepted, s$v))
      cg <- if (is.null(m)) NA_real_ else concurvity_new(m, s$v, accepted)
      # The gate fails closed: a candidate whose concurvity cannot be computed
      # (full-data fit or concurvity() failed) is rejected, not waved through.
      if (!is.finite(cg)) ungated <- union(ungated, s$v)
      if (!is.finite(cg) || cg > concurvity_max) {
        path[[length(path) + 1]] <- data.frame(step = step, var = s$v, action = "reject-concurvity",
                                                concurvity = cg, cv_gain = s$mean, cv_se = s$se)
        # Out for good: worst-case concurvity depends only on the two smooths'
        # bases, so against a growing accepted set it can only rise. Leaving
        # it in would re-fit it on every fold at every later step.
        remaining <- setdiff(remaining, s$v)
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

  if (length(low_folds))
    warning(sprintf("forward selection: %s could be scored on fewer than half the CV folds (fit/prediction failures); their CV gain is unreliable.",
                    paste(low_folds, collapse = ", ")), call. = FALSE)

  if (length(ungated))
    warning(sprintf(paste0("forward selection: concurvity could not be computed for %s (the ",
                           "full-data fit or concurvity() failed); rejected, since the ",
                           "concurvity gate could not be applied."),
                    paste(ungated, collapse = ", ")), call. = FALSE)

  list(accepted = accepted,
       path = if (length(path)) do.call(rbind, path) else NULL,
       scores = if (length(last)) do.call(rbind, last) else NULL)
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

#' Automated GAM covariate validation and selection for mgcv
#'
#' Given a response and a pool of candidate covariates, this:
#' 1. profiles the response and auto-suggests a family (count / proportion /
#'    positive / gaussian branches), ranked by AIC alone;
#' 2. forward-selects a decorrelated, prediction-relevant covariate subset,
#'    scored by cross-validated held-out deviance (spatial-block or
#'    response-stratified folds), with a concurvity gate and a 1-SE stopping
#'    rule;
#' 3. runs a backward elimination re-check (BIC), a final `select = TRUE`
#'    shrinkage fit, and adequacy diagnostics.
#'
#' Every threshold is a logged, overridable argument; nothing is hardcoded to
#' a particular dataset. Structural (spatial/temporal) terms are forced in and
#' not screened -- covariates must earn their place over and above them.
#' Random-effect structural terms are excluded from held-out CV prediction,
#' matching a deployment `predict(..., exclude = ...)`. That covers every
#' `bs = "re"` form (random intercepts, random slopes such as
#' `s(vessel, sst, bs = "re")`, `by =` variables), factor-smooth interactions
#' (`bs = "fs"`), and tensor terms with an `"re"` margin. Each needs a factor
#' grouping variable.
#'
#' Search fits use `mgcv::bam(discrete = TRUE)` with the family shape
#' parameter (negbin theta, Tweedie p, beta precision, scaled-t df and scale)
#' fixed at its structural-only estimate, whether the family was auto-suggested
#' or supplied. This keeps held-out deviance comparable between candidate
#' models and is much faster. Other families with an estimated shape parameter
#' (e.g. `ocat()`, `ziP()`) cannot be fixed this way and give a warning. The
#' final fit is a free `mgcv::gam(method = "REML")`.
#'
#' @param data A data.frame with the response, coordinates, candidate
#'   covariates, and any variables used by `structural` / `offset`.
#' @param response Name of the response column.
#' @param candidates Character vector of candidate covariate names. Each enters
#'   the model as a `s(<name>)` smooth.
#' @param structural Character vector of mgcv terms forced into every model
#'   and never screened, e.g. `c("s(year, bs='re')", "te(lon, lat)")`. May be
#'   `character()` for an intercept-only baseline.
#' @param coords Length-2 character vector naming the longitude and latitude
#'   columns, used for spatial folds and complete-case filtering.
#' @param family `"auto"` (default) profiles the response and takes the
#'   lowest-AIC candidate that fits on the structural-only model. Candidates
#'   are limited to families whose AIC is comparable on the response: Poisson
#'   and negative binomial for counts, squeezed beta alone for proportions
#'   with exact 0/1. Tweedie is not auto-suggested for either (its density
#'   cannot be ranked against count probabilities or a squeezed response);
#'   pass `family = tw` to use it. The reported
#'   `dispersion` and `resid_ks_p` columns are diagnostics to inspect, not
#'   inputs to the choice. Both are measured before covariate selection, so a
#'   low `resid_ks_p` reflects not-yet-modelled covariates and structure as much
#'   as the family -- re-check it on the final fit before acting on it.
#'   Alternatively a family object (e.g. `nb()`) or a zero-argument family
#'   generator (e.g. `nb`) to skip the suggestion step. For a family object
#'   with an estimated shape parameter, theta is reset to its supplied value
#'   before every fit, and the object is modified in place by the fits; pass
#'   the constructor instead to leave your object untouched.
#' @param concurvity_max Reject a candidate whose worst pairwise concurvity
#'   against the already-accepted covariate smooths exceeds this (structural
#'   terms excluded from the gate). A candidate whose concurvity cannot be
#'   computed is rejected with a warning.
#' @param spearman_pre `|Spearman rho|` threshold for the correlation clusters
#'   reported in the decision table (annotation only; the concurvity gate does
#'   the actual removal).
#' @param na_max Drop candidates whose NA/non-finite fraction exceeds this in
#'   the univariate screen. The search runs on the rows complete for every
#'   candidate that passes, so each CV model sees the same data; a sparse
#'   candidate therefore costs the search rows even if it is not selected. The
#'   final model is refitted on all rows complete for the kept covariates.
#' @param cv_k Number of CV folds.
#' @param cv_scheme Fold construction: `"spatial"` (SPCV, the default; see
#'   [make_spatial_folds()]) or `"stratified"` (non-spatial k-fold stratified
#'   on the response; see [make_stratified_folds()]). Use `"stratified"` for
#'   low-prevalence (e.g. binomial presence/absence) targets where spatial
#'   blocks would leave folds with too few positives.
#' @param cv_metric SPCV distance: `"haversine"` (the default; great-circle on
#'   lon/lat degrees, correct anywhere with no projection) or `"euclidean"`
#'   (planar; faster on large n, but see `cv_crs`).
#' @param cv_crs Only used when `cv_metric = "euclidean"`. `NULL` (default)
#'   would compute euclidean distance on raw lon/lat degrees, which is only
#'   valid near the equator -- so for euclidean pass an sf/PROJ CRS (proj
#'   string or EPSG code) suited to your study area to reproject first (needs
#'   the sf package). No study-area projection is assumed; `"haversine"` avoids
#'   the question entirely.
#' @param cv_threshold SPCV AHC block threshold (`NULL` = 10th percentile of
#'   pairwise distances).
#' @param cv_linkage SPCV AHC linkage (`"ward.D2"` default; Ward variants are
#'   downgraded to `"average"` for haversine). See [make_spatial_folds()].
#' @param knots mgcv per-smooth knot positions, passed to every `gam()`/`bam()`
#'   fit. Mainly the cyclic boundary for `s(month, bs='cc')`: pass
#'   `list(month = c(0.5, 12.5))` so the 12-month cycle wraps (Dec to Jan)
#'   instead of collapsing months 1 and 12 onto the same knot. `NULL` lets
#'   mgcv place knots automatically; entries for absent variables are ignored.
#' @param offset Name of an effort column; adds `offset(log(<offset>))` to
#'   every fit (including family suggestion, so AIC stays comparable).
#' @param outdir Directory for the decision-table CSV, diagnostic PNGs, and the
#'   run log (created if absent). Defaults to the session's temporary
#'   directory, which R deletes when the session ends: pass a path to keep
#'   the files.
#' @param logfile Filename (within `outdir`) for a run log that captures all
#'   the console output, plus warnings, messages and the error that stopped a
#'   failed run (each distinct text once, where it first occurred); `NULL`
#'   disables it. Only written when `verbose = TRUE` (there is nothing to
#'   capture otherwise). The sink is unwound on exit, so an error mid-run
#'   still restores the console.
#' @param seed Random seed for fold construction.
#' @param verbose Print progress, the family table, the selection path, and
#'   the final model summary.
#'
#' @return (Invisibly) a list: `model` (the final mgcv gam), `kept` (selected
#'   covariate names), `family` (chosen family label), `decision` (per-candidate
#'   decision table, also written to `covariate_decision_table.csv`; see
#'   below), `path`
#'   (forward-selection steps), `formula` (final formula as a string),
#'   `dev_expl` (deviance explained), `folds` (integer fold vector, one per row used in the search), and
#'   `logfile` (path to the run log, or `NULL` if none was written).
#'
#'   The `decision` table has one row per candidate. `cv_gain` and `cv_se`
#'   are the mean paired per-fold drop in held-out deviance from adding the
#'   candidate, and its standard error, at the step that decided its fate:
#'   the step it was accepted or rejected for concurvity, otherwise the last
#'   step it was scored. `reason` is one of `"selected"`; `"concurvity"`;
#'   `"concurvity check failed"` (concurvity could not be computed, so the
#'   candidate was rejected; a warning names it);
#'   `"backward-BIC"`; `"no CV improvement"` (gain <= 0); `"CV gain below 1
#'   SE"` (positive gain smaller than its standard error); `"CV scoring
#'   failed"` (the fits or held-out predictions failed on too many folds); or
#'   the univariate-screen reason (`"absent"`, `"non-numeric"`, `"NA <n>%"`,
#'   `"zero variance"`), for which the CV columns are `NA`.
#'
#'   Side effects in `outdir`: `covariate_decision_table.csv`,
#'   `selection_diagnostics.png` (gam.check panels), `partial_effects.png`,
#'   and (unless disabled) the run log named by `logfile`.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n  <- 400
#' df <- data.frame(lon = runif(n, -10, 0), lat = runif(n, 35, 45),
#'                  x1 = rnorm(n), x2 = rnorm(n))
#' df$y <- sin(df$x1) + 0.1 * df$lon + rnorm(n, sd = 0.3)
#' res <- select_gam_covariates(
#'   data       = df,
#'   response   = "y",
#'   candidates = c("x1", "x2"),
#'   structural = "te(lon, lat)",
#'   cv_k       = 3,
#'   cv_scheme  = "stratified",
#'   outdir     = tempdir()
#' )
#' res$kept
#' }
#'
#' @export
select_gam_covariates <- function(
    data, response, candidates,
    structural     = c("s(month, bs='cc')", "s(year, bs='re')", "te(lon, lat)"),
    coords         = c("lon", "lat"),
    family         = "auto",
    concurvity_max = 0.8,
    spearman_pre   = 0.7,
    na_max         = 0.30,
    cv_k           = 5,
    cv_scheme      = "spatial",
    cv_metric      = "haversine",
    cv_crs         = NULL,
    cv_threshold   = NULL,
    cv_linkage     = "ward.D2",
    knots          = NULL,
    offset         = NULL,
    outdir         = tempdir(),
    logfile        = "covariate_selection_log.txt",
    seed           = 1,
    verbose        = TRUE) {

  run <- function() select_covariates_impl(
    data, response, candidates, structural, coords, family, concurvity_max,
    spearman_pre, na_max, cv_k, cv_scheme, cv_metric, cv_crs, cv_threshold,
    cv_linkage, knots, offset, outdir, logfile, seed, verbose)
  # All progress output is verbose-gated: nothing to log unless verbose.
  if (is.null(logfile) || !isTRUE(verbose)) return(invisible(run()))

  # Run log: tee console output to outdir/logfile (split = TRUE keeps the
  # console echo). sink() only sees stdout, so warnings, messages and a
  # terminating error are written to the log by calling handlers, each
  # distinct text once. on.exit restores the console even if the run errors.
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  path <- file.path(outdir, logfile)
  file.create(path)                        # start empty: `con` appends
  con <- file(path, "at")
  sink(con, split = TRUE)
  on.exit({ sink(); close(con) }, add = TRUE)
  seen <- character()
  note <- function(prefix) function(cond) {
    txt <- paste0(prefix, trimws(conditionMessage(cond), "right"))
    if (txt %in% seen) return()
    seen <<- c(seen, txt)
    # Append to the file by path: writing to `con` itself would be echoed to
    # the console by the split sink, on top of R's own report of the
    # condition. `con` is in append mode, so its later writes land after this.
    flush(con)
    cat(txt, "
", sep = "", file = path, append = TRUE)
  }
  invisible(withCallingHandlers(run(), warning = note("Warning: "),
                                message = note(""), error = note("Error: ")))
}

select_covariates_impl <- function(
    data, response, candidates, structural, coords, family, concurvity_max,
    spearman_pre, na_max, cv_k, cv_scheme, cv_metric, cv_crs, cv_threshold,
    cv_linkage, knots, offset, outdir, logfile, seed, verbose) {

  if (!is.data.frame(data)) stop("`data` must be a data.frame.")
  .need_cols(data, response, "response")
  cv_scheme <- match.arg(cv_scheme, c("spatial", "stratified"))
  if (cv_scheme == "spatial") {            # fail before the family fits, not after
    if (length(coords) != 2) stop("`coords` must name exactly two columns.")
    .need_cols(data, coords, "coords")
  }
  if (response %in% candidates)
    stop(sprintf("the response '%s' is also listed in `candidates`.", response))
  if (anyDuplicated(candidates)) {
    warning(sprintf("duplicate candidates dropped: %s.",
                    paste(unique(candidates[duplicated(candidates)]), collapse = ", ")),
            call. = FALSE)
    candidates <- unique(candidates)
  }
  if (!length(structural)) structural <- "1"   # nothing forced in: intercept-only baseline
  if (!is.null(offset)) {
    if (!offset %in% names(data)) stop(sprintf("offset column '%s' not found in data.", offset))
    if (any(data[[offset]] <= 0, na.rm = TRUE))
      stop(sprintf("offset column '%s' must be > 0 (it enters as log(%s)).", offset, offset))
  }

  # Random-effect structural terms are excluded from held-out CV prediction
  # (see forward_select). They need a grouping factor: on numeric variables
  # alone, s(x, bs='re') is a penalised slope, not random intercepts.
  for (term in structural) {
    spec <- tryCatch(eval(parse(text = term)[[1]], asNamespace("mgcv")),
                     error = function(e) NULL)
    if (!is_re_smooth(spec)) next
    v <- spec$term
    if (!all(v %in% names(data)))
      stop(sprintf("random-effect term %s: column '%s' not found in data.",
                   term, setdiff(v, names(data))[1]))
    if (!any(vapply(data[v], is.factor, logical(1))))
      stop(sprintf(paste0("random-effect term %s needs %s to be a factor ",
                          "(numeric gives a penalised slope, not random intercepts); ",
                          "convert it first, e.g. data$%s <- factor(data$%s)."),
                   term, paste0("'", v, "'", collapse = " or "), v[1], v[1]))
  }
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  do_log <- !is.null(logfile) && isTRUE(verbose)   # the caller holds the sink

  # -- family --------------------------------------------------------------
  if (identical(family, "auto")) {
    fs <- suggest_family(data, response, structural, knots, offset = offset)
    family_gen <- fs$family_gen; squeeze <- fs$squeeze
    if (verbose) { cat("\n== Family suggestion ==\n"); print(fs$table)
      cat(sprintf("-> chosen: %s%s\n", fs$chosen, if (squeeze) " (SV squeeze)" else "")) }
  } else {
    family_gen <- if (is.function(family)) function() family() else family_object_gen(family)
    squeeze <- FALSE; fs <- list(chosen = "user", table = NULL, fit = NULL)
    # A free shape parameter needs a structural-only fit to fix it for the search.
    if (isTRUE(family_gen()$n.theta > 0)) {
      off_term <- if (!is.null(offset)) sprintf("offset(log(%s))", offset) else NULL
      form0 <- stats::as.formula(paste(response, "~", paste(c(structural, off_term), collapse = " + ")))
      fs$fit <- tryCatch(gam(form0, data = data, family = family_gen(), method = "REML",
                             knots = knots), error = function(e) NULL)
    }
  }
  # fixed-shape family for the search; free family for the final fit
  select_gen <- make_select_gen(family_gen, fs$fit)

  # -- univariate screen ---------------------------------------------------
  scr <- univariate_screen(data, candidates, na_max = na_max)
  if (verbose) cat(sprintf("\n== Univariate screen ==\n  kept %d, dropped %d\n",
                           length(scr$keep), length(scr$dropped)))

  # NA fraction on the data as supplied (before the complete-case filter below,
  # which would otherwise make it 0 for every screened-in candidate)
  na_frac <- vapply(candidates, function(v)
    if (v %in% names(data)) mean(is.na(data[[v]])) else NA_real_, numeric(1))

  # complete-case on the comparison columns so every CV model sees the same rows.
  # The offset is included: an NA offset gives NA held-out predictions, which
  # would silently void every CV fold.
  # The supplied data is kept: the final fit only needs the kept covariates, so
  # it gets back the rows that a sparse, unselected candidate cost the search.
  usable <- function(vars) {
    vars <- intersect(vars, names(full_data))
    finite_ok <- vapply(full_data[vars], function(x) !is.numeric(x) | is.finite(x),
                        logical(nrow(full_data)))            # complete.cases keeps Inf
    stats::complete.cases(full_data[, vars]) & rowSums(!finite_ok) == 0
  }
  full_data  <- data
  model_vars <- c(response, offset,
                  all.vars(stats::as.formula(paste("~", paste(structural, collapse = "+")))))
  data <- full_data[usable(c(model_vars, coords, scr$keep)), ]
  if (verbose) cat(sprintf("  search rows: %d of %d (complete for every screened-in candidate)\n",
                           nrow(data), nrow(full_data)))

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
                       metric = cv_metric, threshold = cv_threshold, crs = cv_crs,
                       linkage = cv_linkage, seed = seed, verbose = verbose)
  }

  # -- forward selection ---------------------------------------------------
  if (verbose) cat("\n== Forward selection (CV-scored, concurvity-gated) ==\n")
  fw <- forward_select(data, search_resp, scr$keep, structural, select_gen,
                       folds, knots, concurvity_max = concurvity_max,
                       offset = offset,
                       verbose = verbose)

  # -- backward re-check ---------------------------------------------------
  if (verbose) cat("\n== Backward re-check (BIC) ==\n")
  kept <- backward_recheck(data, search_resp, fw$accepted, structural, select_gen,
                           knots, offset = offset, verbose = verbose)

  # -- final fit + adequacy ------------------------------------------------
  final_data <- full_data[usable(c(model_vars, kept)), ]
  if (verbose && nrow(final_data) != nrow(data))
    cat(sprintf("\n== Final fit ==\n  %d rows (complete for the kept covariates; the search used %d)\n",
                nrow(final_data), nrow(data)))
  model <- final_fit(final_data, response, kept, structural, family_gen, knots,
                     offset = offset, squeeze = squeeze)

  # -- decision table ------------------------------------------------------
  edf <- {
    s <- summary(model)$s.table                     # NULL when the model has no smooths
    if (is.null(s)) numeric() else stats::setNames(s[, "edf"], gsub("^s\\(|\\)$", "", rownames(s)))
  }
  decision <- data.frame(
    covariate       = candidates,
    na_frac         = na_frac,
    spearman_cluster= clusters[candidates],
    entry_step      = NA_integer_,
    concurvity      = NA_real_,
    cv_gain         = NA_real_,
    cv_se           = NA_real_,
    status          = NA_character_,
    reason          = NA_character_,
    final_edf       = NA_real_,
    row.names = NULL, stringsAsFactors = FALSE
  )
  # screened-out
  decision$status[match(scr$dropped, decision$covariate)] <- "dropped"
  decision$reason[match(scr$dropped, decision$covariate)] <- scr$reason[scr$dropped]
  # CV score behind each decision: the last time the candidate was scored,
  # overridden below by the step at which it was accepted or concurvity-rejected
  if (!is.null(fw$scores)) {
    ms <- match(fw$scores$var, decision$covariate)
    decision$cv_gain[ms] <- fw$scores$cv_gain; decision$cv_se[ms] <- fw$scores$cv_se
  }
  # forward path
  if (!is.null(fw$path)) {
    acc <- fw$path[fw$path$action == "accept", ]
    rej <- fw$path[fw$path$action == "reject-concurvity", ]
    mi <- match(acc$var, decision$covariate)
    decision$entry_step[mi] <- acc$step; decision$concurvity[mi] <- acc$concurvity
    decision$cv_gain[mi] <- acc$cv_gain; decision$cv_se[mi] <- acc$cv_se
    mr <- match(rej$var, decision$covariate)
    decision$status[mr] <- "dropped"
    decision$reason[mr] <- ifelse(is.na(rej$concurvity), "concurvity check failed", "concurvity")
    decision$concurvity[mr] <- rej$concurvity
    decision$cv_gain[mr] <- rej$cv_gain; decision$cv_se[mr] <- rej$cv_se
  }
  # backward drops vs final kept
  fwd_acc <- if (!is.null(fw$path)) fw$path$var[fw$path$action == "accept"] else character()
  back_dropped <- setdiff(fwd_acc, kept)
  decision$status[match(back_dropped, decision$covariate)] <- "dropped"
  decision$reason[match(back_dropped, decision$covariate)] <- "backward-BIC"
  ki <- match(kept, decision$covariate)
  decision$status[ki] <- "kept"; decision$reason[ki] <- "selected"
  decision$final_edf[ki] <- edf[kept]
  # anything still unlabelled = screened-in but never entered; say why, in the
  # order the forward loop tests it
  un <- is.na(decision$status)
  g  <- decision$cv_gain[un]; se <- decision$cv_se[un]
  decision$reason[un] <- ifelse(is.finite(g) & g <= 0, "no CV improvement",
                         ifelse(is.finite(g) & is.finite(se), "CV gain below 1 SE",
                                "CV scoring failed"))
  decision$status[un] <- "dropped"

  utils::write.csv(decision, file.path(outdir, "covariate_decision_table.csv"), row.names = FALSE)

  # -- diagnostics plots ---------------------------------------------------
  grDevices::png(file.path(outdir, "selection_diagnostics.png"), width = 1400, height = 1100, res = 130)
  # gam.check() also prints convergence info and the k-index table; keep that
  # in the console/log when verbose, swallow it otherwise.
  tryCatch({ par(mfrow = c(2, 2))
             if (verbose) gam.check(model) else utils::capture.output(gam.check(model)) },
           finally = grDevices::dev.off())

  grDevices::png(file.path(outdir, "partial_effects.png"), width = 1400, height = 1100, res = 130)
  # all.terms = TRUE hands parametric terms to termplot(), which evaluates
  # model$call$data (the symbol `data`) in the formula's environment -- the
  # global env, where `data` is utils::data. Give a plot-only copy the frame.
  pmodel <- model; pmodel$call$data <- final_data
  # A plot failure (e.g. an intercept-only model has no terms to plot) must not
  # cost the caller the finished selection: put the message in the PNG instead.
  tryCatch(plot(pmodel, pages = 1, scale = 0, all.terms = TRUE, shade = TRUE),
           error = function(e) { graphics::plot.new(); graphics::title(conditionMessage(e)) },
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
                 dev_expl = dev_expl, folds = folds,
                 logfile = if (do_log) file.path(outdir, logfile) else NULL))
}
