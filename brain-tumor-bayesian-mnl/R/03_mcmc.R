# =============================================================================
# 03_mcmc.R  --  the four posterior samplers
#   AMH     : Adaptive Metropolis-Hastings                    (BayesMultiLogit)
#   DA+AMH  : Gamma augmentation + Adaptive MH                (BayesMultiLogit)
#   DA+eSS  : Gamma augmentation + Elliptical Slice Sampling  (BayesMultiLogit)
#   PG      : Polya-Gamma Gibbs sampler, stick-breaking form  (implemented here)
# Output (global): results   (cached in cache/mcmc_results.rds)
# =============================================================================

# ---- Polya-Gamma: binary logistic regression Gibbs sampler -------------------
fit_binary_PG <- function(y, X, n_sample, n_burn, prior_mean, prior_var) {
  n <- nrow(X); P <- ncol(X)
  B_inv  <- solve(prior_var)
  Binv_b <- as.numeric(B_inv %*% prior_mean)
  kappa  <- y - 0.5
  beta   <- rep(0, P)
  n_iter <- n_sample + n_burn
  draws  <- matrix(NA_real_, nrow = P, ncol = n_sample)

  for (iter in seq_len(n_iter)) {
    psi   <- as.numeric(X %*% beta)
    omega <- BayesLogit::rpg(n, h = 1, z = psi)

    Xw    <- X * omega
    V_inv <- crossprod(X, Xw) + B_inv
    V     <- chol2inv(chol(V_inv))
    m     <- as.numeric(V %*% (crossprod(X, kappa) + Binv_b))
    beta  <- as.numeric(m + t(chol(V)) %*% rnorm(P))

    if (iter > n_burn) draws[, iter - n_burn] <- beta
  }
  draws
}

# ---- Polya-Gamma: multinomial via stick-breaking -----------------------------
# Class order: all non-reference classes, then the reference class last.
# Step k models  P(class = order[k] | class in order[k:C])  as a binary logit.
multilogit_PG_stickbreak <- function(Y, X, n_sample, n_burn, prior_mean, prior_var,
                                     reference_cat, progress = TRUE) {
  C <- ncol(Y); P <- ncol(X)
  cls_int   <- max.col(Y)
  order_vec <- c(setdiff(seq_len(C), reference_cat), reference_cat)
  K <- C - 1L
  posterior_coef <- array(0, dim = c(P, C, n_sample))

  for (k in seq_len(K)) {
    target_cls  <- order_vec[k]
    active_mask <- cls_int %in% order_vec[k:C]
    y_k <- as.numeric(cls_int[active_mask] == target_cls)
    X_k <- X[active_mask, , drop = FALSE]

    if (progress) {
      cat(sprintf("  [PG] step %d/%d - target class '%s' | active samples = %d\n",
                  k, K, CLASS_NAMES[target_cls], nrow(X_k)))
    }

    draws_k <- fit_binary_PG(y_k, X_k, n_sample, n_burn, prior_mean, prior_var)
    posterior_coef[, target_cls, ] <- draws_k
    rm(draws_k, X_k, y_k); invisible(gc(FALSE))
  }
  list(posterior_coef = posterior_coef, method = "stickbreak", order = order_vec)
}

# ---- timing wrapper -----------------------------------------------------------
run_timed <- function(name, fun) {
  cat(sprintf("\n=== %s ===\n", name))
  set.seed(SEED)                       # same seed for every sampler
  tic <- proc.time()
  out <- fun()
  rt  <- (proc.time() - tic)[["elapsed"]]
  cat(sprintf("  done: %.1f s (%.1f min)\n", rt, rt / 60))
  invisible(gc(FALSE))
  list(out = out, runtime = rt, name = name)
}

mcmc_cache <- file.path(CACHE_DIR, "mcmc_results.rds")

if (file.exists(mcmc_cache) && !FORCE_MCMC) {
  cat("Loading cached MCMC results from", mcmc_cache,
      "(set FORCE_MCMC=1 to re-run)\n")
  results <- readRDS(mcmc_cache)
} else {
  results <- list()

  # NOTE: see "Known issues" in the README about which BayesMultiLogit function
  # the first entry should call. The calls below are exactly those used for the
  # thesis results.
  results[["AMH"]] <- run_timed("AMH - Adaptive Metropolis-Hastings", function()
    multilogit(Y_train, X_train, n_sample = n_sample, n_burn = n_burn,
               n_sigma_check = 25L, step_size = 0.05,
               prior = "normal", prior_mean = prior_mu, prior_var = prior_cov,
               reference_cat = ref_cat, probs = FALSE, progress = TRUE))

  results[["DA+AMH"]] <- run_timed("DA+AMH - Gamma Augmentation + Adaptive MH", function()
    multilogit(Y_train, X_train, n_sample = n_sample, n_burn = n_burn,
               n_sigma_check = 25L, step_size = 0.1,
               prior = "normal", prior_mean = prior_mu, prior_var = prior_cov,
               reference_cat = ref_cat, probs = FALSE, progress = TRUE))

  results[["DA+eSS"]] <- run_timed("DA+eSS - Gamma Augmentation + Elliptical Slice Sampling", function()
    multilogit_ESS(Y_train, X_train, n_sample = n_sample, n_burn = n_burn,
                   prior_mean = prior_mu, prior_var = prior_cov,
                   reference_cat = ref_cat, probs = FALSE, progress = TRUE))

  results[["PG"]] <- run_timed("PG - Polya-Gamma Bayesian Logit (Polson et al. 2013)", function()
    multilogit_PG_stickbreak(Y_train, X_train, n_sample = n_sample, n_burn = n_burn,
                             prior_mean = prior_mu, prior_var = prior_cov,
                             reference_cat = ref_cat, progress = TRUE))

  # thinning (runtime above is measured BEFORE thinning)
  for (nm in names(results)) {
    obj <- results[[nm]]$out
    if (is.list(obj) && !is.null(obj$posterior_coef)) {
      n_draws <- dim(obj$posterior_coef)[3]
      keep <- seq(1, n_draws, by = THIN)
      obj$posterior_coef <- obj$posterior_coef[, , keep, drop = FALSE]
      results[[nm]]$out <- obj
      cat(sprintf("  [thin] %s: %d -> %d draws\n", nm, n_draws, length(keep)))
    }
  }
  saveRDS(results, mcmc_cache)
}
