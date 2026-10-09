# =============================================================================
# 04_evaluate.R  --  ESS/ESR, convergence diagnostics, accuracy, figures, tables
# Needs (global): results, X_*, true_*, ref_cat, C, CLASS_NAMES, mcols, OUT_DIR
# =============================================================================

is_mcmc_method <- function(res) {
  is.list(res$out) && !is.null(res$out$posterior_coef)
}

# ---- ESS / ESR per sampler ----------------------------------------------------
get_mcmc_perf <- function(res) {
  rt <- res$runtime
  pc <- res$out$posterior_coef
  active <- setdiff(seq_len(C), ref_cat)
  ess_v <- unlist(lapply(active, function(j)
    lapply(seq_len(dim(pc)[1]), function(p)
      tryCatch(as.numeric(coda::effectiveSize(coda::mcmc(pc[p, j, ]))),
               error = function(e) NA_real_))))
  ess_v <- ess_v[!is.na(ess_v) & is.finite(ess_v)]
  if (length(ess_v) == 0) {
    return(data.frame(Method = res$name, Runtime_s = round(rt, 1),
                      Min_ESS = NA_real_, Med_ESS = NA_real_,
                      Min_ESR = NA_real_, Med_ESR = NA_real_))
  }
  data.frame(Method = res$name, Runtime_s = round(rt, 1),
             Min_ESS = round(min(ess_v), 1), Med_ESS = round(median(ess_v), 1),
             Min_ESR = round(min(ess_v) / rt, 3), Med_ESR = round(median(ess_v) / rt, 3))
}

mcmc_results <- results[sapply(results, is_mcmc_method)]
mcmc_perf <- do.call(rbind, lapply(mcmc_results, get_mcmc_perf))
rownames(mcmc_perf) <- NULL
print(mcmc_perf)
write.csv(mcmc_perf, file.path(OUT_DIR, "table_mcmc_performance.csv"), row.names = FALSE)

posterior_methods <- names(results)[sapply(results, is_mcmc_method)]
method_labels <- setNames(posterior_methods, posterior_methods)

# ---- convergence diagnostics for one target class ------------------------------
target_class <- 1L                      # glioma (reference class = pituitary)
P_check   <- min(6L, dim(results[[posterior_methods[1]]]$out$posterior_coef)[1])
coef_idx  <- seq_len(P_check)
coef_names <- c("intercept", paste0("PC", seq_len(P_check - 1)))
cat(sprintf("Convergence check for target class '%s' | coefficients: %s\n",
            CLASS_NAMES[target_class], paste(coef_names, collapse = ", ")))

# trace plots
pdf(file.path(OUT_DIR, "diag_trace_plots.pdf"), width = 9, height = 9)
for (m in posterior_methods) {
  pc <- results[[m]]$out$posterior_coef
  n_row <- ceiling(P_check / 2)
  par(mfrow = c(n_row, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
  for (k in seq_along(coef_idx)) {
    ch <- pc[coef_idx[k], target_class, ]
    plot(ch, type = "l", col = unname(mcols[m]), lwd = 0.6, main = coef_names[k],
         xlab = "Iteration", ylab = expression(beta))
    abline(h = mean(ch), col = "red", lty = 2, lwd = 1.3)
    if (length(ch) > 50) lines(cumsum(ch) / seq_along(ch), col = "orange", lwd = 1.2)
  }
  mtext(sprintf("Trace - %s (class: %s)", results[[m]]$name, CLASS_NAMES[target_class]),
        outer = TRUE, cex = 1.0, font = 2)
  par(mfrow = c(1, 1))
}
dev.off()

# ACF plots
pdf(file.path(OUT_DIR, "diag_acf_plots.pdf"), width = 9, height = 9)
for (m in posterior_methods) {
  pc <- results[[m]]$out$posterior_coef
  n_row <- ceiling(P_check / 2)
  par(mfrow = c(n_row, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
  for (k in seq_along(coef_idx)) {
    acf(pc[coef_idx[k], target_class, ], main = coef_names[k],
        col = unname(mcols[m]), lag.max = 50)
  }
  mtext(sprintf("ACF - %s (class: %s)", results[[m]]$name, CLASS_NAMES[target_class]),
        outer = TRUE, cex = 1.0, font = 2)
  par(mfrow = c(1, 1))
}
dev.off()

diag_table <- do.call(rbind, lapply(posterior_methods, function(m) {
  pc <- results[[m]]$out$posterior_coef
  n_iter <- dim(pc)[3]
  do.call(rbind, lapply(seq_along(coef_idx), function(k) {
    ch <- pc[coef_idx[k], target_class, ]
    mc <- coda::mcmc(ch)
    acf1 <- tryCatch(acf(ch, plot = FALSE, lag.max = 1)$acf[2, 1, 1],
                     error = function(e) NA_real_)
    ess  <- tryCatch(as.numeric(coda::effectiveSize(mc)), error = function(e) NA_real_)
    gz   <- tryCatch(as.numeric(coda::geweke.diag(mc)$z), error = function(e) NA_real_)
    data.frame(Method = m, Coef = coef_names[k], ACF_lag1 = round(acf1, 3),
               ESS = round(ess, 1), ESS_ratio = round(ess / n_iter, 3),
               Geweke_z = round(gz, 2),
               Converged = ifelse(!is.na(gz) && abs(gz) < 2, "yes", "check"))
  }))
}))
print(diag_table, row.names = FALSE)
write.csv(diag_table, file.path(OUT_DIR, "table_convergence_diagnostics.csv"), row.names = FALSE)

summary_by_method <- aggregate(cbind(ACF_lag1, ESS_ratio, Geweke_z) ~ Method,
                               data = diag_table,
                               FUN = function(x) mean(abs(x), na.rm = TRUE))
names(summary_by_method) <- c("Method", "Mean_|ACF_lag1|", "Mean_ESS_ratio", "Mean_|Geweke_z|")
print(summary_by_method, row.names = FALSE)
write.csv(summary_by_method, file.path(OUT_DIR, "table_convergence_summary.csv"), row.names = FALSE)

# ---- prediction (posterior-mean coefficients) ------------------------------------
softmax <- function(eta) {
  eta <- sweep(eta, 1, apply(eta, 1, max), "-")
  exp_eta <- exp(eta)
  exp_eta / rowSums(exp_eta)
}

predict_cls_posterior_coef <- function(pc, X) {
  bm  <- apply(pc, c(1, 2), mean)
  eta <- X %*% bm
  eta[, ref_cat] <- 0
  max.col(softmax(eta))
}

predict_cls_stickbreak <- function(pc, X, order_vec, ref_cat, C) {
  bm   <- apply(pc, c(1, 2), mean)
  n    <- nrow(X)
  prob <- matrix(0, n, C)
  remaining <- rep(1, n)
  K <- C - 1L
  for (k in seq_len(K)) {
    cls   <- order_vec[k]
    eta_k <- as.numeric(X %*% bm[, cls])
    p_k   <- 1 / (1 + exp(-eta_k))
    prob[, cls] <- remaining * p_k
    remaining   <- remaining * (1 - p_k)
  }
  prob[, order_vec[C]] <- remaining
  max.col(prob)
}

predict_cls_any <- function(res, X) {
  obj <- res$out
  if (is.list(obj) && !is.null(obj$posterior_coef)) {
    if (!is.null(obj$method) && obj$method == "stickbreak") {
      return(predict_cls_stickbreak(obj$posterior_coef, X, obj$order, ref_cat, C))
    }
    return(predict_cls_posterior_coef(obj$posterior_coef, X))
  }
  stop("Prediction not implemented for: ", res$name)
}

compute_acc <- function(res, X, true_cls) {
  tryCatch({
    pred <- predict_cls_any(res, X)
    mean(pred == true_cls) * 100
  }, error = function(e) {
    message("Prediction failed for: ", res$name, " | ", conditionMessage(e))
    NA_real_
  })
}

overfit_df <- do.call(rbind, lapply(names(results), function(nm) {
  data.frame(
    Method = results[[nm]]$name, Key = nm,
    Train = round(compute_acc(results[[nm]], X_train, true_train), 2),
    Val   = round(compute_acc(results[[nm]], X_val,   true_val),   2),
    Test  = round(compute_acc(results[[nm]], X_test,  true_test),  2)
  )
}))
overfit_df$Gap_Train_Test <- round(overfit_df$Train - overfit_df$Test, 2)
overfit_df$Gap_Train_Val  <- round(overfit_df$Train - overfit_df$Val, 2)
rownames(overfit_df) <- NULL
cat("\nAccuracy table (Train / Val / Test) and train-test gaps:\n")
print(overfit_df)
write.csv(overfit_df, file.path(OUT_DIR, "table_accuracy.csv"), row.names = FALSE)

# ---- confusion matrices of the best (test-accuracy) method ----------------------
best_m <- overfit_df$Key[which.max(overfit_df$Test)]
pred_test_best <- predict_cls_any(results[[best_m]], X_test)
pred_val_best  <- predict_cls_any(results[[best_m]], X_val)

make_cm <- function(pred, true_cls) {
  table(Predicted = factor(CLASS_NAMES[pred],     levels = CLASS_NAMES),
        Actual    = factor(CLASS_NAMES[true_cls], levels = CLASS_NAMES))
}
cm_test <- make_cm(pred_test_best, true_test)
cm_val  <- make_cm(pred_val_best,  true_val)
cat(sprintf("\nConfusion matrix (Test) - %s:\n", results[[best_m]]$name)); print(cm_test)
cat(sprintf("\nConfusion matrix (Validation) - %s:\n", results[[best_m]]$name)); print(cm_val)
write.csv(as.data.frame.matrix(cm_test), file.path(OUT_DIR, "table2_confusion_test_best.csv"))
write.csv(as.data.frame.matrix(cm_val),  file.path(OUT_DIR, "confusion_val_best.csv"))

# ---- Figure 2: accuracy by split ------------------------------------------------
overfit_long <- overfit_df %>%
  dplyr::select(Key, Train, Val, Test) %>%
  tidyr::pivot_longer(cols = c(Train, Val, Test), names_to = "Split", values_to = "Accuracy") %>%
  dplyr::mutate(Split = factor(Split, levels = c("Train", "Val", "Test")))

p_overfit <- ggplot(overfit_long, aes(x = Key, y = Accuracy, fill = Split)) +
  geom_col(position = position_dodge(0.75), width = 0.65) +
  geom_text(aes(label = paste0(round(Accuracy, 1), "%")),
            position = position_dodge(0.75), vjust = -0.4, size = 3.3) +
  scale_fill_manual(values = c("Train" = "#95A5A6", "Val" = "#F39C12", "Test" = "#2980B9")) +
  coord_cartesian(ylim = c(0, 105)) +
  labs(title = "Accuracy: Train / Validation / Test - overfitting check",
       subtitle = "A large gap between Train and Val/Test indicates overfitting",
       x = NULL, y = "Accuracy (%)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
        legend.position = "bottom")
ggsave(file.path(OUT_DIR, "fig2_accuracy_by_split.png"), p_overfit, width = 8, height = 5, dpi = 200)

# ---- Figure 3: 2x2 comparison panel ------------------------------------------------
plot_perf <- mcmc_perf
plot_perf$Key <- names(mcmc_results)
mcols_use <- mcols[names(mcols) %in% unique(plot_perf$Key)]

th <- theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
        legend.position = "bottom", axis.text.x = element_text(angle = 25, hjust = 1))

p1 <- plot_perf %>%
  dplyr::select(Key, Min_ESS, Med_ESS) %>%
  tidyr::pivot_longer(c(Min_ESS, Med_ESS), names_to = "Type", values_to = "ESS") %>%
  dplyr::mutate(Type = dplyr::recode(Type, "Min_ESS" = "Minimum", "Med_ESS" = "Median")) %>%
  ggplot(aes(x = Key, y = ESS, fill = Key, alpha = Type)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_x_discrete(labels = method_labels) + scale_fill_manual(values = mcols_use) +
  scale_alpha_manual(values = c("Minimum" = 0.5, "Median" = 1.0)) +
  facet_wrap(~Type, scales = "free_y") +
  labs(title = "Effective Sample Size (ESS)",
       subtitle = sprintf("Brain MRI | C=%d | N=%d | P=%d", C, N, P), x = NULL, y = "ESS") +
  th + guides(fill = "none", alpha = "none")

p2 <- plot_perf %>%
  dplyr::select(Key, Min_ESR, Med_ESR) %>%
  tidyr::pivot_longer(c(Min_ESR, Med_ESR), names_to = "Type", values_to = "ESR") %>%
  dplyr::mutate(Type = dplyr::recode(Type, "Min_ESR" = "Minimum", "Med_ESR" = "Median"),
                lESR = log10(pmax(ESR, 1e-5))) %>%
  ggplot(aes(x = Key, y = lESR, fill = Key, alpha = Type)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_x_discrete(labels = method_labels) + scale_fill_manual(values = mcols_use) +
  scale_alpha_manual(values = c("Minimum" = 0.5, "Median" = 1.0)) +
  facet_wrap(~Type) +
  labs(title = "log10(ESR) - Higher = Better Efficiency",
       subtitle = "ESR = Effective Sample Size per second", x = NULL, y = "log10(ESS/s)") +
  th + guides(fill = "none", alpha = "none")

p3 <- ggplot(plot_perf, aes(x = Key, y = Runtime_s, fill = Key)) +
  geom_col(width = 0.6, alpha = 0.9) +
  geom_text(aes(label = paste0(Runtime_s, "s")), vjust = -0.3, size = 4) +
  scale_x_discrete(labels = method_labels) + scale_fill_manual(values = mcols_use) +
  labs(title = "Runtime", subtitle = "Lower is faster", x = NULL, y = "Seconds") +
  th + guides(fill = "none")

p4 <- ggplot(overfit_df, aes(x = Key, y = Test, fill = Key)) +
  geom_col(width = 0.6, alpha = 0.9, na.rm = TRUE) +
  geom_text(aes(label = paste0(Test, "%")), vjust = -0.3, size = 4, na.rm = TRUE) +
  scale_x_discrete(labels = method_labels) + scale_fill_manual(values = mcols_use) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(title = "Test Set Accuracy", subtitle = "Higher is better", x = NULL, y = "Accuracy (%)") +
  th + guides(fill = "none")

png(file.path(OUT_DIR, "fig3_method_comparison.png"), width = 2400, height = 2000, res = 200)
gridExtra::grid.arrange(p1, p2, p3, p4, nrow = 2,
  top = "Comparison of Bayesian Multinomial Logit Methods - Brain Tumor MRI")
dev.off()

# ---- Figure 4: posterior 95% credible intervals (class glioma) ------------------------
n_feat_plot <- min(7, P)
ci_list <- lapply(posterior_methods, function(m) {
  cm2 <- results[[m]]$out$posterior_coef[1:n_feat_plot, 1, , drop = FALSE]
  cm2_mat <- matrix(cm2, nrow = n_feat_plot, ncol = dim(results[[m]]$out$posterior_coef)[3])
  feature_names <- c("intercept", paste0("PC", 1:(n_feat_plot - 1)))
  data.frame(Key = m, Feature = feature_names, Mean = rowMeans(cm2_mat),
             Lo = apply(cm2_mat, 1, quantile, probs = 0.025),
             Hi = apply(cm2_mat, 1, quantile, probs = 0.975))
})
ci_df <- do.call(rbind, ci_list)

p6 <- ggplot(ci_df, aes(x = Feature, y = Mean, color = Key, ymin = Lo, ymax = Hi)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  geom_errorbar(position = position_dodge(0.5), width = 0.3, linewidth = 0.8) +
  geom_point(position = position_dodge(0.5), size = 2.5) +
  scale_color_manual(values = mcols_use, labels = method_labels) +
  coord_flip() +
  labs(title = "Posterior 95% Credible Intervals",
       subtitle = "Class: glioma | Baseline: pituitary", x = NULL, y = expression(beta)) +
  th
ggsave(file.path(OUT_DIR, "fig4_credible_intervals.png"), p6, width = 8, height = 6, dpi = 200)

# ---- final summary ---------------------------------------------------------------------
final_table <- merge(mcmc_perf %>% dplyr::mutate(Key = names(mcmc_results)),
                     overfit_df[, c("Key", "Train", "Val", "Test", "Gap_Train_Test")],
                     by = "Key", all.x = TRUE, sort = FALSE)
final_table <- final_table[match(names(mcmc_results), final_table$Key), ]
write.csv(final_table, file.path(OUT_DIR, "table_final_summary.csv"), row.names = FALSE)
print(final_table, row.names = FALSE)

idx_fast    <- which.min(final_table$Runtime_s)
idx_med_ess <- which.max(final_table$Med_ESS)
idx_min_ess <- which.max(final_table$Min_ESS)
idx_med_esr <- which.max(final_table$Med_ESR)
idx_test    <- which.max(final_table$Test)
idx_gap     <- which.min(abs(final_table$Gap_Train_Test))
cat(sprintf(" * Fastest method:        %s (%.1f s)\n",  final_table$Key[idx_fast],    final_table$Runtime_s[idx_fast]))
cat(sprintf(" * Best median ESS:       %s (%.1f)\n",    final_table$Key[idx_med_ess], final_table$Med_ESS[idx_med_ess]))
cat(sprintf(" * Best minimum ESS:      %s (%.1f)\n",    final_table$Key[idx_min_ess], final_table$Min_ESS[idx_min_ess]))
cat(sprintf(" * Best median ESR:       %s (%.3f)\n",    final_table$Key[idx_med_esr], final_table$Med_ESR[idx_med_esr]))
cat(sprintf(" * Best test accuracy:    %s (%.2f%%)\n",  final_table$Key[idx_test],    final_table$Test[idx_test]))
cat(sprintf(" * Smallest train-test gap: %s (%.2f%%)\n", final_table$Key[idx_gap],    final_table$Gap_Train_Test[idx_gap]))
cat("\nAll tables and figures were written to ", OUT_DIR, "/\n", sep = "")
