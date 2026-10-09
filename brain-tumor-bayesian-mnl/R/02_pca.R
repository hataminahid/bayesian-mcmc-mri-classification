# =============================================================================
# 02_pca.R  --  PCA, choice of number of PCs, standardised design matrices
# Output (global): X_train, X_val, X_test, Y_train, Y_val, Y_test, N, P, ref_cat,
#                  n_pc, prior_mu, prior_cov, true_train/val/test
# =============================================================================

feat_cols <- grep("^f[0-9]+$", names(df_train), value = TRUE)
# drop constant columns (prcomp with scale. = TRUE cannot handle them)
col_sd <- vapply(df_train[, feat_cols], sd, numeric(1))
feat_cols <- feat_cols[is.finite(col_sd) & col_sd > 0]

pca_model <- prcomp(df_train[, feat_cols], center = TRUE, scale. = TRUE)
cum_var <- cumsum(pca_model$sdev^2) / sum(pca_model$sdev^2)

pca_plot_df <- data.frame(
  PC = seq_along(pca_model$sdev),
  Variance = (pca_model$sdev^2) / sum(pca_model$sdev^2),
  CumulativeVariance = cum_var
)

X_train_full_pca <- predict(pca_model, df_train[, feat_cols])
X_val_full_pca   <- predict(pca_model, df_val[, feat_cols])
X_test_full_pca  <- predict(pca_model, df_test[, feat_cols])

# ---- choose the number of PCs by validation accuracy of a plain multinomial GLM
max_pc  <- ncol(X_train_full_pca)
pc_grid <- unique(pmin(c(10, 20, 30, 50, 80, 120, 160), max_pc))

val_acc_by_pc <- sapply(pc_grid, function(k) {
  Xtr_k <- scale(X_train_full_pca[, 1:k, drop = FALSE])
  ctr <- attr(Xtr_k, "scaled:center"); scl <- attr(Xtr_k, "scaled:scale")
  Xv_k <- scale(X_val_full_pca[, 1:k, drop = FALSE], center = ctr, scale = scl)

  fit <- tryCatch(
    nnet::multinom(y ~ ., data = data.frame(y = df_train$label_name, Xtr_k),
                   trace = FALSE, MaxNWts = 20000, maxit = 300),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NA_real_)
  pred <- predict(fit, newdata = data.frame(Xv_k))
  mean(pred == df_val$label_name) * 100
})
names(val_acc_by_pc) <- pc_grid
print(round(val_acc_by_pc, 2))

n_pc <- as.integer(pc_grid[which.max(val_acc_by_pc)])
cat(sprintf("Selected PCs (by validation accuracy): %d (Val acc = %.2f%%, cumulative variance = %.1f%%)\n",
            n_pc, max(val_acc_by_pc, na.rm = TRUE), cum_var[n_pc] * 100))
write.csv(data.frame(n_pc = pc_grid, val_acc = as.numeric(val_acc_by_pc)),
          file.path(OUT_DIR, "pca_validation_accuracy.csv"), row.names = FALSE)

# ---- Figure 1: scree plot and cumulative variance ----------------------------
n_show <- min(nrow(pca_plot_df), max(n_pc * 3, 40))

p_scree <- ggplot(pca_plot_df[1:n_show, ], aes(x = PC, y = Variance)) +
  geom_line(color = "#4C78A8", linewidth = 1.1) +
  geom_point(color = "#4C78A8", size = 2.2) +
  geom_vline(xintercept = n_pc, linetype = "dashed", color = "red", linewidth = 0.9) +
  annotate("text", x = n_pc, y = pca_plot_df$Variance[n_pc],
           label = paste0("Selected PC = ", n_pc), hjust = -0.1, vjust = -0.7,
           size = 4, color = "red") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(title = "Scree Plot of PCA", x = "Principal Component", y = "Explained Variance") +
  theme_minimal(base_size = 13)

p_cumvar <- ggplot(pca_plot_df[1:n_show, ], aes(x = PC, y = CumulativeVariance)) +
  geom_line(color = "#2ca02c", linewidth = 1.1) +
  geom_point(color = "#2ca02c", size = 2.5) +
  geom_vline(xintercept = n_pc, linetype = "dashed", color = "darkgreen", linewidth = 1) +
  geom_hline(yintercept = cum_var[n_pc], linetype = "dashed", color = "gray50") +
  annotate("point", x = n_pc, y = cum_var[n_pc], color = "darkgreen", size = 3.5) +
  annotate("text", x = n_pc + 1, y = cum_var[n_pc] - 0.04,
           label = paste0("Selected PC = ", n_pc, "\nCum. Var = ",
                          round(cum_var[n_pc] * 100, 1), "%"),
           color = "darkgreen", hjust = 0, fontface = "bold") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(title = "Cumulative Explained Variance by PCA (zoom)",
       x = "Principal Component", y = "Cumulative Explained Variance") +
  theme_minimal(base_size = 13) + theme(panel.grid.minor = element_blank())

ggsave(file.path(OUT_DIR, "fig1_pca_cumvar_scree.png"),
       gridExtra::arrangeGrob(p_cumvar, p_scree, ncol = 2),
       width = 12, height = 5, dpi = 200)

# ---- standardised design matrices (intercept + n_pc scaled PCs) --------------
X_tr_pca  <- X_train_full_pca[, 1:n_pc, drop = FALSE]
X_val_pca <- X_val_full_pca[,   1:n_pc, drop = FALSE]
X_te_pca  <- X_test_full_pca[,  1:n_pc, drop = FALSE]
rm(X_train_full_pca, X_val_full_pca, X_test_full_pca); invisible(gc(FALSE))

X_sc  <- scale(X_tr_pca)
sc_c2 <- attr(X_sc, "scaled:center")
sc_s2 <- attr(X_sc, "scaled:scale")
X_train <- cbind(1, X_sc)
X_val   <- cbind(1, scale(X_val_pca, center = sc_c2, scale = sc_s2))
X_test  <- cbind(1, scale(X_te_pca,  center = sc_c2, scale = sc_s2))
colnames(X_train)[1] <- colnames(X_val)[1] <- colnames(X_test)[1] <- "intercept"

make_Y <- function(df) {
  Y <- matrix(0L, nrow(df), C, dimnames = list(NULL, CLASS_NAMES))
  for (j in seq_len(C)) Y[, j] <- as.integer(df$label_name == CLASS_NAMES[j])
  Y
}
Y_train <- make_Y(df_train)
Y_val   <- make_Y(df_val)
Y_test  <- make_Y(df_test)

true_train <- apply(Y_train, 1, which.max)
true_val   <- apply(Y_val,   1, which.max)
true_test  <- apply(Y_test,  1, which.max)

N <- nrow(X_train); P <- ncol(X_train); ref_cat <- C   # reference class = pituitary
cat(sprintf("N=%d | P=%d (intercept+%d PC) | C=%d\n", N, P, n_pc, C))
cat(sprintf("Total parameters: %d x %d = %d\n", P, C - 1, P * (C - 1)))

rm(X_tr_pca, X_val_pca, X_te_pca); invisible(gc(FALSE))

# prior: beta_j ~ N(0, I)
prior_mu  <- rep(0, P)
prior_cov <- diag(1, P)
cat(sprintf("\nMCMC: n_sample=%d | n_burn=%d | thin=%d\n", n_sample, n_burn, THIN))
