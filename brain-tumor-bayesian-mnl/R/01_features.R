# =============================================================================
# 01_features.R  --  image loading, HOG + LBP + intensity features, data split
# Output (global): df_train, df_val, df_test   (cached in cache/features.rds)
# =============================================================================

read_gray <- function(path) {
  ext <- tolower(tools::file_ext(path))
  img <- tryCatch({
    if (ext == "jpg" || ext == "jpeg") {
      jpeg::readJPEG(path)
    } else if (ext == "png") {
      png::readPNG(path)
    } else {
      stop("unknown image format")
    }
  }, error = function(e) NULL)
  if (is.null(img)) return(NULL)

  if (length(dim(img)) == 3) {
    if (dim(img)[3] >= 3) {
      gray <- 0.299 * img[, , 1] + 0.587 * img[, , 2] + 0.114 * img[, , 3]
    } else {
      gray <- img[, , 1]
    }
  } else {
    gray <- img
  }
  gray
}

resize_mat <- function(mat, target = 96L) {
  nr <- nrow(mat); nc <- ncol(mat)
  row_idx <- round(seq(1, nr, length.out = target))
  col_idx <- round(seq(1, nc, length.out = target))
  mat[row_idx, col_idx]
}

# Histogram of Oriented Gradients: 8x8 cells, 9 unsigned orientation bins
hog_features <- function(gray, cell_size = 8L, n_bins = 9L) {
  sz <- nrow(gray)
  n_cells <- sz %/% cell_size

  Gx <- cbind(gray[, -1], gray[, ncol(gray)]) -
        cbind(gray[, 1], gray[, -ncol(gray)])
  Gy <- rbind(gray[-1, ], gray[nrow(gray), ]) -
        rbind(gray[1, ], gray[-nrow(gray), ])

  mag   <- sqrt(Gx^2 + Gy^2)
  angle <- (atan2(Gy, Gx) * 180 / pi) %% 180

  feats <- c()
  for (i in seq_len(n_cells)) {
    for (j in seq_len(n_cells)) {
      r  <- ((i - 1) * cell_size + 1):(i * cell_size)
      cc <- ((j - 1) * cell_size + 1):(j * cell_size)
      r  <- r[r <= sz]; cc <- cc[cc <= sz]
      cell_mag <- mag[r, cc]
      cell_ang <- angle[r, cc]
      hist_vals <- numeric(n_bins)
      bin_width <- 180 / n_bins
      for (b in seq_len(n_bins)) {
        lo <- (b - 1) * bin_width; hi <- b * bin_width
        mask <- cell_ang >= lo & cell_ang < hi
        hist_vals[b] <- sum(cell_mag[mask])
      }
      norm  <- sqrt(sum(hist_vals^2) + 1e-6)
      feats <- c(feats, hist_vals / norm)
    }
  }
  feats
}

# Local Binary Patterns (8 neighbours); the 256-code histogram is aggregated
# into n_bins groups of consecutive codes
lbp_features <- function(gray, n_bins = 16L) {
  nr <- nrow(gray); nc <- ncol(gray)
  lbp_img <- matrix(0L, nr, nc)
  offsets <- list(c(-1, -1), c(-1, 0), c(-1, 1), c(0, 1),
                  c(1, 1),   c(1, 0),  c(1, -1), c(0, -1))
  for (k in seq_along(offsets)) {
    dr <- offsets[[k]][1]; dc <- offsets[[k]][2]
    r_s  <- pmax(1, 1 + dr);  r_e  <- pmin(nr, nr + dr)
    c_s  <- pmax(1, 1 + dc);  c_e  <- pmin(nc, nc + dc)
    r_t  <- pmax(1, 1 - dr);  r_te <- pmin(nr, nr - dr)
    c_t  <- pmax(1, 1 - dc);  c_te <- pmin(nc, nc - dc)
    neighbor <- matrix(0, nr, nc)
    neighbor[r_t:r_te, c_t:c_te] <- gray[r_s:r_e, c_s:c_e]
    lbp_img <- lbp_img + (neighbor >= gray) * (2L^(k - 1L))
  }
  full_hist  <- tabulate(lbp_img + 1L, nbins = 256L)
  group_size <- 256L %/% n_bins
  h <- vapply(seq_len(n_bins), function(b) {
    lo <- (b - 1L) * group_size + 1L
    hi <- if (b == n_bins) 256L else b * group_size
    sum(full_hist[lo:hi])
  }, numeric(1))
  h / (sum(h) + 1e-9)
}

# 8 first-order intensity statistics
intensity_stats <- function(gray) {
  v <- as.vector(gray)
  c(mean     = mean(v),
    sd       = sd(v),
    q25      = quantile(v, 0.25, names = FALSE),
    q75      = quantile(v, 0.75, names = FALSE),
    q90      = quantile(v, 0.90, names = FALSE),
    contrast = max(v) - min(v),
    energy   = mean(v^2),
    entropy  = {
      h <- hist(v, breaks = 32, plot = FALSE)$counts
      p <- h / sum(h); p <- p[p > 0]
      -sum(p * log(p))
    })
}

# One image -> feature vector (1296 HOG + 16 LBP + 8 intensity = 1320)
extract_one <- function(path, img_size = 96L) {
  gray <- read_gray(path)
  if (is.null(gray)) return(NULL)
  gray <- resize_mat(gray, img_size)
  gray <- (gray - min(gray)) / (max(gray) - min(gray) + 1e-9)

  c(
    hog_features(gray, cell_size = 8L, n_bins = 9L),
    lbp_features(gray, n_bins = 16L),
    intensity_stats(gray)
  )
}

# Load every image of one split ("Training" / "Testing") into a data.frame
load_split <- function(split_name, n_per_class) {
  records <- list()
  for (cls in CLASS_NAMES) {
    folder <- file.path(DATA_ROOT, split_name, cls)
    if (!dir.exists(folder)) stop("Folder not found: ", folder)

    imgs <- list.files(folder, pattern = "\\.(jpg|jpeg|png)$",
                       full.names = TRUE, ignore.case = TRUE)
    if (length(imgs) == 0) stop("No images in: ", folder)

    sel <- if (is.infinite(n_per_class) || length(imgs) <= n_per_class) {
      imgs
    } else {
      sample(imgs, n_per_class)
    }

    cat(sprintf("  [%s/%s] processing %d images...\n", split_name, cls, length(sel)))

    ok <- 0L
    for (p in sel) {
      f <- extract_one(p)
      if (!is.null(f)) {
        records[[length(records) + 1]] <- c(label = which(CLASS_NAMES == cls),
                                            label_name = cls, f)
        ok <- ok + 1L
      }
    }
    cat(sprintf("    OK %d/%d\n", ok, length(sel)))
  }

  df <- as.data.frame(do.call(rbind, records))
  df$label      <- as.integer(df$label)
  df$label_name <- as.character(df$label_name)
  feat_cols <- paste0("f", seq_len(ncol(df) - 2))
  names(df)[3:ncol(df)] <- feat_cols
  for (fc in feat_cols) df[[fc]] <- as.numeric(df[[fc]])
  df
}

# Stratified train / validation split of the Training folder
split_train_val <- function(df, val_frac = VAL_FRAC, seed = SEED) {
  set.seed(seed)
  idx_val <- unlist(lapply(split(seq_len(nrow(df)), df$label_name), function(idx) {
    n_val <- max(1L, round(length(idx) * val_frac))
    sample(idx, n_val)
  }))
  list(train = df[-idx_val, ], val = df[idx_val, ])
}

# ---- run (with cache) --------------------------------------------------------
features_cache <- file.path(CACHE_DIR, "features.rds")

if (file.exists(features_cache) && !FORCE_FEATURES) {
  cat("Loading cached features from", features_cache, "\n")
  tmp <- readRDS(features_cache)
  df_train <- tmp$train; df_val <- tmp$val; df_test <- tmp$test
  rm(tmp)
} else {
  cat("Extracting features (this takes a while)...\n")
  df_train_full <- load_split("Training", n_per_class = N_PER_CLASS)
  df_test       <- load_split("Testing",  n_per_class = N_TEST_CLASS)
  sp <- split_train_val(df_train_full)
  df_train <- sp$train
  df_val   <- sp$val
  rm(sp, df_train_full); invisible(gc(FALSE))
  saveRDS(list(train = df_train, val = df_val, test = df_test), features_cache)
}

cat(sprintf("\nTrain: %d | Val: %d | Test: %d samples\n",
            nrow(df_train), nrow(df_val), nrow(df_test)))
cat("Class distribution (Train):\n"); print(table(df_train$label_name))
cat("Class distribution (Val):\n");   print(table(df_val$label_name))
cat("Class distribution (Test):\n");  print(table(df_test$label_name))

# Table 1 of the report
split_table <- rbind(Train = table(df_train$label_name),
                     Val   = table(df_val$label_name),
                     Test  = table(df_test$label_name))
split_table <- cbind(split_table, Total = rowSums(split_table))
write.csv(split_table, file.path(OUT_DIR, "table1_split_counts.csv"))
