# =============================================================================
# 00_config.R  --  packages, paths, constants, MCMC settings
# =============================================================================

required_pkgs <- c("BayesMultiLogit", "BayesLogit", "jpeg", "png", "coda",
                   "ggplot2", "dplyr", "tidyr", "gridExtra", "scales", "nnet")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace,
                                      logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing packages: ", paste(missing_pkgs, collapse = ", "),
       "\nRun  Rscript install_packages.R  first.", call. = FALSE)
}
invisible(lapply(required_pkgs, library, character.only = TRUE))

# ---- reproducibility --------------------------------------------------------
SEED <- 1402L
set.seed(SEED)

# ---- paths ------------------------------------------------------------------
# Folder that contains Training/ and Testing/ (Kaggle "Brain Tumor MRI Dataset").
# Override with an environment variable, e.g.  BRAIN_MRI_DIR=/path/to/archive
DATA_ROOT <- Sys.getenv("BRAIN_MRI_DIR", unset = "data/brain-tumor-mri-dataset")
OUT_DIR   <- "results"
CACHE_DIR <- "cache"
dir.create(OUT_DIR,   showWarnings = FALSE, recursive = TRUE)
dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- data constants ---------------------------------------------------------
CLASS_NAMES <- c("glioma", "meningioma", "notumor", "pituitary")
C <- length(CLASS_NAMES)
N_PER_CLASS  <- Inf    # max training images per class (Inf = use all)
N_TEST_CLASS <- Inf    # max test images per class      (Inf = use all)
VAL_FRAC     <- 0.25   # fraction of Training/ held out for validation

# ---- MCMC settings (identical for all four samplers) ------------------------
# For a quick smoke test:  QUICK=1 Rscript run_all.R
if (Sys.getenv("QUICK", unset = "0") == "1") {
  n_sample <- 300L; n_burn <- 300L
} else {
  n_sample <- 6000L; n_burn <- 6000L
}
THIN <- 5L

# Re-run the (very slow) MCMC even if a cached result exists:  FORCE_MCMC=1
FORCE_MCMC <- Sys.getenv("FORCE_MCMC", unset = "0") == "1"
# Re-extract image features even if a cache exists:            FORCE_FEATURES=1
FORCE_FEATURES <- Sys.getenv("FORCE_FEATURES", unset = "0") == "1"

# Colours / labels used in all figures
mcols <- c("AMH" = "#2C3E50", "DA+AMH" = "#E67E22",
           "DA+eSS" = "#2980B9", "PG" = "#27AE60")
