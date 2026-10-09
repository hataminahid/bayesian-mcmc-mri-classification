# Installs every package needed by the brain-tumour pipeline.
cran <- c("remotes", "jpeg", "png", "coda", "ggplot2", "dplyr", "tidyr",
          "gridExtra", "scales", "nnet")
need <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(need)) install.packages(need)

# BayesMultiLogit (Fisher & McEvoy) is only on GitHub
if (!requireNamespace("BayesMultiLogit", quietly = TRUE))
  remotes::install_github("kylemcevoy/BayesMultiLogit")

# BayesLogit (Polya-Gamma sampler): try CRAN first, then the CRAN archive
if (!requireNamespace("BayesLogit", quietly = TRUE)) {
  ok <- tryCatch({ install.packages("BayesLogit"); requireNamespace("BayesLogit", quietly = TRUE) },
                 error = function(e) FALSE, warning = function(w) FALSE)
  if (!ok) remotes::install_version("BayesLogit", version = "2.1")
}
cat("All packages installed.\n")
