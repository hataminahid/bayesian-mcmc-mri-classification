# Run the whole brain-tumour pipeline.  Execute from THIS folder:
#   Rscript run_all.R
# Environment variables (all optional):
#   BRAIN_MRI_DIR   path to the folder containing Training/ and Testing/
#   QUICK=1         short MCMC (300 + 300 iterations) for a smoke test
#   FORCE_FEATURES=1, FORCE_MCMC=1   ignore caches in cache/
t0 <- Sys.time()
for (f in c("R/00_config.R", "R/01_features.R", "R/02_pca.R",
            "R/03_mcmc.R",   "R/04_evaluate.R")) {
  cat("\n########## ", f, " ##########\n", sep = "")
  source(f, echo = FALSE)
}
cat("\nTotal time:", format(round(difftime(Sys.time(), t0, units = "mins"), 1)), "\n")
sink(file.path("results", "sessionInfo.txt")); print(sessionInfo()); sink()
