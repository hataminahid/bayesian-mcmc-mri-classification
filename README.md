# Bayesian multinomial logistic regression for brain-tumour MRI – sampler comparison

Four posterior samplers for Bayesian multinomial logistic regression are compared
on 4-class brain-MRI data (glioma, meningioma, no tumour, pituitary):

| Key | Sampler | Source |
|---|---|---|
| `AMH` | Adaptive Metropolis–Hastings | `BayesMultiLogit::multilogit` (see [NOTES](../NOTES.md) #1) |
| `DA+AMH` | Gamma data augmentation + adaptive MH | `BayesMultiLogit::multilogit` |
| `DA+eSS` | Gamma data augmentation + elliptical slice sampling | `BayesMultiLogit::multilogit_ESS` |
| `PG` | Pólya–Gamma Gibbs sampler, stick-breaking decomposition | implemented in `R/03_mcmc.R` (uses `BayesLogit::rpg`) |

Pipeline: HOG (1296) + LBP (16) + intensity statistics (8) → PCA (80 PCs chosen on
a validation set) → MCMC (6000 burn-in, 6000 kept draws, thinning 5, N(0, I) prior)
→ ESS, ESS/second, Geweke, ACF, accuracy, confusion matrix.

## 1. Data

Download the *Brain Tumor MRI Dataset* (Masoud Nickparvar) from Kaggle:
<https://www.kaggle.com/datasets/masoudnickparvar/brain-tumor-mri-dataset>

```bash
pip install kaggle        # needs ~/.kaggle/kaggle.json
kaggle datasets download -d masoudnickparvar/brain-tumor-mri-dataset -p data --unzip
```

Expected layout (folder names are case-sensitive on Linux):

```
data/brain-tumor-mri-dataset/     (or any path given in BRAIN_MRI_DIR)
├── Training/{glioma,meningioma,notumor,pituitary}/*.jpg
└── Testing/{glioma,meningioma,notumor,pituitary}/*.jpg
```

If the archive unzips to a different top-level name, simply point
`BRAIN_MRI_DIR` at the folder that contains `Training/` and `Testing/`.

## 2. Install

```bash
Rscript install_packages.R
```
(`BayesMultiLogit` comes from GitHub: `kylemcevoy/BayesMultiLogit`;
`BayesLogit` is installed from CRAN or the CRAN archive.)
A C++ compiler (Rtools on Windows) is required.

## 3. Run

```bash
cd brain-tumor-bayesian-mnl
BRAIN_MRI_DIR=/path/to/dataset Rscript run_all.R
```

| Variable | Effect |
|---|---|
| `BRAIN_MRI_DIR` | folder containing `Training/` and `Testing/` |
| `QUICK=1` | 300 + 300 MCMC iterations (smoke test, minutes) |
| `FORCE_FEATURES=1` | re-extract image features (otherwise `cache/features.rds` is reused) |
| `FORCE_MCMC=1` | re-run the samplers (otherwise `cache/mcmc_results.rds` is reused) |

The full run takes **several hours** (in the report: AMH 4750 s, DA+AMH 5170 s,
DA+eSS 413 s, PG 1258 s for the sampling alone) plus feature extraction.

## 4. Scripts

| File | Purpose |
|---|---|
| `R/00_config.R` | packages, paths, seeds, MCMC settings |
| `R/01_features.R` | image loading, HOG / LBP / intensity features, train/val/test split |
| `R/02_pca.R` | PCA, choice of #PCs, standardised design matrices, Fig. 1 |
| `R/03_mcmc.R` | the four samplers, timing, thinning |
| `R/04_evaluate.R` | ESS/ESR, convergence diagnostics, accuracy, confusion matrix, Figs. 2–4 |
| `run_all.R` | runs everything in order |

## 5. Outputs (`results/`)

| File | Report item |
|---|---|
| `table1_split_counts.csv` | Table 1 |
| `fig1_pca_cumvar_scree.png` | Fig. 1 |
| `table_accuracy.csv`, `fig2_accuracy_by_split.png` | Fig. 2 |
| `table_mcmc_performance.csv`, `fig3_method_comparison.png` | Fig. 3 (ESS, ESR, runtime, test accuracy) |
| `fig4_credible_intervals.png` | Fig. 4 |
| `table2_confusion_test_best.csv` | Table 2 |
| `diag_trace_plots.pdf`, `diag_acf_plots.pdf`, `table_convergence_*.csv` | convergence diagnostics |
| `table_final_summary.csv`, `sessionInfo.txt` | overall summary, R session |

## 6. Troubleshooting

* `Folder not found: …/Training/glioma` → wrong `BRAIN_MRI_DIR`.
* `there is no package called 'BayesLogit'` → run `install_packages.R`; on recent R
  versions install from the CRAN archive (`remotes::install_version("BayesLogit", "2.1")`).
* Out of memory during PG: the PG sampler builds an `N × P` matrix every iteration
  (N ≈ 4000, P = 81) – 8 GB RAM is enough.

---
