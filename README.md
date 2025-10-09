# Non-linear Causal Inference in Observational Data using Additive Models and Kernel-based Independence Testing

## Overview
This repository provides a multi-file R project that:
1. Generates synthetic data with configurable nonlinearity, sparsity, and noise.
2. Implements two causal discovery methods:
   - Proposed Method (additive models + HSIC threshold)
   - LiNGAM (via `fastICA`)
3. Calculates metrics (Graph Accuracy, SHD, Misoriented Edges, Recall, F1-Score, MSE, and Execution Time) to compare performance.
4. Runs experiments in parallel across a grid of parameters (number of variables, sample sizes).


.
├── src/
│   ├── data_generation.R
│   ├── methods/
│   │   ├── proposed_method.R
│   │   └── lingam.R
│   ├── metrics.R
│   └── plotting.R
├── main_causal_analysis_simulated.R
├── main_wine_red_causal_analysis.R
├── results/
│   └── .gitkeep
├── plots/
│   └── .gitkeep
├── README.md
├── .gitignore
└── requirements.txt




# Causal Discovery Simulation Experiment

This repository provides code for nonlinear causal inference using Generalized Additive Models (GAM) and the Hilbert-Schmidt Independence Criterion (HSIC), supporting both simulation and real data analysis on the UCI Wine Quality datasets.

## Quick Start

Copy and run the following code in your R console to **clone the repository**, **install all required packages**, and execute the analysis.


# Clone repository, install all packages, and run any analysis (simulated and red wine) in one go

```r
if (!requireNamespace("git2r", quietly = TRUE)) install.packages("git2r")
library(git2r)
repo_url <- "https://github.com/ash141886/NonLinear-CausalInference.git"
local_path <- "NonLinear-CausalInference"
if (!dir.exists(local_path)) {
  git2r::clone(url = repo_url, local_path = local_path)
}
setwd(local_path)

pkgs_needed <- c(
  "dplyr",
  "tidyr",
  "ggplot2",
  "foreach",
  "doParallel",
  "parallel",
  "gridExtra",
  "mgcv",
  "fastICA",
   "progress"
)
for (pkg in pkgs_needed) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# --- Run any analysis script as needed ---
# For simulated data:
source("main_causal_analysis_simulated.R")

# For red wine data:
source("main_wine_red_causal_analysis.R")

