# Exploring Non-linear Causal Inference using GAM and HSIC

## Overview
This repository provides a multi-file R project that:
1. Generates synthetic data with configurable nonlinearity, sparsity, and noise.
2. Implements two causal discovery methods:
   - Proposed Method (GAM + HSIC threshold)
   - LiNGAM (via `fastICA`)
3. Calculates metrics (Graph Accuracy, SHD, Misoriented Edges, Recall, F1-Score, MSE, and Execution Time) to compare performance.
4. Runs experiments in parallel across a grid of parameters (number of variables, sample sizes).
5. Saves results and produces plots in PDF and PNG formats.


.
├── src/
│ ├── data_generation.R
│ ├── methods/
│ │ ├── proposed_method.R
│ │ └── lingam.R
│ ├── metrics.R
│ └── plotting.R
├── scripts/
│ └── run_experiment.R
├── results/
│ └── .gitkeep
├── plots/
│ └── .gitkeep
├── README.md
├── .gitignore
└── requirements.txt




# Causal Discovery Simulation Experiment

This repository provides code for nonlinear causal inference using Generalized Additive Models (GAM) and the Hilbert-Schmidt Independence Criterion (HSIC), supporting both simulation and real data analysis.

## Quick Start

Copy and run the following code in your R console to clone the repository, install all required packages, and execute the analysis.  
**For simulation, run the `run_main_simulated.R` script. For real data, run `run_main_real.R`.**

```r
if (!requireNamespace("git2r", quietly = TRUE)) install.packages("git2r")
library(git2r)
repo_url <- "https://github.com/ash141886/Exploring-Non-linear-Causal-Inference-using-GAM-and-HSIC.git"
local_path <- "Exploring-Non-linear-Causal-Inference-using-GAM-and-HSIC"
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
  "fastICA"
)
for (pkg in pkgs_needed) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# For simulation study
source("run_main_simulated.R")

# For real data analysis (uncomment the following line)
# source("run_main_real.R")
