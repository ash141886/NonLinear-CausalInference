# Non-linear-Causal-Inference-with-GAM-and-HSIC

Exploring Non-linear Causal Inference using GAM and HSIC

## Overview
This repository provides a multi-file R project that:
1. Generates synthetic data with configurable nonlinearity, sparsity, and noise.
2. Implements two causal discovery methods:
   - Proposed Method (GAM + HSIC threshold)
   - LiNGAM (via `fastICA`)
3. Calculates metrics (Graph Accuracy, SHD, Misoriented Edges, Recall, F1-Score, MSE, and Execution Time) to compare performance.
4. Runs experiments in parallel across a grid of parameters (number of variables, sample sizes).
5. Saves results and produces plots in PDF and PNG formats.

