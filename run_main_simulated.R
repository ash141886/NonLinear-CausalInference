# ==========================================================
# Causal Discovery Simulation Framework
# - Benchmarking "Proposed Method" vs. LiNGAM
# - Computes F1, Accuracy, SHD, MSE, MAE, Runtime, etc.
# - Publication-ready version
# ==========================================================

# Load Required Libraries
library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)
library(foreach)
library(doParallel)

# Source Modular Scripts (Provide full paths or setwd as needed)
source("src/data_generation.R")
source("src/methods/proposed_method.R")
source("src/methods/lingam.R")
source("src/metrics.R")
source("src/analysis.R")


# -------------------------------------------
# Experiment Grid Parameters
# -------------------------------------------
n_vars_range <- c(4, 5, 6, 7)
n_samples_range <- c(300, 400, 500, 560)
nonlinearity <- 0.5
sparsity <- 0.3
noise_level <- 0.1
gam_sigma <- 0.5
threshold_percentile <- 0.1
n_reps <- 10

# -------------------------------------------
# Parallel Experiment Execution
# -------------------------------------------
num_cores <- detectCores() - 1
registerDoParallel(cores = num_cores)
results <- foreach(rep = 1:n_reps, .combine = rbind) %:%
    foreach(n_vars = n_vars_range, .combine = rbind) %:%
    foreach(n_samples = n_samples_range, .combine = rbind) %dopar% {
        analyze_causal_structure(n_vars, n_samples, nonlinearity, sparsity, 
                                 noise_level, gam_sigma, threshold_percentile)
    }
stopImplicitCluster()

# -------------------------------------------
# Save & Summarize Results
# -------------------------------------------
if (!dir.exists("results")) dir.create("results", recursive = TRUE)
write.csv(results, "results/experiment_results.csv", row.names = FALSE)

summary_stats <- results %>%
    group_by(Method, Variables, Samples) %>%
    summarise(
        across(c(F1, Accuracy, SHD, MSE, Time, Misoriented), mean, na.rm = TRUE),
        .groups = "drop"
    )

# Performance Summary (for manuscript or SI)
performance_summary <- summary_stats %>%
    group_by(Method) %>%
    summarise(
        Avg_F1 = round(mean(F1, na.rm = TRUE), 3),
        Avg_Accuracy = round(mean(Accuracy, na.rm = TRUE), 3),
        Avg_Misoriented = round(mean(Misoriented, na.rm = TRUE), 1),
        Avg_Time = round(mean(Time, na.rm = TRUE), 3),
        .groups = "drop"
    )
print(performance_summary)

# Informative message
if (performance_summary$Avg_F1[performance_summary$Method == "Proposed Method"] > 
    performance_summary$Avg_F1[performance_summary$Method == "LiNGAM"]) {
    cat("Proposed Method outperforms LiNGAM on nonlinear data (Avg F1: ",
        performance_summary$Avg_F1[performance_summary$Method == "Proposed Method"], ")\n")
}

# -------------------------------------------
# Generate Result Plots
# -------------------------------------------
plot_combined_results(summary_stats)
cat("Experiment completed! See 'plots/' for figures.\n")
