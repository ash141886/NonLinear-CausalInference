# =============================================================================
# Causal Discovery Simulation Experiment
# =============================================================================

# --- 1. Load required packages ---
library(dplyr)
library(foreach)
library(doParallel)

# Source your existing files
source("src/data_generation.R")
source("src/methods/proposed_method.R")
source("src/methods/lingam.R")
source("src/metrics.R")
source("src/analysis.R")

# --- 2. Experimental Parameters ---
n_vars_range <- c(4, 5, 6, 7)
n_samples_range <- c(300, 400, 500, 560)
nonlinearity <- 0.5
sparsity <- 0.3
noise_level <- 0.1
gam_sigma <- 0.5
threshold_percentile <- 0.1
n_reps <- 3

# --- 3. Parallel Processing Setup ---
num_cores <- max(1, detectCores() - 1)
registerDoParallel(cores = num_cores)

# --- 4. Main Experiment Loop ---
results <- foreach(rep = 1:n_reps, .combine = rbind) %:%
    foreach(n_vars = n_vars_range, .combine = rbind) %:%
    foreach(n_samples = n_samples_range, .combine = rbind) %dopar% {
        analyze_causal_structure(
            n_vars      = n_vars,
            n_samples   = n_samples,
            nonlinearity= nonlinearity,
            sparsity    = sparsity,
            noise_level = noise_level,
            gam_sigma   = gam_sigma,
            threshold_percentile = threshold_percentile
        )
    }

stopImplicitCluster()

# --- 5. Results Saving ---
if (!dir.exists("results")) dir.create("results", recursive = TRUE)
write.csv(results, "results/experiment_results.csv", row.names = FALSE)

# --- 6. Summary Statistics ---
summary_stats <- results %>%
    group_by(Method, Variables, Samples) %>%
    summarise(
        F1_Score_dir = mean(F1_Score_dir, na.rm = TRUE),
        Graph_Accuracy = mean(Graph_Accuracy, na.rm = TRUE),
        SHD = mean(SHD, na.rm = TRUE),
        MSE = mean(MSE, na.rm = TRUE),
        Time = mean(Time, na.rm = TRUE),
        Misoriented = mean(Misoriented, na.rm = TRUE),
        .groups = "drop"
    )

# --- 7. Visualization ---
plot_combined_results(summary_stats)

# --- End of script ---
