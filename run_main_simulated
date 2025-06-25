source("src/data_generation.R")
source("src/methods/proposed_method.R")
source("src/methods/lingam.R")
source("src/metrics.R")
source("src/analysis.R")

library(foreach)
library(doParallel)

# Experiment parameters
n_vars_range <- c(4, 8, 12, 15)
n_samples_range <- c(400, 800, 1200, 1500)
nonlinearity <- 0.5
sparsity <- 0.3
noise_level <- 0.1
gam_sigma <- 0.5
threshold_percentile <- 0.3
n_reps <- 100

# Set up parallel processing
num_cores <- detectCores() - 1
registerDoParallel(cores = num_cores)

# Run experiment
results <- foreach(rep = 1:n_reps, .combine = rbind) %:%
    foreach(n_vars = n_vars_range, .combine = rbind) %:%
    foreach(n_samples = n_samples_range, .combine = rbind) %dopar% {
        analyze_causal_structure(n_vars, n_samples, nonlinearity, sparsity, 
                                 noise_level, gam_sigma, threshold_percentile)
    }

# Plot results and display plots
metrics <- c("Graph_Accuracy", "SHD", "F1_Score_dir", "Time")
for (metric in metrics) {
    p <- plot_results(results, metric)
    print(p)
}
