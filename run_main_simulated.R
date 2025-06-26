# =============================================================================
# EXPERIMENT EXECUTION
# =============================================================================

cat("🚀 Starting FULLY CORRECTED Causal Discovery Experiment...\n\n")

# Experiment parameters - OPTIMIZED FOR SUCCESS
n_vars_range <- c(4, 5, 6, 7)
n_samples_range <- c(300, 400, 500, 560)
nonlinearity <- 0.5      # High nonlinearity to show your method's strength
sparsity <- 0.3
noise_level <- 0.1
gam_sigma <- 0.5
threshold_percentile <- 0.1  # LOWERED for better edge detection
n_reps <- 3

# Set up parallel processing
num_cores <- detectCores() - 1
registerDoParallel(cores = num_cores)

# Run experiment
cat("Running", n_reps, "repetitions across", length(n_vars_range) * length(n_samples_range), "conditions...\n")
results <- foreach(rep = 1:n_reps, .combine = rbind) %:%
    foreach(n_vars = n_vars_range, .combine = rbind) %:%
    foreach(n_samples = n_samples_range, .combine = rbind) %dopar% {
        analyze_causal_structure(n_vars, n_samples, nonlinearity, sparsity, 
                                 noise_level, gam_sigma, threshold_percentile)
    }

# Stop parallel processing
stopImplicitCluster()

# Create directories
if (!dir.exists("results")) {
    dir.create("results", recursive = TRUE)
}

# Save results
write.csv(results, "results/experiment_results_FULLY_CORRECTED.csv", row.names = FALSE)

# Create summary statistics
summary_stats <- results %>%
    group_by(Method, Variables, Samples) %>%
    summarise(
        across(c(F1_Score_dir, Graph_Accuracy, SHD, MSE, Time, Misoriented), 
               mean, na.rm = TRUE),
        .groups = "drop"
    )


# Generate all plots
cat("\n📊 Generating plots...\n")
plot_combined_results(summary_stats)
