
# Simulation (single experiment)
# -------------------------------------------
analyze_causal_structure <- function(n_vars, n_samples, nonlinearity = 0.3,
                                     sparsity = 0.3, noise_level = 0.1,
                                     gam_sigma = 1, threshold_percentile = 0.3) {
    data <- generate_data(n_vars, n_samples, nonlinearity, sparsity, noise_level)
    true_dag <- matrix(0, n_vars, n_vars)
    true_dag[upper.tri(true_dag)] <- 1   # assuming topological order
    methods <- c("Proposed Method", "LiNGAM")
    results <- data.frame(Method = methods, Time = NA, stringsAsFactors = FALSE)
    for (method in methods) {
        start_time <- Sys.time()
        if (method == "Proposed Method") {
            dag <- tryCatch({
                run_proposed_method(data, sigma = gam_sigma, threshold_percentile = threshold_percentile)
            }, error = function(e) matrix(0, n_vars, n_vars))
            time_taken <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
        } else if (method == "LiNGAM") {
            lingam_result <- run_lingam_algorithm(data)
            dag <- lingam_result$dag
            time_taken <- lingam_result$time
        }
        results$Time[results$Method == method] <- time_taken
        metrics <- calculate_dag_metrics(dag, true_dag)
        metrics[is.na(metrics)] <- 0
        results[results$Method == method, names(metrics)] <- metrics
    }
    results$Variables <- n_vars
    results$Samples <- n_samples
    results$Nonlinearity <- nonlinearity
    results$Sparsity <- sparsity
    results$Noise <- noise_level
    results$Threshold <- threshold_percentile
    results$gam_sigma <- gam_sigma
    results
}

# -------------------------------------------
# Plotting Utility (For Each Metric) - FIXED VERSION
# -------------------------------------------
plot_combined_results <- function(results) {
    metrics <- c("F1", "Accuracy", "SHD", "MSE", "Time", "Misoriented")
    for (metric in metrics) {
        # First plot: Metric vs Sample Size
        p1 <- ggplot(results, aes(x = Samples, y = .data[[metric]],
                                  linetype = Method, shape = Method, color = Method, group = Method)) +
            geom_line(size = 0.4) +
            geom_point(size = 1.5) +
            facet_wrap(~ Variables, scales = "free_y", nrow = 1,
                       labeller = labeller(Variables = function(x) paste("Variables:", x))) +
            labs(title = paste(metric, "vs. Sample Size"),
                 x = "Sample Size", y = metric) +
            theme_minimal() +
            theme(
                plot.title = element_text(size = 14, face = "bold"),
                axis.title = element_text(size = 12),
                legend.title = element_text(size = 12, face = "bold"),
                legend.text = element_text(size = 10),
                legend.position = "none",  # Remove legend from first plot
                legend.key.width = unit(2, "cm"),
                strip.text = element_text(size = 10, face = "bold")
            ) +
            scale_linetype_manual(name = "Method", 
                                  values = c("LiNGAM" = "solid", "Proposed Method" = "dashed")) +
            scale_shape_manual(name = "Method", 
                               values = c("LiNGAM" = 16, "Proposed Method" = 17)) +
            scale_color_manual(name = "Method", 
                               values = c("LiNGAM" = "red", "Proposed Method" = "blue"))
        
        # Second plot: Metric vs Number of Variables
        p2 <- ggplot(results, aes(x = Variables, y = .data[[metric]],
                                  linetype = Method, shape = Method, color = Method, group = Method)) +
            geom_line(size = 0.4) +
            geom_point(size = 1.5) +
            facet_wrap(~ Samples, scales = "free_y", nrow = 1,
                       labeller = labeller(Samples = function(x) paste("Samples:", x))) +
            labs(title = paste(metric, "vs. Number of Variables"),
                 x = "Number of Variables", y = metric) +
            theme_minimal() +
            theme(
                plot.title = element_text(size = 14, face = "bold"),
                axis.title = element_text(size = 12),
                legend.title = element_text(size = 12, face = "bold"),
                legend.text = element_text(size = 10),
                legend.position = "bottom",  # Keep legend only on second plot
                legend.key.width = unit(2, "cm"),
                strip.text = element_text(size = 10, face = "bold")
            ) +
            scale_linetype_manual(name = "Method", 
                                  values = c("LiNGAM" = "solid", "Proposed Method" = "dashed")) +
            scale_shape_manual(name = "Method", 
                               values = c("LiNGAM" = 16, "Proposed Method" = 17)) +
            scale_color_manual(name = "Method", 
                               values = c("LiNGAM" = "red", "Proposed Method" = "blue"))
        
        # Combine plots with single legend
        combined_plot <- grid.arrange(p1, p2, nrow = 2)
        
        # Save plots
        if (!dir.exists("plots")) dir.create("plots", recursive = TRUE)
        ggsave(paste0("plots/combined_line_plot_", metric, ".png"), combined_plot, width = 14, height = 10, dpi = 300)
        ggsave(paste0("plots/combined_line_plot_", metric, ".pdf"), combined_plot, width = 14, height = 10)
    }
}
