# Your existing analyze_causal_structure function
analyze_causal_structure <- function(n_vars, n_samples, nonlinearity = 0.3,
                                     sparsity = 0.3, noise_level = 0.1,
                                     gam_sigma = 1, threshold_percentile = 0.3) {
    cat("Running analysis for", n_vars, "variables and", n_samples, "samples...\n")
    data <- generate_data(n_vars, n_samples, nonlinearity, sparsity, noise_level)
    true_dag <- matrix(0, n_vars, n_vars)
    true_dag[lower.tri(true_dag)] <- 1
    methods <- c("Proposed Method", "LiNGAM")
    results <- data.frame(Method = methods, Time = NA, stringsAsFactors = FALSE)
    for (method in methods) {
        cat("Running", method, "method...\n")
        start_time <- Sys.time()
        if (method == "Proposed Method") {
            dag <- tryCatch({
                run_proposed_method(data, gam_sigma, threshold_percentile)
            }, error = function(e) {
                cat("Error in Proposed Method:", e$message, "\n")
                matrix(0, n_vars, n_vars)
            })
            end_time <- Sys.time()
            time_taken <- as.numeric(difftime(end_time, start_time, units = "secs"))
        } else if (method == "LiNGAM") {
            lingam_result <- run_lingam_algorithm(data)
            dag <- lingam_result$dag
            time_taken <- lingam_result$time
        }
        results$Time[results$Method == method] <- time_taken
        metrics <- calculate_dag_metrics_improved(dag, true_dag)
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

# --- MODIFIED PLOTTING FUNCTION ---
plot_combined_results <- function(results) {
    metrics <- c("F1_Score_dir", "Graph_Accuracy", "SHD", "MSE", "Time", "Misoriented")
    
    for (metric in metrics) {
        # Plot vs Sample Size (Fixed Variable Counts)
        p1 <- ggplot(results, aes(x = Samples, y = .data[[metric]],
                                  linetype = Method, shape = Method, color = Method, group = Method)) +
            geom_line(size = 0.7) +
            geom_point(size = 2.5) +
            facet_wrap(~ Variables, scales = "free_y", nrow = 1,
                       labeller = labeller(Variables = function(x) paste(x))) +
            labs(title = paste(metric, "vs. Sample Size (Fixed Variable Counts)"),
                 x = "Sample Size", y = metric) +
            theme_minimal() +
            theme(
                plot.title = element_text(size = 14, face = "bold"),
                axis.title = element_text(size = 12),
                legend.position = "none",  # Remove legend from individual plots
                strip.text = element_text(size = 10, face = "bold")
            ) +
            scale_linetype_manual(values = c("LiNGAM" = "solid", "Proposed Method" = "dashed")) +
            scale_shape_manual(values = c("LiNGAM" = 16, "Proposed Method" = 17)) +  # Circle and triangle
            scale_color_manual(values = c("LiNGAM" = "red", "Proposed Method" = "blue"))

        # Plot vs Number of Variables (Fixed Sample Sizes)
        p2 <- ggplot(results, aes(x = Variables, y = .data[[metric]],
                                  linetype = Method, shape = Method, color = Method, group = Method)) +
            geom_line(size = 0.7) +
            geom_point(size = 2.5) +
            facet_wrap(~ Samples, scales = "free_y", nrow = 1,
                       labeller = labeller(Samples = function(x) paste(x))) +
            labs(title = paste(metric, "vs. Number of Variables (Fixed Sample Sizes)"),
                 x = "Number of Variables", y = metric) +
            theme_minimal() +
            theme(
                plot.title = element_text(size = 14, face = "bold"),
                axis.title = element_text(size = 12),
                legend.position = "none",  # Remove legend from individual plots
                strip.text = element_text(size = 10, face = "bold")
            ) +
            scale_linetype_manual(values = c("LiNGAM" = "solid", "Proposed Method" = "dashed")) +
            scale_shape_manual(values = c("LiNGAM" = 16, "Proposed Method" = 17)) +  # Circle and triangle
            scale_color_manual(values = c("LiNGAM" = "red", "Proposed Method" = "blue"))
        
        # Extract legend from one of the plots
        legend <- get_legend(
            ggplot(results, aes(x = Samples, y = .data[[metric]],
                               linetype = Method, shape = Method, color = Method, group = Method)) +
                geom_line(size = 0.7) +
                geom_point(size = 2.5) +
                scale_linetype_manual(values = c("LiNGAM" = "solid", "Proposed Method" = "dashed")) +
                scale_shape_manual(values = c("LiNGAM" = 16, "Proposed Method" = 17)) +
                scale_color_manual(values = c("LiNGAM" = "red", "Proposed Method" = "blue")) +
                theme_minimal() +
                theme(
                    legend.title = element_text(size = 12, face = "bold"),
                    legend.text = element_text(size = 10),
                    legend.position = "bottom"
                ) +
                guides(linetype = guide_legend(title = "Method"),
                       shape = guide_legend(title = "Method"),
                       color = guide_legend(title = "Method"))
        )
        
        # Combine plots with shared legend at bottom
        combined_plot <- grid.arrange(
            arrangeGrob(p1, p2, nrow = 2),
            legend,
            nrow = 2,
            heights = c(10, 1)  # Adjust heights: plots take most space, legend takes small space
        )
        
        # Save plots
        if (!dir.exists("plots")) {
            dir.create("plots", recursive = TRUE)
        }
        
        ggsave(paste0("plots/combined_line_plot_", metric, ".png"), 
               combined_plot, width = 14, height = 10, dpi = 300)
        ggsave(paste0("plots/combined_line_plot_", metric, ".pdf"), 
               combined_plot, width = 14, height = 10)
    }
}
