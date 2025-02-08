library(dplyr)
library(ggplot2)
library(gridExtra)

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

plot_results <- function(results, metric) {
    p <- ggplot(results, aes(x = Samples, y = .data[[metric]], color = Method, group = Method)) +
        geom_line() +
        geom_point() +
        facet_wrap(~ Variables, scales = "free") +
        labs(title = paste(metric, "vs Sample Size"),
             x = "Sample Size", y = metric) +
        theme_minimal()
    print(p)
    ggsave(paste0("plots/", metric, "_plot.png"), p, width = 10, height = 6)
}
