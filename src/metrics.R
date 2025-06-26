# =============================================================================
# Causal Discovery Project: Metrics and Analysis Functions
# =============================================================================

# ----------------------------------------------------------------------------- 
# Metrics Calculation Function
# -----------------------------------------------------------------------------
calculate_dag_metrics_improved <- function(estimated_dag, true_dag) {
    n <- nrow(true_dag)
    TP_dir <- 0; FP_dir <- 0; FN_dir <- 0; misoriented <- 0; TN_dir <- 0
    for (i in 1:n) {
        for (j in 1:n) {
            if (i == j) next
            if (true_dag[i, j] == 1) {
                if (estimated_dag[i, j] == 1) {
                    TP_dir <- TP_dir + 1
                } else if (estimated_dag[j, i] == 1) {
                    misoriented <- misoriented + 1
                    FN_dir <- FN_dir + 1
                } else {
                    FN_dir <- FN_dir + 1
                }
            } else {
                if (estimated_dag[i, j] == 1) {
                    FP_dir <- FP_dir + 1
                } else {
                    TN_dir <- TN_dir + 1
                }
            }
        }
    }
    Precision_dir <- ifelse((TP_dir + FP_dir) > 0, TP_dir / (TP_dir + FP_dir), 0)
    Recall_dir <- ifelse((TP_dir + FN_dir) > 0, TP_dir / (TP_dir + FN_dir), 0)
    F1_Score_dir <- ifelse((Precision_dir + Recall_dir) > 0,
                           2 * Precision_dir * Recall_dir / (Precision_dir + Recall_dir), 0)
    Graph_Accuracy <- (TP_dir + TN_dir) / (n * (n - 1))
    SHD <- sum(abs(estimated_dag - true_dag))
    MSE <- mean((true_dag - estimated_dag)[row(true_dag) != col(true_dag)]^2)
    MAE <- mean(abs(true_dag - estimated_dag)[row(true_dag) != col(true_dag)])
    c(Precision_dir = Precision_dir, Recall_dir = Recall_dir, F1_Score_dir = F1_Score_dir,
      Graph_Accuracy = Graph_Accuracy, Misoriented = misoriented, SHD = SHD,
      MSE = MSE, MAE = MAE)
}

# ----------------------------------------------------------------------------- 
# Main Experiment Analysis Function
# -----------------------------------------------------------------------------
analyze_causal_structure <- function(n_vars, n_samples, nonlinearity = 0.3,
                                     sparsity = 0.3, noise_level = 0.1,
                                     gam_sigma = 1, threshold_percentile = 0.3) {
    cat("Running analysis for", n_vars, "variables and", n_samples, "samples...\n")
    data <- generate_data(n_vars, n_samples, nonlinearity, sparsity, noise_level)
    true_dag <- matrix(0, n_vars, n_vars)
    true_dag[upper.tri(true_dag)] <- 1
    methods <- c("Proposed Method", "LiNGAM")
    results <- data.frame(Method = methods, Time = NA, stringsAsFactors = FALSE)
    for (method in methods) {
        cat("Running", method, "method...\n")
        start_time <- Sys.time()
        if (method == "Proposed Method") {
            dag <- tryCatch({
                run_proposed_method(data, sigma = gam_sigma, threshold_percentile = threshold_percentile)
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
