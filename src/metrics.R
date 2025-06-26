# -------------------------------------------
# DAG Metric Calculation
# -------------------------------------------
calculate_dag_metrics <- function(estimated_dag, true_dag) {
    n <- nrow(true_dag)
    TP <- FP <- FN <- misoriented <- TN <- 0
    for (i in 1:n) {
        for (j in 1:n) {
            if (i == j) next
            if (true_dag[i, j] == 1) {
                if (estimated_dag[i, j] == 1) {
                    TP <- TP + 1
                } else if (estimated_dag[j, i] == 1) {
                    misoriented <- misoriented + 1
                    FN <- FN + 1
                } else {
                    FN <- FN + 1
                }
            } else {
                if (estimated_dag[i, j] == 1) {
                    FP <- FP + 1
                } else {
                    TN <- TN + 1
                }
            }
        }
    }
    Precision <- ifelse((TP + FP) > 0, TP / (TP + FP), 0)
    Recall <- ifelse((TP + FN) > 0, TP / (TP + FN), 0)
    F1 <- ifelse((Precision + Recall) > 0, 2 * Precision * Recall / (Precision + Recall), 0)
    Accuracy <- (TP + TN) / (n * (n - 1))
    SHD <- sum(abs(estimated_dag - true_dag))
    MSE <- mean((true_dag - estimated_dag)[row(true_dag) != col(true_dag)]^2)
    MAE <- mean(abs(true_dag - estimated_dag)[row(true_dag) != col(true_dag)])
    c(
        Precision = Precision, Recall = Recall, F1 = F1, 
        Accuracy = Accuracy, Misoriented = misoriented,
        SHD = SHD, MSE = MSE, MAE = MAE
    )
}

# -------------------------------------------
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
