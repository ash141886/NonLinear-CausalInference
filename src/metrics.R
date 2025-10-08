# =============================================================================
# SECTION 5: PERFORMANCE EVALUATION
# =============================================================================

evaluate_performance <- function(estimated, truth) {
  n <- nrow(truth)
  true_positives <- sum(estimated == 1 & truth == 1)
  false_positives <- sum(estimated == 1 & truth == 0)
  false_negatives <- sum(estimated == 0 & truth == 1)
  true_negatives <- sum(estimated == 0 & truth == 0) - n
  misoriented <- 0
  for (i in 1:n) for (j in 1:n) if (i < j) {
    if (truth[i, j] == 1 && estimated[j, i] == 1 && estimated[i, j] == 0) misoriented <- misoriented + 1
    else if (truth[j, i] == 1 && estimated[i, j] == 1 && estimated[j, i] == 0) misoriented <- misoriented + 1
  }
  precision <- ifelse((true_positives + false_positives) > 0, true_positives / (true_positives + false_positives), 0)
  recall    <- ifelse((true_positives + false_negatives) > 0, true_positives / (true_positives + false_negatives), 0)
  f1_score  <- ifelse((precision + recall) > 0, 2 * precision * recall / (precision + recall), 0)
  accuracy  <- (true_positives + true_negatives) / (n * (n - 1))
  structural_hamming <- sum(abs(estimated - truth))
  mse <- mean((estimated - truth)^2)
  c(Precision_dir = precision, Recall_dir = recall, F1_Score_dir = f1_score,
    Graph_Accuracy = accuracy, Misoriented = misoriented, SHD = structural_hamming, MSE = mse)
}

# =============================================================================
# SECTION 6: SIMULATION FRAMEWORK
# =============================================================================

conduct_simulation_study <- function(variable_sizes, sample_sizes, replications = 50,
                                     nonlinearity = 0.5, sparsity = 0.3, noise_level = 0.5) {
  complete_results <- data.frame()
  for (n_vars in variable_sizes) {
    for (n_samples in sample_sizes) {
      for (rep in 1:replications) {
        data_generation <- generate_nonlinear_sem_data(
          n_vars, n_samples, nonlinearity, sparsity, noise_level,
          seed = rep * 1000 + n_vars * 10 + n_samples
        )
        data <- data_generation$data
        true_graph <- data_generation$true_dag
        if (sum(true_graph) == 0) next
        time_start <- Sys.time()
        estimated_proposed <- tryCatch(
          discover_causal_structure(data, sigma = NULL, threshold_percentile = 85),
          error = function(e) matrix(0, n_vars, n_vars)
        )
        time_proposed <- as.numeric(difftime(Sys.time(), time_start, units = "secs"))
        metrics_proposed <- evaluate_performance(estimated_proposed, true_graph)
        lingam_output <- lingam_discovery(data)
        metrics_lingam <- evaluate_performance(lingam_output$dag, true_graph)
        iteration_results <- rbind(
          data.frame(Method = "Proposed Method", Variables = n_vars, Samples = n_samples,
                     Replication = rep, Time = time_proposed, t(metrics_proposed)),
          data.frame(Method = "LiNGAM", Variables = n_vars, Samples = n_samples,
                     Replication = rep, Time = lingam_output$time, t(metrics_lingam))
        )
        complete_results <- rbind(complete_results, iteration_results)
      }
    }
  }
  aggregated_results <- complete_results %>%
    group_by(Method, Variables, Samples) %>%
    summarise(across(c(Precision_dir, Recall_dir, F1_Score_dir,
                       Graph_Accuracy, Misoriented, SHD, MSE, Time),
                     mean, na.rm = TRUE), .groups = 'drop')
  list(detailed = complete_results, summary = aggregated_results)
}
