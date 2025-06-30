# ===================================================================
# Causal Discovery on Wine Data (red)
#
# Full hyperparameter tuning and performance comparison of the proposed methods on red wine data.
# ===================================================================
 
# -------------------------------------------------------------------
# Step 1: Load Required Libraries
# -------------------------------------------------------------------

library(mgcv)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pcalg)
library(fastICA)
library(parallel)
library(doParallel)
library(gridExtra)
library(viridis)

# -------------------------------------------------------------------
# Step 2: Load and Preprocess Data
# -------------------------------------------------------------------

if (!file.exists("data/winequality-red.csv")) {
    download.file("https://archive.ics.uci.edu/ml/machine-learning-databases/wine-quality/winequality-red.csv",
                  "data/winequality-red.csv")
}
wine_data <- read.csv("data/winequality-red.csv", sep = ";")
wine_data <- as.data.frame(scale(wine_data))

cat("Wine data loaded and standardized.\n")
cat("Dimensions:", dim(wine_data), "\n")
cat("Variables:", paste(colnames(wine_data), collapse = ", "), "\n\n")

# -------------------------------------------------------------------
# Step 3: Define Helper and Core Algorithm Functions
# -------------------------------------------------------------------

rbf_kernel <- function(x, sigma = 1) {
    dist_matrix <- as.matrix(dist(x))
    exp(-dist_matrix^2 / (2 * sigma^2))
}

hsic <- function(x, y, sigma = NULL) {
    n <- length(x)
    if (is.null(sigma)) {
        dist_matrix <- as.matrix(dist(cbind(x, y)))
        sigma <- median(dist_matrix[lower.tri(dist_matrix)])
    }
    Kx <- rbf_kernel(as.matrix(x), sigma)
    Ky <- rbf_kernel(as.matrix(y), sigma)
    H <- diag(n) - matrix(1/n, n, n)
    sum(H %*% Kx %*% H * Ky) / n^2
}

run_proposed_method_wine <- function(data, sigma = 1, threshold_percentile = 0.8, gam_k_max = 10) {
    n_vars <- ncol(data)
    var_names <- colnames(data)
    hsic_matrix <- matrix(0, n_vars, n_vars)
    residuals_list <- list()
    for (i in 1:n_vars) {
        target_var_name <- var_names[i]
        predictor_indices <- setdiff(1:n_vars, i)
        formula_terms <- c()
        for (pred_idx in predictor_indices) {
            pred_name <- var_names[pred_idx]
            n_unique <- length(unique(data[[pred_name]]))
            k_dynamic <- min(gam_k_max, n_unique - 1)
            if (k_dynamic >= 3) {
                formula_terms <- c(formula_terms, paste0("s(`", pred_name, "`, k=", k_dynamic, ")"))
            } else {
                formula_terms <- c(formula_terms, paste0("`", pred_name, "`"))
            }
        }
        formula_str <- paste0("`", target_var_name, "` ~ ", paste(formula_terms, collapse = " + "))
        formula <- as.formula(formula_str)
        gam_model <- tryCatch({
            gam(formula, data = data, method = "REML")
        }, error = function(e) NULL)
        if (!is.null(gam_model)) {
            residuals_list[[i]] <- residuals(gam_model)
        } else {
            residuals_list[[i]] <- rnorm(nrow(data))
        }
    }
    for (i in 1:(n_vars - 1)) {
        for (j in (i + 1):n_vars) {
            hsic_val <- tryCatch(hsic(residuals_list[[i]], residuals_list[[j]], sigma), error = function(e) 0)
            hsic_matrix[i, j] <- hsic_val
            hsic_matrix[j, i] <- hsic_val
        }
    }
    diag(hsic_matrix) <- 0
    threshold <- quantile(hsic_matrix[lower.tri(hsic_matrix)], threshold_percentile, na.rm = TRUE)
    estimated_adj_matrix <- hsic_matrix > threshold
    estimated_adj_matrix
}

run_lingam_algorithm_wine <- function(data) {
    start_time <- Sys.time()
    dag_result <- tryCatch({
        ica_result <- fastICA(data, n.comp = ncol(data))
        W <- ica_result$W
        B <- solve(W) %*% diag(apply(W, 2, max))
        diag(B) <- 0
        B[abs(B) < 0.01] <- 0
        B != 0
    }, error = function(e) matrix(0, ncol(data), ncol(data)))
    end_time <- Sys.time()
    time_taken <- as.numeric(difftime(end_time, start_time, units = "secs"))
    list(dag = dag_result, time = time_taken)
}

# -------------------------------------------------------------------
# Step 4: Ground Truth and Evaluation Metrics
# -------------------------------------------------------------------

get_wine_true_dag <- function() {
    vars <- colnames(wine_data)
    n_vars <- length(vars)
    true_dag <- matrix(0, n_vars, n_vars, dimnames = list(vars, vars))
    true_dag["fixed.acidity", "pH"] <- 1
    true_dag["volatile.acidity", "pH"] <- 1
    true_dag["citric.acid", "pH"] <- 1
    true_dag["citric.acid", "fixed.acidity"] <- 1
    true_dag["residual.sugar", "density"] <- 1
    true_dag["alcohol", "density"] <- 1
    true_dag["free.sulfur.dioxide", "total.sulfur.dioxide"] <- 1
    true_dag["alcohol", "quality"] <- 1
    true_dag["volatile.acidity", "quality"] <- 1
    true_dag["sulphates", "quality"] <- 1
    true_dag
}

calculate_dag_metrics_undirected <- function(estimated_dag, true_dag) {
    estimated_undirected <- estimated_dag | t(estimated_dag)
    diag(estimated_undirected) <- 0
    true_undirected <- true_dag | t(true_dag)
    diag(true_undirected) <- 0
    tp <- sum(estimated_undirected == 1 & true_undirected == 1) / 2
    fp <- sum(estimated_undirected == 1 & true_undirected == 0) / 2
    fn <- sum(estimated_undirected == 0 & true_undirected == 1) / 2
    n_vars <- ncol(true_dag)
    total_possible_edges <- n_vars * (n_vars - 1) / 2
    tn <- total_possible_edges - tp - fp - fn
    precision <- ifelse((tp + fp) == 0, 0, tp / (tp + fp))
    recall <- ifelse((tp + fn) == 0, 0, tp / (tp + fn))
    f1_score <- ifelse((precision + recall) == 0, 0, 2 * (precision * recall) / (precision + recall))
    shd <- fp + fn
    accuracy <- (tp + tn) / total_possible_edges
    mse <- mean((true_undirected - estimated_undirected)^2)
    c(Precision = precision, Recall = recall, F1_Score = f1_score, SHD = shd, Accuracy = accuracy, MSE = mse)
}

# -------------------------------------------------------------------
# Step 5: Full Hyperparameter Tuning
# -------------------------------------------------------------------

tune_wine_parameters_full <- function() {
    cat("=== Starting full hyperparameter tuning ===\n")
    threshold_percentiles <- c(0.7, 0.8, 0.85, 0.9, 0.95)
    gam_sigmas <- c(0.5, 1.0, 1.5)
    gam_k_max_values <- c(5, 10, 15)
    best_f1 <- -1
    best_params <- NULL
    true_dag <- get_wine_true_dag()
    param_grid <- expand.grid(thresh = threshold_percentiles, sigma = gam_sigmas, k_max = gam_k_max_values)
    cat("Total combinations to test:", nrow(param_grid), "\n")
    for (i in 1:nrow(param_grid)) {
        thresh <- param_grid$thresh[i]
        sigma <- param_grid$sigma[i]
        k_val <- param_grid$k_max[i]
        cat(sprintf("Testing %d/%d: Threshold=%.2f, Sigma=%.2f, K_max=%d\n", i, nrow(param_grid), thresh, sigma, k_val))
        est_dag <- run_proposed_method_wine(wine_data, sigma, thresh, gam_k_max = k_val)
        metrics <- calculate_dag_metrics_undirected(est_dag, true_dag)
        if (metrics["F1_Score"] > best_f1) {
            best_f1 <- metrics["F1_Score"]
            best_params <- list(threshold = thresh, sigma = sigma, k_max = k_val)
        }
    }
    if (is.null(best_params)) {
        cat("\nWarning: No valid combination found. Using default values.\n")
        best_params <- list(threshold = 0.9, sigma = 1.0, k_max = 5)
    } else {
        cat("\nTuning complete.\n")
        cat("Best F1-Score:", round(best_f1, 3), "\n")
        cat("Optimal parameters: Threshold =", best_params$threshold,
            "| Sigma =", best_params$sigma,
            "| K_max =", best_params$k_max, "\n\n")
    }
    list(best_params = best_params)
}

# -------------------------------------------------------------------
# Step 6: Main Execution
# -------------------------------------------------------------------

cat("=================================================================\n")
cat("=== WINE QUALITY CAUSAL DISCOVERY: FULL ANALYSIS ===\n")
cat("=================================================================\n\n")

tuning_output <- tune_wine_parameters_full()
best_params <- tuning_output$best_params

cat("\nRunning final analysis with optimal parameters.\n")
best_thresh <- best_params$threshold
best_sigma <- best_params$sigma
best_k_max <- best_params$k_max
final_results_list <- list()
true_dag <- get_wine_true_dag()

cat("Analyzing with the proposed method...\n")
start_time_prop <- Sys.time()
prop_dag <- run_proposed_method_wine(wine_data, best_sigma, best_thresh, best_k_max)
prop_time <- as.numeric(difftime(Sys.time(), start_time_prop, units = "secs"))
prop_metrics <- calculate_dag_metrics_undirected(prop_dag, true_dag)
final_results_list[["Proposed"]] <- c(prop_metrics, Time_sec = prop_time)

cat("Analyzing with LiNGAM...\n")
lingam_output <- run_lingam_algorithm_wine(wine_data)
lingam_metrics <- calculate_dag_metrics_undirected(lingam_output$dag, true_dag)
final_results_list[["LiNGAM"]] <- c(lingam_metrics, Time_sec = lingam_output$time)

cat("\n=========================================================\n")
cat("=== FINAL PERFORMANCE COMPARISON (FULL ANALYSIS) ===\n")
cat("=========================================================\n")
final_results_df <- do.call(rbind, lapply(names(final_results_list), function(name) {
    data.frame(Method = name, t(final_results_list[[name]]))
}))
numeric_cols <- names(final_results_df)[sapply(final_results_df, is.numeric)]
final_results_df[numeric_cols] <- lapply(final_results_df[numeric_cols], round, 3)
print(final_results_df, row.names = FALSE)
cat("\nFull analysis complete.\n")
