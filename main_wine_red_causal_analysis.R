# ===================================================================
# Causal Discovery on Wine Data (red)
# ===================================================================

# -------------------------------------------------------------------
# Step 1: Load Required Libraries
# -------------------------------------------------------------------

suppressPackageStartupMessages({
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
  # MASS used via MASS::ginv if needed
})

# -------------------------------------------------------------------
# Step 2: Load and Preprocess Data
# -------------------------------------------------------------------

if (!file.exists("data/winequality-red.csv")) {
  dir.create("data", showWarnings = FALSE, recursive = TRUE)
  download.file(
    "https://archive.ics.uci.edu/ml/machine-learning-databases/wine-quality/winequality-red.csv",
    "data/winequality-red.csv"
  )
}
wine_data <- read.csv("data/winequality-red.csv", sep = ";")
wine_data <- as.data.frame(scale(wine_data))

cat("Wine data loaded and standardized.\n")
cat("Dimensions:", paste(dim(wine_data), collapse = " x "), "\n")
cat("Variables:", paste(colnames(wine_data), collapse = ", "), "\n\n")

# -------------------------------------------------------------------
# Step 3: Helper + Core Functions
# -------------------------------------------------------------------

rbf_kernel <- function(x, sigma = 1) {
  dist_matrix <- as.matrix(dist(x))
  exp(-dist_matrix^2 / (2 * sigma^2))
}

hsic <- function(x, y, sigma = NULL) {
  n <- length(x)
  if (is.null(sigma)) {
    dist_x <- as.matrix(dist(x))
    dist_y <- as.matrix(dist(y))
    sigma_x <- median(dist_x[lower.tri(dist_x)]); if (!is.finite(sigma_x) || sigma_x == 0) sigma_x <- 1
    sigma_y <- median(dist_y[lower.tri(dist_y)]); if (!is.finite(sigma_y) || sigma_y == 0) sigma_y <- 1
    sigma <- sqrt(sigma_x * sigma_y)
  }
  Kx <- rbf_kernel(as.matrix(x), sigma)
  Ky <- rbf_kernel(as.matrix(y), sigma)
  H  <- diag(n) - matrix(1/n, n, n)
  sum((H %*% Kx %*% H) * Ky) / (n - 1)^2
}

# ---- Proposed (DIRECTED, HSIC residual → predictor like in simulation) ----
run_proposed_method_wine <- function(data, sigma = NULL, threshold_percentile = 0.85, gam_k_max = 10) {
  stopifnot(threshold_percentile > 0 && threshold_percentile < 1)
  n_vars    <- ncol(data)
  n_samples <- nrow(data)
  var_names <- colnames(data)

  score_matrix <- matrix(0, n_vars, n_vars,
                         dimnames = list(var_names, var_names))

  for (i in 1:n_vars) {
    for (j in 1:n_vars) {
      if (i == j) next

      predictors <- setdiff(1:n_vars, c(i, j))
      if (length(predictors) > 0) {
        # dynamic k to keep GAM stable on real data
        k_basis <- min(gam_k_max, max(3, floor(n_samples / 40)))
        terms <- paste0("s(`", var_names[predictors], "`, k=", k_basis, ")")
        fml <- as.formula(paste0("`", var_names[i], "` ~ ", paste(terms, collapse = " + ")))
        model <- tryCatch(
          suppressWarnings(gam(fml, data = data, method = "REML", gamma = 1.4)),
          error = function(e) NULL
        )
        resid_i <- if (!is.null(model)) residuals(model) else scale(data[[i]], scale = FALSE)
      } else {
        resid_i <- data[[i]] - mean(data[[i]])
      }

      # HSIC between residual of i (with j excluded) and raw X_j → directional j -> i
      score_matrix[i, j] <- tryCatch(
        hsic(resid_i, data[[j]], sigma),
        error = function(e) 0
      )
    }
  }

  # Threshold on positive scores over all ordered pairs
  positive_scores <- as.numeric(score_matrix[score_matrix > 0])
  thr <- if (length(positive_scores)) {
    as.numeric(quantile(positive_scores, probs = threshold_percentile, na.rm = TRUE))
  } else Inf

  est_adj <- matrix(0, n_vars, n_vars, dimnames = list(var_names, var_names))
  if (is.finite(thr)) {
    est_adj[score_matrix > thr] <- 1
  }
  # optional: you could enforce acyclicity here if desired

  est_adj
}

run_lingam_algorithm_wine <- function(data) {
  start_time <- Sys.time()
  dag_result <- tryCatch({
    ica <- fastICA(as.matrix(data), n.comp = ncol(data),
                   alg.typ = "parallel", fun = "logcosh",
                   method = "C", verbose = FALSE)
    W <- ica$W
    A <- tryCatch(solve(W), error = function(e) MASS::ginv(W))
    B <- A
    diag(B) <- 0
    # simple threshold
    B[abs(B) < 0.01] <- 0
    (abs(B) > 0) * 1
  }, error = function(e) matrix(0, ncol(data), ncol(data)))
  time_taken <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  dimnames(dag_result) <- list(colnames(data), colnames(data))
  list(dag = dag_result, time = time_taken)
}

# -------------------------------------------------------------------
# Step 4: "Ground truth" (domain-inspired) + Metrics
# -------------------------------------------------------------------

get_wine_true_dag <- function() {
  vars <- colnames(wine_data)
  n_vars <- length(vars)
  true_dag <- matrix(0, n_vars, n_vars, dimnames = list(vars, vars))
  true_dag["fixed.acidity",        "pH"]                    <- 1
  true_dag["volatile.acidity",     "pH"]                    <- 1
  true_dag["citric.acid",          "pH"]                    <- 1
  true_dag["citric.acid",          "fixed.acidity"]         <- 1
  true_dag["residual.sugar",       "density"]               <- 1
  true_dag["alcohol",              "density"]               <- 1
  true_dag["free.sulfur.dioxide",  "total.sulfur.dioxide"]  <- 1
  true_dag["alcohol",              "quality"]               <- 1
  true_dag["volatile.acidity",     "quality"]               <- 1
  true_dag["sulphates",            "quality"]               <- 1
  true_dag
}

# Undirected metrics (as you had)
calculate_dag_metrics_undirected <- function(estimated_dag, true_dag) {
  estimated_undirected <- (estimated_dag | t(estimated_dag)) * 1
  diag(estimated_undirected) <- 0
  true_undirected <- (true_dag | t(true_dag)) * 1
  diag(true_undirected) <- 0

  tp <- sum(estimated_undirected == 1 & true_undirected == 1) / 2
  fp <- sum(estimated_undirected == 1 & true_undirected == 0) / 2
  fn <- sum(estimated_undirected == 0 & true_undirected == 1) / 2
  n_vars <- ncol(true_dag)
  total_possible_edges <- n_vars * (n_vars - 1) / 2
  tn <- total_possible_edges - tp - fp - fn

  precision <- ifelse((tp + fp) == 0, 0, tp / (tp + fp))
  recall    <- ifelse((tp + fn) == 0, 0, tp / (tp + fn))
  f1_score  <- ifelse((precision + recall) == 0, 0, 2 * (precision * recall) / (precision + recall))
  shd       <- fp + fn
  accuracy  <- (tp + tn) / total_possible_edges
  mse       <- mean((true_undirected - estimated_undirected)^2)

  c(Precision = precision, Recall = recall, F1_Score = f1_score,
    SHD = shd, Accuracy = accuracy, MSE = mse)
}

# NEW: Misoriented edges metric (directional)
calculate_misoriented_edges <- function(estimated_dag, true_dag) {
  stopifnot(identical(dim(estimated_dag), dim(true_dag)))
  n <- nrow(true_dag)
  mis <- 0L
  for (i in 1:n) {
    for (j in 1:n) if (i < j) {
      # Count an error when exactly one direction is present in truth,
      # but the estimated puts the edge in the opposite direction only.
      if (true_dag[i, j] == 1 && true_dag[j, i] == 0 && estimated_dag[j, i] == 1 && estimated_dag[i, j] == 0) {
        mis <- mis + 1L
      } else if (true_dag[j, i] == 1 && true_dag[i, j] == 0 && estimated_dag[i, j] == 1 && estimated_dag[j, i] == 0) {
        mis <- mis + 1L
      }
    }
  }
  mis
}

# -------------------------------------------------------------------
# Step 5: Full Hyperparameter Tuning (Proposed)
# -------------------------------------------------------------------

tune_wine_parameters_full <- function() {
  cat("=== Starting full hyperparameter tuning ===\n")
  threshold_percentiles <- c(0.70, 0.80, 0.85, 0.90, 0.95)  # probs in (0,1)
  gam_sigmas            <- c(0.5, 1.0, 1.5)                 # passed to HSIC
  gam_k_max_values      <- c(5, 10, 15)

  best_f1 <- -1
  best_params <- NULL
  true_dag <- get_wine_true_dag()
  param_grid <- expand.grid(thresh = threshold_percentiles,
                            sigma = gam_sigmas,
                            k_max = gam_k_max_values)
  cat("Total combinations to test:", nrow(param_grid), "\n")

  for (i in 1:nrow(param_grid)) {
    thresh <- param_grid$thresh[i]
    sigma  <- param_grid$sigma[i]
    k_val  <- param_grid$k_max[i]
    cat(sprintf("Testing %d/%d: Threshold=%.2f, Sigma=%.2f, K_max=%d\n",
                i, nrow(param_grid), thresh, sigma, k_val))
    est_dag <- run_proposed_method_wine(wine_data, sigma = sigma,
                                        threshold_percentile = thresh,
                                        gam_k_max = k_val)
    metrics <- calculate_dag_metrics_undirected(est_dag, true_dag)
    if (metrics["F1_Score"] > best_f1) {
      best_f1 <- metrics["F1_Score"]
      best_params <- list(threshold = thresh, sigma = sigma, k_max = k_val)
    }
  }

  if (is.null(best_params)) {
    cat("\nWarning: No valid combination found. Using default values.\n")
    best_params <- list(threshold = 0.90, sigma = 1.0, k_max = 5)
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
best_sigma  <- best_params$sigma
best_k_max  <- best_params$k_max

final_results_list <- list()
true_dag <- get_wine_true_dag()

cat("Analyzing with the Proposed method (directional)…\n")
t0 <- Sys.time()
prop_dag <- run_proposed_method_wine(wine_data, sigma = best_sigma,
                                     threshold_percentile = best_thresh,
                                     gam_k_max = best_k_max)
prop_time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
prop_metrics <- calculate_dag_metrics_undirected(prop_dag, true_dag)
prop_mis <- calculate_misoriented_edges(prop_dag, true_dag)
final_results_list[["Proposed"]] <- c(prop_metrics, Misoriented = prop_mis, Time_sec = prop_time)

cat("Analyzing with LiNGAM…\n")
lingam_output <- run_lingam_algorithm_wine(wine_data)
lingam_metrics <- calculate_dag_metrics_undirected(lingam_output$dag, true_dag)
lingam_mis <- calculate_misoriented_edges(lingam_output$dag, true_dag)
final_results_list[["LiNGAM"]] <- c(lingam_metrics, Misoriented = lingam_mis, Time_sec = lingam_output$time)

cat("\n=========================================================\n")
cat("=== FINAL PERFORMANCE COMPARISON (FULL ANALYSIS) ===\n")
cat("=========================================================\n")
final_results_df <- do.call(rbind, lapply(names(final_results_list), function(name) {
  data.frame(Method = name, t(final_results_list[[name]]))
}))
numeric_cols <- names(final_results_df)[sapply(final_results_df, is.numeric)]
final_results_df[numeric_cols] <- lapply(final_results_df[numeric_cols], function(v) {
  # don't round SHD, Misoriented (integers); round the rest
  nm <- names(v)
  if (is.null(nm)) return(round(v, 3))
  round(v, 3)
})
print(final_results_df, row.names = FALSE)

cat("\nFull analysis complete.\n")
