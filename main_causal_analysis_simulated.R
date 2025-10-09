# =============================================================================
# Nonlinear Causal Discovery: Additive Models with Kernel-based Independence Testing 
# =============================================================================

suppressPackageStartupMessages({
  library(mgcv)
  library(fastICA)
  library(ggplot2)
  library(gridExtra)
  library(grid)
  library(dplyr)
 
})

# ---- Time formatter (for final runtime message) ------------------------------
.format_secs <- function(s) {
  s <- as.numeric(s)
  if (is.na(s) || !is.finite(s)) return("NA")
  h <- floor(s / 3600); m <- floor((s %% 3600) / 60); sec <- round(s %% 60)
  paste0(if (h > 0) sprintf("%dh ", h) else "",
         if (m > 0) sprintf("%dm ", m) else "",
         sprintf("%ds", sec))
}

# =============================================================================
# SECTION 1: KERNEL-BASED INDEPENDENCE TESTING
# =============================================================================

rbf_kernel <- function(x, sigma = 1) {
  # x can be a vector or a matrix; we compute pairwise Euclidean distances
  dist_matrix <- as.matrix(dist(x))
  exp(-dist_matrix^2 / (2 * sigma^2))
}

hsic <- function(x, y, sigma = NULL) {
  n <- length(x)

  if (is.null(sigma)) {
    dist_x <- as.matrix(dist(x))
    dist_y <- as.matrix(dist(y))
    sigma_x <- median(dist_x[lower.tri(dist_x)])
    sigma_y <- median(dist_y[lower.tri(dist_y)])
    if (!is.finite(sigma_x) || sigma_x == 0) sigma_x <- 1
    if (!is.finite(sigma_y) || sigma_y == 0) sigma_y <- 1
    sigma <- sqrt(sigma_x * sigma_y)
  }

  Kx <- rbf_kernel(as.matrix(x), sigma)
  Ky <- rbf_kernel(as.matrix(y), sigma)
  H  <- diag(n) - matrix(1 / n, n, n)

  sum((H %*% Kx %*% H) * Ky) / (n - 1)^2
}

# =============================================================================
# SECTION 2: DATA GENERATION MECHANISM
# =============================================================================

generate_nonlinear_sem_data <- function(
  n_vars, n_samples = 1000, nonlinearity = 0.5, sparsity = 0.3,
  noise_level = 0.5, seed = 123
) {
  set.seed(seed)

  true_dag <- matrix(0, n_vars, n_vars)
  data     <- matrix(0, nrow = n_samples, ncol = n_vars)

  # Non-Gaussian root
  data[, 1] <- rexp(n_samples, rate = 1) - 1

  for (i in 2:n_vars) {
    potential_parents <- 1:(i - 1)

    if (length(potential_parents) == 1) {
      selected_parents <- potential_parents
    } else {
      n_parents <- max(1, rbinom(1, length(potential_parents), 1 - sparsity))
      selected_parents <- sample(potential_parents, n_parents)
    }

    for (parent in selected_parents) true_dag[parent, i] <- 1

    effects <- numeric(n_samples)
    for (parent in selected_parents) {
      if (runif(1) < nonlinearity) {
        transform_func <- sample(list(
          function(x) sin(2 * x),
          function(x) cos(2 * x),
          function(x) x^2,
          function(x) tanh(2 * x),
          function(x) exp(-abs(x))
        ), 1)[[1]]
        effects <- effects + transform_func(data[, parent])
      } else {
        coefficient <- runif(1, 0.5, 1.5) * sample(c(-1, 1), 1)
        effects <- effects + coefficient * data[, parent]
      }
    }

    # Non-Gaussian noise types
    nd <- sample(1:3, 1)
    noise <- if (nd == 1) {
      rexp(n_samples, rate = 1 / noise_level) - noise_level
    } else if (nd == 2) {
      rt(n_samples, df = 5) * noise_level
    } else {
      runif(n_samples, -1, 1) * noise_level
    }

    data[, i] <- effects + noise
  }

  colnames(data) <- paste0("V", 1:n_vars)
  list(data = as.data.frame(scale(data)), true_dag = true_dag)
}

# =============================================================================
# SECTION 3: CAUSAL DISCOVERY ALGORITHM
# =============================================================================

discover_causal_structure <- function(data, sigma = NULL, threshold_percentile = 85) {
  n_vars    <- ncol(data)
  n_samples <- nrow(data)
  score_matrix <- matrix(0, n_vars, n_vars)

  for (i in 1:n_vars) {
    for (j in 1:n_vars) {
      if (i == j) next

      predictors <- setdiff(1:n_vars, c(i, j))
      if (length(predictors) > 0) {
        k_basis <- min(5, floor(n_samples / 40))
        k_basis <- max(k_basis, 3)  # ensure >= 3 to keep GAM stable
        formula_str <- paste0(
          "V", i, " ~ ",
          paste0("s(V", predictors, ", k = ", k_basis, ")", collapse = " + ")
        )

        model <- tryCatch(
          suppressWarnings(
            gam(as.formula(formula_str), data = data, method = "REML", gamma = 1.4)
          ),
          error = function(e) {
            lm(
              as.formula(paste0("V", i, " ~ ", paste0("V", predictors, collapse = " + "))),
              data = data
            )
          }
        )
        residuals_i <- residuals(model)
      } else {
        residuals_i <- data[[paste0("V", i)]] - mean(data[[paste0("V", i)]])
      }

      score_matrix[i, j] <- hsic(residuals_i, data[[paste0("V", j)]], sigma)
    }
  }

  positive_scores <- score_matrix[score_matrix > 0]
  threshold <- if (length(positive_scores) > 0) {
    suppressWarnings(quantile(positive_scores, threshold_percentile / 100, na.rm = TRUE))
  } else {
    Inf  # no positive scores -> no edges selected
  }

  adjacency <- matrix(0, n_vars, n_vars)
  for (i in 1:n_vars) {
    for (j in 1:n_vars) {
      if (i != j && is.finite(threshold) && score_matrix[i, j] > threshold) {
        # edge j -> i (i depends on j)
        adjacency[j, i] <- 1
      }
    }
  }

  enforce_acyclicity(adjacency, score_matrix)
}

enforce_acyclicity <- function(adjacency, weights) {
  n <- nrow(adjacency)
  iteration <- 0
  max_iterations <- n * n

  repeat {
    cycle <- detect_cycle(adjacency)
    if (is.null(cycle) || iteration >= max_iterations) break

    min_weight <- Inf
    weakest_edge <- NULL
    for (k in seq_along(cycle)) {
      from_node <- cycle[k]
      to_node   <- cycle[ifelse(k == length(cycle), 1, k + 1)]
      if (adjacency[from_node, to_node] == 1) {
        edge_weight <- weights[to_node, from_node]  # weight of to->from in score_matrix
        if (edge_weight < min_weight) {
          min_weight  <- edge_weight
          weakest_edge <- c(from_node, to_node)
        }
      }
    }

    if (!is.null(weakest_edge)) {
      adjacency[weakest_edge[1], weakest_edge[2]] <- 0
    }
    iteration <- iteration + 1
  }
  adjacency
}

detect_cycle <- function(adjacency) {
  n <- nrow(adjacency)
  color <- rep("white", n)
  parent <- rep(NA_integer_, n)

  dfs <- function(v) {
    color[v] <<- "gray"
    children <- which(adjacency[v, ] == 1)
    for (u in children) {
      if (color[u] == "gray") {
        cycle_path <- c(u)
        current <- v
        while (current != u && !is.na(parent[current])) {
          cycle_path <- c(cycle_path, current)
          current <- parent[current]
        }
        return(cycle_path)
      }
      if (color[u] == "white") {
        parent[u] <<- v
        cycle <- dfs(u)
        if (!is.null(cycle)) return(cycle)
      }
    }
    color[v] <<- "black"
    NULL
  }

  for (v in 1:n) {
    if (color[v] == "white") {
      cycle <- dfs(v)
      if (!is.null(cycle)) return(cycle)
    }
  }
  NULL
}

# =============================================================================
# SECTION 4: BASELINE METHOD - LINEAR NON-GAUSSIAN ACYCLIC MODEL (LiNGAM-ish)
# =============================================================================

lingam_discovery <- function(data) {
  start_time <- Sys.time()
  result <- tryCatch({
    ica_result <- fastICA(
      as.matrix(data),
      n.comp  = ncol(data),
      alg.typ = "parallel",
      fun     = "logcosh",
      method  = "C",
      verbose = FALSE
    )

    W <- ica_result$W
    A <- tryCatch(solve(W), error = function(e) MASS::ginv(W))
    n <- ncol(data)

    # crude permutation scoring to approximate ordering
    B <- A
    diag(B) <- 0
    perm <- integer(n)
    Btmp <- B

    for (i in 1:n) {
      scores <- colSums(abs(Btmp))
      next_var <- which.min(scores)
      perm[i] <- next_var
      Btmp[next_var, ] <- 0
      Btmp[, next_var] <- 0
    }

    B <- tryCatch(solve(W), error = function(e) MASS::ginv(W))
    B <- B[perm, perm]
    diag(B) <- 0
    B[abs(B) < 0.1] <- 0

    dag <- (abs(B) > 0) * 1
    dag[upper.tri(dag)] <- 0

    inv_perm <- order(perm)
    dag <- dag[inv_perm, inv_perm]

    list(
      dag  = dag,
      time = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    )
  }, error = function(e) {
    n <- ncol(data)
    list(
      dag  = matrix(0, n, n),
      time = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    )
  })
  result
}

# =============================================================================
# SECTION 5: PERFORMANCE EVALUATION
# =============================================================================

evaluate_performance <- function(estimated, truth) {
  n <- nrow(truth)

  true_positives  <- sum(estimated == 1 & truth == 1)
  false_positives <- sum(estimated == 1 & truth == 0)
  false_negatives <- sum(estimated == 0 & truth == 1)
  true_negatives  <- sum(estimated == 0 & truth == 0) - n  # exclude diagonals

  # Misoriented edges (count edges pointing the wrong way where only one direction is in truth)
  misoriented <- 0
  for (i in 1:n) {
    for (j in 1:n) if (i < j) {
      if (truth[i, j] == 1 && estimated[j, i] == 1 && estimated[i, j] == 0) {
        misoriented <- misoriented + 1
      } else if (truth[j, i] == 1 && estimated[i, j] == 1 && estimated[j, i] == 0) {
        misoriented <- misoriented + 1
      }
    }
  }

  precision <- ifelse((true_positives + false_positives) > 0,
                      true_positives / (true_positives + false_positives), 0)
  recall    <- ifelse((true_positives + false_negatives) > 0,
                      true_positives / (true_positives + false_negatives), 0)
  f1_score  <- ifelse((precision + recall) > 0,
                      2 * precision * recall / (precision + recall), 0)
  accuracy  <- (true_positives + true_negatives) / (n * (n - 1))
  structural_hamming <- sum(abs(estimated - truth))
  mse <- mean((estimated - truth)^2)

  c(
    Precision_dir  = precision,
    Recall_dir     = recall,
    F1_Score_dir   = f1_score,
    Graph_Accuracy = accuracy,
    Misoriented    = misoriented,
    SHD            = structural_hamming,
    MSE            = mse
  )
}

# =============================================================================
# SECTION 6: SIMULATION FRAMEWORK (simple console counter, no packages)
# =============================================================================

conduct_simulation_study <- function(
  variable_sizes, sample_sizes, replications = 50,
  nonlinearity = 0.5, sparsity = 0.3, noise_level = 0.5
) {
  complete_results <- data.frame()
  t0 <- Sys.time()

  total_iters <- length(variable_sizes) * length(sample_sizes) * replications
  iter <- 0

  for (n_vars in variable_sizes) {
    for (n_samples in sample_sizes) {
      for (rep in 1:replications) {
        iter <- iter + 1
        cat(sprintf("%d/%d Vars=%d, Samples=%d\n", iter, total_iters, n_vars, n_samples))

        data_generation <- generate_nonlinear_sem_data(
          n_vars, n_samples, nonlinearity, sparsity, noise_level,
          seed = rep * 1000 + n_vars * 10 + n_samples
        )
        data <- data_generation$data
        true_graph <- data_generation$true_dag

        # If no true edges, skip metrics collection for this iteration
        if (sum(true_graph) == 0) next

        # Proposed method
        time_start <- Sys.time()
        estimated_proposed <- tryCatch(
          discover_causal_structure(data, sigma = NULL, threshold_percentile = 85),
          error = function(e) matrix(0, n_vars, n_vars)
        )
        time_proposed <- as.numeric(difftime(Sys.time(), time_start, units = "secs"))
        metrics_proposed <- evaluate_performance(estimated_proposed, true_graph)

        # LiNGAM baseline
        lingam_output  <- lingam_discovery(data)
        metrics_lingam <- evaluate_performance(lingam_output$dag, true_graph)

        iteration_results <- rbind(
          data.frame(
            Method      = "Proposed Method",
            Variables   = n_vars,
            Samples     = n_samples,
            Replication = rep,
            Time        = time_proposed,
            t(metrics_proposed)
          ),
          data.frame(
            Method      = "LiNGAM",
            Variables   = n_vars,
            Samples     = n_samples,
            Replication = rep,
            Time        = lingam_output$time,
            t(metrics_lingam)
          )
        )

        complete_results <- rbind(complete_results, iteration_results)
      }
    }
  }

  cat(sprintf(
    "Completed %d iterations in %s.\n",
    total_iters, .format_secs(difftime(Sys.time(), t0, units = "secs"))
  ))

  aggregated_results <- complete_results %>%
    dplyr::group_by(Method, Variables, Samples) %>%
    dplyr::summarise(
      dplyr::across(
        c(Precision_dir, Recall_dir, F1_Score_dir,
          Graph_Accuracy, Misoriented, SHD, MSE, Time),
        \(x) mean(x, na.rm = TRUE)
      ),
      .groups = "drop"
    )

  list(detailed = complete_results, summary = aggregated_results)
}

# =============================================================================
# SECTION 7: VISUALIZATION (force display in RStudio + save files)
# =============================================================================

generate_figures <- function(
  results,
  base_size = 9,
  title_size = 10,
  axis_title_size = 9,
  axis_text_size = 8,
  strip_size = 9,
  legend_title_size = 9,
  legend_text_size = 8
) {
  if (!dir.exists("figures")) dir.create("figures", recursive = TRUE)

  method_colors <- c("LiNGAM" = "#E31A1C", "Proposed Method" = "#1F78B4")
  method_shapes <- c("LiNGAM" = 16, "Proposed Method" = 17)
  method_lines  <- c("LiNGAM" = "solid", "Proposed Method" = "dashed")

  metrics_to_plot <- c(
    "Precision_dir", "F1_Score_dir", "Graph_Accuracy",
    "Misoriented", "SHD", "MSE", "Time"
  )

  plain_theme <- ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      text             = element_text(face = "plain"),
      plot.title       = element_text(size = title_size, face = "plain", hjust = 0.5),
      axis.title       = element_text(size = axis_title_size, face = "plain"),
      axis.text        = element_text(size = axis_text_size, face = "plain"),
      legend.title     = element_text(size = legend_title_size, face = "plain"),
      legend.text      = element_text(size = legend_text_size, face = "plain"),
      legend.position  = "bottom",
      legend.key.width = unit(1.2, "cm"),
      strip.text       = element_text(size = strip_size, face = "plain"),
      panel.grid.minor = element_blank()
    )

  for (metric in metrics_to_plot) {
    p_samples <- ggplot(
      results$summary,
      aes(x = Samples, y = .data[[metric]],
          color = Method, shape = Method, linetype = Method)
    ) +
      geom_line(linewidth = 0.5) +
      geom_point(size = 1.3) +
      facet_wrap(~ Variables, scales = "free_y", nrow = 1) +
      labs(title = paste(metric, "vs. Sample Size"),
           x = "Sample Size", y = metric) +
      plain_theme +
      theme(legend.position = "none") +
      scale_color_manual(values = method_colors) +
      scale_shape_manual(values = method_shapes) +
      scale_linetype_manual(values = method_lines)

    p_variables <- ggplot(
      results$summary,
      aes(x = Variables, y = .data[[metric]],
          color = Method, shape = Method, linetype = Method)
    ) +
      geom_line(linewidth = 0.5) +
      geom_point(size = 1.3) +
      facet_wrap(~ Samples, scales = "free_y", nrow = 1) +
      labs(title = paste(metric, "vs. Number of Variables"),
           x = "Number of Variables", y = metric) +
      plain_theme +
      scale_color_manual(name = "Method", values = method_colors) +
      scale_shape_manual(name = "Method", values = method_shapes) +
      scale_linetype_manual(name = "Method", values = method_lines) +
      scale_x_continuous(breaks = unique(results$summary$Variables))

    if (metric %in% c("Precision_dir", "F1_Score_dir", "Graph_Accuracy")) {
      p_samples   <- p_samples   + scale_y_continuous(limits = c(0, 1))
      p_variables <- p_variables + scale_y_continuous(limits = c(0, 1))
    }

    g <- gridExtra::grid.arrange(p_samples, p_variables, nrow = 2)

    ggsave(filename = file.path("figures", paste0(metric, ".pdf")), plot = g,
           width = 12, height = 8, dpi = 300)
    ggsave(filename = file.path("figures", paste0(metric, ".png")), plot = g,
           width = 12, height = 8, dpi = 300)
  }
}

# =============================================================================
# SECTION 8: MAIN EXECUTION
# =============================================================================

simulation_results <- conduct_simulation_study(
  variable_sizes = c(4, 5, 6, 7),
  sample_sizes   = c(200, 250, 300, 350),
  replications   = 1,
  nonlinearity   = 0.5,
  sparsity       = 0.3,
  noise_level    = 0.5
)

generate_figures(simulation_results)

saveRDS(simulation_results, file = "causal_discovery_results.rds")
write.csv(simulation_results$summary, file = "summary_statistics.csv", row.names = FALSE)
