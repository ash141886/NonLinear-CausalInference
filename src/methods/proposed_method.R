# =============================================================================
# SECTION 3: CAUSAL DISCOVERY ALGORITHM
# =============================================================================

discover_causal_structure <- function(data, sigma = NULL, threshold_percentile = 85) {
  n_vars <- ncol(data); n_samples <- nrow(data)
  score_matrix <- matrix(0, n_vars, n_vars)
  for (i in 1:n_vars) {
    for (j in 1:n_vars) {
      if (i == j) next
      predictors <- setdiff(1:n_vars, c(i, j))
      if (length(predictors) > 0) {
        k_basis <- min(5, floor(n_samples/40))
        formula_str <- paste0("V", i, " ~ ",
                              paste0("s(V", predictors, ", k=", k_basis, ")", collapse = " + "))
        model <- tryCatch(
          suppressWarnings(gam(as.formula(formula_str), data = data, method = "REML", gamma = 1.4)),
          error = function(e) lm(as.formula(paste0("V", i, " ~ ", paste0("V", predictors, collapse = " + "))), data = data)
        )
        residuals_i <- residuals(model)
      } else {
        residuals_i <- data[[paste0("V", i)]] - mean(data[[paste0("V", i)]])
      }
      score_matrix[i, j] <- hsic(residuals_i, data[[paste0("V", j)]], sigma)
    }
  }
  positive_scores <- score_matrix[score_matrix > 0]
  threshold <- suppressWarnings(quantile(positive_scores, threshold_percentile/100, na.rm = TRUE))
  adjacency <- matrix(0, n_vars, n_vars)
  for (i in 1:n_vars) for (j in 1:n_vars)
    if (i != j && score_matrix[i, j] > threshold) adjacency[j, i] <- 1
  enforce_acyclicity(adjacency, score_matrix)
}

enforce_acyclicity <- function(adjacency, weights) {
  n <- nrow(adjacency); iteration <- 0; max_iterations <- n * n
  repeat {
    cycle <- detect_cycle(adjacency)
    if (is.null(cycle) || iteration >= max_iterations) break
    min_weight <- Inf; weakest_edge <- NULL
    for (k in seq_along(cycle)) {
      from_node <- cycle[k]
      to_node <- cycle[ifelse(k == length(cycle), 1, k + 1)]
      if (adjacency[from_node, to_node] == 1) {
        edge_weight <- weights[to_node, from_node]
        if (edge_weight < min_weight) { min_weight <- edge_weight; weakest_edge <- c(from_node, to_node) }
      }
    }
    if (!is.null(weakest_edge)) adjacency[weakest_edge[1], weakest_edge[2]] <- 0
    iteration <- iteration + 1
  }
  adjacency
}

detect_cycle <- function(adjacency) {
  n <- nrow(adjacency)
  color <- rep("white", n); parent <- rep(NA, n)
  dfs <- function(v) {
    color[v] <<- "gray"
    children <- which(adjacency[v, ] == 1)
    for (u in children) {
      if (color[u] == "gray") {
        cycle_path <- c(u); current <- v
        while (current != u && !is.na(parent[current])) { cycle_path <- c(cycle_path, current); current <- parent[current] }
        return(cycle_path)
      }
      if (color[u] == "white") {
        parent[u] <<- v
        cycle <- dfs(u); if (!is.null(cycle)) return(cycle)
      }
    }
    color[v] <<- "black"; NULL
  }
  for (v in 1:n) if (color[v] == "white") { cycle <- dfs(v); if (!is.null(cycle)) return(cycle) }
  NULL
}
