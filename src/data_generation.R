# =============================================================================
# SECTION 2: DATA GENERATION MECHANISM
# =============================================================================

generate_nonlinear_sem_data <- function(n_vars, n_samples = 1000, 
                                        nonlinearity = 0.5, sparsity = 0.3, 
                                        noise_level = 0.5, seed = 123) {
  set.seed(seed)
  true_dag <- matrix(0, n_vars, n_vars)
  data <- matrix(0, nrow = n_samples, ncol = n_vars)
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
    effects <- rep(0, n_samples)
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
    nd <- sample(1:3, 1)
    noise <- if (nd == 1) rexp(n_samples, rate = 1/noise_level) - noise_level
      else if (nd == 2) rt(n_samples, df = 5) * noise_level
      else runif(n_samples, -1, 1) * noise_level
    data[, i] <- effects + noise
  }
  colnames(data) <- paste0("V", 1:n_vars)
  list(data = as.data.frame(scale(data)), true_dag = true_dag)
}
