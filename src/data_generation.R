# =============================================================================
# Causal Discovery Project: Data Generation Functions
# =============================================================================

generate_data <- function(n_vars, n_samples = 1000, nonlinearity = 0.3, 
                          sparsity = 0.3, noise_level = 0.1) {
    set.seed(123)
    data <- matrix(0, nrow = n_samples, ncol = n_vars)
    data[, 1] <- rnorm(n_samples)
    for (i in 2:n_vars) {
        parents <- 1:(i - 1)
        parent_subset <- parents[rbinom(length(parents), 1, 1 - sparsity) == 1]
        if (length(parent_subset) == 0) {
            parent_subset <- sample(parents, 1)
        }
        if (length(parent_subset) > 0) {
            parent_contribution <- rep(0, n_samples)
            for (p in parent_subset) {
                if (runif(1) < nonlinearity) {
                    f <- sample(c(sin, cos, function(x) x^2,
                                  function(x) log(abs(x) + 1), 
                                  function(x) tanh(x)), 1)[[1]]
                    transformed_values <- tryCatch({
                        f(data[, p])
                    }, error = function(e) {
                        data[, p]
                    })
                    if (any(!is.finite(transformed_values))) {
                        transformed_values <- data[, p]
                    }
                    parent_contribution <- parent_contribution + transformed_values
                } else {
                    parent_contribution <- parent_contribution + data[, p]
                }
            }
            data[, i] <- parent_contribution
        } else {
            data[, i] <- rnorm(n_samples)
        }
        noise_scale <- noise_level * (1 + 0.1 * abs(data[, i]))
        data[, i] <- data[, i] + rnorm(n_samples, 0, noise_scale)
    }
    colnames(data) <- paste0("V", 1:n_vars)
    return(as.data.frame(scale(data)))
}
