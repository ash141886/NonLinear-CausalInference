generate_data <- function(n_vars, n_samples = 1000, nonlinearity = 0.3, 
                          sparsity = 0.3, noise_level = 0.1) {
    set.seed(123)
    data <- matrix(0, nrow = n_samples, ncol = n_vars)
    data[, 1] <- rnorm(n_samples)
    for (i in 2:n_vars) {
        parents <- 1:(i - 1)
        parent_subset <- parents[rbinom(length(parents), 1, 1 - sparsity) == 1]
        if (length(parent_subset) > 0) {
            data[, i] <- rowSums(sapply(parent_subset, function(p) {
                if (runif(1) < nonlinearity) {
                    f <- sample(c(sin, cos, exp, function(x) x^2,
                                  function(x) log(abs(x) + 1), tanh), 1)[[1]]
                    f(data[, p])
                } else {
                    data[, p]
                }
            }))
        }
        data[, i] <- data[, i] + rnorm(n_samples, 0, noise_level * (1 + abs(data[, i])))
    }
    return(as.data.frame(scale(data)))
}
