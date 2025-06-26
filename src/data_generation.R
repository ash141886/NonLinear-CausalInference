generate_data <- function(n_vars, n_samples = 1000, nonlinearity = 0.3, 
                          sparsity = 0.3, noise_level = 0.1) {
    set.seed(123)
    data <- matrix(0, nrow = n_samples, ncol = n_vars)
    
    # First variable is independent
    data[, 1] <- rnorm(n_samples)
    
    # Generate subsequent variables with causal dependencies
    for (i in 2:n_vars) {
        parents <- 1:(i - 1)
        
        # Select parent subset based on sparsity (ensure at least one parent)
        parent_subset <- parents[rbinom(length(parents), 1, 1 - sparsity) == 1]
        
        # If no parents selected due to sparsity, force at least one parent
        if (length(parent_subset) == 0) {
            parent_subset <- sample(parents, 1)
        }
        
        # Generate values based on parents
        if (length(parent_subset) > 0) {
            # Initialize contribution from parents
            parent_contribution <- rep(0, n_samples)
            
            for (p in parent_subset) {
                if (runif(1) < nonlinearity) {
                    # Apply nonlinear function
                    f <- sample(c(sin, cos, function(x) x^2,
                                  function(x) log(abs(x) + 1), 
                                  function(x) tanh(x)), 1)[[1]]
                    
                    # Handle potential numerical issues
                    transformed_values <- tryCatch({
                        f(data[, p])
                    }, error = function(e) {
                        # Fallback to linear if transformation fails
                        data[, p]
                    })
                    
                    # Check for infinite or NaN values
                    if (any(!is.finite(transformed_values))) {
                        transformed_values <- data[, p]  # Fallback to linear
                    }
                    
                    parent_contribution <- parent_contribution + transformed_values
                } else {
                    # Linear relationship
                    parent_contribution <- parent_contribution + data[, p]
                }
            }
            
            data[, i] <- parent_contribution
        } else {
            # No parents (shouldn't happen now, but safety check)
            data[, i] <- rnorm(n_samples)
        }
        
        # Add noise (with adaptive scaling)
        noise_scale <- noise_level * (1 + 0.1 * abs(data[, i]))
        data[, i] <- data[, i] + rnorm(n_samples, 0, noise_scale)
    }
    
    # Convert to data frame with proper column names
    colnames(data) <- paste0("V", 1:n_vars)
    return(as.data.frame(scale(data)))
}
