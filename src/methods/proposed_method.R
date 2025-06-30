# =============================================================================
# CORE FUNCTIONS
# =============================================================================

# RBF Kernel function
rbf_kernel <- function(x, sigma = 1) {
    dist_matrix <- as.matrix(dist(x))
    exp(-dist_matrix^2 / (2 * sigma^2))
}

# HSIC function
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

# Transitive closure function
transitive_closure <- function(graph) {
    n <- nrow(graph)
    closure <- graph
    for (k in 1:n) {
        for (i in 1:n) {
            for (j in 1:n) {
                closure[i, j] <- closure[i, j] | (closure[i, k] & closure[k, j])
            }
        }
    }
    closure
}



# Proposed Method function
run_proposed_method <- function(data, sigma = 1, threshold_percentile = 0.3) {
    n_vars <- ncol(data)
    hsic_matrix <- matrix(0, n_vars, n_vars)
    residuals <- vector("list", n_vars)
    
    # Step 1: Fit GAM for each variable and extract residuals
    for (i in 1:n_vars) {
        formula <- as.formula(paste0("V", i, " ~ ", 
                                     paste0("s(V", (1:n_vars)[-i], ", k=5)", collapse = " + ")))
        gam_model <- tryCatch({
            gam(formula, data = data, method = "REML")
        }, error = function(e) {
            cat("Error fitting GAM for variable", i, ":", e$message, "\n")
            return(NULL)
        })
        if (!is.null(gam_model)) {
            residuals[[i]] <- residuals(gam_model)
        } else {
            residuals[[i]] <- rnorm(nrow(data))
        }
    }
    
    # Step 2: Calculate HSIC between residuals
    for (i in 1:(n_vars - 1)) {
        for (j in (i + 1):n_vars) {
            hsic_matrix[i, j] <- tryCatch({
                hsic(residuals[[i]], residuals[[j]], sigma)
            }, error = function(e) {
                cat("Error in HSIC for variables", i, "and", j, ":", e$message, "\n")
                return(0)
            })
        }
    }
    
    # Step 3: Make symmetric and apply threshold
    hsic_matrix <- hsic_matrix + t(hsic_matrix)
    diag(hsic_matrix) <- 1
    threshold <- quantile(hsic_matrix[lower.tri(hsic_matrix)], threshold_percentile, na.rm = TRUE)
    proposed_method_dag <- hsic_matrix > threshold
    
    proposed_method_dag[lower.tri(proposed_method_dag)] <- 0
    
    diag(proposed_method_dag) <- 0
    
    # Convert to numeric matrix 
    result_dag <- matrix(as.numeric(proposed_method_dag), nrow = n_vars, ncol = n_vars)
    transitive_closure(result_dag)
}
