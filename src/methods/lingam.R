
# REAL LiNGAM implementation using fastICA
run_lingam_algorithm <- function(data) {
    start_time <- Sys.time()
    res <- tryCatch({
        ica_result <- fastICA(data, n.comp = ncol(data))
        W <- ica_result$W
        B <- solve(W) %*% diag(apply(W, 2, max))
        diag(B) <- 0
        B[abs(B) < 0.01] <- 0
        end_time <- Sys.time()
        time_taken <- as.numeric(difftime(end_time, start_time, units = "secs"))
        
        # Convert to 0/1 matrix (not logical) and ensure proper dimensions
        dag_matrix <- matrix(as.numeric(B != 0), nrow = nrow(B), ncol = ncol(B))
        
        list(dag = dag_matrix, time = time_taken)
    }, error = function(e) {
        cat("Error in LiNGAM:", e$message, "\n")
        end_time <- Sys.time()
        time_taken <- as.numeric(difftime(end_time, start_time, units = "secs"))
        list(dag = matrix(0, ncol(data), ncol(data)), time = time_taken)
    })
    res
}
