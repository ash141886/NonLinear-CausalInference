# =============================================================================
# SECTION 4: BASELINE METHOD - LINEAR NON-GAUSSIAN ACYCLIC MODEL
# =============================================================================

lingam_discovery <- function(data) {
  start_time <- Sys.time()
  result <- tryCatch({
    ica_result <- fastICA(as.matrix(data), n.comp = ncol(data),
                          alg.typ = "parallel", fun = "logcosh", method = "C", verbose = FALSE)
    W <- ica_result$W
    A <- solve(W)
    n <- ncol(data)
    B <- A; diag(B) <- 0
    perm <- 1:n
    for (i in 1:n) {
      scores <- colSums(abs(B))
      next_var <- which.min(scores)
      perm[i] <- next_var
      B[next_var, ] <- 0; B[, next_var] <- 0
    }
    B <- solve(W); B <- B[perm, perm]; diag(B) <- 0; B[abs(B) < 0.1] <- 0
    dag <- (abs(B) > 0) * 1; dag[upper.tri(dag)] <- 0
    inv_perm <- order(perm); dag <- dag[inv_perm, inv_perm]
    list(dag = dag, time = as.numeric(difftime(Sys.time(), start_time, units = "secs")))
  }, error = function(e) {
    n <- ncol(data)
    list(dag = matrix(0, n, n), time = as.numeric(difftime(Sys.time(), start_time, units = "secs")))
  })
  result
}
