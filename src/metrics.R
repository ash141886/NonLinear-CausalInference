# -------------------------------------------
# DAG Metric Calculation
# -------------------------------------------
calculate_dag_metrics <- function(estimated_dag, true_dag) {
    n <- nrow(true_dag)
    TP <- FP <- FN <- misoriented <- TN <- 0
    for (i in 1:n) {
        for (j in 1:n) {
            if (i == j) next
            if (true_dag[i, j] == 1) {
                if (estimated_dag[i, j] == 1) {
                    TP <- TP + 1
                } else if (estimated_dag[j, i] == 1) {
                    misoriented <- misoriented + 1
                    FN <- FN + 1
                } else {
                    FN <- FN + 1
                }
            } else {
                if (estimated_dag[i, j] == 1) {
                    FP <- FP + 1
                } else {
                    TN <- TN + 1
                }
            }
        }
    }
    Precision <- ifelse((TP + FP) > 0, TP / (TP + FP), 0)
    Recall <- ifelse((TP + FN) > 0, TP / (TP + FN), 0)
    F1 <- ifelse((Precision + Recall) > 0, 2 * Precision * Recall / (Precision + Recall), 0)
    Accuracy <- (TP + TN) / (n * (n - 1))
    SHD <- sum(abs(estimated_dag - true_dag))
    MSE <- mean((true_dag - estimated_dag)[row(true_dag) != col(true_dag)]^2)
    MAE <- mean(abs(true_dag - estimated_dag)[row(true_dag) != col(true_dag)])
    c(
        Precision = Precision, Recall = Recall, F1 = F1, 
        Accuracy = Accuracy, Misoriented = misoriented,
        SHD = SHD, MSE = MSE, MAE = MAE
    )
}
