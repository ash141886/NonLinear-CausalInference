# Your existing calculate_dag_metrics_improved function
calculate_dag_metrics_improved <- function(estimated_dag, true_dag) {
    n <- nrow(true_dag)
    TP_dir <- 0; FP_dir <- 0; FN_dir <- 0; misoriented <- 0; TN_dir <- 0
    for (i in 1:n) {
        for (j in 1:n) {
            if (i == j) next
            if (true_dag[i, j] == 1) {
                if (estimated_dag[i, j] == 1) {
                    TP_dir <- TP_dir + 1
                } else if (estimated_dag[j, i] == 1) {
                    misoriented <- misoriented + 1
                    FN_dir <- FN_dir + 1
                } else {
                    FN_dir <- FN_dir + 1
                }
            } else {
                if (estimated_dag[i, j] == 1) {
                    FP_dir <- FP_dir + 1
                } else {
                    TN_dir <- TN_dir + 1
                }
            }
        }
    }
    Precision_dir <- ifelse((TP_dir + FP_dir) > 0, TP_dir / (TP_dir + FP_dir), 0)
    Recall_dir <- ifelse((TP_dir + FN_dir) > 0, TP_dir / (TP_dir + FN_dir), 0)
    F1_Score_dir <- ifelse((Precision_dir + Recall_dir) > 0,
                           2 * Precision_dir * Recall_dir / (Precision_dir + Recall_dir), 0)
    Graph_Accuracy <- (TP_dir + TN_dir) / (n * (n - 1))
    SHD <- sum(abs(estimated_dag - true_dag))
    MSE <- mean((true_dag - estimated_dag)[row(true_dag) != col(true_dag)]^2)
    MAE <- mean(abs(true_dag - estimated_dag)[row(true_dag) != col(true_dag)])
    c(Precision_dir = Precision_dir, Recall_dir = Recall_dir, F1_Score_dir = F1_Score_dir,
      Graph_Accuracy = Graph_Accuracy, Misoriented = misoriented, SHD = SHD,
      MSE = MSE, MAE = MAE)
}
