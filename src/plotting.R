# =============================================================================
# Causal Discovery Project: Plotting Functions
# =============================================================================

library(ggplot2)
library(gridExtra)

plot_combined_results <- function(results) {
    metrics <- c("F1_Score_dir", "Graph_Accuracy", "SHD", "MSE", "Time", "Misoriented")

    for (metric in metrics) {
        p1 <- ggplot(results, aes(x = Samples, y = .data[[metric]],
                                  linetype = Method, shape = Method, color = Method, group = Method)) +
            geom_line(linewidth = 0.5) +
            geom_point(size = 1.5) +
            facet_wrap(~ Variables, scales = "free_y", nrow = 1,
                       labeller = labeller(Variables = function(x) paste("Variables:", x))) +
            labs(title = paste(metric, "vs. Sample Size"),
                 x = "Sample Size", y = metric) +
            theme_minimal() +
            theme(
                plot.title = element_text(size = 14, face = "bold"),
                axis.title = element_text(size = 12),
                legend.position = "none",
                strip.text = element_text(size = 10, face = "bold")
            ) +
            scale_linetype_manual(name = "Method",
                                 values = c("LiNGAM" = "solid", "Proposed Method" = "dashed")) +
            scale_shape_manual(name = "Method",
                              values = c("LiNGAM" = 16, "Proposed Method" = 17)) +
            scale_color_manual(name = "Method",
                              values = c("LiNGAM" = "#E31A1C", "Proposed Method" = "#1F78B4"))

        p2 <- ggplot(results, aes(x = Variables, y = .data[[metric]],
                                  linetype = Method, shape = Method, color = Method, group = Method)) +
            geom_line(linewidth = 0.5) +
            geom_point(size = 1.5) +
            facet_wrap(~ Samples, scales = "free_y", nrow = 1,
                       labeller = labeller(Samples = function(x) paste("Samples:", x))) +
            labs(title = paste(metric, "vs. Number of Variables"),
                 x = "Number of Variables", y = metric) +
            theme_minimal() +
            theme(
                plot.title = element_text(size = 14, face = "bold"),
                axis.title = element_text(size = 12),
                legend.title = element_text(size = 12, face = "bold"),
                legend.text = element_text(size = 10),
                legend.position = "bottom",
                legend.key.width = unit(2, "cm"),
                strip.text = element_text(size = 10, face = "bold")
            ) +
            scale_linetype_manual(name = "Method",
                                 values = c("LiNGAM" = "solid", "Proposed Method" = "dashed")) +
            scale_shape_manual(name = "Method",
                              values = c("LiNGAM" = 16, "Proposed Method" = 17)) +
            scale_color_manual(name = "Method",
                              values = c("LiNGAM" = "#E31A1C", "Proposed Method" = "#1F78B4"))

        # Force y-axis to start at 0 for F1_Score_dir and Graph_Accuracy
        if (metric %in% c("F1_Score_dir", "Graph_Accuracy")) {
            p1 <- p1 + scale_y_continuous(limits = c(0, NA), expand = c(0, 0))
            p2 <- p2 + scale_y_continuous(limits = c(0, NA), expand = c(0, 0))
        }

        combined_plot <- grid.arrange(p1, p2, nrow = 2)

        if (!dir.exists("plots")) {
            dir.create("plots", recursive = TRUE)
        }

        ggsave(paste0("plots/combined_line_plot_", metric, ".png"),
                combined_plot, width = 14, height = 10, dpi = 300)

        ggsave(paste0("plots/combined_line_plot_", metric, ".pdf"),
                combined_plot, width = 14, height = 10)
    }
}
