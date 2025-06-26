# =============================================================================
# Causal Discovery Project: Plotting Functions (Single Shared Legend, Clean Console)
# =============================================================================

library(ggplot2)
library(cowplot)

plot_combined_results <- function(results) {
    metrics <- c("F1_Score_dir", "Graph_Accuracy", "SHD", "MSE", "Time", "Misoriented")
    for (metric in metrics) {
        p1 <- ggplot(results, aes(x = Samples, y = .data[[metric]],
                                  linetype = Method, shape = Method, color = Method, group = Method)) +
            geom_line(linewidth = 0.4) +
            geom_point(size = 2) +
            facet_wrap(~ Variables, scales = "free_y", nrow = 1,
                       labeller = labeller(Variables = function(x) paste("Variables:", x))) +
            labs(title = paste(metric, "vs. Sample Size"),
                 x = "Sample Size", y = metric) +
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
            scale_linetype_manual(name = "Method", values = c("LiNGAM" = "solid", "Proposed Method" = "dashed")) +
            scale_shape_manual(name = "Method", values = c("LiNGAM" = 16, "Proposed Method" = 17)) +
            scale_color_manual(name = "Method", values = c("LiNGAM" = "#E31A1C", "Proposed Method" = "#1F78B4"))

        p2 <- ggplot(results, aes(x = Variables, y = .data[[metric]],
                                  linetype = Method, shape = Method, color = Method, group = Method)) +
            geom_line(linewidth = 0.4) +
            geom_point(size = 2) +
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
            scale_linetype_manual(name = "Method", values = c("LiNGAM" = "solid", "Proposed Method" = "dashed")) +
            scale_shape_manual(name = "Method", values = c("LiNGAM" = 16, "Proposed Method" = 17)) +
            scale_color_manual(name = "Method", values = c("LiNGAM" = "#E31A1C", "Proposed Method" = "#1F78B4"))

        # Extract legend (suppress warning about multiple guides)
        legend <- suppressWarnings(cowplot::get_legend(p1))

        # Remove legends from both plots
        p1_clean <- p1 + theme(legend.position = "none")
        p2_clean <- p2 + theme(legend.position = "none")

        # Combine plots and add legend once at bottom
        combined_plot <- cowplot::plot_grid(
            p1_clean, p2_clean, ncol = 1, align = "v", rel_heights = c(1, 1)
        )
        final_plot <- cowplot::plot_grid(combined_plot, legend, ncol = 1, rel_heights = c(1, 0.08))

        if (!dir.exists("plots")) dir.create("plots", recursive = TRUE)
        ggsave(paste0("plots/combined_line_plot_", metric, ".png"),
               final_plot, width = 14, height = 10, dpi = 300)
        ggsave(paste0("plots/combined_line_plot_", metric, ".pdf"),
               final_plot, width = 14, height = 10)
    }
}
