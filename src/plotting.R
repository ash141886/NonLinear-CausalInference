# Ensure you have the necessary libraries installed and loaded
# install.packages(c("ggplot2", "gridExtra"))
library(ggplot2)
library(gridExtra)

# Enhanced plotting function
plot_combined_results <- function(results) {
    metrics <- c("F1_Score_dir", "Graph_Accuracy", "SHD", "MSE", "Time", "Misoriented")

    for (metric in metrics) {
        # Plot vs Sample Size (Fixed Variable Counts)
        # Legend is removed from this plot to avoid duplication in the final combined image.
        p1 <- ggplot(results, aes(x = Samples, y = .data[[metric]],
                                  linetype = Method, shape = Method, color = Method, group = Method)) +
            geom_line(size = 0.4) + # Made line thicker
            geom_point(size = 1.5) +
            facet_wrap(~ Variables, scales = "free_y", nrow = 1,
                       labeller = labeller(Variables = function(x) paste("Variables:", x))) +
            labs(title = paste(metric, "vs. Sample Size"),
                 x = "Sample Size", y = metric) +
            theme_minimal() +
            theme(
                plot.title = element_text(size = 14, face = "bold"),
                axis.title = element_text(size = 12),
                legend.position = "none", # Removed legend from the first plot
                strip.text = element_text(size = 10, face = "bold")
            ) +
            scale_linetype_manual(name = "Method",
                                 values = c("LiNGAM" = "solid", "Proposed Method" = "dashed")) +
            scale_shape_manual(name = "Method",
                              values = c("LiNGAM" = 16, "Proposed Method" = 17)) +
            scale_color_manual(name = "Method",
                              values = c("LiNGAM" = "#E31A1C", "Proposed Method" = "#1F78B4"))

        # Plot vs Number of Variables (Fixed Sample Sizes)
        # This plot retains the legend, which will be displayed at the bottom of the combined plot.
        p2 <- ggplot(results, aes(x = Variables, y = .data[[metric]],
                                  linetype = Method, shape = Method, color = Method, group = Method)) +
            geom_line(size = 0.4) + # Made line thicker
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

        # Combine both plots. Since p1 has no legend, grid.arrange will produce a clean layout.
        combined_plot <- grid.arrange(p1, p2, nrow = 2)

        # Save plots
        if (!dir.exists("plots")) {
            dir.create("plots", recursive = TRUE)
        }

        ggsave(paste0("plots/combined_line_plot_", metric, ".png"),
                combined_plot, width = 14, height = 10, dpi = 300)

        ggsave(paste0("plots/combined_line_plot_", metric, ".pdf"),
                combined_plot, width = 14, height = 10)
    }
}
