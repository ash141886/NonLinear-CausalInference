# =============================================================================
# SECTION 7: VISUALIZATION 
# =============================================================================

generate_figures <- function(results,
                             base_size = 9,
                             title_size = 10,
                             axis_title_size = 9,
                             axis_text_size = 8,
                             strip_size = 9,
                             legend_title_size = 9,
                             legend_text_size = 8) {
  if (!dir.exists("figures")) dir.create("figures", recursive = TRUE)
  method_colors <- c("LiNGAM" = "#E31A1C", "Proposed Method" = "#1F78B4")
  method_shapes <- c("LiNGAM" = 16, "Proposed Method" = 17)
  method_lines  <- c("LiNGAM" = "solid", "Proposed Method" = "dashed")
  metrics_to_plot <- c("Precision_dir", "F1_Score_dir", "Graph_Accuracy", "Misoriented", "SHD", "MSE", "Time")
  plain_theme <- theme_minimal(base_size = base_size) +
    theme(
      text            = element_text(face = "plain"),
      plot.title      = element_text(size = title_size, face = "plain", hjust = 0.5),
      axis.title      = element_text(size = axis_title_size, face = "plain"),
      axis.text       = element_text(size = axis_text_size, face = "plain"),
      legend.title    = element_text(size = legend_title_size, face = "plain"),
      legend.text     = element_text(size = legend_text_size, face = "plain"),
      legend.position = "bottom",
      legend.key.width= unit(1.2, "cm"),
      strip.text      = element_text(size = strip_size, face = "plain"),
      panel.grid.minor= element_blank()
    )
  for (metric in metrics_to_plot) {
    p_samples <- ggplot(results$summary,
                        aes(x = Samples, y = .data[[metric]],
                            color = Method, shape = Method, linetype = Method)) +
      geom_line(linewidth = 0.5) +
      geom_point(size = 1.3) +
      facet_wrap(~ Variables, scales = "free_y", nrow = 1) +
      labs(title = paste(metric, "vs. Sample Size"), x = "Sample Size", y = metric) +
      plain_theme +
      theme(legend.position = "none") +
      scale_color_manual(values = method_colors) +
      scale_shape_manual(values = method_shapes) +
      scale_linetype_manual(values = method_lines)

    p_variables <- ggplot(results$summary,
                          aes(x = Variables, y = .data[[metric]],
                              color = Method, shape = Method, linetype = Method)) +
      geom_line(linewidth = 0.5) +
      geom_point(size = 1.3) +
      facet_wrap(~ Samples, scales = "free_y", nrow = 1) +
      labs(title = paste(metric, "vs. Number of Variables"), x = "Number of Variables", y = metric) +
      plain_theme +
      scale_color_manual(name = "Method", values = method_colors) +
      scale_shape_manual(name = "Method", values = method_shapes) +
      scale_linetype_manual(name = "Method", values = method_lines) +
      scale_x_continuous(breaks = unique(results$summary$Variables))

    if (metric %in% c("Precision_dir", "F1_Score_dir", "Graph_Accuracy")) {
      p_samples   <- p_samples   + scale_y_continuous(limits = c(0, 1))
      p_variables <- p_variables + scale_y_continuous(limits = c(0, 1))
    }

    g <- gridExtra::arrangeGrob(p_samples, p_variables, nrow = 2)
    print(g)  # display on-screen in interactive sessions
    ggsave(paste0("figures/", metric, ".pdf"), g, width = 12, height = 8, dpi = 300)
    ggsave(paste0("figures/", metric, ".png"), g, width = 12, height = 8, dpi = 300)
  }
}
