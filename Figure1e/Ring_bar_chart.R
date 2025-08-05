gd <- read.csv("/code/Figure1/Figure 1e/Dementia-gradient_importance_Gradually_Decreasing_Deep_ordinal_regression.csv", stringsAsFactors = FALSE)
gi <- read.csv("/code/Figure1/Figure 1e/Dementia-gradient_importance_dementia_Gradually_Increasing_Deep_ordinal_regression.csv", stringsAsFactors = FALSE)
fi <- read.csv("/code/Figure1/Figure 1e/Dementia-gradient_importance_dementia_Fluctuating_Increasing_Deep_ordinal_regression.csv", stringsAsFactors = FALSE)
fd <- read.csv("/code/Figure1/Figure 1e/Dementia-gradient_importance_dementia_Fluctuating_Decreasing_Deep_ordinal_regression.csv", stringsAsFactors = FALSE)

# Visualization
library(dplyr)
library(ggplot2)
library(scales)
library(patchwork)

plot_single_ring <- function(df, source_name, fill_colors) {
  # Compute the maximum value for the group, round to nearest integer or 0.5
  max_value <- max(df$Gradient.Importance)
  axis_limit <- ifelse(
    max_value %% 1 == 0, 
    max_value,  # If already an integer, keep it
    ceiling(max_value * 2) / 2  # Otherwise round up to nearest 0.5
  )
  
  # Set axis height based on the rounded max value
  axis_height <- axis_limit
  
  # Compute bar height, position and angle for top 10 proteins
  df_top <- df %>% 
    slice_max(Gradient.Importance, n = 10) %>%
    mutate(
      height_ratio = Gradient.Importance / max_value,
      ymin = 10,
      ymax = ymin + height_ratio * axis_height,
      id = row_number(),
      angle = 360 * (id - 1)/10,
      label_angle = ifelse(angle > 180, angle - 180, angle),
      hjust = ifelse(angle > 180, 1, 0)
    )
  
  axis_angle <- 360
  
  # Create axis line and tick marks
  axis_line <- data.frame(
    x = axis_angle,
    xend = axis_angle,
    y = 10,
    yend = 10 + axis_height
  )
  
  tick_height <- axis_height + 10
  ticks <- data.frame(
    x = axis_angle,
    xend = axis_angle,
    y = tick_height - 0.2,
    yend = tick_height + 0.2,
    label = round(axis_limit, 1)
  )
  
  # Generate the ring plot
  p=ggplot(df_top) +
    geom_rect(
      aes(xmin = angle - 18, xmax = angle + 18, ymin = ymin, ymax = ymax, fill = Gradient.Importance),
      color = "white", linewidth = 0.3
    ) +
    scale_fill_gradientn(
      colours = fill_colors,
      name = paste(source_name, "\nImportance")
    ) +
    geom_segment(data = axis_line,
                 aes(x = x, xend = xend, y = y, yend = yend),
                 inherit.aes = FALSE, color = "gray20", linewidth = 0.5) +
    geom_segment(data = ticks,
                 aes(x = x, xend = xend, y = y, yend = yend),
                 inherit.aes = FALSE, color = "gray20", linewidth = 0.5) +
    geom_text(data = ticks,
              aes(x = x, y = yend + 0.2, label = label),
              inherit.aes = FALSE, hjust = 0.5, vjust = 0, size = 2.5, color = "gray20") +
    geom_text(
      aes(x = angle, y = ymax + 0.3, label = Protein),
      size = 2.5, angle = 0, hjust = 0.5
    ) +
    coord_polar(theta = "x", start = pi/2) +
    theme_void() +
    ggtitle(paste("Ring Plot -", source_name)) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
      plot.margin = margin(5, 5, 5, 5),
      legend.position = "right",
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7)
    )
}

# 2. Color settings for each group
# Purple tones with reduced saturation (more gray/pink for softness)
fd_colors <- c("#d4cad4", "#c39cc2", "#9a6f98", "#7c4d72")

# Cyan tones with reduced saturation (less sharpness)
fi_colors <- c("#c4f0f0", "#6cc5c5", "#459999", "#1a6666")

# Orange tones with reduced saturation (less fiery)
gd_colors <- c("#f5dede", "#f5b4a7", "#d97d61", "#b3522a")

# Blue tones with reduced saturation (more matte texture)
gi_colors <- c("#dde8f2", "#8cbfe8", "#5a88b4", "#2a5680")

# Generate plots
p1 <- plot_single_ring(fd, "fd", fd_colors)
p2 <- plot_single_ring(fi, "fi", fi_colors)
p3 <- plot_single_ring(gd, "gd", gd_colors)
p4 <- plot_single_ring(gi, "gi", gi_colors)

# Combine plots with a title
final_plot <- (p1 | p2) / (p3 | p4) +
  plot_annotation(title = "Comparative Ring Plots - Gradient Enhancement Version",
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold")))

# Save output

ggsave('/results/Figure1e.pdf', plot = final_plot,  height = 8, width = 12, units = "in")
