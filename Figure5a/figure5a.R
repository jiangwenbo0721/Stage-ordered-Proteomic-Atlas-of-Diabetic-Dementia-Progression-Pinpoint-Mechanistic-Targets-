library(ggplot2)
library(dplyr)
library(patchwork)
library(ggrepel)



DEG <- read.csv("/code/Figure5/Figure5a/Metabolomics_Volcano_Plot_dementia.csv", header = TRUE)%>%
  filter(p_fdr != 0)

# Setting the Theme
mytheme <- theme_classic() +
  theme(
    legend.key = element_rect(fill = 'transparent'),
    legend.background = element_rect(fill = 'transparent'),
    legend.position = "bottom",
    axis.text.x = element_text(hjust = 0.5, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.line = element_line(),
    plot.title = element_text(size = 16, hjust = 0.5)
  )

# Add group information
DEG <- DEG %>%
  mutate(Group = case_when(
    p_fdr < 0.05 & beta_std > 0 ~ "Up",
    p_fdr < 0.05 & beta_std < 0 ~ "Down",
    TRUE ~ "Stable"
  )) %>%
  mutate(Group = factor(Group, levels = c("Up", "Stable", "Down")))

# Define volcano plot drawing function (with labeling)
draw_volcano <- function(data, protein_name) {
  # Extract top 10 significantly up-regulated & down-regulated results (for labeling)
  label_data <- data %>%
    filter(Group != "Stable") %>%
    group_by(Group) %>%
    slice_max(order_by = -log10(p_fdr), n = 5, with_ties = FALSE) %>%
    ungroup()
  
  # Count the number of Up and Down
  up_count <- sum(data$p_fdr < 0.05 & data$beta_std > 0)
  down_count <- sum(data$p_fdr < 0.05 & data$beta_std < 0)
  
  ggplot(data, aes(x = beta_std, y = -log10(p_fdr))) +
    geom_point(aes(color = beta_std), size = 2, alpha = 0.8)+
    geom_text_repel(
      data = label_data,
      aes(label = Abbreviation),
      color = "black",
      size = 4,
      min.segment.length = 0,
      segment.linetype = 1,
      force = 5,
      force_pull = 4,
      max.overlaps = Inf,
      box.padding = unit(0.6, "lines"),
      point.padding = unit(0.5, "lines"),
      segment.color = "black",
      segment.size = 0.4
    ) +
    # Add quantity annotation
    annotate("text", x = 0.23, y = max(-log10(data$p_fdr), na.rm = TRUE), 
             label = paste0("Up: ", up_count), hjust = 1, linewidth = 4.5, color = "#b81f25") +
    annotate("text", x = -0.23, y = max(-log10(data$p_fdr), na.rm = TRUE), 
             label = paste0("Down: ", down_count), hjust = 0, linewidth = 4.5, color = "#39489f") +
    scale_color_gradientn(colors = c("#39489f", "#39bbec", "#f9ed36", "#f38466", "#b81f25"),
                          limits = c(-0.25, 0.25)) +
    geom_hline(yintercept = -log10(0.05), size = 0.7, color = "black", lty = "dashed") +
    labs(
      x = "Standardized β",
      y = bquote(~-Log[10]~"(P-FDR)"),
      color = "Standardized β",
      title = protein_name
    ) +
    mytheme
}


# Extract protein name list
protein_list <- unique(DEG$exposure)

# Batch drawing
plots <- lapply(protein_list, function(prot) {
  df_sub <- DEG %>% filter(exposure == prot)
  draw_volcano(df_sub, prot)
})

# Merged display
combined_plot <- wrap_plots(plots, ncol = 3)

# Save Image
ggsave(
  plot = combined_plot,
  filename = "/results/Figure5a.pdf",
  height = 8,
  width = 12,
  device = "pdf",
  dpi = 600
)

