# Load necessary packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(readxl)
library(patchwork)
library(cowplot) 

# 1. Read and preprocess data
df <- read_excel("/code/Figure3/Figure3a/brain_data.xlsx", sheet=1)

# Define a function to extract the base brain region name (remove left/right identifiers)
get_base_region <- function(outcome_name) {
  # Handle IDP T1 FAST ROIs type (L/R prefix)
  if (str_detect(outcome_name, "^IDP.T1.FAST.ROIs.[LR].")) {
    return(str_replace(outcome_name, "^IDP.T1.FAST.ROIs.[LR].", ""))
  }
  # Handle dMRI type (.left./.right.)
  else if (str_detect(outcome_name, "\\.left\\.|\\.right\\.")) {
    return(str_replace_all(outcome_name, "\\.left\\.|\\.right\\.", "\\."))
  }
  # Return original name for other cases
  else {
    return(outcome_name)
  }
}

# Add a base_region column for all rows
df <- df %>% 
  mutate(base_region = sapply(outcome, get_base_region))

# Find all base regions with at least one significant result (p.value < 0.05)
significant_regions <- df %>%
  filter(FDR < 0.05) %>%
  distinct(base_region) %>%
  pull(base_region)

# Filter all rows corresponding to these brain regions (including non-significant results)
result <- df %>%
  filter(base_region %in% significant_regions) %>%
  select(-base_region)  # Remove the temporarily added column

# View the result
print(result)

# 1. Read and preprocess data

aa <- result
aa$outcome <- make.names(aa$outcome)

# Data transformation (retain Region information)
transformed_data <- aa %>%
  select(term, standardized_estimate, FDR, outcome, Region) %>%
  pivot_wider(
    names_from = term,
    values_from = c(standardized_estimate, FDR),
    names_glue = "{term}_{.value}"
  ) %>%
  rename(Metric = outcome) %>%
  select(
    Metric,
    Region,
    GDF15_standardized_estimate,
    LGALS4_standardized_estimate,
    GDF15_FDR,
    LGALS4_FDR
  ) %>%
  rename(
    GDF15_status = GDF15_standardized_estimate,
    LGALS4_status = LGALS4_standardized_estimate,
    GDF15_FDR = GDF15_FDR,
    LGALS4_FDR = LGALS4_FDR
  )

aa <- transformed_data

# 2. Clean up metric names (remove left/right markers)
aa <- aa %>%
  mutate(
    Metric = gsub("^IDP\\.T1\\.FAST\\.ROIs\\.", "", Metric),
    Metric = gsub("^Weighted\\.mean\\.", "", Metric),
    Metric = gsub("\\.orientation\\.dispersion\\.index\\.\\.in\\.tract\\.", "", Metric),
    Metric = gsub("\\.diffusion\\.tensor\\.mode\\.\\.in\\.tract\\.", "", Metric),
    Metric = gsub("\\.mean\\.diffusivity\\.\\.in\\.tract\\.", "", Metric),
    Metric = gsub("\\.isotropic\\.or\\.free\\.water\\.volume\\.fraction\\.\\.in\\.tract\\.", "", Metric),
    Metric = gsub("\\.intra\\.cellular\\.volume\\.fraction\\.\\.in\\.tract\\.", "", Metric),
    Metric = gsub("\\.fractional\\.anisotropy\\.\\.in\\.tract\\.", "", Metric),
    Metric = gsub("\\.\\.from\\.dMRI\\.data\\.$", "", Metric),
    # Remove left/right markers (retain main region name)
    Metric_modified = str_remove(Metric, "^(FA|MD|MO|OD|ICVF|ISOVF)\\.(left|right)\\."),
    Metric_modified = str_remove(Metric_modified, "^(L|R)\\.")
  )

# 3. Add Category classification
aa <- aa %>%
  mutate(
    Category = case_when(
      grepl("^L|^R", Metric) ~ "GMV",
      grepl("^FA", Metric) ~ "FA",
      grepl("^MD", Metric) ~ "MD",
      grepl("^MO", Metric) ~ "MO",
      grepl("^OD", Metric) ~ "OD",
      grepl("^ICVF", Metric) ~ "ICVF",
      grepl("^ISOVF", Metric) ~ "ISOVF",
      TRUE ~ "Other"
    )
  )

# 4. Calculate -log10(FDR)
aa <- aa %>%
  mutate(
    GDF15_log_FDR = -log10(GDF15_FDR),
    LGALS4_log_FDR = -log10(LGALS4_FDR)
  )

# 5. Sort (by Category and maximum FDR value, GMV first)
aa <- aa %>% 
  arrange(factor(Category, levels = c("GMV", "FA", "MD", "MO", "OD", "ICVF", "ISOVF")), 
          desc(pmax(GDF15_log_FDR, LGALS4_log_FDR))) %>%
  mutate(Metric = factor(Metric, levels = unique(Metric)))

# 6. Create plotting data (handle NA issues)
plot_data <- bind_rows(
  # GDF15 Left
  aa %>% filter(Region == "L") %>% 
    select(Metric, Metric_modified, Category, status = GDF15_status, FDR = GDF15_FDR) %>%
    mutate(model_region = "GDF15_L"),
  # GDF15 Right
  aa %>% filter(Region == "R") %>% 
    select(Metric, Metric_modified, Category, status = GDF15_status, FDR = GDF15_FDR) %>%
    mutate(model_region = "GDF15_R"),
  # LGALS4 Left
  aa %>% filter(Region == "L") %>% 
    select(Metric, Metric_modified, Category, status = LGALS4_status, FDR = LGALS4_FDR) %>%
    mutate(model_region = "LGALS4_L"),
  # LGALS4 Right
  aa %>% filter(Region == "R") %>% 
    select(Metric, Metric_modified, Category, status = LGALS4_status, FDR = LGALS4_FDR) %>%
    mutate(model_region = "LGALS4_R")
) %>%
  mutate(
    model_region = factor(model_region,
                          levels = c("GDF15_L", "GDF15_R", "LGALS4_L", "LGALS4_R"),
                          labels = c("GDF15 (L)", "GDF15 (R)", "LGALS4 (L)", "LGALS4 (R)"))
  )

# 7. Define color mapping
category_colors <- c(
  GMV = "#458B74",
  FA = "#4682B4",
  MD = "#1D1D8F",
  MO = "#DAA520",
  OD = "#CD5B45",
  ICVF = "#C71585",
  ISOVF = "#91D1C2",
  Other = "gray"
)

# 8. Create bubble plot (hide legend for now, add later)
bubble_plot <- ggplot(plot_data, aes(x = model_region, y = Metric_modified)) +
  geom_point(aes(size = abs(status)), color = "black", shape = 1, stroke = 0.8) +
  geom_point(aes(size = abs(status), color = status), alpha = 0.8) +
  scale_size_continuous(range = c(3, 10), name = "Effect size") +
  scale_color_gradient2(
    low = "#08519C",
    mid = "white",
    high = "#A50F15",
    midpoint = 0,
    name = "Effect direction"
  ) +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",  # Hide legend for now
    panel.grid.major.y = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 0, face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 10),
    strip.text = element_text(face = "bold", size = 12),
    panel.spacing = unit(0.5, "lines"),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5)
  ) +
  facet_grid(Category ~ ., scales = "free_y", space = "free_y")

# 9. Create FDR bar plot (distinguish left/right)
# Prepare data
fdr_data <- plot_data %>%
  mutate(
    model = str_extract(model_region, "GDF15|LGALS4"),
    region = str_extract(model_region, "\\(L\\)|\\(R\\)") %>% 
      str_remove_all("[()]"),
    log_FDR = -log10(FDR)
  ) %>%
  unite("model_region", model, region, sep = "_") %>%
  mutate(model_region = factor(model_region, 
                               levels = c("GDF15_L", "GDF15_R", "LGALS4_L", "LGALS4_R"),
                               labels = c("GDF15 (L)", "GDF15 (R)", "LGALS4 (L)", "LGALS4 (R)")))

fdr_plot <- ggplot(fdr_data, aes(y = Metric_modified, x = log_FDR, 
                                 fill = Category, alpha = model_region)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, 
           color = "black", linewidth = 0.3) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red") +
  scale_fill_manual(values = category_colors, name = "Category") +
  scale_alpha_manual(
    values = c("GDF15 (L)" = 0.1, "GDF15 (R)" = 0.4, "LGALS4 (L)" = 0.7, "LGALS4 (R)" = 1),
    name = "Model & Region"
  ) +
  labs(x = "-log10(FDR)", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
    legend.position = "right",
    legend.box = "vertical",
    legend.spacing.y = unit(0.2, "cm"),
    legend.key.size = unit(0.8, "cm"),
    plot.margin = margin(5, 15, 5, 5),
    strip.text = element_blank()
  ) +
  facet_grid(Category ~ ., scales = "free_y", space = "free_y")

# 10. Combine plots
# Extract legends
legend_bubble <- cowplot::get_legend(
  bubble_plot + 
    guides(size = guide_legend(order = 1), color = guide_legend(order = 2)) +
    theme(legend.position = "right",
          legend.box = "vertical",
          legend.spacing.y = unit(0.5, "cm"))
)

legend_fdr <- cowplot::get_legend(
  fdr_plot + 
    guides(fill = guide_legend(order = 1), alpha = guide_legend(order = 2)) +
    theme(legend.position = "right",
          legend.box = "vertical",
          legend.spacing.y = unit(0.5, "cm"))
)

# Combine main plots (without legends)
main_plot <- bubble_plot + fdr_plot + 
  plot_layout(widths = c(4, 1.5)) & 
  theme(legend.position = "none")

# Add legends to main plot
final_plot <- main_plot +
  inset_element(legend_bubble, left = 0.9, bottom = 0.7, right = 1, top = 1) +
  inset_element(legend_fdr, left = 0.9, bottom = 0.4, right = 1, top = 0.65) +
  plot_annotation(theme = theme(plot.margin = margin(1, 5, 1, 1, "cm")))

# 1. Create bubble plot (with legend)
bubble_plot_with_legend <- bubble_plot +
  theme(legend.position = "right",
        legend.box = "vertical",
        legend.spacing.y = unit(0.5, "cm"))

# 2. Create FDR bar plot (with legend)
fdr_plot_with_legend <- fdr_plot +
  theme(legend.position = "right",
        legend.box = "vertical",
        legend.spacing.y = unit(0.5, "cm"))

# 3. Extract legends (optional, if you need separate control over legend position)
legend_bubble <- get_legend(bubble_plot_with_legend)
legend_fdr <- get_legend(fdr_plot_with_legend)

# 4. Create main plots without legends
bubble_plot_no_legend <- bubble_plot + theme(legend.position = "none")
fdr_plot_no_legend <- fdr_plot + theme(legend.position = "none")

# 5. Combine plots using plot_grid
combined_plot <- plot_grid(
  bubble_plot_no_legend, 
  fdr_plot_no_legend,
  nrow = 1,
  align = "h",
  axis = "tb",
  rel_widths = c(4, 1.5)
)

# 6. Add legends (two options)

# Option 1: Place legends on the right
final_plot <- plot_grid(
  combined_plot,
  plot_grid(legend_bubble, legend_fdr, ncol = 1),
  nrow = 1,
  rel_widths = c(8, 2)
)

# 7. Save the plot
ggsave("/results/figure3a.pdf",
       final_plot, 
       width = 8,
       height = 4 + length(unique(plot_data$Metric))*0.15,
       device = cairo_pdf,
       limitsize = FALSE)

# Display the plot
print(final_plot)
