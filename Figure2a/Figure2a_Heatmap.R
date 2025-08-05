# Load necessary libraries
library(ggplot2)
library(dplyr)
library(readr)
library(RColorBrewer)
library(readxl)

# Read data
df <- read_excel("/code/Figure2/Figure 2a/Figure2a_Heatmap.xlsx")

# Remove unwanted proteins
df <- df %>%
  filter(!Protein %in% c("PAEP", "CD58", "CTRC", "MEPE", "SELE", "CA14", "FGFBP1", "PRCP", "CDON",
                         "NCLN", "ACAN", "CTSO", "MICB_MICA", "FST",
                         "BCAN", "PTPRN2", "CREG1", "FASLG", "FOLR3", "RRM2" ))

# Set protein order
protein_order <- c("MATN3", "IGFBP1", "GHRL", "UMOD", "PON3",
                   "APOM", "GDF15", "CD276", "PLXNB2", "RTN4R", "NFASC", "LGALS4",
                   "SCARA5", "DSG2", "IGFBP2", 
                   "EGFR", "EPS8L2", "REN", "FABP1",
                   "RHOC" )

# Data preprocessing
df <- df %>%
  mutate(
    direction = ifelse(HR >= 1, "positive", "negative"),
    Outcome_group = paste0(Outcome, "_", group),
    # Set size for p-values: P <= 0.01 fixed at 0.05; linear mapping for others; NA if P > 0.05
    size_value = case_when(
      p_fdr <= 0.01 ~ 0.05,
      p_fdr <= 0.05 ~ 0.05 - p_fdr,
      TRUE ~ NA_real_
    ),
    # Cap HR values for color mapping
    HR_display = ifelse(HR > 2.5, 2.5, HR),
    Protein = factor(Protein, levels = protein_order)
  )

# Set factor levels for Outcome and Outcome_group
df <- df %>%
  mutate(
    Outcome = factor(Outcome, levels = unique(Outcome)),
    Outcome_group = paste0(Outcome, "_", group)
  ) %>%
  arrange(Outcome, group) %>%
  mutate(
    Outcome_group = factor(Outcome_group, levels = unique(Outcome_group)),
    Protein = factor(Protein, levels = protein_order)
  )

# Prepare rectangles for group annotations
df_rect <- df %>%
  distinct(Outcome, Outcome_group) %>%
  arrange(Outcome_group) %>%
  group_by(Outcome) %>%
  summarise(
    y_min = which(levels(df$Outcome_group) == Outcome_group[1]) - 0.5,
    y_max = which(levels(df$Outcome_group) == Outcome_group[2]) + 0.5,
    y_mid = mean(c(y_min, y_max)),
    .groups = "drop"
  )

# Custom color palette
custom_colors <- c(
  "#1f78b4", "#33a0cc", "#66c2a5", "#abdda4", "#e6f598",
  "#fee08b", "#fdae61", "#f46d43", "#d53e4f", "#b2182b"
)

# Breakpoints for HR legend
min_hr <- min(df$HR, na.rm = TRUE)
legend_breaks <- unique(sort(c(seq(min_hr, 2.5, length.out = 5), 2.5)))
legend_labels <- ifelse(legend_breaks == 2.5, ">2.5", sprintf("%.1f", legend_breaks))

# Plot
p <- ggplot(df, aes(x = Protein, y = Outcome_group)) +
  geom_point(
    data = df %>% filter(!is.na(size_value)),
    aes(size = size_value, fill = HR_display),
    shape = 21, color = "black", stroke = 0.4
  ) +
  annotate("text", x = 0.3, y = df_rect$y_mid, label = df_rect$Outcome, hjust = 1, size = 4.2) +
  scale_size_continuous(
    name = "FDR",
    range = c(4, 10),
    breaks = c(0.05, 0.03, 0.02, 0.01),
    labels = c("0.01", "0.02", "0.03", "0.04")
  ) +
  scale_fill_gradientn(
    colours = custom_colors,
    name = "HR",
    limits = c(min_hr, 2),
    breaks = legend_breaks,
    labels = legend_labels,
    oob = scales::squish
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "grey90", size = 0.4),
    plot.margin = margin(5.5, 10, 5.5, 80)
  ) +
  labs(
    x = "Protein",
    y = "Outcome",
    title = "Association Heatmap (HR color & P-value size)"
  )

# Display plot
print(p)

# Save plot
ggsave("/results/Figure2a_Heatmap.pdf", plot = p, width = 13, height = 5.5, dpi = 300)

