library(readxl)
library(ggplot2)
library(dplyr)

df <- read_excel("/code/Figure2/Figure 2a/Figure2a_ForestPlot.xlsx")

# Convert data types
df <- df %>%
  mutate(
    OR = as.numeric(OR),
    CI_low = as.numeric(CI_low),
    CI_high = as.numeric(CI_high),
    ymin_rect = OR - (CI_high - CI_low) / 4,
    ymax_rect = OR + (CI_high - CI_low) / 4,
    x_numeric = as.numeric(factor(Outcome, levels = Outcome))  # For plotting
  )

# Color palette
color_palette <- c(
  "Fluctuating Decreasing" = "#59355E",
  "Fluctuating Increasing" = "#118A8C",
  "Gradually Decreasing"  = "#EDA598",
  "Gradually Increasing"  = "#64A5DA"
)

p <- ggplot(df) +
  # Draw confidence interval lines
  geom_errorbar(aes(x = x_numeric, ymin = CI_low, ymax = CI_high, color = Group),
                width = 0.2, linewidth = 1.2) +
  
  # Use rectangles to represent OR values (height is half of CI range)
  geom_rect(aes(xmin = x_numeric - 0.15, xmax = x_numeric + 0.15,
                ymin = ymin_rect, ymax = ymax_rect, fill = Group),
            color = NA, alpha = 1) +
  
  # Reference line at OR = 1
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.8) +
  
  # Axis settings
  scale_y_continuous(breaks = seq(0.8, 1.2, 0.05)) +
  scale_x_continuous(breaks = df$x_numeric, labels = df$Outcome, expand = expansion(mult = 0.1)) +
  scale_color_manual(values = color_palette) +
  scale_fill_manual(values = color_palette) +
  
  # Theme settings
  labs(x = "Outcome", y = "OR (95% CI)", color = "Protein Group", fill = "Protein Group") +
  theme_classic() +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold")
  )


# Save as PDF
ggsave("/results/Figure2a_ForestPlot.pdf", height = 4, width = 12)
