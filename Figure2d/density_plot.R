# Load necessary libraries
library(ggplot2)
library(ggridges)
library(dplyr)
library(readxl)
library(tidyr)

# Read the data from Excel file
df <- read_excel("/code/Figure2/Figure2d/density_plot_data.xlsx")

# Data preprocessing
df_filtered <- df %>% 
  filter(group %in% c("Q1", "Q10")) %>%  # Filter data for groups Q1 and Q10
  mutate(
    group = factor(group, levels = c("Q1", "Q10")),  # Set the order of the group factor
    type = factor(type, levels = c("Pairs matching accuracy", "Numeric memory", 
                                   "Trail making numeric path", "Matrix pattern completion", 
                                   "Fluid intelligence", "Reaction time", 
                                   "Paired learning", "Picture vocabulary")),  # Set the order of the 'type' factor
    predicted_cognition_norm = as.numeric(predicted_cognition_norm)  # Ensure numeric type for the predicted cognition
  ) %>%
  drop_na(predicted_cognition_norm)  # Remove rows with NA values in the 'predicted_cognition_norm' column

# Calculate mean values for each group, type, and protein
mean_values <- df_filtered %>%
  group_by(protein, type, group) %>%
  summarise(
    mean_val = mean(predicted_cognition_norm, na.rm = TRUE),  # Calculate mean of 'predicted_cognition_norm'
    .groups = "drop"  # Drop the grouping after summarizing to avoid grouping warnings
  )

# Set x-axis range (example range, adjust according to your data)
x_range <- c(0, 1.2)  # Alternatively, you can use range(df_filtered$predicted_cognition_norm, na.rm = TRUE)

# Read pfdr data
pp <- read_excel("/code/Figure2/Figure2d/pfdr_data.xlsx")

# Merge df_filtered and pp data frames by 'protein' and 'type'
df_with_text <- df_filtered %>%
  left_join(pp, by = c("protein", "type"))  # Merge df_filtered and pp data frames by 'protein' and 'type'

# Create a new column to flag the first row for each protein and type combination
df_with_text <- df_with_text %>%
  group_by(protein, type) %>%
  mutate(show_pfdr = row_number() == 1) %>%  # Flag only the first row in each protein-type group
  ungroup()

# Create the plot
p=ggplot(df_with_text, aes(x = predicted_cognition_norm, y = group, fill = group)) +
  # Density ridge plot
  geom_density_ridges(
    alpha = 0.7, 
    scale = 0.9, 
    size = 0.3,  # Set the line thickness
    panel_scaling = TRUE  # Ensure scaling of each panel
  ) +
  # Add mean line
  geom_vline(
    data = mean_values,
    aes(xintercept = mean_val, color = group),
    linetype = "dashed",
    size = 0.6,
    show.legend = FALSE
  ) +
  # Add text (pfdr) to the plot, only display pfdr at the first position for each protein-type combination
  geom_text(
    aes(label = ifelse(show_pfdr, as.character(pfdr), "")),  # Display pfdr value only at the first position
    position = position_jitter(width = 0.05, height = 0),  # Adjust text position to avoid overlap
    size = 3,  # Set text size
    color = "black",  # Set text color
    show.legend = FALSE  # Do not show a legend for text
  ) +
  # Facet settings
  facet_grid(
    protein ~ type, 
    scales = "free_y", 
    space = "free_y"  # Allow free space for y-axis in each facet
  ) +
  # Color settings
  scale_fill_manual(values = c("Q1" = "#b8d9bc", "Q10" = "#f4bf9f")) +
  scale_color_manual(values = c("Q1" = "#b8d9bc", "Q10" = "#f4bf9f")) +
  # x-axis settings
  scale_x_continuous(
    limits = x_range,
    breaks = seq(x_range[1], x_range[2], by = 0.2),
    expand = c(0, 0)
  ) +
  # Labels and title
  labs(
    x = "Predicted Cognition (Normalized)",
    y = "Group",
    title = "Q1/Q10 Distribution by Protein"
  ) +
  # Theme settings
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",  # No legend for the plot
    strip.text = element_text(face = "bold", size = 10),  # Bold text for facet labels
    strip.text.y = element_text(angle = 0, hjust = 0),  # Vertical facet label alignment
    panel.spacing = unit(0.2, "cm"),  # Space between panels
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),  # Rotate x-axis labels for better visibility
    axis.text.y = element_text(size = 10),  # Set y-axis label size
    panel.grid.major.x = element_blank(),  # Remove major grid lines for x-axis
    panel.grid.major.y = element_line(color = "gray90", size = 0.2),  # Add light gray grid lines for y-axis
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    plot.title = element_text(hjust = 0.5, face = "bold")  # Center and bold the plot title
  )
ggsave('/results/Figure2d.pdf', plot = p,  height = 6, width = 13, units = "in")