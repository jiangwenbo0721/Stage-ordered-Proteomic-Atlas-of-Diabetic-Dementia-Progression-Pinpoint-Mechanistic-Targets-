# Load required packages
library(ggplot2)    # For data visualization
library(dplyr)      # For data manipulation
library(readxl)     # For reading Excel files
library(forcats)    # For factor manipulation
library(scales)     # For percentage formatting (used in scale_y_continuous)

# Read data from Excel file
A <- read_excel("/code/Figure1/Figure1d/figure1_d_Stacked Bar Plot.xlsx")
head(A)  # Preview first few rows of data

# Define custom color palette for categories
my_colors <- c(
  "#4E79A7", "#F28E2B", "#E15759", "lightgrey", "#59A14F", 
  "#EDC948", "#B07AA1", "#FF9DA7", "#bcbd22",  "#fdae61", 
  "#86BCB6","#C70039", "#FABFD2",  "#17becf", "#8A3B1B", 
  "#e377c2", "#7E6B3C", "#F0E6A2", "#6B4F8B", "#D35400", 
  "#B9D1E3", "#BB8D6F", "#D3B8E5", "#97C5E2"
)

# Calculate percentage composition by group
A_percent <- A %>%
  group_by(Group) %>%                          # Group data by the 'Group' column
  mutate(Percent = Count / sum(Count) * 100) %>% # Calculate percentage of each category within group
  arrange(Group, Percent) %>%                  # Sort by Group and Percent
  ungroup() %>%                                # Remove grouping structure
  mutate(Category = factor(Category))          # Convert Category to factor

# Create stacked bar plot with percentages
p=ggplot(A_percent, 
       aes(x = Group,                          # X-axis: Group variable
           y = Percent,                        # Y-axis: Percentage values
           fill = fct_reorder2(Category, Group, Percent))) +  # Order categories by percentage within groups
  geom_bar(stat = "identity",                  # Create bar plot using exact y-values
           width = 0.7,                        # Set bar width
           color = "black") +                  # Add black borders to bars
  geom_text(aes(label = paste0(round(Percent, 1), "%")),  # Add percentage labels
            position = position_stack(vjust = 0.5),  # Center labels in each stack
            size = 3,                           # Label text size
            color = "white") +                  # White text for contrast
  scale_fill_manual(values = my_colors) +       # Apply custom color palette
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +  # Format y-axis as percentages
  labs(x = "Protein",                          # X-axis label
       y = "Pathway Percentage(%)",            # Y-axis label
       fill = "Category") +                    # Legend title
  theme_minimal() +                            # Use minimal theme as base
  theme(
    axis.text.x = element_text(size = 12,      # X-axis text size
                               angle = 45,              # Tilt x-axis labels 45 degrees
                               hjust = 1),              # Adjust horizontal justification
    axis.text.y = element_text(size = 12),      # Y-axis text size
    axis.title = element_text(size = 14),       # Axis title size
    legend.title = element_text(size = 12),     # Legend title size
    legend.text = element_text(size = 10)       # Legend item text size
  )
ggsave('/results/Figure1d.pdf', plot = p,  height = 6, width = 13, units = "in")
