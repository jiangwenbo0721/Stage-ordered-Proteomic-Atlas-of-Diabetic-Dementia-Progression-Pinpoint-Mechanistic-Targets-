# Load required packages
library(tidyverse)            # For data manipulation and piping (%>%)
library(funkyheatmap)         # For creating customized heatmaps

# Load the R data file containing all necessary objects
load('/code/Figure4/Figure4c/Hdata.Rdata')

# Create the funky heatmap visualization
p02 = funky_heatmap(
  data          = j01_CSEA,           # Main data matrix for the heatmap
  column_info   = j02_column_info,    # Column metadata/annotations
  palettes      = j02_palettes,       # Color palettes for visualization
  # position_args = position_arguments(
  #   col_annot_offset = 3
  # ), # Uncomment if column annotations need more space
  scale_column = F # Important: Disable column scaling (keep raw values)
)
ggsave('/results/figure4c.pdf', plot = p02, height = 13, width = 10, units = "in")

