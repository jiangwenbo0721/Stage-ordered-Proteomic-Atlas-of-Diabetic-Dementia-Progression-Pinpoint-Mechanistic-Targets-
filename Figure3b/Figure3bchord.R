
library(tidyr)
library(dplyr)
library(stringr)
library(ComplexHeatmap)
library(combinat)
library(dplyr)
library(circlize)
library(circlize)
library(dplyr)
library(purrr)

load("/code/Figure3/Figure3b/cirden_Test/res_206_anno.Rdata")
# Function definition -----------------------------------------------------

head(res_206_anno)
res_206_anno$region2_short[grep("^R-", res_206_anno$region1)]=res_206_anno$region2_short[grep("^R-", res_206_anno$region1)] %>% 
  gsub("Accumbens","R-Accumbens",.) %>% 
  gsub("Amygdala", "R-Amygdala",.) %>% 
  gsub("Caudate", "R-Caudate",.) %>% 
  gsub("Hippocampus","R-Hippocampus",.) %>% 
  gsub("Pallidum", "R-Pallidum",.) %>% 
  gsub("Putamen", "R-Putamen",.) %>% 
  gsub("Thalamus","R-Thalamus",.)
res_206_anno$region2_short[grep("^L-", res_206_anno$region1)]=res_206_anno$region2_short[grep("^L-", res_206_anno$region1)] %>% 
  gsub("Accumbens","L-Accumbens",.) %>% 
  gsub("Amygdala", "L-Amygdala",.) %>% 
  gsub("Caudate", "L-Caudate",.) %>% 
  gsub("Hippocampus","L-Hippocampus",.) %>% 
  gsub("Pallidum", "L-Pallidum",.) %>% 
  gsub("Putamen", "L-Putamen",.) %>% 
  gsub("Thalamus","L-Thalamus",.)


res_to_plot=rio::import("/code/Figure3/Figure3b/data.csv")
## Columns extracted for FDR plotting
res_to_plot_narrow=res_to_plot %>% select(c("id.exposure","outcome","b","method","pval","or"))

res_to_plot_narrow=res_to_plot_narrow %>% separate(col="outcome",into="outcome",sep=" \\|\\| ")
res_to_plot_narrow_ivm=res_to_plot_narrow[which(res_to_plot_narrow$method=="Inverse variance weighted"),]
all_res_plot=inter_method_narrow_final=res_206_anno %>% left_join(res_to_plot_narrow_ivm)


col_fun = colorRamp2(c(-0.02,-0.01,-0.005,0,0.005,0.01,0.02), 
                     c("#E85450", "#F08E88","#F8C3BC", 
                       "#FEF7EF",
                       "#C8D8ED","#94BBEC", "#5C9BEB"), 
                     transparency = 0)

sorted_list <- c("L-Cont", "L-Default", "L-DorsAttn", "L-Limbic", "L-SalVentAttn", 
                 "L-SomMot", "L-Vis",
                 paste0("L-",c("Accumbens", "Amygdala", "Caudate", "Hippocampus", "Pallidum", "Putamen", "Thalamus")),
                 "R-Cont", "R-Default", "R-DorsAttn", "R-Limbic", "R-SalVentAttn", 
                 "R-SomMot", "R-Vis",
                 paste0("R-",c("Accumbens", "Amygdala", "Caudate", "Hippocampus", "Pallidum", "Putamen", "Thalamus")))

# Plotting section --------------------------------------------------------

sorted_list <- c(paste0("L-",c("Accumbens", "Amygdala", "Caudate", "Hippocampus", "Pallidum", "Putamen", "Thalamus")),
                 "L-Cont", "L-Default", "L-DorsAttn", "L-Limbic", "L-SalVentAttn", 
                 "L-SomMot", "L-Vis",
                 
                 rev(c(paste0("R-",c("Accumbens", "Amygdala", "Caudate", "Hippocampus", "Pallidum", "Putamen", "Thalamus")),
                       "R-Cont", "R-Default", "R-DorsAttn", "R-Limbic", "R-SalVentAttn", 
                       "R-SomMot", "R-Vis"
                 )
                 )
)

# Data preprocessing function ---------------------------------------------
prepare_data <- function(df) {
  # Define subregions that follow the hemisphere of region1
  follow_region1 <- c("Accumbens", "Amygdala", "Caudate", "Hippocampus", 
                      "Pallidum", "Putamen", "Thalamus")
  
  df %>%
    filter(class == "normal") %>%  # Exclude global data
    mutate(
      # Determine the hemisphere of 'from'
      from_hemisphere = case_when(
        grepl("^L-", region1) ~ "Left",
        grepl("^R-", region1) ~ "Right",
        TRUE ~ NA_character_
      ),
      from_region = gsub("^[LR]-", "", region1),
      from_region=region1,
      # Determine the hemisphere of 'to' (special regions follow region1, others follow region2)
      to_hemisphere = case_when(
        region2_short %in% follow_region1 ~ from_hemisphere,
        grepl("^L-", region2_short) ~ "Left",
        grepl("^R-", region2_short) ~ "Right",
        TRUE ~ NA_character_
      ),
      to_region = ifelse(region2_short == "NANA", "Global",
                         region2_short)
    ) %>%
    filter(!is.na(from_hemisphere) & !is.na(to_hemisphere))  # Ensure valid hemisphere information
}

# Visualization function --------------------------------------------------

p=plot_brain_chord <- function(df, exposure_id,pcut=0.05) {
  # Filter data for a specific gene
  plot_data <- df %>% 
    filter(id.exposure == exposure_id)
  
  # Create connection data frame
  chord_df <- plot_data %>% 
    transmute(
      # from = paste(from_hemisphere, from_region, sep = "_"),
      # to = paste(to_hemisphere, to_region, sep = "_"),
      from = from_region,
      to =to_region,
      value = abs(log(or)),
      p_value = pval,
      direction = ifelse(or > 1, "positive", "negative"),
      
      # Connection type: intra-hemisphere/inter-hemisphere
      connection_type = ifelse(from_hemisphere == to_hemisphere, 
                               "intra-hemisphere", "inter-hemisphere")
    )
  
  # Get all nodes
  all_nodes <- unique(c(chord_df$from, chord_df$to))
  all_nodes=sorted_list
  # Create grouping system
  group_vec <- setNames(
    ifelse(grepl("^L-", all_nodes), "Left", "Right"),
    all_nodes
  )
  
  # Color settings modification 4
  color_palette <- c(
    Left = "#E64B35FF",       # Left hemisphere - warm red (from Nature colors)
    Right = "#357EBD",      # Right hemisphere - classic dark blue
    intra_positive = "#D73027",   # Intra-hemisphere positive - bright cyan blue
    intra_negative = "#56B4E9",   # Intra-hemisphere negative - dark blue gray
    inter_positive = "#F4A261",   # Inter-hemisphere positive - natural green
    inter_negative = "#00A087FF",    # Inter-hemisphere negative - complementary orange-red
    Nosig = "#D3D3D3"       # Non-significant - light gray
  )
  
  # Create connection colors
  chord_df <- chord_df %>%
    mutate(
      link_color = case_when(
        connection_type == "intra-hemisphere" & direction == "positive" ~ color_palette["intra_positive"],
        connection_type == "intra-hemisphere" & direction == "negative" ~ color_palette["intra_negative"],
        connection_type == "inter-hemisphere" & direction == "positive" ~ color_palette["inter_positive"],
        connection_type == "inter-hemisphere" & direction == "negative" ~ color_palette["inter_negative"],
      )
    )
  chord_df$link_color[which(chord_df$p_value >= pcut)]=color_palette["Nosig"]
  # Set graphic parameters
  circos.clear()
  circos.par(
    gap.after = c(
      rep(2, sum(group_vec == "Left")), 
      20,  # Gap between left and right hemispheres
      rep(2, sum(group_vec == "Right"))
    ),
    start.degree = 78+180,
    track.margin = c(0.01, 0.01))
  
  # Draw the base chord diagram
   chordDiagram(
    
    x = chord_df[, c("from", "to", "value")],
    group = group_vec,
    grid.col = setNames(
      ifelse(grepl("^L", all_nodes), 
             color_palette["Left"], 
             color_palette["Right"]),
      all_nodes
    ),
    col = chord_df$link_color,
    transparency = 0.3,
    directional = 1,
    direction.type = "arrows",
    link.arr.type = "big.arrow",
    link.arr.width = 0.15, 
    annotationTrack = "grid",
    preAllocateTracks = list(
      list(track.height = 0.15),  # Main label track
      list(track.height = 0.1)    # P-value track
    ),
    order=sorted_list
  )
  
  # Add region labels (improved version)
  circos.trackPlotRegion(
    track.index = 1,
    bg.border = NA,
    panel.fun = function(x, y) {
      sector <- CELL_META$sector.index
      parts <- strsplit(sector, "-")[[1]]
      region_name <- gsub(" ", "\n", parts[2])  # Automatically handle names with spaces
      
      # Display brain region name
      circos.text(
        CELL_META$xcenter,
        CELL_META$ylim[1]+0.1,
        label = sector,
        col = ifelse(parts[1] == "L", color_palette["Left"], color_palette["Right"]),
        cex = 0.65,
        facing = "clockwise",
        adj = c(0.5, 0),
        niceFacing = TRUE,
        font = 2
      )
    }
  )
  
  # Add legend
  legend("bottomleft", 
         legend = c("Left Hemisphere", "Right Hemisphere"),
         fill = color_palette[c("Left", "Right")],
         border = NA,
         bty = "n",
         cex = 0.8)
  
  legend("bottomright",
         legend = c("Intra-hemisphere +", "Intra-hemisphere -",
                    "Inter-hemisphere +", "Inter-hemisphere -",
                    "Nosig"),
         fill = color_palette[c("intra_positive", "intra_negative",
                                "inter_positive", "inter_negative",
                                "Nosig")],
         border = NA,
         bty = "n",
         cex = 0.8)
  
  # Add title
  title(paste("Gene:", exposure_id), cex.main = 1.2)
}

# Main program ------------------------------------------------------------
# Prepare data
processed_data <- prepare_data(all_res_plot)

# Get all gene types
exposure_ids <- unique(processed_data$id.exposure)

# Batch generate charts
walk(exposure_ids, ~{
  png(filename = paste0("/results/Brain_Connectivity_", gsub("[^A-Za-z0-9]", "_", .x), ".png"), 
      width = 1800, height = 1800, res = 300)
  plot_brain_chord(processed_data, .x)
  dev.off()
})


walk(exposure_ids, ~{
  pdf(paste0("/results/Brain_Connectivity_", gsub("[^A-Za-z0-9]", "_", .x), ".pdf"), 
      width = 7, height = 7)
  plot_brain_chord(processed_data, .x)
  dev.off()
})

walk(exposure_ids, ~{
  png(filename = paste0("/results/Brain_Connectivity_", gsub("[^A-Za-z0-9]", "_", .x), "001.png"), 
      width = 1800, height = 1800, res = 300)
  plot_brain_chord(processed_data, .x,0.01)
  dev.off()
})
walk(exposure_ids, ~{
  pdf(paste0("/results/Brain_Connectivity_", gsub("[^A-Za-z0-9]", "_", .x), "001.pdf"), 
      width = 7, height = 7)
  plot_brain_chord(processed_data, .x,0.01)
  dev.off()
})

