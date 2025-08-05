.libPaths(c("/code/brainconn", .libPaths()))
library(brainconn)
library(rio)
library(tidyr)
library(dplyr)
library(combinat)
library(dplyr)
library(circlize)
library(grid)
library(cowplot)

# Define sorted list of subcortical regions
sorted_list <- c("Accumbens", "Amygdala", "Caudate", "Hippocampus", "Pallidum", "Putamen", "Thalamus")
sub_cor <- sorted_list %>% paste0(c("L-","R-"), .)

# List available brain atlases
brainconn::list_atlases()

# Process Schaefer 100-parcel atlas with 7 networks
atlas1 <- schaefer100_n7
atlas1$X <- paste0(atlas1$hemi,"-",atlas1$network) %>% 
  gsub("LH-","L-", .) %>% 
  gsub("RH-","R-", .) %>% 
  gsub(" control","Cont", .) %>% 
  gsub(" somatomotor","SomMot", .) %>% 
  gsub(" default mode","Default", .) %>% 
  gsub(" dorsal attention","DorsAttn", .) %>% 
  gsub(" limbic","Limbic", .) %>% 
  gsub(" visual","Vis", .) %>% 
  gsub(" salience/ventral attention","SalVentAttn", .) %>% 
  gsub(" network","", .)

# Process DK82 atlas for subcortical regions
atlas2 <- dk82_aspree[which(dk82_aspree$X %in% sorted_list),]
atlas2$X <- paste0(atlas2$hemi,"-",atlas2$X)

# Combine both atlases
self_atlas <- rbind(atlas1, atlas2[,which(colnames(atlas2) %in% colnames(atlas1))])
self_atlas$index <- 1:nrow(self_atlas)
self_atlas$x.mni <- as.integer(self_atlas$x.mni)
self_atlas$y.mni <- as.integer(self_atlas$y.mni)
self_atlas$z.mni <- as.integer(self_atlas$z.mni)
rownames(self_atlas) <- 1:nrow(self_atlas)

# Validate the combined atlas
check_atlas(self_atlas)

# Load annotation data
load("/res_206_anno.Rdata")
res_206_anno_new <- res_206_anno
res_206_anno_new <- res_206_anno_new %>% separate("region1", c("region1_hemi"), sep="-", remove = F)

# Define subcortical regions
subcor <- c("Accumbens", "Amygdala", "Caudate", "Hippocampus", "Pallidum", "Putamen", "Thalamus")
n_to_sub <- which(res_206_anno_new$region2_short %in% subcor)

# Add hemisphere prefix to subcortical regions based on region1's hemisphere
res_206_anno_new$region2_short[n_to_sub] <- res_206_anno_new$region2_short[n_to_sub] %>% 
  paste0(res_206_anno_new$region1_hemi[n_to_sub], "-", .)

# Function to generate brain connectivity visualization with annotations
plot_brain_anno <- function(plot2_chor_plot, ref) {
  # Prepare annotation data for region1 and region2
  self_atlas_anno1 <- self_atlas %>% select(c("X", "index"))
  self_atlas_anno2 <- self_atlas %>% select(c("X", "index"))
  colnames(self_atlas_anno1) <- c("region1", "index1")
  colnames(self_atlas_anno2) <- c("region2_short", "index2")
  
  # Filter significant connections
  plot2_chor_plot_sig <- plot2_chor_plot[which(plot2_chor_plot$if_sig == "True"),] %>% 
    select(c("region1", "region2_short", "or"))
  
  # Map regions to indices and create connectivity matrix
  plot2_chor_plot_inter <- plot2_chor_plot_sig %>% 
    left_join(self_atlas_anno1) %>% 
    left_join(self_atlas_anno2) %>% 
    select(c("index1", "index2", "or"))
  
  # Initialize connectivity matrix
  cor_matrix <- matrix(0, nrow = 214, ncol = 214)
  
  # Populate the connectivity matrix (symmetric)
  for (i in 1:nrow(plot2_chor_plot_inter)) {
    row_index <- plot2_chor_plot_inter$index1[i]
    col_index <- plot2_chor_plot_inter$index2[i]
    cor_value <- plot2_chor_plot_inter$or[i]
    cor_matrix[row_index, col_index] <- cor_value
    cor_matrix[col_index, row_index] <- cor_value # Assume symmetric matrix
  }
  
  # Generate brain connectivity plot
  p1 <<- brainconn2(self_atlas, conmat = cor_matrix, node.size = 3, view = "ortho", edge.color.weighted = T)
  
  return(p1)
}

# Load and preprocess results data
res_to_plot <- rio::import("data.csv")
# Extract columns for FDR plotting
res_to_plot_narrow <- res_to_plot %>% select(c("id.exposure", "outcome", "method", "pval", "or", "p_fdr"))

# Split outcome column
res_to_plot_narrow <- res_to_plot_narrow %>% separate(col = "outcome", into = "outcome", sep = " \\|\\| ")

# Get unique exposures
all_exopsure <- unique(res_to_plot_narrow$id.exposure)

# Calculate maximum scale value for visualization
max_scale_value <- max(abs(log(res_to_plot_narrow$or[which(res_to_plot_narrow$pval < 0.05 & res_to_plot_narrow$method == "Inverse variance weighted")])))

# Loop through each exposure and generate plots
for (i in 1:length(all_exopsure)) {
  inter_exposure <- all_exopsure[i]
  inter_res_narrow <- res_to_plot_narrow[which(res_to_plot_narrow$id.exposure == inter_exposure),]
  
  # Get unique methods for current exposure
  all_methods <- inter_res_narrow$method %>% unique()
  inter_method <- "Inverse variance weighted"
  inter_method_narrow <- inter_res_narrow[which(inter_res_narrow$method == inter_method),]
  
  # Join with annotation data
  inter_method_narrow_n <- inter_method_narrow %>% select(c("outcome", "or", "pval", "p_fdr"))
  inter_method_narrow_final <- res_206_anno_new %>% left_join(inter_method_narrow_n)
  
  # Flag significant results (p < 0.05)
  inter_method_narrow_final$if_sig <- case_when(inter_method_narrow_final$pval < 0.05 ~ "True",
                                                .default = "False")
  
  # Calculate FDR-adjusted p-values
  inter_method_narrow_final$p_fdr_new <- inter_method_narrow_final$pval %>% p.adjust("fdr")
  
  # Print progress information
  print(all_exopsure[i])
  print(paste0(length(which(inter_method_narrow_final$pval < 0.05)), "/", length(which(inter_method_narrow_final$p_fdr_new > 0))))
  
  # Generate plots if there are both significant and non-significant results
  if (length(unique(inter_method_narrow_final$if_sig)) > 1) {
    inter_method_narrow_final_f <- inter_method_narrow_final[which(inter_method_narrow_final$if_sig == "True"),]
    plot1_bar_plot <- inter_method_narrow_final[which(inter_method_narrow_final$class == "global"),]
    plot2_chor_plot <- inter_method_narrow_final[which(inter_method_narrow_final$class == "normal"),]
    
    # Generate reference name and save plot
    ref <- paste0(inter_exposure, "_", inter_method) %>% gsub("finngen_R11_", "", .)
    p2 <- plot_brain_anno(plot2_chor_plot, ref)
    ggsave(p2, file = paste0(ref, '_brainconn.pdf'), width = 8, height = 8)
  }
}