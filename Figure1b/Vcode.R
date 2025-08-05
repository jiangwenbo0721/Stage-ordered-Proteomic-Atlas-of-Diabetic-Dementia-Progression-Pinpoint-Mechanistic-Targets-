# Load required R packages
library(tidyverse)            # For data manipulation and visualization
library(ggrepel)              # For adding labels with repulsion
library(cowplot)              # For combining plots (alternative to patchwork)

# Load the data file
load('/code/Figure1/Figure1b/Vdata.Rdata') # Load the Rdata file containing necessary data and variables

# Select top 10 genes by absolute logFC for each cluster
clabel = mresult %>%
  group_by(cluster) %>%
  arrange(desc(abs(logFC))) %>% # Sort by absolute logFC value
  dplyr::slice(1:10) # Take top 10 rows per cluster

# Initialize list for storing plots
spic = list()

# Create volcano plots for each cluster
for(f01 in sname) {
  sdata = mresult %>% filter(jcluster == f01)
  snum = which(sname == f01)
  
  p01 <- ggplot() +
    # Points for non-significant genes
    geom_point(data = sdata[sdata$group == 'not',], alpha = 0.4, fill = '#CECECE', 
               shape = 21, stroke = .3,
               aes(y = logFC, x = -log10(adj.P.Val), size = abs(logFC))) +
    
    # Points for significant genes
    geom_point(data = subset(sdata, group != 'not'), alpha = 0.8, 
               shape = 21, stroke = .3,
               aes(y = logFC, x = -log10(adj.P.Val), fill = abs(logFC), size = abs(logFC))) +
    
    # Size scaling
    scale_size_continuous(limits = c(0, 1.4), range = c(.2, 6.5)) +
    
    # Color gradient
    scale_fill_gradientn(colors = c('white', wcolor[snum]),
                         limits = c(0, 1.4)) +
    scale_color_manual(values = wcolor[snum]) +
    
    # Reference lines
    geom_hline(yintercept = c(-.1,0.1), lty = 4, col = "black", lwd = 0.5) +
    
    # Axis labels and title
    labs(x = "-log10(P.adj)", y = "log2FC", size = '|log2FC|', 
         title = f01) +
    
    # Coordinate limits
    coord_cartesian(ylim = c(-1.3, 1.5), xlim = c(0, 321)) +
    
    # Legend adjustments
    guides(fill = "none", color = 'none', size = 'none') +
    
    # Theme settings
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.background = element_rect(
            color="black", fill=wcolor[snum]
          )) +
    
    # Gene labels
    geom_text_repel(
      data = clabel %>% filter(jcluster == f01),           
      aes(y = logFC, x = -log10(adj.P.Val), label = gene),
      color = 'grey50',
      box.padding = 0.3,
      size = 4
    )
  
  spic[[f01]] = p01
}

# Combine plots into a grid
vsplot = plot_grid(spic[[1]], spic[[2]], spic[[3]], spic[[4]],
                   ncol = 4, align = 'v', axis = 'r')

# Save the final plot
ggsave('/results/Figure1b.pdf', plot = vsplot,  height = 6, width = 13, units = "in")
