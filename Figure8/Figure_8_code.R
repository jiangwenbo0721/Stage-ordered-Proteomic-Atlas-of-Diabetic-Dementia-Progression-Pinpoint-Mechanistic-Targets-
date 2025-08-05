# Add to R package--------------------------------------------------------------------
# Figure 8 a
library(tidyverse)
library(ggrepel)
library(ggplot2)
# data
d <- read_csv('/code/Figure8/Figure_8_a_data.csv')
colnames(d)
head(d)

d <- d %>%
    group_by(PhecodeCategory) %>%
    mutate(
        pfdr = p.adjust(pvalue, method = "fdr"),
    ) %>%
    ungroup()
# Standardization + Group Processing
d_scaled <- d %>%
    group_by(PhecodeCategory) %>%
    mutate(coef = scale(coef)) %>%
    ungroup() %>%
    mutate(
        point_type = case_when(
            pfdr >= 0.05 ~ "ns",
            coef > 0 ~ "pos",
            TRUE ~ "neg"
        )
    )
d_scaled$PhecodeString
# 1. Replace the underscores with spaces.
d_scaled$PhecodeString <- gsub("_0", " ", d_scaled$PhecodeString)

# 2.Remove the trailing ".0" (if it exists).
d_scaled$PhecodeString <- sub("_", "", d_scaled$PhecodeString)


label_data <- d_scaled %>%
    filter(pfdr < 0.05) %>%
    group_by(PhecodeCategory) %>%
    slice_max(coef, n = 3, with_ties = FALSE) %>% 
    bind_rows(
        d_scaled %>%
            filter(pfdr < 0.05) %>% 
            group_by(PhecodeCategory) %>%
            slice_min(coef, n = 3, with_ties = FALSE)
    ) %>%
    ungroup() %>%
    distinct()  


label_data$coef <- as.numeric(label_data[[3]])


shape_mapping <- c(
    "APOM"   = 15,
    "CD276"  = 16,
    "EGFR"   = 17,
    "GDF15"  = 18,
    "GHRL"   = 21,
    "LGALS4" = 22
)

# Draw a picture
p=ggplot(d_scaled, aes(
    x = PhecodeCategory,
    y = coef,
    size = abs(coef),
    alpha = abs(coef), 
    color = point_type, 
    shape = Row_Name  
)) +
    geom_jitter(width = 0.45, height = 0) +
    geom_text_repel(
        data = label_data,
        aes(label = PhecodeString),
        size = 3.5,
        max.overlaps = 100,
        segment.color = "grey50",
        show.legend = FALSE
    ) +
    scale_color_manual(
        values = c("ns" = "lightgrey", "pos" = '#ed3e2e', "neg" = "#1196b5")
    ) +
    scale_shape_manual(values = shape_mapping)+
    scale_size_continuous(range = c(1, 10)) + 
    scale_alpha_continuous(range = c(0.7, 3)) +  
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
    theme_classic() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),  
        legend.position = "right"
    )
ggsave('/results/figure8a.pdf', plot = p, height = 13, width = 10, units = "in")

# Figure 8 b --------------------------------------------------------------

library(readxl)
# Read MAEA data from Excel file
df_APOM <- read_excel("/code/Figure8/figure_8_b_data.xlsx", sheet = "APOM")
df_CD276 <- read_excel("/code/Figure8/figure_8_b_data.xlsx", sheet = "CD276")
df_EGFR <- read_excel("/code/Figure8/figure_8_b_data.xlsx", sheet = "EGFR")
df_GDF15 <- read_excel("/code/Figure8/figure_8_b_data.xlsx", sheet = "GDF15")
df_GHRL <- read_excel("/code/Figure8/figure_8_b_data.xlsx", sheet = "GHRL")
df_LGALS4 <- read_excel("/code/Figure8/figure_8_b_data.xlsx", sheet = "LGALS4")

table(df_APOM$Category)

# Define the order of categories for visualization
levels_order <- c(
    "Biochemical indicators", "Body exam",'Condition of illness', "Food intake", "Food liking",
    "Psychosocial factors", "Lifestyles"
)

# Process CEP164 data:
# 1. Convert Category to factor with specified order
# 2. Arrange data by Category and outcome1
# 3. Convert outcome1 to factor with unique levels
df_APOM <- df_APOM %>%
    mutate(
        Category = factor(Category, levels = levels_order)
    ) %>%
    arrange(Category, outcome1) %>%
    mutate(
        outcome1 = factor(outcome1, levels = unique(outcome1))
    )


# Create bubble plot for CEP164 data
APOM_P<- ggplot(df_APOM, aes(x =exposure , y = outcome1,size = abs(beta_std),fill=beta_std))+#global
    # Draw bubbles with diamond shape (shape 23)
    geom_point(shape = 23 )+
    # Set bubble size range
    scale_size_continuous(range = c(2, 7))+
    # Set color gradient from light orange to dark red
    scale_fill_gradient(low = "#1196b5", high = "#ed3e2e",
                        guide = guide_colorbar(title.position = "top",
                                               order = 1)) + 
    # Position x-axis labels at top
    scale_x_discrete(position = "top") + 
    
    # Set plot labels and theme
    labs(
        x = ""
    ) + theme_bw() +
    # Customize theme settings
    theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5),
          # Remove panel border
          panel.border = element_blank())

# Display CEP164 plot
ggsave('/results/figure 8b_APOM_P.pdf', plot = APOM_P, height = 13, width = 10, units = "in")
# Process CEP164 data:
# 1. Convert Category to factor with specified order
# 2. Arrange data by Category and outcome1
# 3. Convert outcome1 to factor with unique levels
df_GHRL <- df_GHRL %>%
    mutate(
        Category = factor(Category, levels = levels_order)
    ) %>%
    arrange(Category, outcome1) %>%
    mutate(
        outcome1 = factor(outcome1, levels = unique(outcome1))
    )


# Create bubble plot for CEP164 data
GHRL_P<- ggplot(df_GHRL, aes(x =exposure , y = outcome1,size = abs(beta_std),fill=beta_std))+#global
    # Draw bubbles with diamond shape (shape 23)
    geom_point(shape = 23 )+
    # Set bubble size range
    scale_size_continuous(range = c(2, 7))+
    # Set color gradient from light orange to dark red
    scale_fill_gradient(low = "#1196b5", high = "#ed3e2e",
                        guide = guide_colorbar(title.position = "top",
                                               order = 1)) + 
    # Position x-axis labels at top
    scale_x_discrete(position = "top") + 
    
    # Set plot labels and theme
    labs(
        x = ""
    ) + theme_bw() +
    # Customize theme settings
    theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5),
          # Remove panel border
          panel.border = element_blank())

# Display CEP164 plot
ggsave('/results/figure 8b_GHRL_P.pdf', plot = GHRL_P, height = 13, width = 10, units = "in")
# Process CEP164 data:
# 1. Convert Category to factor with specified order
# 2. Arrange data by Category and outcome1
# 3. Convert outcome1 to factor with unique levels
df_GDF15 <- df_GDF15 %>%
    mutate(
        Category = factor(Category, levels = levels_order)
    ) %>%
    arrange(Category, outcome1) %>%
    mutate(
        outcome1 = factor(outcome1, levels = unique(outcome1))
    )


# Create bubble plot for CEP164 data
GDF15_P<- ggplot(df_GDF15, aes(x =exposure , y = outcome1,size = abs(beta_std),fill=beta_std))+#global
    # Draw bubbles with diamond shape (shape 23)
    geom_point(shape = 23 )+
    # Set bubble size range
    scale_size_continuous(range = c(2, 7))+
    # Set color gradient from light orange to dark red
    scale_fill_gradient(low = "#1196b5", high = "#ed3e2e",
                        guide = guide_colorbar(title.position = "top",
                                               order = 1)) + 
    # Position x-axis labels at top
    scale_x_discrete(position = "top") + 
    
    # Set plot labels and theme
    labs(
        x = ""
    ) + theme_bw() +
    # Customize theme settings
    theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5),
          # Remove panel border
          panel.border = element_blank())

# Display GDF15_P plot
ggsave('/results/figure 8b_GDF15_P.pdf', plot = GDF15_P, height = 13, width = 10, units = "in")
# Process CEP164 data:
# 1. Convert Category to factor with specified order
# 2. Arrange data by Category and outcome1
# 3. Convert outcome1 to factor with unique levels
df_EGFR <- df_EGFR %>%
    mutate(
        Category = factor(Category, levels = levels_order)
    ) %>%
    arrange(Category, outcome1) %>%
    mutate(
        outcome1 = factor(outcome1, levels = unique(outcome1))
    )


# Create bubble plot for CEP164 data
EGFR_P<- ggplot(df_EGFR, aes(x =exposure , y = outcome1,size = abs(beta_std),fill=beta_std))+#global
    # Draw bubbles with diamond shape (shape 23)
    geom_point(shape = 23 )+
    # Set bubble size range
    scale_size_continuous(range = c(2, 7))+
    # Set color gradient from light orange to dark red
    scale_fill_gradient(low = "#1196b5", high = "#ed3e2e",
                        guide = guide_colorbar(title.position = "top",
                                               order = 1)) + 
    # Position x-axis labels at top
    scale_x_discrete(position = "top") + 
    
    # Set plot labels and theme
    labs(
        x = ""
    ) + theme_bw() +
    # Customize theme settings
    theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5),
          # Remove panel border
          panel.border = element_blank())

# Display EGFR_P plot
ggsave('/results/figure 8b_EGFR_P.pdf', plot = EGFR_P, height = 13, width = 10, units = "in")
# Process CEP164 data:
# 1. Convert Category to factor with specified order
# 2. Arrange data by Category and outcome1
# 3. Convert outcome1 to factor with unique levels
df_APOM <- df_APOM %>%
    mutate(
        Category = factor(Category, levels = levels_order)
    ) %>%
    arrange(Category, outcome1) %>%
    mutate(
        outcome1 = factor(outcome1, levels = unique(outcome1))
    )


# Create bubble plot for CEP164 data
APOM_P<- ggplot(df_APOM, aes(x =exposure , y = outcome1,size = abs(beta_std),fill=beta_std))+#global
    # Draw bubbles with diamond shape (shape 23)
    geom_point(shape = 23 )+
    # Set bubble size range
    scale_size_continuous(range = c(2, 7))+
    # Set color gradient from light orange to dark red
    scale_fill_gradient(low = "#1196b5", high = "#ed3e2e",
                        guide = guide_colorbar(title.position = "top",
                                               order = 1)) + 
    # Position x-axis labels at top
    scale_x_discrete(position = "top") + 
    
    # Set plot labels and theme
    labs(
        x = ""
    ) + theme_bw() +
    # Customize theme settings
    theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5),
          # Remove panel border
          panel.border = element_blank())

# Display CEP164 plot
ggsave('/results/figure 8b_APOM_P.pdf', plot = APOM_P, height = 13, width = 10, units = "in")

df_CD276 <- df_CD276 %>%
    mutate(
        Category = factor(Category, levels = levels_order)
    ) %>%
    arrange(Category, outcome1) %>%
    mutate(
        outcome1 = factor(outcome1, levels = unique(outcome1))
    )


# Create bubble plot for CEP164 data
CD276_P<- ggplot(df_CD276, aes(x =exposure , y = outcome1,size = abs(beta_std),fill=beta_std))+#global
    # Draw bubbles with diamond shape (shape 23)
    geom_point(shape = 23 )+
    # Set bubble size range
    scale_size_continuous(range = c(2, 7))+
    # Set color gradient from light orange to dark red
    scale_fill_gradient(low = "#1196b5", high = "#ed3e2e",
                        guide = guide_colorbar(title.position = "top",
                                               order = 1)) + 
    # Position x-axis labels at top
    scale_x_discrete(position = "top") + 
    
    # Set plot labels and theme
    labs(
        x = ""
    ) + theme_bw() +
    # Customize theme settings
    theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5),
          # Remove panel border
          panel.border = element_blank())

# Display CD276_P plot
ggsave('/results/figure 8b_CD276_P.pdf', plot = CD276_P, height = 13, width = 10, units = "in")
