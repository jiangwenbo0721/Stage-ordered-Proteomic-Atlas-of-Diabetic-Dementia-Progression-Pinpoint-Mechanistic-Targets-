# Add to R package
library(tibble)
library(dplyr)
library(ggplot2)
library(forcats)
library(readxl)
d <- read_xlsx('/code/Figure9/Figure9b/Figure_9_data.xlsx')

library(dplyr)
library(tidyr)
library(stringr)

# Extract numerical values
d_clean <- d %>%
    mutate(
        HR = if_else(str_detect(`HR (95%CI)`, "Reference"), NA_character_, `HR (95%CI)`),
        HR = str_remove_all(HR, " "), 
        HR_value = as.numeric(str_extract(HR, "^[0-9.]+")),
        lower_CI = as.numeric(str_match(HR, "\\(([^,]+),")[, 2]),
        upper_CI = as.numeric(str_match(HR, ",([^)]+)\\)")[, 2])
    )

d_clean <- d_clean %>%
    mutate(
        group = case_when(
            database %in% c("Q1", "Q2", "Q3", "Q4") ~ NA_character_,
            TRUE ~ database
        )
    ) %>%
    fill(group) %>%
    mutate(
        Quartile = if_else(database %in% c("Q1", "Q2", "Q3", "Q4"), database, "Overall"),
        label = paste(group, Quartile, sep = "-")
    )

library(ggplot2)
library(forcats)

# Set factor order reversal
d_clean <- d_clean %>%
    mutate(label = fct_rev(factor(label, levels = unique(label))))

# Draw a forest map
p=ggplot(d_clean %>% filter(!is.na(HR_value)), aes(x = HR_value, y = label)) +
    geom_point(shape = 18, size = 3, color = "steelblue") +
    geom_errorbarh(aes(xmin = lower_CI, xmax = upper_CI), height = 0.2, color = "steelblue") +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
    scale_x_continuous(trans = "log10", limits = c(0.2, 1.5)) +
    labs(x = "Hazard Ratio (log scale)", y = NULL, title = "Forest Plot of CVD/T2D Risk by Quartile") +
    theme_minimal(base_size = 14)
ggsave("/results/figure9b.pdf", plot = p, height = 6, width = 13, units = "in")
