library(survival)
library(rms)
library(ggplot2)
library(dplyr)
df_filtered<-read.csv("rcs_data_charls_0727.csv")
df_filtered$sex<-as.factor(df_filtered$sex)
df_filtered$age_quartile<-as.factor(df_filtered$age_quartile)
# 
surv_obj <- Surv(time = df_filtered$time_diff_days_new, event = df_filtered$dementia_status)



df_model <- df_filtered %>% 
    select(total_score, time_diff_days_new, dementia_status, sex, age_quartile, bmi)

dd <- datadist(df_model)
options(datadist = "dd")


fit <- cph(
    surv_obj ~ rcs(total_score, 4) + sex + age_quartile + bmi,
    data = df_model,
    x = TRUE,
    y = TRUE
)

#  p for overall 和 p for linear
anova_fit <- anova(fit)
p_overall   <- signif(anova_fit[1, "P"], 1)
p_nonlinear <- signif(anova_fit[2, "P"], 2)

#  HR 
pred <- Predict(fit, total_score, ref.zero = TRUE, fun = exp)
plot_df <- as.data.frame(pred)

# plot
ukb_fig <- ggplot(plot_df, aes(x = total_score, y = yhat)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "lightblue", alpha = 0.4) +
    geom_line(color = "blue", size = 1.2) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    labs(
        x = "Total Score",
        y = "Hazard Ratio (HR) for Dementia",
        title = "Restricted Cubic Spline of Total Score (Cox Model)"
    ) +
    annotate("text", 
             x = min(plot_df$total_score) + 0.1 * diff(range(plot_df$total_score)), 
             y = max(plot_df$yhat) * 0.95,
             label = paste0("P for overall = ", p_overall, 
                            "\nP for nonlinear = ", p_nonlinear),
             hjust = 0, size = 5, fontface = "italic") +
    theme_minimal(base_size = 14) +
    theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title = element_text(face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        panel.border = element_blank()
    )

# 打印图形
ggsave("RCS_TotalScore_HR_charls_0727.pdf", plot = ukb_fig, width = 2.66, height = 2.99, units = "in")
