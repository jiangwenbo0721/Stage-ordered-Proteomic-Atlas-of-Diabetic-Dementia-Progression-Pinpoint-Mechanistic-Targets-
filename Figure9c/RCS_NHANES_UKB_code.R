# RCS --- NHANES
library(survival)
library(rms)
library(ggplot2)
library(dplyr)
d <- read_rds('NHANES_50_data_2.rds')
dd <- datadist(d)
options(datadist = "dd")
fit <- cph(Surv(permth_int, Death_of_AD) ~ rcs(total_score, 3) + +ageQ+sex+race+BMI,
           data = d,
           x = TRUE, y = TRUE, surv = TRUE)
# Predicting object, set ref = 0 for HR curve (relative risk)
pred <- Predict(fit, total_score, ref.zero = TRUE, fun = exp)  # fun=exp 得到 HR
# changing dataframe
plot_df <- as.data.frame(pred)

NHANES_fig <- ggplot(plot_df, aes(x = total_score, y = yhat)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "lightblue", alpha = 0.4) +
    geom_line(color = "blue", size = 1.2) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    labs(
        x = "BrainGuard scores",
        y = "Hazard Ratio (HR) for AD Death",
        title = "Restricted Cubic Spline Curve of Total Score (Cox Model)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title = element_text(face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        panel.border = element_blank()
    )
ggsave(filename = 'NHANES_50_RCS_0731.pdf',NHANES_fig,width = 6,height = 5.8,dpi = 300)

#UKB --- RCS
d <- read_rds('ukb_age50_data.rds')

# RCS
# d$total_score <- d$total_score* -1
dd <- datadist(d)
options(datadist = "dd")
# model
fit <- cph(Surv(Dementia_interval, Dementia_type) ~ rcs(total_score1, 4) + 
               age + sex + race + edu + TDI + APOEe4_carrier,
           data = d,
           x = TRUE, y = TRUE, surv = TRUE)
pred <- Predict(fit, total_score1, ref.zero = TRUE, fun = exp)  # HR curve
plot_df <- as.data.frame(pred)
ukb_fig<- ggplot(plot_df, aes(x = total_score1, y = yhat)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "lightblue", alpha = 0.4) +
    geom_line(color = "blue", size = 1.2) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    labs(
        x = "BrainGuard scores",
        y = "Hazard Ratio (HR) for Dementia",
        title = "Restricted Cubic Spline of Total Score (Cox Model)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title = element_text(face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        panel.border = element_blank()
    )
ggsave(filename = 'UKB_RCS.pdf',ukb_fig,width = 6,height = 5.8,dpi = 300)
