require(tidyverse)
require(rms)
require(survival)
dementia_with_coef_for_cox_df <- read_tsv("mhas_rcs_data.tsv.gz")
dd <- datadist(dementia_with_coef_for_cox_df)
options(datadist = "dd")
rcs_fit <- cph(
  Surv(dementia_with_coef_for_cox_df$dementia_interval, dementia_with_coef_for_cox_df$dementia_status) ~
    rcs(comprehensive_score, 3) + bmi + age + gender,
  data = dementia_with_coef_for_cox_df
)
anova(rcs_fit)
# 查看hr预测表
rcs_hr <- Predict(rcs_fit, comprehensive_score, fun = exp)
# 绘图
ggplot() +
  geom_line(data = rcs_hr, aes(comprehensive_score, yhat),
            linetype = "solid", size = 1, alpha = 0.7, colour = "#0070b9") +
  geom_ribbon(data = rcs_hr,
              aes(comprehensive_score, ymin = lower, ymax = upper),
              alpha = 0.1, fill = "#0070b9") +
  theme_classic() +
  geom_hline(yintercept = 1, linetype = 2, size = 1) +
  labs(title = "dementia risk", x = "comprehensive_score", y = "HR (95%CI)")



