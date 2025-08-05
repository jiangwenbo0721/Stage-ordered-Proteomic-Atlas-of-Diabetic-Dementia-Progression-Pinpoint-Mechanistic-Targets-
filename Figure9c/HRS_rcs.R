setwd("/formal/a1 HRS/res")

library(dplyr)
library(tidyr)
library(survival)
library(zoo)
library(lubridate)

# process_raw_data --------------------------------------------------------

# Load raw HRS data
load("/formal/a1 HRS/data/raw/hrs_raw.Rdata")

# part1_score -------------------------------------------------------------
# Custom function: Find first non-NA value in subsequent waves
get_next_non_na <- function(values) {
  non_na_index <- which(!is.na(values))
  if (length(non_na_index) == 0) return(NA)
  return(values[non_na_index[1]])
}

# Process behavioral data
behavior_base <- hrs %>%
  select(hhidpn, wave, diabe, rxdiabo, hba1c, cesd, lblonely3, met) %>%
  arrange(hhidpn, wave)

# Step 1: Identify baseline and count total waves per participant
behavior_step1 <- behavior_base %>%
  group_by(hhidpn) %>%
  mutate(
    is_baseline = ifelse(wave == min(wave), 1, 0),
    total_waves = n()
  ) %>%
  ungroup()

# Step 2: Process glycemic indicators (hba1c/diabe)
behavior_step2 <- behavior_step1 %>%
  group_by(hhidpn) %>%
  mutate(
    # Find first non-NA hba1c value in subsequent waves (original logic preserved)
    next_hba1c = ifelse(is_baseline == 1 & is.na(hba1c),
                        hba1c[which(!is.na(hba1c) & wave > wave[is_baseline == 1])[1]],
                        hba1c),
    
    # New glycemic scoring system
    glycemic_score = case_when(
      # No diabetes (hba1c<6.5 and diabe=0)
      (!is.na(next_hba1c) & next_hba1c < 6.5) | 
        (!is.na(diabe) & diabe == 0) ~ 4,
      
      # With diabetes (hba1c≥6.5 or diabe=1)
      (!is.na(next_hba1c) & next_hba1c >= 6.5) | 
        (!is.na(diabe) & diabe == 1) ~ 
        ifelse(!is.na(rxdiabo) & rxdiabo == 1, 2, 0),
      
      # Other cases (both indicators missing)
      TRUE ~ NA_real_
    ),
    
    # Add glycemic status flag (optional)
    glycemic_status = case_when(
      (!is.na(next_hba1c) & next_hba1c >= 6.5) | 
        (!is.na(diabe) & diabe == 1) ~ "Diabetic",
      TRUE ~ "Non-diabetic"
    )
  ) %>%
  ungroup()

# Step 3: Process CESD scores (depression scale)
# Step 3.1: Calculate quartile thresholds by participant (only cesd>0 observations)
cesd_quantiles <- behavior_step2 %>%
  filter(cesd > 0) %>%
  group_by(hhidpn) %>%
  summarise(
    q25 = quantile(cesd, probs = 0.25, na.rm = TRUE),
    q50 = quantile(cesd, probs = 0.50, na.rm = TRUE),
    q75 = quantile(cesd, probs = 0.75, na.rm = TRUE)
  )

# Step 3.2: Merge back and calculate mood scores
behavior_step3 <- behavior_step2 %>%
  # Find first non-NA cesd value in subsequent waves (original imputation logic preserved)
  group_by(hhidpn) %>%
  mutate(
    next_cesd = ifelse(is_baseline == 1 & is.na(cesd),
                       cesd[which(!is.na(cesd) & wave > wave[is_baseline == 1])[1]],
                       cesd)
  ) %>%
  ungroup() %>%
  # Merge quartile thresholds
  left_join(cesd_quantiles, by = "hhidpn") %>%
  # Calculate new mood scores
  mutate(
    mood_score = case_when(
      is.na(next_cesd) ~ NA_real_,
      next_cesd == 0 ~ 0,  # 0 points gets 0
      next_cesd > 0 & next_cesd <= q25 ~ -1,
      next_cesd > 0 & next_cesd <= q50 ~ -2,
      next_cesd > 0 & next_cesd <= q75 ~ -3,
      next_cesd > 0 & next_cesd > q75 ~ -4
    ),
    # Remove temporary columns (optional)
    across(c(q25, q50, q75), ~NULL)
  )

# Step 4: Process loneliness scores
behavior_step4 <- behavior_step3 %>%
  group_by(hhidpn) %>%
  mutate(
    # Original next-wave imputation logic preserved
    next_lonely = ifelse(is_baseline == 1 & is.na(lblonely3),
                         lblonely3[which(!is.na(lblonely3) & wave > wave[is_baseline == 1])[1]],
                         lblonely3),
    
    # Mathematical linear transformation: 3 - next_lonely
    social_score = ifelse(!is.na(next_lonely), 
                          round(3 - next_lonely, digits = 2),  # Keep 2 decimal places
                          NA_real_),
    
    # Ensure score is within 0-3 range (handle edge cases)
    social_score = pmax(0, pmin(3, social_score)),  # Clamp to [0,3] range
    
    social_source = "Loneliness_scale"
  ) %>%
  ungroup()

# Step 5: Process physical activity scores
behavior_step5 <- behavior_step4 %>%
  group_by(hhidpn) %>%
  mutate(
    # Original next-wave imputation logic preserved
    next_met = ifelse(is_baseline == 1 & is.na(met),
                      met[which(!is.na(met) & wave > wave[is_baseline == 1])[1]],
                      met)
  ) %>%
  ungroup() %>%
  # Calculate global quartiles (excluding NA)
  mutate(
    met_q25 = quantile(met, probs = 0.25, na.rm = TRUE),
    met_q50 = quantile(met, probs = 0.50, na.rm = TRUE),
    met_q75 = quantile(met, probs = 0.75, na.rm = TRUE)
  ) %>%
  # Assign quartile scores
  mutate(
    activity_score = case_when(
      is.na(next_met) ~ NA_real_,
      next_met <= met_q25 ~ 1,
      next_met <= met_q50 ~ 2,
      next_met <= met_q75 ~ 3,
      next_met > met_q75 ~ 4
    ),
    # Add activity level labels (optional)
    activity_level = case_when(
      activity_score == 1 ~ "Q1 (Lowest)",
      activity_score == 2 ~ "Q2",
      activity_score == 3 ~ "Q3",
      activity_score == 4 ~ "Q4 (Highest)",
      TRUE ~ "Missing"
    ),
    # Remove temporary quartile columns
    across(starts_with("met_q"), ~NULL)
  )

# Final behavior score calculation
behavior_final <- behavior_step5 %>%
  filter(is_baseline == 1) %>%
  mutate(
    # Count valid components
    valid_components = 4 - is.na(glycemic_score) - is.na(mood_score) - 
      is.na(social_score) - is.na(activity_score),
    
    # Calculate total score (requires at least 3 valid components)
    behavior_score = ifelse(valid_components >= 3,
                            rowSums(select(., glycemic_score, mood_score, 
                                           social_score, activity_score),
                                    na.rm = FALSE),  # Preserve NA
                            NA_real_),
    
    # Group labels
    behavior_group = case_when(
      behavior_score == 4 ~ "Optimal",
      behavior_score >= 2 ~ "Moderate",
      behavior_score < 2 ~ "Poor",
      TRUE ~ "Missing"
    )
  ) %>%
  select(hhidpn, wave, 
         glycemic_score, mood_score, social_score, activity_score,
         valid_components, behavior_score, behavior_group)

# Save final behavior data
save(behavior_final, file="/formal/a1 HRS/data/processed/behavior_final.Rdata")
load("/formal/a1 HRS/data/processed/behavior_final.Rdata")

# part2_dementia ----------------------------------------------------------
# Step 1: Define functional impairment (ADLs) (unchanged)
hrs_clean <- hrs %>%
  mutate(
    across(c(dressa, batha, eata, beda, toilta, walkra), ~ ifelse(. == 1, 1, 0)),
    adl_total = rowSums(select(., dressa, batha, eata, beda, toilta, walkra), na.rm = TRUE),
    adl_impair = ifelse(adl_total >= 1, 1, 0)
  )

# Step 2: Build education-specific cognitive reference values (grouped by raeducl)
cognitive_norms <- hrs_clean %>%
  filter(ragey_b >= 50) %>%  # Only use population aged 50+ to establish norms
  group_by(raeducl) %>%      # Group by original education codes
  summarise(
    memory_mean = mean(rowMeans(cbind(imrc, dlrc), na.rm = TRUE), na.rm = TRUE),
    memory_sd = sd(rowMeans(cbind(imrc, dlrc), na.rm = TRUE), na.rm = TRUE),
    exec_mean = mean(rowMeans(cbind(ser7, bwc20), na.rm = TRUE), na.rm = TRUE),
    exec_sd = sd(rowMeans(cbind(ser7, bwc20), na.rm = TRUE), na.rm = TRUE),
    orient_mean = mean(orient, na.rm = TRUE),
    orient_sd = sd(orient, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    # Cognitive impairment threshold (mean - 1.5SD)
    memory_threshold = memory_mean - 1.5 * memory_sd,
    exec_threshold = exec_mean - 1.5 * exec_sd,
    orient_threshold = orient_mean - 1.5 * orient_sd,
    
    # Cognitive improvement threshold (1SD change)
    memory_improve = memory_sd,
    exec_improve = exec_sd,
    orient_improve = orient_sd
  )

# Step 3: Apply cognitive impairment definition
hrs_clean <- hrs_clean %>%
  left_join(cognitive_norms, by = "raeducl") %>%
  mutate(
    # Calculate current cognitive scores
    memory_score = rowMeans(cbind(imrc, dlrc), na.rm = TRUE),
    exec_score = rowMeans(cbind(ser7, bwc20), na.rm = TRUE),
    orient_score = orient,
    
    # Flag impaired domains
    memory_impaired = ifelse(!is.na(memory_score), memory_score < memory_threshold, FALSE),
    exec_impaired = ifelse(!is.na(exec_score), exec_score < exec_threshold, FALSE),
    orient_impaired = ifelse(!is.na(orient_score), orient_score < orient_threshold, FALSE),
    
    # Cognitive impairment definition (≥2 impaired domains)
    cog_impair = ifelse(
      (memory_impaired + exec_impaired + orient_impaired) >= 2, 1, 0)
  )

# Step 4: Process physician diagnosis
hrs_clean <- hrs_clean %>%
  mutate(
    doctor_diagnosis = ifelse(demen == 1 | alzhe == 1 | demene == 1 | alzhee == 1, 1, 0)
  )

# Step 5: Define initial dementia cases
hrs_clean <- hrs_clean %>%
  mutate(
    dementia_initial = ifelse(
      (adl_impair == 1 & cog_impair == 1) | doctor_diagnosis == 1,
      1, 0)
  )

# Step 6: Dynamic exclusion rules (education-stratified version)
hrs_clean <- hrs_clean %>%
  group_by(hhidpn) %>%
  arrange(hhidpn, wave) %>%
  mutate(
    # ADL recovery: no ADL impairment in next two waves
    adl_recovery = ifelse(
      dementia_initial == 1 & doctor_diagnosis == 0,
      (lead(adl_impair, default = 0) == 0 & lead(adl_impair, 2, default = 0) == 0),
      FALSE
    ),
    
    # Cognitive improvement: ≥2 domains improve by >1SD
    cog_improvement = ifelse(
      dementia_initial == 1 & doctor_diagnosis == 0,
      {
        mem_diff = lead(memory_score, default = NA) - memory_score
        exe_diff = lead(exec_score, default = NA) - exec_score
        ori_diff = lead(orient_score, default = NA) - orient_score
        
        sum(
          c(!is.na(mem_diff) & mem_diff > memory_improve,
            !is.na(exe_diff) & exe_diff > exec_improve,
            !is.na(ori_diff) & ori_diff > orient_improve),
          na.rm = TRUE
        ) >= 2
      },
      FALSE
    ),
    
    # Final dementia definition
    dementia_final = case_when(
      doctor_diagnosis == 1 ~ 1,
      adl_recovery | cog_improvement ~ 0,
      dementia_initial == 1 ~ 1,
      TRUE ~ 0
    )
  ) %>%
  ungroup() %>%
  select(-matches("_mean|_sd|_threshold|_improve|_score"))  # Clean intermediate variables

# Save cleaned dementia data
save(hrs_clean, file="/formal/a1 HRS/data/processed/dementia_final2.Rdata")
load("/formal/a1 HRS/data/processed/dementia_final2.Rdata")

# part3_combine ----------------------------------------------------------

# Extract dementia definition
dementia_def <- hrs_clean %>%
  select(hhidpn, wave, dementia_final)

# Combine with behavior data
final_cohort <- hrs %>%
  left_join(dementia_def) %>% # Note: only joining by hhidpn
  left_join(behavior_final)

# Check date completeness
date_check <- final_cohort %>%
  summarise(
    iwendy_NA = mean(is.na(iwendy)),
    iwendm_NA = mean(is.na(iwendm)),
    radyear_NA = mean(is.na(radyear)),
    radmonth_NA = mean(is.na(radmonth))
  )

# Correct dates and calculate follow-up time
final_cohort_corrected <- final_cohort %>%
  # Ensure valid dates
  mutate(
    valid_date_flag = !is.na(iwendy) & !is.na(iwendm) & between(iwendm, 1, 12)
  ) %>%
  
  # Process by participant
  group_by(hhidpn) %>%
  arrange(hhidpn, wave) %>%
  
  # Flag baseline records (first wave for each participant)
  mutate(
    is_baseline = (wave == min(wave))
  ) %>%
  
  # Calculate key dates for each participant
  mutate(
    # Baseline date (only using baseline records)
    baseline_date = if_else(
      is_baseline & valid_date_flag,
      ymd(paste(iwendy, iwendm, "15", sep = "-")),
      NA_Date_
    ),
    
    # Dementia date (using diagnosis wave dates)
    dementia_date = if_else(
      dementia_final == 1 & valid_date_flag,
      ymd(paste(iwendy, iwendm, "15", sep = "-")),
      NA_Date_
    ),
    
    # Last follow-up date (last valid record date for each participant)
    last_followup_date = if_else(
      wave == max(wave[valid_date_flag]),
      ymd(paste(iwendy, iwendm, "15", sep = "-")),
      NA_Date_
    )
  ) %>%
  
  # Fill baseline dates for all records
  fill(baseline_date, last_followup_date, .direction = "downup")

# Final cohort processing
final_cohort_corrected2 <- final_cohort_corrected %>%
  group_by(hhidpn) %>%
  mutate(
    # Standardize dementia status
    dementia_final = as.integer(any(dementia_final == 1, na.rm = TRUE)),
    
    # Safely get first dementia date (avoid empty vector issue)
    first_dementia_date = if (any(!is.na(dementia_date))) {
      min(dementia_date[!is.na(dementia_date)])
    } else {
      NA_Date_
    }
  ) %>%
  ungroup() %>%
  
  # Determine final outcome date
  mutate(
    outcome_date = coalesce(first_dementia_date, last_followup_date),
    
    # Validate date sequence
    date_valid = !is.na(baseline_date) & !is.na(outcome_date) & (outcome_date >= baseline_date)
  ) %>%
  
  # Keep only baseline records
  filter(is_baseline) %>%
  
  # Calculate follow-up time
  mutate(
    followup_years = if_else(
      date_valid,
      as.numeric(outcome_date - baseline_date) / 365.25,
      NA_real_
    ),
    
    # Create outcome indicator
    status = if_else(
      !is.na(first_dementia_date) & dementia_final == 1, 
      1L,  # Dementia event
      0L   # Censored
    ),
    
    # Flag unusual follow-up times
    followup_issue = case_when(
      is.na(date_valid) ~ "Invalid dates",
      followup_years < 0 ~ "Negative time",
      followup_years > 30 ~ "Over 30 years",
      dementia_final == 1 & is.na(first_dementia_date) ~ "Dementia without date",
      TRUE ~ "Valid"
    )
  ) %>%
  
  # Ensure one record per participant
  distinct(hhidpn, .keep_all = TRUE)

# Save corrected cohort
save(final_cohort_corrected, file="/formal/a1 HRS/data/processed/final_cohort_corrected.Rdata")
load("/formal/a1 HRS/data/processed/final_cohort_corrected.Rdata")

# part4_analysis ----------------------------------------------------------
# Prepare complete analysis dataset
final_cohort_complete <- final_cohort_corrected2 %>%
  # Add baseline BMI data
  left_join(
    hrs %>%
      select(hhidpn, wave, bmi, pmbmi) %>%
      group_by(hhidpn) %>%
      filter(wave == min(wave)) %>%  # Get baseline BMI
      mutate(
        final_bmi = coalesce(pmbmi, bmi),  # Prefer measured BMI
        bmi_source = ifelse(!is.na(pmbmi), "measured", "self-reported")
      ) %>%
      select(hhidpn, final_bmi, bmi_source),
    by = "hhidpn"
  ) %>%
  # Exclude records with missing BMI
  filter(!is.na(final_bmi)) %>%
  # Prepare analysis variables
  mutate(
    age = ragey_b,
    sex = factor(ragender),
    bmi_category = cut(final_bmi,
                       breaks = c(0, 18.5, 25, 30, Inf),
                       labels = c("Underweight", "Normal", "Overweight", "Obese"))
  )

# Standardize dementia status
final_cohort_complete <- final_cohort_complete %>%
  group_by(hhidpn) %>%
  mutate(
    dementia_final = as.integer(any(dementia_final == 1)) # Standardize to maximum (diagnosis) value
  ) %>%
  ungroup()

# Extract baseline data (first record per participant)
analysis_data <- final_cohort_complete %>%
  group_by(hhidpn) %>%
  arrange(wave) %>%
  filter(row_number() == 1) %>%  # Take baseline record
  ungroup() %>%
  select(
    hhidpn,
    baseline_behavior_score = behavior_score, # Baseline behavior score
    age = ragey_b,                           # Baseline age
    sex = ragender,                          # Sex
    final_bmi,                               # Baseline BMI
    dementia_final,                          # Standardized dementia outcome
    followup_years                           # Follow-up time
  ) %>%
  filter(!(followup_years == 0 & dementia_final == 0)) %>% 
  filter(
    !is.na(baseline_behavior_score),         # Ensure behavior score not missing
    !is.na(followup_years),                  # Ensure follow-up time not missing
    followup_years > 0                       # Ensure positive follow-up time
  ) %>% na.omit()

# Prepare final analysis dataset
analysis_data2 <- analysis_data %>% 
  select(hhidpn, dementia_final, followup_years,
         baseline_behavior_score,
         sex, age, final_bmi) %>% 
  na.omit()

# Cox proportional hazards model results table function
cox_result_table <- function(cox_model, var_name = "baseline_behavior_score") {
  # Check input
  if (!inherits(cox_model, "coxph")) {
    stop("Input must be a coxph model object")
  }
  
  # Get model summary
  model_summary <- summary(cox_model)
  
  # Check if specified variable is in model
  if (!var_name %in% names(cox_model$coefficients)) {
    stop(sprintf("Variable not found in model: %s", var_name))
  }
  
  # Extract results for specified variable
  coef_idx <- which(names(cox_model$coefficients) == var_name)
  hr <- exp(cox_model$coefficients[coef_idx])
  ci_lower <- model_summary$conf.int[coef_idx, "lower .95"]
  ci_upper <- model_summary$conf.int[coef_idx, "upper .95"]
  p_value <- model_summary$coefficients[coef_idx, "Pr(>|z|)"]
  
  # Create single-row data frame
  result_table <- data.frame(
    Variable = var_name,
    HR = sprintf("%.2f", hr),
    CI_95 = sprintf("%.2f-%.2f", ci_lower, ci_upper),
    P_value = format(p_value, scientific = TRUE, digits = 2),
    Cases = model_summary$nevent,
    N = model_summary$n,
    stringsAsFactors = FALSE
  )
  
  # Prettify column names
  colnames(result_table) <- c("Variable", "HR", "95% CI", "P value", "Cases", "Total N")
  
  return(result_table)
}

# Run Cox model
cox_model <- coxph(
  Surv(followup_years, dementia_final) ~ 
    baseline_behavior_score + 
    age + 
    sex + 
    final_bmi,  # Account for within-individual correlation
  data = analysis_data2
)

# Save results
summary(cox_model)
int_res <- cox_result_table(cox_model)
rio::export(int_res, file="/formal/a1 HRS/res/int_res2.xlsx")
save(final_cohort_complete, file="/formal/a1 HRS/res/final_cohort_complete.Rdata")
save(analysis_data, file="/formal/a1 HRS/res/analysis_data2.Rdata")
save(cox_model, file="/formal/a1 HRS/res/cox_model2.Rdata")
load(file="/formal/a1 HRS/res/analysis_data2.Rdata")
load(file="/formal/a1 HRS/res/cox_model2.Rdata")

# Quartile analysis -------------------------------------------------------
# 1. Divide baseline_behavior_score into quartiles
quartiles <- quantile(analysis_data2$baseline_behavior_score, 
                      probs = c(0, 0.25, 0.5, 0.75, 1), 
                      na.rm = TRUE)

analysis_data2$behavior_score_quartile <- factor(
  cut(
    analysis_data2$baseline_behavior_score,
    breaks = quartiles,
    labels = c("Q1", "Q2", "Q3", "Q4"),
    include.lowest = TRUE
  ),
  levels = c("Q1", "Q2", "Q3", "Q4")  # Ensure Q1 is reference
)

# 2. Run Cox model with quartile as factor variable
quartile_model <- coxph(
  Surv(followup_years, dementia_final) ~ 
    behavior_score_quartile +  # As factor variable, Q1 automatically reference
    age + 
    sex + 
    final_bmi,
  data = analysis_data2
)

# 3. Quartile results extraction function
extract_quartile_results <- function(cox_model, data, quartile_var = "behavior_score_quartile") {
  # Get model summary
  model_summary <- summary(cox_model)
  
  # Get all comparison groups
  coef_names <- names(cox_model$coefficients)
  quartile_coefs <- coef_names[grepl(paste0("^", quartile_var), coef_names)]
  
  # Initialize results data frame
  results <- data.frame()
  
  # Add reference group info
  ref_level <- levels(data[[quartile_var]])[1]
  ref_cases <- sum(data$dementia_final[data[[quartile_var]] == ref_level])
  ref_total <- sum(data[[quartile_var]] == ref_level)
  
  results <- rbind(results, data.frame(
    Comparison = paste0(ref_level, " (Ref)"),
    HR = "1.00",
    CI_95 = "Reference",
    P_value = "",
    Cases = ref_cases,
    Total_N = ref_total,
    stringsAsFactors = FALSE
  ))
  
  # Extract results for each comparison group
  for (coef in quartile_coefs) {
    # Extract comparison group name
    comparison <- gsub(quartile_var, "", coef)
    
    # Extract HR and CI
    hr <- exp(cox_model$coefficients[coef])
    ci_lower <- model_summary$conf.int[coef, "lower .95"]
    ci_upper <- model_summary$conf.int[coef, "upper .95"]
    p_value <- model_summary$coefficients[coef, "Pr(>|z|)"]
    
    # Get case counts and totals for this group
    cases <- sum(data$dementia_final[data[[quartile_var]] == comparison])
    total <- sum(data[[quartile_var]] == comparison)
    
    results <- rbind(results, data.frame(
      Comparison = comparison,
      HR = sprintf("%.2f", hr),
      CI_95 = sprintf("%.2f-%.2f", ci_lower, ci_upper),
      P_value = format(p_value, scientific = TRUE, digits = 2),
      Cases = cases,
      Total_N = total,
      stringsAsFactors = FALSE
    ))
  }
  
  colnames(results) <- c("Comparison", "HR", "95% CI", "P value", "Cases", "Total N")
  return(results)
}

# 4. Extract quartile comparison results
quartile_results <- extract_quartile_results(quartile_model, analysis_data2)

# 5. Add score range information
score_ranges <- c(
  paste(round(quartiles[1], 2), "-", round(quartiles[2], 2)),
  paste(round(quartiles[2], 2), "-", round(quartiles[3], 2)),
  paste(round(quartiles[3], 2), "-", round(quartiles[4], 2)),
  paste(round(quartiles[4], 2), "-", round(quartiles[5], 2))
)

quartile_results$Score_Range <- c(score_ranges[1], score_ranges[2:4])
quartile_results$Variable <- "baseline_behavior_score"

# 6. Combine continuous and quartile results
int_res$Analysis_Type <- "Continuous"
int_res$Comparison <- "Continuous"
int_res$Score_Range <- "Full range"

quartile_results$Analysis_Type <- "Quartile comparison"

# Merge results
combined_results <- rbind(
  int_res[, c("Analysis_Type", "Comparison", "Score_Range", "Variable", "HR", "95% CI", "P value", "Cases", "Total N")],
  quartile_results[, c("Analysis_Type", "Comparison", "Score_Range", "Variable", "HR", "95% CI", "P value", "Cases", "Total N")]
)

# 7. Format final table
combined_results$HR_CI <- ifelse(combined_results$`95% CI` == "Reference",
                                 "Reference",
                                 paste0(combined_results$HR, " (", combined_results$`95% CI`, ")"))

final_table <- combined_results[, c("Analysis_Type", "Comparison", "Score_Range", "HR_CI", "P value", "Cases", "Total N")]
colnames(final_table) <- c("Analysis Type", "Comparison", "Score Range", "HR (95% CI)", "P Value", "Cases", "Total N")

# Export results
rio::export(final_table, file="/formal/a1 HRS/res/baseline_behavior_score_quartile_comparison.xlsx")

# Restricted cubic splines analysis ---------------------------------------
analysis_data2 <- analysis_data %>% 
  select(hhidpn, dementia_final, followup_years,
         baseline_behavior_score,
         sex, age, final_bmi) %>% 
  na.omit()

library(plotRCS)
library(ggplot2)

# Rename for plotting
colnames(analysis_data2)[4] <- "behavior_score"

# Run RCS analysis
rcs_res <- rcsplot(data = analysis_data2,
                   outcome = "dementia_final",
                   time = "followup_years",
                   exposure = "behavior_score",
                   covariates = c("sex", "final_bmi", "age"),
                   linesize = 1,
                   linecolor = "blue")

# Save plots
ggsave(plot = rcs_res, file = "/formal/a1 HRS/res/rcs_plot3.pdf", width = 2.5, height = 2.5)
ggsave(plot = rcs_res, file = "/formal/a1 HRS/res/rcs_plot3.png", width = 2.5, height = 2.5)