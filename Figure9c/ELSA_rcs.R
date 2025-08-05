# Load required libraries
library(tidyverse)   # For data manipulation and visualization
library(lubridate)   # For date handling
library(survival)    # For survival analysis
library(broom)       # For tidying model outputs
library(Hmisc)       # For statistical functions
library(haven)       # For importing SAS/Stata files
library(dplyr)       # For data manipulation
library(tidyr)       # For data tidying


# Load raw ELSA dataset
load("/code/Figure9/Figure 9b/elsa.Rdata")

# Define function to fill forward missing values
fill_forward <- function(data, id_var, time_var, fill_vars) {
  data %>%
    group_by({{id_var}}) %>%
    arrange({{time_var}}) %>%
    mutate(across(all_of(fill_vars), 
                  ~ifelse(is.na(.) & row_number() == 1,  # If baseline and missing
                          .[which(!is.na(.) & row_number() > 1)[1]],  # Fill with next non-missing
                          .))) %>%
    ungroup()
}

# Initialize clean dataset
elsa_clean0 = elsa

# PART 1: BEHAVIOR SCORE CONSTRUCTION -------------------------------------

# Step 1: Build glycemic score (prioritizing HbA1c, then diabetes diagnosis)
score_clean <- elsa_clean0 %>%
  # Fill forward relevant variables
  fill_forward(idauniqc, wave, c("hba1c", "diabe", "rxdiabi", "rxdiabo")) %>%
  mutate(
    glycemic_score = case_when(
      # Criteria for diabetes
      (!is.na(hba1c) & hba1c > 6.5) | 
        (diabe == 1) | 
        (rxdiabi == 1) ~ 
        case_when(
          rxdiabo == 1 ~ 2,   # On medication: 2 points
          TRUE ~ 0            # No medication: 0 points
        ),
      
      # No diabetes criteria
      (!is.na(hba1c) & hba1c <= 6.5) | 
        (diabe == 0) | 
        (is.na(diabe) & !(rxdiabi == 1 | rxdiabo == 1)) ~ 4,  # No diabetes: 4 points
      
      TRUE ~ NA_real_
    ),
    glycemic_source = case_when(
      !is.na(hba1c) ~ "HbA1c",
      !is.na(diabe) | !is.na(rxdiabi) | !is.na(rxdiabo) ~ "Diabetes_record",
      TRUE ~ "Missing"
    )
  )

# Step 2: Build physical activity score
score_clean <- score_clean %>%
  # Fill forward activity variables
  fill_forward(idauniqc, wave, c("vgactx_e", "mdactx_e", "ltactx_e")) %>%
  group_by(wave) %>%  # Calculate activity levels by wave
  mutate(
    # Convert activity frequency to scores
    vg_score = case_when(
      vgactx_e == 5 ~ 0,
      vgactx_e == 4 ~ 1,
      vgactx_e == 3 ~ 2,
      vgactx_e == 2 ~ 3,
      TRUE ~ NA_real_
    ),
    md_score = case_when(
      mdactx_e == 5 ~ 0,
      mdactx_e == 4 ~ 1,
      mdactx_e == 3 ~ 2,
      mdactx_e == 2 ~ 3,
      TRUE ~ NA_real_
    ),
    lt_score = case_when(
      ltactx_e == 5 ~ 0,
      ltactx_e == 4 ~ 1,
      ltactx_e == 3 ~ 2,
      ltactx_e == 2 ~ 3,
      TRUE ~ NA_real_
    ),
    
    # Calculate weighted total score
    activity_score = 
      (ifelse(is.na(vg_score), 0, vg_score) * 4) + 
      (ifelse(is.na(md_score), 0, md_score) * 2) + 
      (ifelse(is.na(lt_score), 0, lt_score) * 1),
  ) %>%
  ungroup()

# Step 3: Build mood and loneliness scores
score_clean <- score_clean %>%
  # Fill forward mood variables
  fill_forward(idauniqc, wave, c("cesd", "lnlys3")) %>%
  # Calculate CESD quartiles
  group_by(wave) %>%
  mutate(
    cesd_q25 = quantile(cesd[cesd > 0], 0.25, na.rm = TRUE),
    cesd_q50 = quantile(cesd[cesd > 0], 0.50, na.rm = TRUE),
    cesd_q75 = quantile(cesd[cesd > 0], 0.75, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    # CES-D score
    mood_score = case_when(
      is.na(cesd) ~ NA_real_,
      cesd == 0 ~ 0,
      cesd > 0 & cesd <= cesd_q25 ~ -1,
      cesd > 0 & cesd <= cesd_q50 ~ -2,
      cesd > 0 & cesd <= cesd_q75 ~ -3,
      cesd > 0 & cesd > cesd_q75 ~ -4
    ),
    
    # Loneliness score
    social_score = ifelse(!is.na(lnlys3), 3 - lnlys3, NA_real_),
    
    # Remove temporary columns
    across(starts_with("cesd_q"), ~NULL),
    
    # Mark data source
    social_source = "UCLA_3item_continuous"
  )

# Step 4: Build complete behavior score
score_clean <- score_clean %>%
  mutate(
    # Count valid components
    valid_components = 4 - is.na(glycemic_score) - is.na(mood_score) - 
      is.na(social_score) - is.na(activity_score),
    
    # Total behavior score (requires ≥3 valid components)
    behavior_score = ifelse(
      valid_components >= 3,
      rowSums(select(., glycemic_score, mood_score, social_score, activity_score), 
              na.rm = FALSE),
      NA_real_
    ),
    
    # Behavior score categories
    behavior_group = case_when(
      behavior_score == 4 ~ "Optimal",
      behavior_score >= 2 ~ "Intermediate",
      behavior_score < 2 ~ "Poor",
      TRUE ~ "Missing"
    )
  )

# Save behavior score components
score_clean_sub = score_clean %>% 
  select(idauniqc, wave, glycemic_score, activity_score, mood_score, social_score, 
         behavior_score, behavior_group)

# Save processed data
save(score_clean, file = "/data/processed/score_clean.Rdata")


# PART 2: DEMENTIA ASSESSMENT --------------------------------------------

# Step 1: Calculate ADL impairment
elsa_clean <- elsa_clean0 %>%
  mutate(
    across(c(dressa, batha, eata, beda, toilta, walkra), ~ ifelse(. == 1, 1, 0)),
    adl_total = rowSums(select(., dressa, batha, eata, beda, toilta, walkra), na.rm = TRUE),
    adl_impair = ifelse(adl_total >= 1, 1, 0)
  )

# Step 2: Build education-specific cognitive norms (based on raeducl groups)
cognitive_norms <- elsa_clean %>%
  filter(agey >= 50) %>%  # Only use age ≥50 for norm building
  group_by(raeducl) %>%   # Group by original education codes
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
    # Cognitive impairment thresholds (mean - 1.5SD)
    memory_threshold = memory_mean - 1.5 * memory_sd,
    exec_threshold = exec_mean - 1.5 * exec_sd,
    orient_threshold = orient_mean - 1.5 * orient_sd,
    
    # Cognitive improvement thresholds (1SD change)
    memory_improve = memory_sd,
    exec_improve = exec_sd,
    orient_improve = orient_sd
  )

# Step 3: Apply cognitive impairment definitions
elsa_clean <- elsa_clean %>%
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

# Step 4: Process physician diagnoses
elsa_clean <- elsa_clean %>%
  mutate(
    doctor_diagnosis = ifelse(demene == 1 | alzhe == 1, 1, 0)
  )

# Step 5: Define initial dementia cases
elsa_clean <- elsa_clean %>%
  mutate(
    dementia_initial = ifelse(
      (adl_impair == 1 & cog_impair == 1) | doctor_diagnosis == 1,
      1, 0)
  )

# Step 6: Dynamic exclusion rules (education-stratified)
elsa_clean <- elsa_clean %>%
  group_by(idauniqc) %>%
  arrange(idauniqc, wave) %>%
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
  select(-matches("_mean|_sd|_threshold|_score"))  # Clean intermediate variables

# Step 7: Process time variables
elsa_clean2 <- elsa_clean %>%
  group_by(idauniqc) %>%
  mutate(
    # Baseline date (using 15th of month)
    baseline_date = if_else(
      !is.na(iwindy) & !is.na(iwindm) & between(iwindm, 1, 12),
      ymd(paste(iwindy, iwindm, "15", sep = "-")),
      NA_Date_
    ),
    
    # Death date (three-level logic)
    death_date = case_when(
      !is.na(radyear) & iwindy == radyear & iwindm > 6 ~ 
        ymd(paste(radyear, iwindm, "15", sep = "-")),
      !is.na(radyear) ~ ymd(paste(radyear, 6, "15", sep = "-")),
      TRUE ~ NA_Date_
    ),
    
    # Last follow-up date
    last_followup_date = {
      valid_waves <- which(!is.na(iwindy) & !is.na(iwindm) & between(iwindm, 1, 12))
      if(length(valid_waves) > 0) {
        last_wave <- max(valid_waves)
        ymd(paste(iwindy[last_wave], iwindm[last_wave], "15", sep = "-"))
      } else {
        NA_Date_
      }
    },
    
    # Date validity checks
    death_date = if_else(
      !is.na(death_date) & !is.na(baseline_date) & death_date < baseline_date,
      baseline_date + days(1), death_date
    ),
    death_date = if_else(
      !is.na(death_date) & !is.na(last_followup_date) & death_date > last_followup_date,
      last_followup_date, death_date
    )
  ) %>%
  ungroup()

# Step 8: Determine outcomes and follow-up time
elsa_clean3 <- elsa_clean2 %>%
  group_by(idauniqc) %>%
  # Calculate dementia date
  mutate(dementia_date = if_else(dementia_initial == 1, baseline_date, NA_Date_)) %>%
  # Calculate outcome date
  mutate(outcome_date = coalesce(dementia_date, death_date, last_followup_date)) %>%
  # Calculate follow-up years
  mutate(followup_years = ifelse(
    !is.na(baseline_date) & !is.na(outcome_date),
    as.numeric(outcome_date - baseline_date) / 365.25,
    NA_real_
  )) %>%
  # Final dementia classification
  mutate(dementia_final = case_when(
    doctor_diagnosis == 1 ~ 1,
    adl_recovery | cog_improvement ~ 0,
    dementia_initial == 1 ~ 1,
    TRUE ~ 0
  )) %>%
  # Flag valid records
  mutate(valid_record = !is.na(baseline_date) & !is.na(outcome_date) & (outcome_date >= baseline_date)) %>%
  ungroup()

# Save processed data
save(elsa_clean3, file = "/data/processed/elsa_clean3.Rdata")


# PART 3: DATA MERGING AND ANALYSIS ---------------------------------------

# Merge behavior scores with dementia outcomes
elsa_clean4 <- elsa_clean3 %>% 
  left_join(score_clean_sub)

# Prepare analysis dataset
elsa_clean4 <- elsa_clean4 %>% 
  arrange(idauniqc, wave) %>%
  group_by(idauniqc) %>%
  # Forward-fill behavior-related scores
  mutate(
    glycemic_score_filled = ifelse(
      is.na(glycemic_score) & wave == min(wave),
      as.numeric(first(na.omit(glycemic_score))),
      as.numeric(glycemic_score)
    ),
    activity_score_filled = ifelse(
      is.na(activity_score) & wave == min(wave),
      as.numeric(first(na.omit(activity_score))),
      as.numeric(activity_score)
    ),
    mood_score_filled = ifelse(
      is.na(mood_score) & wave == min(wave),
      as.numeric(first(na.omit(mood_score))),
      as.numeric(mood_score)
    ),
    social_score_filled = ifelse(
      is.na(social_score) & wave == min(wave),
      as.numeric(first(na.omit(social_score))),
      as.numeric(social_score)
    )
  ) %>%
  # Recalculate behavior score and groups
  mutate(
    valid_components = 4 - is.na(glycemic_score_filled) - is.na(mood_score_filled) - 
      is.na(social_score_filled) - is.na(activity_score_filled),
    behavior_score = ifelse(
      valid_components >= 3,
      rowSums(
        cbind(glycemic_score_filled, mood_score_filled, social_score_filled, activity_score_filled),
        na.rm = FALSE
      ),
      NA_real_
    ),
    behavior_group = case_when(
      behavior_score == 4 ~ "Optimal",
      behavior_score >= 2 ~ "Intermediate",
      behavior_score < 2 ~ "Poor",
      TRUE ~ "Missing"
    )
  ) %>%
  ungroup()

# Create final analysis dataset
analysis_data <- elsa_clean4 %>%
  group_by(idauniqc) %>%
  # Forward-fill baseline BMI
  mutate(
    bmi_filled = ifelse(
      is.na(mbmi) & wave == min(wave),
      first(na.omit(mbmi)),
      mbmi
    )
  ) %>%
  # Flag dementia events across all waves
  mutate(
    any_dementia = max(dementia_final, na.rm = TRUE),  # Any dementia event
    any_doc = max(as.numeric(doctor_diagnosis), na.rm = TRUE)  # Any physician diagnosis
  ) %>%
  # Keep baseline records only
  filter(wave == min(wave)) %>%
  # Update outcome variable (using all follow-up dementia status)
  mutate(
    dementia_event = ifelse(any_dementia == 1, 1, 0)
  ) %>%
  # Select final variables
  select(
    idauniqc, wave,
    age = agey, sex = ragender, bmi = bmi_filled,
    behavior_score, behavior_group,
    glycemic_score_filled, activity_score_filled, mood_score_filled, social_score_filled,
    dementia_event,  # All follow-up events
    followup_years   # Time from baseline to outcome
  ) %>%
  # Exclude records with zero follow-up and no dementia
  filter(!(followup_years == 0 & dementia_event == 0)) %>% 
  # Exclude records with missing or zero follow-up time
  filter(followup_years > 0,
         !is.na(followup_years)) %>% 
  na.omit()

# Function to create Cox results table
cox_result_table <- function(cox_model, var_name = "behavior_score") {
  # Input validation
  if (!inherits(cox_model, "coxph")) {
    stop("Input must be a coxph model object")
  }
  
  # Get model summary
  model_summary <- summary(cox_model)
  
  # Check if variable exists in model
  if (!var_name %in% names(cox_model$coefficients)) {
    stop(sprintf("Variable not in model: %s", var_name))
  }
  
  # Extract variable results
  coef_idx <- which(names(cox_model$coefficients) == var_name)
  hr <- exp(cox_model$coefficients[coef_idx])
  ci_lower <- model_summary$conf.int[coef_idx, "lower .95"]
  ci_upper <- model_summary$conf.int[coef_idx, "upper .95"]
  p_value <- model_summary$coefficients[coef_idx, "Pr(>|z|)"]
  
  # Create single-row dataframe
  result_table <- data.frame(
    Variable = var_name,
    HR = sprintf("%.2f", hr),
    CI_95 = sprintf("%.2f-%.2f", ci_lower, ci_upper),
    P_value = format(p_value, scientific = TRUE, digits = 2),
    Cases = model_summary$nevent,
    N = model_summary$n,
    stringsAsFactors = FALSE
  )
  
  # Format column names
  colnames(result_table) <- c("Variable", "HR", "95% CI", "P value", "Cases", "Total N")
  
  return(result_table)
}

# Prepare analysis dataset with complete cases
analysis_data2 = analysis_data %>% 
  select(idauniqc, dementia_event, followup_years,
         behavior_score,
         glycemic_score_filled, activity_score_filled, mood_score_filled, social_score_filled,
         sex, age, bmi) %>% 
  na.omit()

# Run Cox model with continuous behavior score
cox_model <- coxph(
  Surv(followup_years, dementia_event) ~ 
    behavior_score + age + sex + bmi,
  data = analysis_data2
)

# Generate results table
int_res = cox_result_table(cox_model)
rio::export(int_res, file = "int_res3.xlsx")

# QUARTILE ANALYSIS -------------------------------------------------------

# Function to extract quartile results
extract_quartile_results <- function(cox_model, data, quartile_var = "behavior_score_quartile") {
  # Get model summary
  model_summary <- summary(cox_model)
  
  # Get all comparison groups
  coef_names <- names(cox_model$coefficients)
  quartile_coefs <- coef_names[grepl(paste0("^", quartile_var), coef_names)]
  
  # Initialize results dataframe
  results <- data.frame()
  
  # Add reference group info
  ref_level <- levels(data[[quartile_var]])[1]
  ref_cases <- sum(data$dementia_event[data[[quartile_var]] == ref_level])
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
    # Get comparison group name
    comparison <- gsub(quartile_var, "", coef)
    
    # Extract HR and CI
    hr <- exp(cox_model$coefficients[coef])
    ci_lower <- model_summary$conf.int[coef, "lower .95"]
    ci_upper <- model_summary$conf.int[coef, "upper .95"]
    p_value <- model_summary$coefficients[coef, "Pr(>|z|)"]
    
    # Get case counts and totals
    cases <- sum(data$dementia_event[data[[quartile_var]] == comparison])
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

# Create quartile groups
quartiles <- quantile(analysis_data2$behavior_score, 
                      probs = c(0, 0.25, 0.5, 0.75, 1), 
                      na.rm = TRUE)

analysis_data2$behavior_score_quartile <- factor(
  cut(
    analysis_data2$behavior_score,
    breaks = quartiles,
    labels = c("Q1", "Q2", "Q3", "Q4"),
    include.lowest = TRUE
  ),
  levels = c("Q1", "Q2", "Q3", "Q4")  # Ensure Q1 is reference
)

# Run Cox model with quartile groups
quartile_model <- coxph(
  Surv(followup_years, dementia_event) ~ 
    behavior_score_quartile +  # As factor variable
    age + 
    sex + 
    bmi,
  data = analysis_data2
)

# Extract quartile results
quartile_results <- extract_quartile_results(quartile_model, analysis_data2)

# Add score range information
score_ranges <- c(
  paste(round(quartiles[1], 2), "-", round(quartiles[2], 2)),
  paste(round(quartiles[2], 2), "-", round(quartiles[3], 2)),
  paste(round(quartiles[3], 2), "-", round(quartiles[4], 2)),
  paste(round(quartiles[4], 2), "-", round(quartiles[5], 2))
)

quartile_results$Score_Range <- c(score_ranges[1], score_ranges[2:4])

# Reorder columns
final_results <- quartile_results[, c("Comparison", "Score_Range", "HR", "95% CI", "P value", "Cases", "Total N")]

# Combine continuous and quartile results
int_res$Analysis_Type <- "Continuous"
int_res$Quartile <- "Overall"
int_res$Score_Range <- "Full range"
int_res$Comparison <- "Continuous"

# Rename columns to match final_results
colnames(int_res) <- c("Variable", "HR", "95% CI", "P value", "Cases", "Total N", 
                       "Analysis_Type", "Quartile", "Score_Range", "Comparison")

final_results$Analysis_Type <- "Quartile comparison"
final_results$Variable <- "behavior_score_quartile"

# Merge results
combined_results <- rbind(
  int_res[, c("Comparison", "Score_Range", "Variable", "HR", "95% CI", "P value", "Cases", "Total N", "Analysis_Type")],
  final_results[, c("Comparison", "Score_Range", "Variable", "HR", "95% CI", "P value", "Cases", "Total N", "Analysis_Type")]
)

# Format final table
combined_results$HR_CI <- ifelse(combined_results$`95% CI` == "Reference",
                                 "Reference",
                                 paste0(combined_results$HR, " (", combined_results$`95% CI`, ")"))

final_table <- combined_results[, c("Analysis_Type", "Comparison", "Score_Range", "HR_CI", "P value", "Cases", "Total N")]
colnames(final_table) <- c("Analysis Type", "Comparison", "Score Range", "HR (95% CI)", "P Value", "Cases", "Total N")

# Export results
rio::export(final_table, file = "/res/behavior_score_quartile_comparison.xlsx")

# RESTRICTED CUBIC SPLINES ANALYSIS ---------------------------------------

rcs_res = rcsplot(data = analysis_data2,
                  outcome = "dementia_event",
                  time = "followup_years",
                  exposure = "behavior_score",
                  covariates = c("sex", "bmi", "age"),
                  linesize = 1,
                  linecolor = "blue")

# Save RCS plots
ggsave(plot = rcs_res, file = "/res/rcs_plot3.pdf", width = 2.5, height = 2.5)
ggsave(plot = rcs_res, file = "/res/rcs_plot3.png", width = 2.5, height = 2.5)