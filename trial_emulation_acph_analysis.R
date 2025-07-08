### Trial emulation: EMPA-REG outcome ###
### Adjusted CPH model ###

gc()
rm(list = ls())

#Import packages
library(dplyr)
library(survival)
library(survey)
library(lubridate)
library(ggplot2)
library(patchwork)
library(mice)
library(mitools)
library(knitr)
library(kableExtra)
library(parallel)

#Import data
data <- read.csv("S:/23_010_CVD_trial_emu/march_analysis/rw_230325.csv")

#Re-format dates
data$last_date <- ymd(data$last_contact_correct)
data$first_date <- ymd(data$cohort_entry_date)

#Define survival object 
data$surv_obj <- with(data, Surv(time = as.numeric(difftime(last_date, first_date, units="days"))/365.25, event = died))
summary(data$surv_obj)

#Add additional variables 
data$egfr_sq <- data$egfr_at_initiation^2

#Set dataformat 
data$empa <- as.factor(data$empa)
data$eth_recoded <- as.factor(data$eth_recoded)
data$cohort_year_recoded <- as.factor(data$cohort_year_recoded)
data$smoking_status <- as.factor(data$smoking_status)
data$imd_quintile <- as.factor(data$imd_quintile)
data$male <- as.factor(data$male)
data$su_before_initiation <- as.factor(data$su_before_initiation)
data$metformin_monotherapy_before_ini <- as.factor(data$metformin_monotherapy_before_ini)
data$lipid_lowering_before_initiation <- as.factor(data$lipid_lowering_before_initiation)
data$anti_htn_before_initiation <- as.factor(data$anti_htn_before_initiation)
data$glp1ra_before_initiation <- as.factor(data$glp1ra_before_initiation)
data$anti_platelet_before_initiation <- as.factor(data$anti_platelet_before_initiation)
data$anticoag_before_initiation <- as.factor(data$anticoag_before_initiation)
data$insulin_before_initiation <- as.factor(data$insulin_before_initiation)
data$hf <- as.factor(data$hf)
data$smi <- as.factor(data$smi)
data$copd <- as.factor(data$copd)
data$ra <- as.factor(data$ra)
data$asthma <- as.factor(data$asthma)
data$epilepsy <- as.factor(data$epilepsy)
data$dementia <- as.factor(data$dementia)
data$ibd <- as.factor(data$ibd)
data$liver <- as.factor(data$liver)
data$imd_quintile <- as.factor(data$imd_quintile)
data$cvd_diagnosis_before_cohort <- as.factor(data$cvd_diagnosis_before_cohort)
data$cancer_diagnosis_5y_before_cohor <- as.factor(data$cancer_diagnosis_5y_before_cohor)

#Unadjusted CPH model CCA  -----------------------------------------------
unadj_cph_cca <- coxph(surv_obj ~ empa, data = data)
summary(unadj_cph_cca)
#exp(coef) exp(-coef) lower .95 upper .95


#Adjusted CPH model CCA  -------------------------------------------------
adj_cph_cca <- coxph(surv_obj ~ empa + cohort_year_recoded + male + 
                       age_at_initiation + eth_recoded + smoking_status + 
                       sbp_at_initiation + bmi_at_initiation + hba1c_at_initiation + 
                       ldl_at_initiation + ldl_sq + hdl_log + egfr_at_initiation + egfr_sq + 
                       su_before_initiation + metformin_monotherapy_before_ini + 
                       lipid_lowering_before_initiation + anti_htn_before_initiation + glp1ra_before_initiation + 
                       anti_platelet_before_initiation + anticoag_before_initiation + insulin_before_initiation + 
                       hf + smi + copd + ra + asthma + cancer_diagnosis_5y_before_cohor + cvd_diagnosis_before_cohort 
                     + epilepsy + dementia + ibd + liver + imd_quintile, data = data)
summary(adj_cph_cca)


# PH diagnostics --------------------------------------------------------------------

# Schoenfeld residual test 
ph_test <- cox.zph(adj_cph_cca) #PH holds for all except SBP, cancer diagnosis, dementia
ph_test # no global evidence of evidence against PH assumption
plot(ph_test, var='empa', resid=FALSE)


#MICE --------------------------------------------------------------------

#Nelson-Aalen estimator of baseline cumulative hazard at each person's time of event/censoring 
cox_na <- coxph(surv_obj ~ 1, data=data)
na_estimates <- basehaz(cox_na, centered=FALSE)
data$na_est <- sapply(data$surv_obj[,1], function(t){
  na_estimates$hazard[which.min(abs(na_estimates$time - t))]
})

#select relevant data
data_selected <- data %>% 
  select(
    empa, cohort_year_recoded, male, age_at_initiation, bmi_at_initiation, eth_recoded, smoking_status, 
    sbp_at_initiation, hba1c_at_initiation, ldl_at_initiation, ldl_sq, hdl_log, egfr_at_initiation, egfr_sq, 
    su_before_initiation, metformin_monotherapy_before_ini, lipid_lowering_before_initiation, 
    anti_htn_before_initiation, glp1ra_before_initiation, anti_platelet_before_initiation, 
    anticoag_before_initiation, insulin_before_initiation, hf, smi, copd, ra, asthma, epilepsy, 
    dementia, ibd, liver, imd_quintile, cvd_diagnosis_before_cohort, cancer_diagnosis_5y_before_cohor, died, na_est, last_contact_correct, person_id, cohort_entry_date
  )

#set up predictor matrix 
vars_to_impute <- c(
  "bmi_at_initiation",
  "sbp_at_initiation", 
  "hba1c_at_initiation", 
  "ldl_at_initiation", 
  "hdl_log", 
  "egfr_at_initiation", 
  "eth_recoded", 
  "smoking_status"
)

meth <- make.method(data_selected)

meth[c("bmi_at_initiation",
       "sbp_at_initiation", 
       "hba1c_at_initiation", 
       "ldl_at_initiation", 
       "hdl_log", 
       "egfr_at_initiation")] <- "pmm"

meth[c("eth_recoded", "smoking_status")] <- "polyreg"

meth["ldl_sq"] <- "~I(ldl_at_initiation^2)"
meth["egfr_sq"] <- "~I(egfr_at_initiation^2)"

meth[!(names(meth) %in% c(vars_to_impute, "ldl_sq", "egfr_sq"))] <- ""

pm <- make.predictorMatrix(data_selected)

pm[c("ldl_sq", "egfr_sq"), ] <- 0 
pm[, c("ldl_sq", "egfr_sq")] <- 1
pm[c('last_contact_correct', 'person_id', 'cohort_entry_date'), ] <- 0 
pm[, c('last_contact_correct', 'person_id', 'cohort_entry_date')] <- 0 

View(pm[vars_to_impute, ])
View(meth)
set.seed(3001)

#imputation process 
imps <- mice(data_selected, 
             m=5, 
             maxit=5, 
             predictorMatrix=pm, 
             method = meth,
             printFlag=TRUE)

#Save files
save(imps, file="S:/23_010_CVD_trial_emu/march_analysis/tte.RData")
#load("S:/23_010_CVD_trial_emu/march_analysis/tte.RData")

#Adjusted CPH model pooled mi  -------------------------------------------
adjusted_cph_mi <- with(imps, coxph(Surv(time=as.numeric(difftime(last_contact_correct, cohort_entry_date, units="days"))/365.25, event = died) ~
                                      empa + cohort_year_recoded + male + 
                                      age_at_initiation + eth_recoded + smoking_status + 
                                      sbp_at_initiation + bmi_at_initiation + hba1c_at_initiation + 
                                      ldl_at_initiation + ldl_sq + hdl_log + egfr_at_initiation + egfr_sq + 
                                      su_before_initiation + metformin_monotherapy_before_ini + 
                                      lipid_lowering_before_initiation + anti_htn_before_initiation + glp1ra_before_initiation + 
                                      anti_platelet_before_initiation + anticoag_before_initiation + insulin_before_initiation + 
                                      hf + smi + copd + ra + asthma + cancer_diagnosis_5y_before_cohor +
                                      epilepsy + dementia + ibd + liver + imd_quintile + cvd_diagnosis_before_cohort))

summary_cph_mi <- summary(pool(adjusted_cph_mi))
summary_cph_mi$exp_coef <- exp(summary_cph_mi$estimate)
print(summary_cph_mi)

summary_cph_mi <- summary_cph_mi %>%
  mutate(
    variable = term,
    HR = exp(estimate),
    lower_ci = exp(estimate - 1.96*std.error), 
    upper_ci = exp(estimate + 1.96*std.error)
  ) 

summary_cph_mi <- summary_cph_mi %>%
  filter(term == "empa1")

print(summary_cph_mi)

#save output for e-value analysis 
write.csv(summary_cph_mi, "S:/23_010_CVD_trial_emu/march_analysis/tte_estimates_mice.csv", row.names=TRUE)

#Interaction analysis  -------------------------------------------

#Import data
data <- read.csv("S:/23_010_CVD_trial_emu/march_analysis/rw_230325.csv")

#Re-format dates
data$last_date <- ymd(data$last_contact_correct)
data$first_date <- ymd(data$cohort_entry_date)

#Add additional variables 
data$egfr_sq <- data$egfr_at_initiation^2
data$rct_eligible <- factor(data$rct_ineligible, levels = c(0,1 ), labels = c(1, 0))

data$empa <- as.numeric(as.character(data$empa))
data$rct_eligible <- as.numeric(as.character(data$rct_eligible))

#create interaction for factor variables - empa/rct_ineligible 
data$empa_rct_interaction <- with(data, empa*rct_eligible)

tabulate(data$empa_rct_interaction)
tabulate(data$empa)

tabulate(data$empa_rct_interaction)/tabulate(data$empa)

#Define survival object 
data$surv_obj <- with(data, Surv(time = as.numeric(difftime(last_date, first_date, units="days"))/365.25, event = died))
summary(data$surv_obj)

#Set dataformat 
data$eth_recoded <- as.factor(data$eth_recoded)
data$cohort_year_recoded <- as.factor(data$cohort_year_recoded)
data$smoking_status <- as.factor(data$smoking_status)
data$imd_quintile <- as.factor(data$imd_quintile)
data$male <- as.factor(data$male)
data$su_before_initiation <- as.factor(data$su_before_initiation)
data$metformin_monotherapy_before_ini <- as.factor(data$metformin_monotherapy_before_ini)
data$lipid_lowering_before_initiation <- as.factor(data$lipid_lowering_before_initiation)
data$anti_htn_before_initiation <- as.factor(data$anti_htn_before_initiation)
data$glp1ra_before_initiation <- as.factor(data$glp1ra_before_initiation)
data$anti_platelet_before_initiation <- as.factor(data$anti_platelet_before_initiation)
data$anticoag_before_initiation <- as.factor(data$anticoag_before_initiation)
data$insulin_before_initiation <- as.factor(data$insulin_before_initiation)
data$hf <- as.factor(data$hf)
data$smi <- as.factor(data$smi)
data$copd <- as.factor(data$copd)
data$ra <- as.factor(data$ra)
data$asthma <- as.factor(data$asthma)
data$epilepsy <- as.factor(data$epilepsy)
data$dementia <- as.factor(data$dementia)
data$ibd <- as.factor(data$ibd)
data$liver <- as.factor(data$liver)
data$imd_quintile <- as.factor(data$imd_quintile)
data$cvd_diagnosis_before_cohort <- as.factor(data$cvd_diagnosis_before_cohort)
data$cancer_diagnosis_5y_before_cohor <- as.factor(data$cancer_diagnosis_5y_before_cohor)

#Nelson-Aalen estimator of baseline cumulative hazard at each person's time of event/censoring 
cox_na <- coxph(surv_obj ~ 1, data=data)
na_estimates <- basehaz(cox_na, centered=FALSE)
data$na_est <- sapply(data$surv_obj[,1], function(t){
  na_estimates$hazard[which.min(abs(na_estimates$time - t))]
})


data_selected <- data %>% 
  select(
    empa, rct_eligible, empa_rct_interaction, cohort_year_recoded, male, age_at_initiation, bmi_at_initiation, eth_recoded, smoking_status, 
    sbp_at_initiation, hba1c_at_initiation, ldl_at_initiation, ldl_sq, hdl_log, egfr_at_initiation, egfr_sq, 
    su_before_initiation, metformin_monotherapy_before_ini, lipid_lowering_before_initiation, 
    anti_htn_before_initiation, glp1ra_before_initiation, anti_platelet_before_initiation, 
    anticoag_before_initiation, insulin_before_initiation, hf, smi, copd, ra, asthma, epilepsy, 
    dementia, ibd, liver, imd_quintile, cvd_diagnosis_before_cohort, cancer_diagnosis_5y_before_cohor, died, na_est, last_contact_correct, person_id, cohort_entry_date)


#arrange the methods - this is passive imputation of derived variables LDL and eGFR square
vars_to_impute <- c(
  "bmi_at_initiation",
  "sbp_at_initiation", 
  "hba1c_at_initiation", 
  "ldl_at_initiation", 
  "hdl_log", 
  "egfr_at_initiation", 
  "eth_recoded", 
  "smoking_status"
)

meth <- make.method(data_selected)

meth[c("bmi_at_initiation",
       "sbp_at_initiation", 
       "hba1c_at_initiation", 
       "ldl_at_initiation", 
       "hdl_log", 
       "egfr_at_initiation")] <- "pmm"

meth[c("eth_recoded", "smoking_status")] <- "polyreg"

meth["ldl_sq"] <- "~I(ldl_at_initiation^2)"
meth["egfr_sq"] <- "~I(egfr_at_initiation^2)"
meth["empa_rct_interaction"] <- "~I(empa*rct_eligible)"

meth[!(names(meth) %in% c(vars_to_impute, "ldl_sq", "egfr_sq", "empa_rct_interaction"))] <- ""

pm <- make.predictorMatrix(data_selected)

pm[c("ldl_sq", "egfr_sq", "empa_rct_interaction"), ] <- 0 
pm[, c("ldl_sq", "egfr_sq", "empa_rct_interaction")] <- 1
pm[c('last_contact_correct', 'person_id', 'cohort_entry_date'), ] <- 0 
pm[, c('last_contact_correct', 'person_id', 'cohort_entry_date')] <- 0 

View(pm[vars_to_impute, ])
View(meth)
set.seed(3001)

#imputation process 
imps <- mice(data_selected, 
             m=5, 
             maxit=5, 
             predictorMatrix=pm, 
             method = meth,
             printFlag=TRUE)

#Save files
save(imps, file="S:/23_010_CVD_trial_emu/march_analysis/tte_int.RData")
#load("S:/23_010_CVD_trial_emu/march_analysis/tte_int.RData")

#Adjusted CPH model pooled mi  -------------------------------------------
adjusted_cph_mi <- with(imps, coxph(Surv(time=as.numeric(difftime(last_contact_correct, cohort_entry_date, units="days"))/365.25, event = died) ~
                                      empa + rct_eligible + empa_rct_interaction + cohort_year_recoded + male + 
                                      age_at_initiation + eth_recoded + smoking_status + 
                                      sbp_at_initiation + bmi_at_initiation + hba1c_at_initiation + 
                                      ldl_at_initiation + ldl_sq + hdl_log + egfr_at_initiation + egfr_sq + 
                                      su_before_initiation + metformin_monotherapy_before_ini + 
                                      lipid_lowering_before_initiation + anti_htn_before_initiation + glp1ra_before_initiation + 
                                      anti_platelet_before_initiation + anticoag_before_initiation + insulin_before_initiation + 
                                      hf + smi + copd + ra + asthma + cancer_diagnosis_5y_before_cohor +
                                      epilepsy + dementia + ibd + liver + imd_quintile + cvd_diagnosis_before_cohort))

summary_cph_mi <- summary(pool(adjusted_cph_mi))
summary_cph_mi$exp_coef <- exp(summary_cph_mi$estimate)
print(summary_cph_mi)

summary_cph_mi <- summary_cph_mi %>%
  mutate(
    variable = term,
    HR = exp(estimate),
    p_value = p.value,
    lower_ci = exp(estimate - 1.96*std.error), 
    upper_ci = exp(estimate + 1.96*std.error)
  ) 

print(summary_cph_mi)

empa_coef <- summary_cph_mi %>%
  filter(term == "empa") %>% 
  pull(estimate)

empa_se <- summary_cph_mi %>%
  filter(term == "empa") %>% 
  pull(std.error)

interaction_coef <- summary_cph_mi %>%
  filter(term == "empa_rct_interaction") %>% 
  pull(estimate)

interaction_se <- summary_cph_mi %>%
  filter(term == "empa_rct_interaction") %>% 
  pull(std.error)

HR_empa_RCT_ineligible <- exp(empa_coef)
HR_empa_RCT_eligible <- exp(empa_coef + interaction_coef)

HR_empa_RCT_eligible_CI <- exp(c((empa_coef + interaction_coef)+1.96*sqrt(empa_se**2 + interaction_se**2), 
                                 (empa_coef + interaction_coef)-1.96*sqrt(empa_se^2 + interaction_se^2)))

HR_empa_RCT_ineligible_CI <- exp(c(empa_coef -1.96*empa_se, empa_coef + 1.96*empa_se))

print(c(HR_empa_RCT_ineligible, HR_empa_RCT_ineligible_CI))
print(c(HR_empa_RCT_eligible, HR_empa_RCT_eligible_CI))



