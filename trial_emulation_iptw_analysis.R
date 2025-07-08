### Trial emulation: EMPA-REG outcome ###
### IPTW CPH model - both as ATE and ATT ###

library(dplyr)
library(tidyr)
library(mice)
library(tidyverse)
library(tidymodels)
library(bootImpute)
library(survival)
library(survminer)
library(lubridate)
library(rms)
library(pec)
library(survey)
library(MatchThem)
library(cobalt)

gc()
rm(list=ls())

# Data organisation -------------------------------------------------------
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
             predictorMatrix=pm, 
             method = meth,
             printFlag=TRUE)

#Save files
save(imps, file="S:/23_010_CVD_trial_emu/march_analysis/tte_2.RData")
load("S:/23_010_CVD_trial_emu/march_analysis/imputation/imp_file_2.RData")

#Re-format dates
data$last_date <- ymd(data$last_contact_correct)
data$first_date <- ymd(data$cohort_entry_date)

#multiple imputation analysis: ATE
ate_imps <- weightthem(
  empa ~cohort_year_recoded+
    male+ 
    age_at_initiation+
    eth_recoded+
    smoking_status +
    sbp_at_initiation+
    bmi_at_initiation+
    hba1c_at_initiation+
    ldl_at_initiation+
    ldl_sq+
    hdl_log+
    egfr_at_initiation+
    egfr_sq+
    su_before_initiation+
    metformin_monotherapy_before_ini+ 
    lipid_lowering_before_initiation+
    anti_htn_before_initiation+
    glp1ra_before_initiation+
    anti_platelet_before_initiation+
    anticoag_before_initiation+
    insulin_before_initiation+
    hf+
    smi+
    copd+
    ra+
    asthma+
    epilepsy+
    dementia+
    ibd+
    liver+
    imd_quintile+
    cvd_diagnosis_before_cohort+
    cancer_diagnosis_5y_before_cohor,
  data = imps, 
  method = "ps",
  approach = "within", 
  estimand = "ATE", 
  include = TRUE)

iptw_ate <- with(ate_imps, coxph(Surv(as.numeric(difftime(last_contact_correct, cohort_entry_date, units="days"))/365.25, 
                                      event = died, 
                                      origin=0) ~ empa, weights = weights, robust = TRUE))

pooled_cox_ate <- pool(iptw_ate)

res_ate <- summary(pooled_cox_ate, conf.int = TRUE, exponentiate = TRUE)
res_ate

#density curves 
df1 <- complete(ate_imps, 3)
df1$weights
ggplot(df1, aes(x=weights, fill = factor(empa)))+ 
  geom_density(alpha = 0.5)+
  #scale_x_log10()+
  coord_cartesian(xlim = c(0, 8))

summary(df1$weights)

#SMD 
bal_ate <- bal.tab(ate_imps, un = TRUE)

smd_data_ate <- bal_ate$Balance.Across.Imputations

smd_data_ate <- smd_data_ate %>%
  tibble::rownames_to_column(var="Variable") %>%
  select(Variable, Mean.Diff.Un, Mean.Diff.Adj )

plot_data_ate <- smd_data_ate %>% 
  select(Variable, Mean.Diff.Un, Mean.Diff.Adj) %>%
  pivot_longer(cols = c(Mean.Diff.Un, Mean.Diff.Adj), 
               names_to = "SMD_Type", 
               values_to = "SMD") %>% 
  mutate(SMD_Type = recode(SMD_Type, 
                           "Mean.Diff.Un" = "Unweighted", 
                           "Mean.Diff.Adj" = "Weighted"))

plot_data_ate <- plot_data_ate %>%
  filter(!Variable %in% c("distance", "prop.score"))

#Rename terms 

list(data$term)

# Define mappings for variables
mappings <- c(
  "eth_recoded_3" = "Other ethnicity",
  "lipid_lowering_before_initiation" = "Lipid lowering agent",
  "ibd" = "Inflammatory bowel disease",
  "ldl_at_initiation" = "LDL cholesterol",
  "metformin_monotherapy_before_ini" = "Metformin",
  "egfr_at_initiation" = "eGFR at initiation",
  "asthma" = "Asthma",
  "smi" = "Significant mental illness",
  "bmi_at_initiation" = "BMI at initiation",
  "sbp_at_initiation" = "SBP at initiation",
  "hba1c_at_initiation" = "HbA1c at initiation",
  "egfr_sq" = "eGFR squared",
  "hdl_log" = "HDL cholesterol",
  "anti_htn_before_initiation" = "Antihypertensive agent",
  "smoking_status_1" = "Current smoker",
  "smoking_status_3" = "Ex-smoker",
  "smoking_status_0" = "Non-smoker", 
  "ldl_sq" = "LDL cholesterol squared",
  "anti_platelet_before_initiation" = "Anti-platelet agent",
  "age_at_initiation" = "Age at initiation",
  "imd_quintile_4th" = "Index of Multiple Deprivation 4 quintile",
  "su_before_initiation" = "Sulfonylurea",
  "epilepsy" = "Epilepsy",
  "eth_recoded_2" = "Asian ethnicity",
  "imd_quintile_3rd" = "Index of Multiple Deprivation 3 quintile",
  "eth_recoded_num" = "Asian ethnicity",
  "anticoag_before_initiation" = "Anticoagulant",
  "glp1ra_before_initiation" = "GLP-1 receptor agonist",
  "liver" = "Liver disease",
  "male" = "Male",
  "hf" = "Heart failure",
  "insulin_before_initiation" = "Insulin",
  "eth_recoded_4" = "Other ethnicity", 
  "dementia"= "Dementia", 
  "eth_recoded_2" = "Asian ethnicity", 
  "cancer_diagnosis_5y_before_cohor" = "Recent cancer diagnosis", 
  "copd" = "COPD", 
  "ra" = "Rheumatoid arthritis", 
  "cohort_year_recoded_2" = "Cohort entry: >= 2020", 
  "imd_quintile_2nd" = "Index of Multiple Deprivation 2 quintile", 
  "cohort_year_recoded_0" = "Cohort entry: 2017 - 2019", 
  "cohort_year_recoded_1" = "Cohort entry: 2014 - 2016", 
  "eth_recoded_0" = "White ethnicity", 
  "imd_quintile_1st" = "Index of Multiple Deprivation 1 quintile", 
  "eth_recoded_1" = "Black ethnicity", 
  "eth_recoded_3" = "Mixed ethnicity", 
  "cvd_diagnosis_before_cohort" = "Cardiovascular disease", 
  "imd_quintile_5th" = "Index of Multiple Deprivation 5 quintile")


plot_data_ate <- plot_data_ate %>%
  mutate(Variable = recode(Variable,!!!mappings))

ggplot(plot_data_ate, aes(x = SMD, y=reorder(Variable, SMD), color = SMD_Type))+
  labs(y="Variable")+
  guides(fill = guide_legend(title = "Weighting status"))+
  geom_point(size = 3)+
  geom_vline(xintercept = 0.1, linetype = 'dashed', color='gray40')+
  geom_vline(xintercept = -0.1, linetype = 'dashed', color='gray40')+
  labs(
    title = "Standardised mean difference: ATE weighting", 
    x = "Standardised mean difference", 
    Y = "Variable", 
    color = "Weighting status"
  )+ theme_minimal()+theme(plot.title=element_text(hjust = 0.5))

### ATT estimation 
att_imps <- weightthem(
  empa ~cohort_year_recoded+
    male+ 
    age_at_initiation+
    eth_recoded+
    smoking_status +
    sbp_at_initiation+
    bmi_at_initiation+
    hba1c_at_initiation+
    ldl_at_initiation+
    ldl_sq+
    hdl_log+
    egfr_at_initiation+
    egfr_sq+
    su_before_initiation+
    metformin_monotherapy_before_ini+ 
    lipid_lowering_before_initiation+
    anti_htn_before_initiation+
    glp1ra_before_initiation+
    anti_platelet_before_initiation+
    anticoag_before_initiation+
    insulin_before_initiation+
    hf+
    smi+
    copd+
    ra+
    asthma+
    epilepsy+
    dementia+
    ibd+
    liver+
    imd_quintile+
    cvd_diagnosis_before_cohort+
    cancer_diagnosis_5y_before_cohor,
  data = imps, 
  method = "ps",
  approach = "within", 
  estimand = "ATT", 
  include = TRUE)

iptw_att <- with(att_imps, coxph(Surv(as.numeric(difftime(last_contact_correct, cohort_entry_date, units="days"))/365.25, 
                                      event = died, 
                                      origin=0) ~ empa, weights = weights, robust = TRUE))

pooled_cox_att <- pool(iptw_att)

res_att <- summary(pooled_cox_att, conf.int = TRUE, exponentiate = TRUE)
res_att

#check att balance 
att_bal <- bal.tab(att_imps, un = TRUE)

att_smd_data <- att_bal$Balance.Across.Imputations

att_smd_data <- att_smd_data %>%
  tibble::rownames_to_column(var="Variable") %>%
  select(Variable, Mean.Diff.Un, Mean.Diff.Adj )

att_plot_data <- att_smd_data %>% 
  select(Variable, Mean.Diff.Un, Mean.Diff.Adj) %>%
  pivot_longer(cols = c(Mean.Diff.Un, Mean.Diff.Adj), 
               names_to = "SMD_Type", 
               values_to = "SMD") %>% 
  mutate(SMD_Type = recode(SMD_Type, 
                           "Mean.Diff.Un" = "Unweighted", 
                           "Mean.Diff.Adj" = "Weighted"))

att_plot_data <- att_plot_data %>%
  filter(!Variable %in% c("distance", "prop.score"))

att_plot_data <- att_plot_data %>%
  mutate(Variable = recode(Variable,!!!mappings))

ggplot(att_plot_data, aes(x = SMD, y=reorder(Variable, SMD), color = SMD_Type))+
  labs(y="Variable")+
  guides(fill = guide_legend(title = "Weighting status"))+
  geom_point(size = 3)+
  geom_vline(xintercept = 0.1, linetype = 'dashed', color='gray40')+
  geom_vline(xintercept = -0.1, linetype = 'dashed', color='gray40')+
  labs(
    title = "Standardised mean difference: ATT weighting", 
    x = "Standardised mean difference", 
    Y = "Variable", 
    color = "Weighting status"
  )+ theme_minimal()+theme(plot.title=element_text(hjust = 0.5))

res_att
res_ate
