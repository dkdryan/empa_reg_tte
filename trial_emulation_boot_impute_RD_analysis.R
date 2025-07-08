# Bootstrap MI procedure for estimating risk difference -------------------
library(dplyr)
library(data.table)
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

# Overall analysis (not stratified) -------------------------------------------------------

#Import data
data <- read.csv("S:/23_010_CVD_trial_emu/march_analysis/rw_230325.csv")
dim(data)

#Add additional variables 
data$egfr_sq <- data$egfr_at_initiation^2

#Re-format dates
data$last_date <- ymd(data$last_contact_correct)
data$first_date <- ymd(data$cohort_entry_date)

#Define survival object 
data$surv_obj <- with(data, Surv(time = as.numeric(difftime(last_date, first_date, units="days"))/365.25, event = died))
summary(data$surv_obj)

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
View(pm)

pm[c("ldl_sq", "egfr_sq"), ] <- 0 
pm[, c("ldl_sq", "egfr_sq")] <- 1
pm[c('last_contact_correct', 'person_id', 'cohort_entry_date'), ] <- 0 
pm[, c('last_contact_correct', 'person_id', 'cohort_entry_date')] <- 0 

View(pm[vars_to_impute, ])
View(meth)

#run the bootstrap imputation 
imp_file <- bootMice(data_selected, nBoot = 20, nImp = 2, predictorMatrix=pm, method=meth, seed = 3001) #maxit = 5 by default 
save(imp_file, file = "S:/23_010_CVD_trial_emu/march_analysis/imputation/imp_file.RData")

#load imps 
#load("S:/23_010_CVD_trial_emu/march_analysis/imp_file.RData")

#set counter 
iter_counter <- 0

#set up datatable to save estimates 
results_file <- "S:/23_010_CVD_trial_emu/march_analysis/imputation/imps_res_3yrd_nnt_ate.csv"

header_df <- data.frame(RiskDiff = numeric(0), 
                        NNT = numeric(0))

write.table(
  x = header_df, 
  file = results_file, 
  sep=",", 
  row.names = FALSE, 
  col.names = TRUE, 
  append=FALSE
)

#function 
anaylseImps_ra <- function(inputData){
  library(data.table)
  library(dplyr)
  library(mice)
  library(tidyverse)
  library(tidymodels)
  library(bootImpute)
  library(survival)
  library(survminer)
  library(lubridate)
  library(rms)
  library(pec)
  library(parallel)
  
  tolerance <- 1e-6
  
  #counter 
  iter_counter <<- iter_counter + 1
  
  #obtain imputation number and dataframe 
  df <- setDT(inputData)
  
  #write individual log file for each imputation 
  log_file <- paste0("S:/23_010_CVD_trial_emu/march_analysis/imputation/imp_log.log")
  cat("", file = log_file, append=TRUE)
  cat(paste0("Dataframe ", iter_counter, " analysis started at ", Sys.time(), "\n"), file = log_file, append=TRUE)
  
  surv_obj <- Surv(time = as.numeric(difftime(df$last_contact_correct, df$cohort_entry_date, units="days"))/365.25, event = df$died, origin=0)
  
  cox_adj <- coxph(surv_obj ~ empa + 
                     cohort_year_recoded+ 
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
                     cancer_diagnosis_5y_before_cohor+
                     cvd_diagnosis_before_cohort,
                   data=df, x=TRUE,
                   method="breslow")
  
  
  #empa df 
  empa_df <- df
  empa_df[, empa:=as.factor(1)]
  
  #dpp4i_df 
  dpp4i_df <- copy(empa_df)
  dpp4i_df[, empa:=as.factor(0)]
  
  #predict survival at 3 years  
  surv.empa <- summary(survfit(cox_adj, empa_df), times = 3, extend=TRUE)$surv
  surv.dpp4i <- summary(survfit(cox_adj, dpp4i_df), times = 3, extend=TRUE)$surv
  
  #find mean survival at time t
  survmean.empa <- mean(surv.empa, na.rm=TRUE)
  survmean.dpp4i <- mean(surv.dpp4i, na.rm=TRUE)
  
  #find risk at time t 
  riskmean.empa <- 1-survmean.empa
  riskmean.dpp4i <- 1-survmean.dpp4i
  
  #risk difference at time t 
  riskdiff = riskmean.empa - riskmean.dpp4i
  
  #NNT 
  nnt = ifelse(abs(riskdiff)>tolerance, abs(1/riskdiff), NA) 
  
  #Show ouput
  res = c(riskdiff, nnt)
  print(c(iter_counter, res))
  
  #Post to log file 
  cat(paste0("Results:\n"), file = log_file, append=TRUE)
  cat(paste(capture.output(print(res)), collapse = "\n"), "\n", file = log_file, append=TRUE)
  cat(paste0("Finished imputation for dataframe: ", iter_counter, " at ", Sys.time(), "\n\n"), file = log_file, append=TRUE)
  res_df <- as.data.table(t(res))
  
  write.table(res_df, 
              file = results_file,  
              row.names = FALSE, 
              append=TRUE, 
              col.names =FALSE,
              sep = ",")
  
  #Save output
  return(res)
}

#run analysis
ests_ra <-bootImputeAnalyse(imp_file, anaylseImps_ra, quiet = FALSE)
write.csv(ests_ra, file = "S:/23_010_CVD_trial_emu/march_analysis/imputation/pooled_3yrd_nnt_main_tte.csv")

#VH pooling: RD  
res <- read.csv("S:/23_010_CVD_trial_emu/march_analysis/imputation/imps_res_3yrd_nnt_ate.csv")

boot_estimates <- res$RiskDiff[4:length(res$RiskDiff)]
point_est <- mean(boot_estimates)
point_est

ci <- quantile(boot_estimates, probs = c(0.025, 0.975))
pooled_result <- list(estimate = point_est, ci.lower = ci[1], ci.upper = ci[2])
pooled_result

#VH pooling: NNT 
res <- read.csv("S:/23_010_CVD_trial_emu/march_analysis/imputation/imps_res_3yrd_nnt_ate.csv")
boot_estimates <- res$NNT[3:length(res$NNT)]
point_est <- mean(boot_estimates)
point_est

ci <- quantile(boot_estimates, probs = c(0.025, 0.975))
pooled_result <- list(estimate = point_est, ci.lower = ci[1], ci.upper = ci[2])
pooled_result

#VH pooling RD SE 
res <- read.csv("S:/23_010_CVD_trial_emu/march_analysis/imputation/imps_res_3yrd_nnt_ate.csv")

boot_estimates <- res$RiskDiff[4:length(res$RiskDiff)]

B <- length(boot_estimates)
D <- 2

boot_estimates

boot_estimates <- res$RiskDiff[3:length(res$RiskDiff)]
theta_matrix <- matrix(boot_estimates, nrow = B, ncol = D, byrow=TRUE)



# RCT eligible population-------------------------------------------------------
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
imp_int_file <- bootMice(data_selected, nBoot = 20, nImp = 2, predictorMatrix=pm, method=meth, seed = 3001) #maxit = 5 by default

#Save files
save(imp_int_file, file="S:/23_010_CVD_trial_emu/march_analysis/tte_int_bsmi.RData")
load("S:/23_010_CVD_trial_emu/march_analysis/tte_int_bsmi.RData")

#set counter 
iter_counter <- 0

#set up datatable to save estimates 
res_rct_eligible <- "S:/23_010_CVD_trial_emu/march_analysis/imputation/imps_res_3yrd_nnt_rct_eligible.csv"

header_df <- data.frame(RiskDiff = numeric(0), 
                        NNT = numeric(0))

write.table(
  x = header_df, 
  file = res_rct_eligible, 
  sep=",", 
  row.names = FALSE, 
  col.names = TRUE, 
  append=FALSE
)

# function rct eligible population-------------------------------------------------------
anaylseImps_ra_rct_eligible <- function(inputData){
  library(data.table)
  library(dplyr)
  library(mice)
  library(tidyverse)
  library(tidymodels)
  library(bootImpute)
  library(survival)
  library(survminer)
  library(lubridate)
  library(rms)
  library(pec)
  library(parallel)
  
  tolerance <- 1e-6
  
  #counter 
  iter_counter <<- iter_counter + 1
  
  #arrange data in datatable format 
  df <- setDT(inputData)
  
  #write individual log file for each imputation 
  log_file <- paste0("S:/23_010_CVD_trial_emu/march_analysis/imputation/imp_log_rct_eligible.log")
  cat("", file = log_file, append=TRUE)
  cat(paste0("Dataframe ", iter_counter, " analysis started at ", Sys.time(), "\n"), file = log_file, append=TRUE)
  
  #set the survival object
  surv_obj <- Surv(time = as.numeric(difftime(df$last_contact_correct, df$cohort_entry_date, units="days"))/365.25, event = df$died, origin=0)
  
  #build the CPH model 
  cox_adj <- coxph(surv_obj ~ empa + rct_eligible+ empa_rct_interaction + 
                     cohort_year_recoded+ 
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
                     cancer_diagnosis_5y_before_cohor+
                     cvd_diagnosis_before_cohort,
                   data=df, x=TRUE,
                   method="breslow")
  
  cat(paste0("Cox model completed ", "\n"), file = log_file, append=TRUE)
  
  #empa and rct eligible df 
  empa_df <- data.table::copy(df)
  dpp4i_df <- data.table::copy(df)
  
  empa_df[, empa := 1] 
  empa_df[, empa := factor(empa, levels = c(0, 1))] 
  
  empa_df[, rct_eligible := 1]
  empa_df[, rct_eligible := factor(rct_eligible, levels = c(0, 1))]
  
  empa_df[, empa_rct_interaction := as.numeric(as.character(empa))*as.numeric(as.character(rct_eligible))]
  empa_df[, empa_rct_interaction := factor(empa_rct_interaction, levels = c(0, 1))]
  
  cat(paste0("Empa: empa: ", empa_df[, .N, by = empa], "\n"), file = log_file, append=TRUE)
  cat(paste0("Empa: rct_eligible: ", empa_df[, .N, by = rct_eligible], "\n"), file = log_file, append=TRUE)
  cat(paste0("Empa: interaction: ", empa_df[, .N, by = empa_rct_interaction], "\n"), file = log_file, append=TRUE)
  
  #dpp4i and rct eligible 
  dpp4i_df[, empa := 0] 
  dpp4i_df[, empa:= factor(empa, levels = c(0, 1))]
  
  dpp4i_df[, rct_eligible := 1]
  dpp4i_df[, rct_eligible:= factor(rct_eligible, levels = c(0, 1))]
  
  dpp4i_df[, empa_rct_interaction := as.numeric(as.character(empa))*as.numeric(as.character(rct_eligible))]
  dpp4i_df[, empa_rct_interaction := factor(empa_rct_interaction, levels = c(0, 1))]
  
  cat(paste0("DPP4i: empa: ", dpp4i_df[, .N, by = empa], "\n"), file = log_file, append=TRUE)
  cat(paste0("DPP4i: rct_eligible: ", dpp4i_df[, .N, by = rct_eligible], "\n"), file = log_file, append=TRUE)
  cat(paste0("DPP4i: interaction: ", dpp4i_df[, .N, by = empa_rct_interaction], "\n"), file = log_file, append=TRUE)
  
  
  #predict survival at 3 years  
  surv.empa <- summary(survfit(cox_adj, empa_df), times = 3, extend=TRUE)$surv
  surv.dpp4i <- summary(survfit(cox_adj, dpp4i_df), times = 3, extend=TRUE)$surv
  
  #find mean survival at time t
  survmean.empa <- mean(surv.empa, na.rm=TRUE)
  survmean.dpp4i <- mean(surv.dpp4i, na.rm=TRUE)
  cat(paste0("Empa: mean survival: ", survmean.empa, "\n"), file = log_file, append=TRUE)
  cat(paste0("DPP4i: mean survival: ", survmean.dpp4i, "\n"), file = log_file, append=TRUE)
  print(c(surv.empa, surv.dpp4i))
  
  #find risk at time t 
  riskmean.empa <- 1-survmean.empa
  riskmean.dpp4i <- 1-survmean.dpp4i
  
  #risk difference at time t 
  riskdiff = riskmean.empa - riskmean.dpp4i
  
  #NNT 
  nnt = ifelse(abs(riskdiff)>tolerance, abs(1/riskdiff), NA) 
  
  #Show ouput
  res = c(riskdiff, nnt)
  print(res)
  
  #Post to log file 
  cat(paste0("Results:\n"), file = log_file, append=TRUE)
  cat(paste(capture.output(print(res)), collapse = "\n"), "\n", file = log_file, append=TRUE)
  cat(paste0("Finished imputation for dataframe: ", iter_counter, " at ", Sys.time(), "\n\n"), file = log_file, append=TRUE)
  res_df <- as.data.table(t(res))
  
  write.table(res_df, 
              file = res_rct_eligible,  
              row.names = FALSE, 
              append=TRUE, 
              col.names =FALSE,
              sep = ",")
  
  #Save output
  return(res)
}
#run analysis
ests_ra_rct_eligible <-bootImputeAnalyse(imp_int_file, anaylseImps_ra_rct_eligible, quiet = FALSE)
write.csv(ests_ra_rct_eligible, file = "S:/23_010_CVD_trial_emu/march_analysis/imputation/pooled_3yrd_nnt_rct_eligible.csv")

# RCT ineligible population-------------------------------------------------------
#set counter 
iter_counter <- 0

#set up datatable to save estimates 
res_rct_ineligible <- "S:/23_010_CVD_trial_emu/march_analysis/imputation/imps_res_3yrd_nnt_rct_ineligible.csv"

header_df <- data.frame(RiskDiff = numeric(0), 
                        NNT = numeric(0))

write.table(
  x = header_df, 
  file = res_rct_ineligible, 
  sep=",", 
  row.names = FALSE, 
  col.names = TRUE, 
  append=FALSE
)

anaylseImps_ra_rct_ineligible <- function(inputData){
  library(data.table)
  library(dplyr)
  library(mice)
  library(tidyverse)
  library(tidymodels)
  library(bootImpute)
  library(survival)
  library(survminer)
  library(lubridate)
  library(rms)
  library(pec)
  library(parallel)
  
  tolerance <- 1e-6
  
  #counter 
  iter_counter <<- iter_counter + 1
  
  #obtain imputation number and dataframe 
  df <- setDT(inputData)
  
  #write individual log file for each imputation 
  log_file <- paste0("S:/23_010_CVD_trial_emu/march_analysis/imputation/imp_rct_ineligible_log.log")
  cat("", file = log_file, append=TRUE)
  cat(paste0("Dataframe ", iter_counter, " rct ineligible analysis started at ", Sys.time(), "\n"), file = log_file, append=TRUE)
  
  surv_obj <- Surv(time = as.numeric(difftime(df$last_contact_correct, df$cohort_entry_date, units="days"))/365.25, event = df$died, origin=0)
  
  cox_adj <- coxph(surv_obj ~ empa + empa_rct_interaction+rct_eligible+
                     cohort_year_recoded+ 
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
                     cancer_diagnosis_5y_before_cohor+
                     cvd_diagnosis_before_cohort,
                   data=df, x=TRUE,
                   method="breslow")
  
  
  #empa and rct ineligible df 
  empa_df <- df
  empa_df[, empa:=as.factor(1)]
  empa_df[, rct_eligible:=as.factor(0)]
  empa_df[, empa_rct_interaction := as.factor(0)]
  
  #dpp4i and rct ineligible 
  dpp4i_df <- df
  dpp4i_df[, empa:=as.factor(0)]
  dpp4i_df[, rct_eligible:=as.factor(0)]
  dpp4i_df[, empa_rct_interaction := as.factor(0)]
  
  #predict survival at 3 years  
  surv.empa <- summary(survfit(cox_adj, empa_df), times = 3, extend=TRUE)$surv
  surv.dpp4i <- summary(survfit(cox_adj, dpp4i_df), times = 3, extend=TRUE)$surv
  
  #find mean survival at time t
  survmean.empa <- mean(surv.empa, na.rm=TRUE)
  survmean.dpp4i <- mean(surv.dpp4i, na.rm=TRUE)
  
  #find risk at time t 
  riskmean.empa <- 1-survmean.empa
  riskmean.dpp4i <- 1-survmean.dpp4i
  
  #risk difference at time t 
  riskdiff = riskmean.empa - riskmean.dpp4i
  
  #NNT 
  nnt = ifelse(abs(riskdiff)>tolerance, abs(1/riskdiff), NA) 
  
  #Show ouput
  res = c(riskdiff, nnt)
  print(c(iter_counter, res))
  
  #Post to log file 
  cat(paste0("Results:\n"), file = log_file, append=TRUE)
  cat(paste(capture.output(print(res)), collapse = "\n"), "\n", file = log_file, append=TRUE)
  cat(paste0("Finished imputation for dataframe: ", iter_counter, " at ", Sys.time(), "\n\n"), file = log_file, append=TRUE)
  res_df <- as.data.table(t(res))
  
  write.table(res_df, 
              file = res_rct_ineligible,  
              row.names = FALSE, 
              append=TRUE, 
              col.names =FALSE,
              sep = ",")
  
  #Save output
  return(res)
}

#run analysis
ests_ra_rct_ineligible <-bootImputeAnalyse(imp_int_file, anaylseImps_ra_rct_ineligible, quiet = FALSE)
write.csv(ests_ra_rct_ineligible, file = "S:/23_010_CVD_trial_emu/march_analysis/imputation/pooled_3yrd_nnt_rct_ineligible.csv")
