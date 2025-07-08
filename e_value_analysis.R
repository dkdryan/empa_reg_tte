### Trial emulation: EMPA-REG outcome ###
### E-value sensitivity analysis by David Ryan ###
### Note E-values are derived from the CPH model from the trial emulation population ###

#Import packages
library(dplyr)
library(survival)
library(survey)
library(lubridate)
library(sandwich)
library(ggplot2)
library(patchwork)
library(mice)
library(mitools)
library(knitr)
library(kableExtra)

#Import data 
data <- read.csv('S:/23_010_CVD_trial_emu/march_analysis/tte_estimates_mice.csv')

#E-value calculator 
calculate_evalue <- function(HR){
  if (HR < 1) {
    return((1/HR)+sqrt((1/HR)*((1/HR)-1)))
  } 
  else {
    return(HR+sqrt(HR*(HR-1)))
  }
  
}

e_value_est <- calculate_evalue(data$HR[data$term=="empa1"])
e_value_u95 <- calculate_evalue(data$lower_ci[data$term=="empa1"])
e_value_l95 <- calculate_evalue(data$upper_ci[data$term=="empa1"])

print(c(e_value_est, e_value_l95, e_value_u95))

#Create the outcome dataframe
e_outcome <- data.frame(
  term = character(), 
  RR = numeric(), 
  lower_95 = numeric(), 
  upper_95 = numeric()
)

#Add the e-values to outcome table 
e_outcome[1, ] <- c("E-value", e_value_est, e_value_l95, e_value_u95)
View(e_outcome)

#Add all the other HR/RR and 95% CI to the outcome table 
other_terms <- data[data$term != "empa1", c("term", "HR", "lower_ci", "upper_ci")]
colnames(other_terms) <- c("term", "RR", "lower_95", "upper_95")
data <- rbind(e_outcome, other_terms)
data$RR <- as.numeric(data$RR)
data <- data[order(-data$RR), ]

#Rename terms 

list(data$term)

# Define mappings for variables
mappings <- c(
  "eth_recoded_num3" = "Other ethnicity",
  "lipid_lowering_before_initiation1" = "Lipid lowering agent",
  "ibd1" = "Inflammatory bowel disease",
  "ldl_at_initiation" = "LDL cholesterol",
  "metformin_monotherapy_before_ini1" = "Metformin",
  "egfr_at_initiation" = "eGFR at initiation",
  "asthma1" = "Asthma",
  "smi1" = "Significant mental illness",
  "bmi_at_initiation" = "BMI at initiation",
  "sbp_at_initiation" = "SBP at initiation",
  "hba1c_at_initiation" = "HbA1c at initiation",
  "egfr_sq" = "eGFR squared",
  "hdl_log" = "HDL cholesterol",
  "anti_htn_before_initiation1" = "Antihypertensive agent",
  "smoking_status1" = "Current smoker",
  "smoking_status3" = "Ex-smoker",
  "ldl_sq" = "LDL cholesterol squared",
  "anti_platelet_before_initiation1" = "Anti-platelet agent",
  "age_at_initiation" = "Age at initiation",
  "cohort_year_recoded" = "Cohort year",
  "imd_quintile4" = "Index of Multiple Deprivation 4",
  "su_before_initiation1" = "Sulfonylurea",
  "epilepsy1" = "Epilepsy",
  "eth_recoded_num2" = "Black ethnicity",
  "imd_quintile3" = "Index of Multiple Deprivation 3",
  "eth_recoded_num1" = "Asian ethnicity",
  "imd_quintile5" = "Index of Multiple Deprivation 5",
  "anticoag_before_initiation1" = "Anticoagulant",
  "imd_quintile2" = "Index of Multiple Deprivation 2",
  "glp1ra_before_initiation1" = "GLP-1 receptor agonist",
  "liver1" = "Liver disease",
  "male1" = "Male",
  "hf1" = "Heart failure",
  "insulin_before_initiation1" = "Insulin",
  "eth_recoded4" = "Other ethnicity", 
  "dementia1"= "Dementia", 
  "eth_recoded2" = "Asian ethnicity", 
  "cancer_diagnosis_5y_before_cohor1" = "Recent cancer diagnosis", 
  "copd1" = "COPD", 
  "ra1" = "Rheumatoid arthritis", 
  "cohort_year_recoded2" = "Cohort entry: >= 2020", 
  "imd_quintile5th" = "5th IMD quintile", 
  "imd_quintile2nd" = "2nd IMD quintile", 
  "cohort_year_recoded1" = "Cohort entry: 2017 - 2019", 
  "imd_quintile4th" = "4th IMD quintile", 
  "imd_quintile3rd" = "3rd IMD quintile", 
  "eth_recoded1" = "Black ethnicity", 
  "eth_recoded3" = "Mixed ethnicity", 
  "cvd_diagnosis_before_cohort1" = "Cardiovascular disease")

data <- data %>%
  mutate(term = recode(term,!!!mappings))

View(data)

#Create forestplot for E-value analysis 
data$lower_95 <- as.numeric(data$lower_95)
data$upper_95 <- as.numeric(data$upper_95)

data <- data %>%
  arrange(-desc(abs(RR))) %>%
  mutate(term = factor(term, levels = term))

data <- data %>%
  mutate(highlight = ifelse(term == "E-value", "E-value", "Other"))

data$highlight <- factor(data$highlight, levels = unique(data$highlight))

fp <- ggplot(data, aes(x = term, y = RR, ymin= lower_95, ymax = upper_95))+ 
  geom_pointrange(data = subset(data, highlight != "E-value"), color = "black", size = 1.1, alpha = 0.6) +
  geom_pointrange(data = subset(data, highlight == "E-value"), color = "#7B1FA2", size = 1.1, alpha = 1)+
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") + 
  coord_flip()+ 
  scale_color_manual(name = "Legend", values = c("E-value" = "red", "Other" = "blue"))+
  theme_minimal()+
  labs(
    title = "E-value analysis: Forest plot of relative risks for confounders", 
    x="Variable", 
    y="Relative risk (RR) and 95% confidence interval")+ 
  theme(
    axis.text.y = element_text(size=12), 
    axis.text.x = element_text(size=12), 
    plot.title = element_text(size=14, face="bold"), 
    legend.position="bottom"
  )

fp

ggsave("S:/23_010_CVD_trial_emu/march_analysis/forest_plot_evalues.png", plot=fp, width = 10, height = 6, dpi = 300 )

