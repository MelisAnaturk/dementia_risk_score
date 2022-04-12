# Written by Melis Anaturk (July 2020)
# This r script calculates the DRS

# STEPS IN SCRIPT
# 1. Data set-up
# 2. Caide caculation

#------ 1. Data set up -----
# 1.1 load libraries
library(dvmisc)
library(dplyr)
library(psych)
library(data.table)
library(tidyverse)

# 1.2 data pathway
data_pathway = "../../raw_data/"
load(file = paste0(data_pathway, "ukbdata_diagnoses_baseline_diseasestatus_baselinemedications_ANUADRI_CAIDE_FRS_recoded.rda"))

# 1.3 subset data as DRS relies on quantiles of townsend deprivation
myvars <- c("Age_when_attended_assesment_centre_0_0","education_years", "Townsend_deprivation_0_0", "BMI_0_0",
            "Sex", "Sleeplesness_insomnia_0_0_bin", "family_history_of_dementia",
            "Diabetes_BIN_FINAL_0_0", "LDL_0_0", "Total_cholesterol_0_0", "HDL_cholesterol_0_0", 
            "depression_BIN_FINAL_0_0","TBI_BIN_FINAL_0_0", "stroke_TIA_BIN_FINAL", "Smoker_bin", "units_combined", 
            "Systolic_BP_bin", "IPAQ_activity_group_0_0", "Hearing_prob", "Sleep_BIN", "Antihypertensive_meds_0_0",
            "Fish_intake_BIN", "Social_engagement_0_2", "Atrial_Fibrillation_BIN_FINAL_0_0",
            "Social_engagement_0_1","dementia_BIN_TOTAL", "APOE_genotype_bin", "NSAIDs_0_0", "HRT_0_0", "statins_0_0")

pMiss <- function(x){sum(is.na(x))/length(x)*100}
apply(df[myvars],2,pMiss)
#n = 500827
# 1.4 subset to complete cases
df <- df[complete.cases(df[myvars]),]
#n=212192 after subsetting for complete cases across all variables

summary(df$dementia_BIN_TOTAL)
#0      1 
#209143   3049

# 1.5 restrict sample to middle aged = 40+
df <- subset(df, Age_when_attended_assesment_centre_0_0>=40)
summary(df$Age_when_attended_assesment_centre_0_0)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   40.00   50.00   57.00   56.23   63.00   73.00 

summary(df$dementia_BIN_TOTAL)
#0      1 
#209138   3049  

# 1.6 recode dementia variable such that only HES and primary care variables contribute to dementia ascertain
myvars <- c("secondary_care_diagnosis_of_Dementia", "death_report_dementia", "primary_care_diagnosis_for_dementia", "primary_care_prescription_for_Dementia")
sapply(df[myvars], class) 
 
 df$dementia_BIN_TOTAL <- apply(df[, myvars], 1, function(x) {
   if(any(x %in% c("1", 1))) { # "1" or 1 so it looks for 1 coded as a factor and 1 as a numeric
     return(1)
   } else {
     return(0)
   }
 })
 
df$dementia_BIN_TOTAL <- as.factor(df$dementia_BIN_TOTAL) 
summary(df$dementia_BIN_TOTAL)
 #0      1 
 #209177   3010 
 
#------ 2. Compute DRS ------
# as our sample is between 40-73 years old, we will use the weights developed within a sample of 60-79 year olds
# 2.1 Age
df$DRS_age_score <- df$Age_when_attended_assesment_centre_0_0
 
# 2.2 BMI
df$DRS_BMI_score <- df$BMI_0_0

# 2.3 Sex
df$DRS_sex_score <- ifelse(df$Sex==0,1,0)

# 2.4 Hypertension
df$DRS_hypertensive_med_score <- ifelse(df$Antihypertensive_meds_0_0==1,1,0)

# 2.5 Townsend deprivation
summary(df$Townsend_deprivation_0_0)
df$Townsend_deprivation_Groups_0_0 <- quant_groups(df$Townsend_deprivation_0_0, groups = 5)
df$Townsend_deprivation_Groups_0_0 <- as.factor(df$Townsend_deprivation_Groups_0_0)
summary(df$Townsend_deprivation_Groups_0_0)

df$Townsend_deprivation_Groups_0_0 <- ifelse(df$Townsend_deprivation_Groups_0_0=="[-6.26,-3.98]", 0,
                                             ifelse(df$Townsend_deprivation_Groups_0_0=="(-3.98,-2.87]", 1,
                                                    ifelse(df$Townsend_deprivation_Groups_0_0=="(-2.87,-1.5]", 2,
                                                           ifelse(df$Townsend_deprivation_Groups_0_0=="(-1.5,0.972]", 3,
                                                                  ifelse(df$Townsend_deprivation_Groups_0_0=="(0.972,10.6]", 4,NA)))))

summary(df$Townsend_deprivation_Groups_0_0)


df$DRS_Townsend_deprivation_0_0_group2 <- ifelse(df$Townsend_deprivation_Groups_0_0==1, 1, 0)
df$DRS_Townsend_deprivation_0_0_group3 <- ifelse(df$Townsend_deprivation_Groups_0_0==2, 1, 0)
df$DRS_Townsend_deprivation_0_0_group4 <- ifelse(df$Townsend_deprivation_Groups_0_0==3, 1, 0)
df$DRS_Townsend_deprivation_0_0_group5 <- ifelse(df$Townsend_deprivation_Groups_0_0==4, 1, 0)

# 2.6 Smoking
df$DRS_current_smoker_score <- ifelse(df$Smoking_Status_0_0==2, 1, ifelse(df$Smoking_Status_0_0==1, 0,
                               ifelse(df$Smoking_Status_0_0==0,0,NA)))

df$DRS_ex_smoker_score <- ifelse(df$Smoking_Status_0_0==2, 0, ifelse(df$Smoking_Status_0_0==1, 1,
                          ifelse(df$Smoking_Status_0_0==0,0,NA)))

# 2.7 Alcohol
df$DRS_alcohol_score <- ifelse(df$Ever_addicted_to_alcohol_0_0==1|df$Sex==0&df$units_combined>=14|df$Sex==1&df$units_combined>=21, 1, 0)
df$DRS_alcohol_score[is.na(df$DRS_alcohol_score)] <- 0

View(df[,c("DRS_alcohol_score", "Sex", "units_combined", "Ever_addicted_to_alcohol_0_0")])

# 2.8 Depression
df$DRS_depression_score <- ifelse(df$depression_BIN_FINAL_0_0==1|df$Antidepressant_meds_0_0==1, 1,0)

# 2.9 Stroke
df$DRS_stroke_score <- ifelse(df$stroke_TIA_BIN_FINAL==1,1,0)

# 2.10 Atrial fibrillation
df$DRS_atrial_fib_score <- ifelse(df$Atrial_Fibrillation_BIN_FINAL_0_0==1,1,0)

# 2.11 Diabetes
df$DRS_diabetes_score <- ifelse(df$Diabetes_II_BIN_FINAL_0_0==1,1,0)
#mar 7 rp tries:
#df$DRS_diabetes_score <- ifelse(df$Diabetes_BIN_FINAL_0_0==1,1,0)
#keeping definition based on Diabetes_II_BIN for now, to confirm.

# 2.12 Calendar year
df$DRS_calendar_year_score <- as.integer(format(df$baseline_date, "%Y"))

# 2.13 Aspirin
df$DRS_aspirin_score <- ifelse(df$Aspirin_0_0==1,1,0)

# 2.14 Calculate total DRS score
df$DRS_total <-  0.20921*(df$Age_when_attended_assesment_centre_0_0 - 65.608) + -0.00339*(df$Age_when_attended_assesment_centre_0_0 - 65.608)*(df$Age_when_attended_assesment_centre_0_0 - 65.608) + 
                      -0.0616*(df$BMI_0_0 - 27.501) + 0.002508*(df$BMI_0_0 - 27.501)*(df$BMI_0_0 - 27.501) + 
                      0.12854*df$DRS_sex_score + -0.13199*df$DRS_hypertensive_med_score + 0.04477*(df$DRS_calendar_year_score - 2003.719) + 
                      0.013371*(df$DRS_Townsend_deprivation_0_0_group2) + 0.117904*(df$DRS_Townsend_deprivation_0_0_group3) + 
                      0.201776*(df$DRS_Townsend_deprivation_0_0_group4) + 0.225529*df$DRS_Townsend_deprivation_0_0_group5 + 
                      -0.06792*df$DRS_ex_smoker_score + -0.08657*df$DRS_current_smoker_score +
                      0.443535*df$DRS_alcohol_score + 0.833612*df$DRS_depression_score + 0.252833*df$DRS_aspirin_score + 
                      0.577207*df$DRS_stroke_score + 0.220728*df$DRS_atrial_fib_score + 0.286701*df$DRS_diabetes_score

# predicted probability
df$DRS_predicted_prob <- 1-0.9969^exp(df$DRS_total)
summary(df$DRS_total)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-7.8012 -3.5008 -1.2630 -1.7712  0.1145  3.4190 
summary(df$DRS_predicted_prob)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#1.270e-06 9.368e-05 8.777e-04 2.647e-03 3.476e-03 9.046e-02 
# examine missingness
vars <- c("DRS_total", "DRS_predicted_prob", "DRS_sex_score", "DRS_hypertensive_med_score", "DRS_calendar_year_score", "DRS_Townsend_deprivation_0_0_group2", "DRS_Townsend_deprivation_0_0_group3", "DRS_Townsend_deprivation_0_0_group4", "DRS_Townsend_deprivation_0_0_group5", "DRS_ex_smoker_score", "DRS_current_smoker_score", "DRS_alcohol_score", "DRS_depression_score", "DRS_aspirin_score", "DRS_stroke_score", "DRS_atrial_fib_score", "DRS_diabetes_score")
apply(df[vars],2,pMiss)

save(df, file = paste0(data_pathway, "ukbdata_diagnoses_baseline_diseasestatus_baselinemedications_ANUADRI_CAIDE_FRS_recoded_DRS.rda"))


