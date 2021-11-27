# Created by Raihaan Patel
# Modified by Melis Anaturk (Feb 2020)
# This r script calculates three versions of the Framingham 10-year CVD risk score: 
#(1) beta-weighted with/ (2) without BMI
# and (3) points based FRS

# STEPS IN SCRIPT
# 1. Data cleaning
# 2. Data recoding
# 3. Caide caculation

#------ 1. DATA SET UP ----
# 1.1 load packages
#library(car)
#library(sjstats)
#library(pROC)
#library(ggplot2)
#library(janitor)
#library(dplyr)

# 1.2 data pathway
data_pathway = "~/Documents/Oxford_DPhil/Biobank_Analyses/Dementia_risk_project/ukb_data/For_Rai/raw_data/"

#------ 2. DATA RECODING ----
# 2.1 Age
df$frs_log_age    <- ifelse(is.na(df$Age_when_attended_assesment_centre_0_0), NA, 
                    log(df$Age_when_attended_assesment_centre_0_0))

# 2.2 Systolic BP
df$frs_log_sbp    <- ifelse(is.na(df$Systolic_BP_auto_mean), NA, log(df$Systolic_BP_auto_mean))

# 2.3 Smoking
df$frs_smok_score <- ifelse(df$Smoking_Status_0_0==2, 1, ifelse(df$Smoking_Status_0_0==1, 0,
                    ifelse(df$Smoking_Status_0_0==0,0,NA)))
summary(as.factor(df$frs_smok_score))

# 2.4 Diabetes
df$frs_diab_score <- ifelse(df$Diabetes_BIN_FINAL_0_0==1, 1, 
                     ifelse(df$Diabetes_BIN_FINAL_0_0==0, 0, NA))
summary(as.factor(df$frs_diab_score))

# 2.5 Cholesterol
#To  convert from mmol/L to mg/dL, you need to multiply by 38.67 
df$HDL_cholesterol_0_0_mgdl <-  df$HDL_cholesterol_0_0*38.67
summary(df$HDL_cholesterol_0_0_mgdl)

df$Total_cholesterol_0_0_mgdl <-  df$Total_cholesterol_0_0*38.67
summary(df$Total_cholesterol_0_0_mgdl)

df$frs_log_HDL_cholesterol_score <-  ifelse(is.na(df$HDL_cholesterol_0_0_mgdl), NA, log(df$HDL_cholesterol_0_0_mgdl))
df$frs_log_total_cholesterol_score <- ifelse(is.na(df$Total_cholesterol_0_0_mgdl), NA, log(df$Total_cholesterol_0_0_mgdl))

# BMI
df$frs_bmi_score <-  ifelse(df$BMI_0_0>0, df$BMI_0_0, NA)

# 2.6 Anti-hypertension medication 
df$antihyp <- ifelse(df$Antihypertensive_meds_0_0==1, 1,0)
summary(as.factor(df$antihyp))
df$frs_lgsbpt <- ifelse(df$antihyp==1, df$frs_log_sbp, 0)
df$frs_lgsbpu <- ifelse(df$antihyp==0, df$frs_log_sbp, 0)

#------ 3. Calculate sex-specific FRS score  ------
# Three different variants calculated:
# 1. beta-weighted frs with BMI 
# 2. beta-weighted frs without BMI 
# 3. point-based frs

# 3.1 BETA-WEIGHTED FRS (WITHOUT BMI)
# Following weights are from Agostino et al. (2008).
df$frs_frm_female_no_bmi <-            df$frs_log_age*2.32888 + df$frs_log_total_cholesterol_score*1.20904 + df$frs_log_HDL_cholesterol_score*-0.70833 + df$frs_lgsbpu*2.76157 +
                                       df$frs_lgsbpt*2.82263+ df$frs_smok_score*0.52873 + df$frs_diab_score*0.69154
df$frs_frm_female_no_bmi_avgrsk <-     exp(df$frs_frm_female_no_bmi - 26.1931) 
df$beta_frs_frm_female_no_bmi_final <-   1 -  0.95012^exp(df$frs_frm_female_no_bmi_avgrsk)

df$frs_frm_male_no_bmi <-              df$frs_log_age*3.06117 + df$frs_log_total_cholesterol_score*1.12370 + df$frs_log_HDL_cholesterol_score*-0.93263 + df$frs_lgsbpu*1.93303 +
                                       df$frs_lgsbpt*1.99881+ df$frs_smok_score*0.65451 + df$frs_diab_score*0.57367
df$frs_frm_male_no_bmi_avgrsk <-       exp(df$frs_frm_male_no_bmi - 23.9802) 
df$beta_frs_frm_male_no_bmi_final   <-  1 -  0.88936^exp(df$frs_frm_male_no_bmi_avgrsk)

summary(df$beta_frs_frm_female_no_bmi_final)
summary(df$beta_frs_frm_male_no_bmi_final)

df$FRS_predicted_prob <- ifelse(df$Sex==0, df$beta_frs_frm_female_no_bmi_final , ifelse(df$Sex==1, df$beta_frs_frm_male_no_bmi_final, NA))
summary(df$FRS_predicted_prob)

# 3.2 BETA-WEIGHTED FRS (WITH BMI)
# Weights are from Unnikrishnan et al. (2016). 
df$frs_frm_female <- df$frs_log_age*2.72107 + df$BMI_0_0*0.51125 + df$frs_lgsbpu*2.81291 + 
                     df$frs_lgsbpt*2.88267 + df$frs_smok_score*0.61868 + 
                     df$frs_diab_score*0.77763 - 26.1931

df$frs_female <-     1 - exp(df$frs_frm_female * log(0.94833))

df$frs_frm_male <-   df$frs_log_age*3.11296 + df$BMI_0_0*0.79299 + df$frs_lgsbpu*1.85505 + 
                     df$frs_lgsbpt*1.92672 + df$frs_smok_score*0.70953 + 
                     df$frs_diab_score*0.53160 - 23.9802

df$frs_male <-        1 - exp(df$frs_frm_male * log(0.88431))

df$FRS_BMI_pred_prob <-        ifelse(df$Sex==0, df$frs_female, ifelse(df$Sex==1, df$frs_male, NA))

summary(df$FRS_BMI_pred_prob)

# check missingness at item level
apply(df[,grep("frs|FRS", names(df))],2,pMiss)

# save file
save(file = paste0(data_pathway,"ukb_data_orig_merged_ANU-ADRI_CAIDE_FRS_oct22.rda"))
