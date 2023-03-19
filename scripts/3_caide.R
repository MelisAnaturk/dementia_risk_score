# Created by Raihaan Patel
# Modified by Melis Anaturk (Feb 2020)
# This r script calculates two versions of the CAIDE score, with APOE and without

library(car)
library(sjstats)
library(pROC)
library(ggplot2)
library(dplyr)
library(psych)
library(tidyverse)

data_pathway = "../../raw_data/"

load(file = paste0(data_pathway, "ukbdata_diagnoses_baseline_diseasestatus_baselinemedications_ANUADRI.rda"))

# STEPS IN SCRIPT
# 1. Data recoding
# 2. Caide caculation

#------ 1. Data recoding --------
# 1.1 Age
df$beta_age_caide_recoded  <-  ifelse(df$Age_when_attended_assesment_centre_0_0<47, 0, 
                               ifelse(df$Age_when_attended_assesment_centre_0_0>=47 & 
                               df$Age_when_attended_assesment_centre_0_0<=53,  1.084,
                               ifelse(df$Age_when_attended_assesment_centre_0_0>53, 1.762,NA)))

summary(df$beta_age_caide_recoded)

df$beta_age_caide_APOE_recoded  <-  ifelse(df$Age_when_attended_assesment_centre_0_0<47, 0, 
                                    ifelse(df$Age_when_attended_assesment_centre_0_0>=47 & 
                                    df$Age_when_attended_assesment_centre_0_0<=53,  1.155,
                                    ifelse(df$Age_when_attended_assesment_centre_0_0>53, 1.874,NA)))

# 1.2 Sex
#Female coded as 0 for score, if male, 1
df$beta_sex_caide_recoded  <-       ifelse(df$Sex==0, 0, ifelse(df$Sex==1, 0.47, NA))
df$beta_sex_caide_APOE_recoded  <-  ifelse(df$Sex==0, 0, ifelse(df$Sex==1, 0.438, NA))

summary(df$beta_sex_caide_recoded)

# 1.3 BMI
df$beta_bmi_caide_recoded  <-  ifelse(df$BMI_0_0<=30, 0, ifelse(df$BMI_0_0>30, 0.831, NA))
df$beta_bmi_caide_APOE_recoded  <-  ifelse(df$BMI_0_0<=30, 0, ifelse(df$BMI_0_0>30, 0.608, NA))

summary(df$beta_bmi_caide_recoded)

# 1.4 Cholesterol 
df$beta_cholesterol_caide_recoded  <-  ifelse(df$Total_cholesterol_0_0<=6.5, 0, ifelse(df$Total_cholesterol_0_0>6.5, 0.631, NA))
df$beta_cholesterol_caide_APOE_recoded  <-  ifelse(df$Total_cholesterol_0_0<=6.5, 0, ifelse(df$Total_cholesterol_0_0>6.5, 0.460, NA))

summary(df$beta_cholesterol_caide_recoded)

# 1.5 Systolic BP
# first calculate average from two measures of BP
df$Systolic_BP_auto_mean <- (df$Systolic_BP_auto_0_0 + df$Systolic_BP_auto_0_1)/2
### if average missing use one of the two assessments
df$beta_sbp_caide_recoded  <-  ifelse(df$Systolic_BP_auto_mean<=140, 0, ifelse(df$Systolic_BP_auto_mean>140, 0.791 ,NA))
df$beta_sbp_caide_APOE_recoded  <-  ifelse(df$Systolic_BP_auto_mean<=140, 0, ifelse(df$Systolic_BP_auto_mean>140, 0.817 ,NA))

summary(df$beta_sbp_caide_recoded)

# 1.6 Physical activity
#Take sum of mod and vig activity per day -> minutes exercising perday
#Multiply by 7 to get wkly total, threshold at <40 (~ 2 days of 20 minute sessions)
#getting weekly total to account for ppl who input eg 15mins/day to represent 1 session of 30 mins every other day
df$Duration_of_moderate_and_vigorous_PA_0_0 <- rowSums(df[,c("Duration_of_moderatePA_0_0","Duration_of_vigorousPA_0_0")], 
                                                       na.rm=TRUE)
df$Duration_of_moderate_and_vigorous_PA_perweek_0_0 <- df$Duration_of_moderate_and_vigorous_PA_0_0*7
summary(df$Duration_of_moderate_and_vigorous_PA_0_0)

df$beta_pa_caide_recoded   <- ifelse(df$Duration_of_moderate_and_vigorous_PA_perweek_0_0<40, 0.527, 
                              ifelse(df$Duration_of_moderate_and_vigorous_PA_perweek_0_0>=40, 0, NA))

df$beta_pa_caide_APOE_recoded   <- ifelse(df$Duration_of_moderate_and_vigorous_PA_perweek_0_0<40, 0.057, 
                                   ifelse(df$Duration_of_moderate_and_vigorous_PA_perweek_0_0>=40, 0, NA))

summary(df$beta_pa_caide_recoded)

# 1.7 Education
#based on a computed years of education 0-6: 3, 7-9: ?
df$beta_edu_caide_recoded  <- ifelse(df$education_years<=6, 1.281,
                              ifelse(df$education_years>=7&df$education_years <=9, 0.910,
                              ifelse(df$education_years>=10, 0, NA)))

df$beta_edu_caide_APOE_recoded  <- ifelse(df$education_years<=6, 1.587,
                              ifelse(df$education_years>=7&df$education_years <=9, 1.149,
                              ifelse(df$education_years>=10, 0, NA)))
                              
summary(df$beta_edu_caide_recoded)

# 1.8 APOE status
df$beta_APOE_caide_recoded  <- ifelse(df$APOE_e4==1, 0.890,
                               ifelse(df$APOE_e4==0,0,NA))


#------ 2. Calculate CAIDE score ----
#DONT do row sums, because if someone is missing a score, you penalize them because they are missing data
#as opposed to them actually having high age or high bmi etc

# 2.1 beta-weighted caide score
df$beta_caide_score <- df$beta_age_caide_recoded + df$beta_sex_caide_recoded + df$beta_bmi_caide_recoded + 
                       df$beta_cholesterol_caide_recoded + df$beta_sbp_caide_recoded + df$beta_pa_caide_recoded + 
                       df$beta_edu_caide_recoded

# 2.2 beta-weighted caide score - WITH APOE
df$beta_caide_APOE_score <- df$beta_age_caide_APOE_recoded + df$beta_sex_caide_APOE_recoded + df$beta_bmi_caide_APOE_recoded + 
                            df$beta_cholesterol_caide_APOE_recoded + df$beta_sbp_caide_APOE_recoded + df$beta_pa_caide_APOE_recoded + 
                            df$beta_edu_caide_APOE_recoded + df$beta_APOE_caide_recoded

hist(df$beta_caide_APOE_score)
summary(df$beta_caide_APOE_score)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   0.00    1.87    2.69    2.52    3.21    6.67  200152 
   
# 2.3 check missingness
pMiss <- function(x){sum(is.na(x))/length(x)*100}

apply(df[,grep("CAIDE|caide", names(df))],2,pMiss)

# 2.4 Calculate Probability(dementia) - CAIDE+APOE
df$CAIDE_APOE_predicted_prob              <- (exp(-8.083 +    1.020 + (0.390*df$beta_caide_APOE_score)))/
                                           (1+(exp(-8.083 + 1.020 + (0.390*df$beta_caide_APOE_score))))
summary(df$CAIDE_APOE_predicted_prob*100)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   0.09    0.18    0.24    0.25    0.30    1.14  200152 
   
# 2.5 Calculate Probability(dementia) - CAIDE without APOE
df$CAIDE_predicted_prob <- (exp(-7.406 + 0.796 + (0.401*df$beta_caide_score)))/
                           (1+(exp(-7.406 + 0.796 + (0.401*df$beta_caide_score))))

summary(df$CAIDE_predicted_prob*100)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   0.13    0.27    0.35    0.38    0.46    1.65   80215 

save(df, file = paste0(data_pathway, "ukbdata_diagnoses_baseline_diseasestatus_baselinemedications_ANUADRI_CAIDE.rda"))
