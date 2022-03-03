# Written by Raihaan Patel
# Modified by Melis Anaturk (Feb 2020)

# This r script calculates two versions of the ANU-ADRI score:
# beta-weighted and points-weighted

#------ 1. Data set up -------
# 1.1 load packages
library(car)
library(sjstats)
library(pROC)
library(ggplot2)
library(dplyr)
library(psych)
library(tidyverse)
# 1.2 data pathway
data_pathway = "../../raw_data/"

# 1.3 read in csv files
load(file = paste0(data_pathway, "ukbdata_diagnoses_baseline_diseasestatus_baselinemedications.rda"))

# 1.4 create df of variable names
col_names <- data.frame(names(df))

#------ 2. Data recoding ------
# 2.1 APOE status 
class(df$rs429358)
summary(as.factor(df$rs429358))

df$APOE_genotype      <-  ifelse(df$rs429358=="C C" & df$rs7412=="T T", "E1/E1",
                          ifelse(df$rs429358=="C T" & df$rs7412=="T T", "E1/E2",
                          ifelse(df$rs429358=="C T" & df$rs7412=="C T", "E1/E3 or E2/E3",
                          ifelse(df$rs429358=="C C" & df$rs7412=="C T", "E1/E4",
                          ifelse(df$rs429358=="T T" & df$rs7412=="T T", "E2/E2",
                          ifelse(df$rs429358=="T T" & df$rs7412=="C T", "E2/E3",
                          ifelse(df$rs429358=="T T" & df$rs7412=="C C", "E3/E3",
                          ifelse(df$rs429358=="C T" & df$rs7412=="C C", "E3/E4",
                          ifelse(df$rs429358=="C C" & df$rs7412=="C C", "E4/E4", NA
                                                          )))))))))

df$APOE_e4 <-            ifelse(df$APOE_genotype=="E1/E2"|df$APOE_genotype=="E2/E2"|df$APOE_genotype=="E3/E3", "0",
                         ifelse(df$APOE_genotype=="E3/E4"|df$APOE_genotype=="E4/E4", "1",NA))


summary(as.factor(df$APOE_genotype))
#E1/E2  E2/E2  E3/E3  E3/E4  E4/E4   NA's 
#     3   2330 242097  97318   9791 149722 
summary(as.factor(df$APOE_e4))
#0      1   NA's 
#244430 107109 149722 

# 2.2 Age - sex specific scores
# N.B. Weights are from Antsey et al (2013).
# Sex coded as 0 == Women, 1 == Men

# 3.1 age weights for men
df$point_age_anu_adri_recoded <- ifelse(df$Sex==1 & df$Age_when_attended_assesment_centre_0_0<65, 0,
                                 ifelse(df$Sex==1 & df$Age_when_attended_assesment_centre_0_0>=65 & 
                                 df$Age_when_attended_assesment_centre_0_0<70, 1,
                                 ifelse(df$Sex==1 & df$Age_when_attended_assesment_centre_0_0>=70 &
                                 df$Age_when_attended_assesment_centre_0_0<75 , 12,    
                                 ifelse(df$Sex==1 & df$Age_when_attended_assesment_centre_0_0>=75 & 
                                 df$Age_when_attended_assesment_centre_0_0<80 , 18,  
                                 ifelse(df$Sex==1 & df$Age_when_attended_assesment_centre_0_0>=80 & 
                                 df$Age_when_attended_assesment_centre_0_0<85 , 26, 
                                 ifelse(df$Sex==1 & df$Age_when_attended_assesment_centre_0_0>=85 & 
                                 df$Age_when_attended_assesment_centre_0_0<=90 , 33, 
                                 ifelse(df$Sex==1 & df$Age_when_attended_assesment_centre_0_0>90, 38, 
                                 ifelse(df$Sex==0 & df$Age_when_attended_assesment_centre_0_0<65, 0,
                                 ifelse(df$Sex==0 & df$Age_when_attended_assesment_centre_0_0>=65 & 
                                 df$Age_when_attended_assesment_centre_0_0<70, 5,
                                 ifelse(df$Sex==0 & df$Age_when_attended_assesment_centre_0_0>=70 &
                                 df$Age_when_attended_assesment_centre_0_0<75 , 14,    
                                 ifelse(df$Sex==0 & df$Age_when_attended_assesment_centre_0_0>=75 & 
                                 df$Age_when_attended_assesment_centre_0_0<80 , 21,  
                                 ifelse(df$Sex==0 & df$Age_when_attended_assesment_centre_0_0>=80 & 
                                 df$Age_when_attended_assesment_centre_0_0<85 , 29, 
                                 ifelse(df$Sex==0 & df$Age_when_attended_assesment_centre_0_0>=85 & 
                                 df$Age_when_attended_assesment_centre_0_0<=90 , 35, 
                                 ifelse(df$Sex==0 & df$Age_when_attended_assesment_centre_0_0>90,41, NA))))))))))))))
# 2.3 age weights for beta-weighted 
df$beta_age_anu_adri_recoded <- ifelse(df$Sex==1 & df$Age_when_attended_assesment_centre_0_0<65, 0,
                                ifelse(df$Sex==1 & df$Age_when_attended_assesment_centre_0_0>=65 & 
                                df$Age_when_attended_assesment_centre_0_0<70, 0.13,
                                ifelse(df$Sex==1 & df$Age_when_attended_assesment_centre_0_0>=70 &
                                df$Age_when_attended_assesment_centre_0_0<75 , 1.57,    
                                ifelse(df$Sex==1 & df$Age_when_attended_assesment_centre_0_0>=75 & 
                                df$Age_when_attended_assesment_centre_0_0<80 , 2.04,  
                                ifelse(df$Sex==1 & df$Age_when_attended_assesment_centre_0_0>=80 & 
                                df$Age_when_attended_assesment_centre_0_0<85 , 3.37, 
                                ifelse(df$Sex==1 & df$Age_when_attended_assesment_centre_0_0>=85 & 
                                df$Age_when_attended_assesment_centre_0_0<=90 , 4.24, 
                                ifelse(df$Sex==1 & df$Age_when_attended_assesment_centre_0_0>90, 4.93,
                                ifelse(df$Sex==0 & df$Age_when_attended_assesment_centre_0_0<65, 0,
                                ifelse(df$Sex==0 & df$Age_when_attended_assesment_centre_0_0>=65 & 
                                df$Age_when_attended_assesment_centre_0_0<70, 0.64,
                                ifelse(df$Sex==0 & df$Age_when_attended_assesment_centre_0_0>=70 &
                                df$Age_when_attended_assesment_centre_0_0<75 , 1.87,    
                                ifelse(df$Sex==0 & df$Age_when_attended_assesment_centre_0_0>=75 & 
                                df$Age_when_attended_assesment_centre_0_0<80 , 2.75,  
                                ifelse(df$Sex==0 & df$Age_when_attended_assesment_centre_0_0>=80 & 
                                df$Age_when_attended_assesment_centre_0_0<85 , 3.71, 
                                ifelse(df$Sex==0 & df$Age_when_attended_assesment_centre_0_0>=85 & 
                                df$Age_when_attended_assesment_centre_0_0<=90 , 4.58, 
                                ifelse(df$Sex==0 & df$Age_when_attended_assesment_centre_0_0>90, 5.28,NA))))))))))))))

summary(df$point_age_anu_adri_recoded)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0000  0.0000  0.0000  0.6251  0.0000 14.0000
summary(df$beta_age_anu_adri_recoded)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00000 0.00000 0.00000 0.08047 0.00000 1.87000 

# quick check
myvars <- c("Sex", "Age_when_attended_assesment_centre_0_0", "point_age_anu_adri_recoded", "beta_age_anu_adri_recoded")
age_check <- df[myvars]
rm(age_check)

# 2.4 Education
summary(df$education_years)
#point weighted
df$point_edu_anu_adri_recoded <- ifelse(df$education_years>11, 0, 
                                 ifelse(df$education_years>=8 & 
                                 df$education_years<=11, 3,
                                 ifelse(df$education_years<8, 6, NA)))
#beta weighted
df$beta_edu_anu_adri_recoded <- ifelse(df$education_years>11, 0, 
                                ifelse(df$education_years>=8 & 
                                df$education_years<=11, 0.42,
                                ifelse(df$education_years<8, 0.80, NA)))

summary(df$point_edu_anu_adri_recoded)
summary(df$beta_edu_anu_adri_recoded)

# quick check
myvars <- c("education_years", "point_edu_anu_adri_recoded", "beta_edu_anu_adri_recoded")
edu_check <- df[myvars]
rm(edu_check)

# 2.5 BMI
df$point_BMI_anu_adri_recoded <-ifelse(df$Age_when_attended_assesment_centre_0_0<60 & df$BMI_0_0>=18.5 & df$BMI_0_0<=24.9, 0, 
                                ifelse(df$Age_when_attended_assesment_centre_0_0<60 & df$BMI_0_0>=25 & df$BMI_0_0<=29.9, 2,
                                ifelse(df$Age_when_attended_assesment_centre_0_0<60 & df$BMI_0_0>=30, 5, NA)))
df$beta_BMI_anu_adri_recoded <- ifelse(df$Age_when_attended_assesment_centre_0_0<60 & df$BMI_0_0>=18.5 & df$BMI_0_0<=24.9, 0, 
                                ifelse(df$Age_when_attended_assesment_centre_0_0<60 & df$BMI_0_0>=25 & df$BMI_0_0<=29.9, 0.3,
                                ifelse(df$Age_when_attended_assesment_centre_0_0<60 & df$BMI_0_0>=30, 0.71, NA)))

# individuals =>60 coded as 0 - to retain in analysis - check with sana 
df$point_BMI_anu_adri_recoded[df$Age_when_attended_assesment_centre_0_0>=60] <- 0
df$beta_BMI_anu_adri_recoded[df$Age_when_attended_assesment_centre_0_0>=60] <- 0

summary(df$point_BMI_anu_adri_recoded)
summary(df$beta_BMI_anu_adri_recoded)

#check
myvars <- c("Age_when_attended_assesment_centre_0_0", "BMI_0_0", "beta_BMI_anu_adri_recoded","point_BMI_anu_adri_recoded")
bmi_check <- df[myvars]
rm(bmi_check)

# 2.6 Diabetes
df$point_diabetes_anu_adri_recoded <- ifelse(df$Diabetes_BIN_FINAL_0_0==1, 3, 
                                      ifelse(df$Diabetes_BIN_FINAL_0_0==0, 0, NA))
df$beta_diabetes_anu_adri_recoded  <- ifelse(df$Diabetes_BIN_FINAL_0_0==1, 0.33, 
                                      ifelse(df$Diabetes_BIN_FINAL_0_0==0, 0, NA))
summary(as.factor(df$point_diabetes_anu_adri_recoded))
summary(as.factor(df$beta_diabetes_anu_adri_recoded))

# 2.7 Current depression
#df$depression_BIN[is.na(df$depression_BIN)] <- 0
#summary(as.factor(df$depression_BIN))
summary(as.factor(df$depression_BIN_FINAL_0_0))

df$point_depression_anu_adri_score <- ifelse(df$depression_BIN_FINAL_0_0==1, 2, 0)
df$beta_depression_anu_adri_score  <-  ifelse(df$depression_BIN_FINAL_0_0==1, 0.29, 0)


# 2.8 Cholesterol 
df$point_cholesterol_anu_adri_recoded <- ifelse(df$Age_when_attended_assesment_centre_0_0<60 & df$Total_cholesterol_0_0<6.2, 0, 
                                         ifelse(df$Age_when_attended_assesment_centre_0_0<60 & df$Total_cholesterol_0_0>=6.2, 3, NA))
df$beta_cholesterol_anu_adri_recoded <-  ifelse(df$Age_when_attended_assesment_centre_0_0<60 & df$Total_cholesterol_0_0<6.2, 0, 
                                         ifelse(df$Age_when_attended_assesment_centre_0_0<60 & df$Total_cholesterol_0_0>=6.2, 0.41, NA))

# individuals older than 60 coded as 0 - to retain in analysis 
df$point_cholesterol_anu_adri_recoded[df$Age_when_attended_assesment_centre_0_0>=60] <- 0
df$beta_cholesterol_anu_adri_recoded[df$Age_when_attended_assesment_centre_0_0>=60] <- 0

summary(df$point_BMI_anu_adri_recoded)
summary(df$beta_BMI_anu_adri_recoded)

# 2.9 TBI 
df$point_tbi_anu_adri_score <- ifelse(df$TBI_BIN_FINAL_0_0==1, 4, 0)
df$beta_tbi_anu_adri_score  <- ifelse(df$TBI_BIN_FINAL_0_0==1, 0.46, 0)
summary(as.factor(df$point_tbi_anu_adri_score))
summary(as.factor(df$beta_tbi_anu_adri_score))

# 2.10 SMOKING
df$point_smok_anu_adri_score <- ifelse(df$Smoking_Status_0_0==2, 0.58, ifelse(df$Smoking_Status_0_0==1, 0.19,
                                                              ifelse(df$Smoking_Status_0_0==0,0,NA)))
df$beta_smok_anu_adri_score <- ifelse(df$Smoking_Status_0_0==2, 4,     ifelse(df$Smoking_Status_0_0==1, 1,
                                                              ifelse(df$Smoking_Status_0_0==0,0,NA)))
summary(as.factor(df$point_smok_anu_adri_score))
summary(as.factor(df$beta_smok_anu_adri_score))

# 2.11 Alcohol
df_alc <- read.csv(file = paste0(data_pathway, "alcohol_oct22_clean.csv"), header=TRUE, sep=",", stringsAsFactors = FALSE)
cols_alc <- read.csv(file = paste0(data_pathway, "alc_new_col_names.csv"), header=TRUE, sep=",", stringsAsFactors = FALSE)
names(df_alc) <- as.vector(cols_alc[,1])

# Below code credit is from Anya Topiwala
names(df_alc)[names(df_alc) == "1568-0.0"] <- "rw" # wkly rw
names(df_alc)[names(df_alc) == "1578-0.0"] <- "ww"
names(df_alc)[names(df_alc) == "1588-0.0"] <- "pint"
names(df_alc)[names(df_alc) == "1598-0.0"] <- "sp"
names(df_alc)[names(df_alc) == "5364-0.0"] <- "other"
names(df_alc)[names(df_alc) == "1608-0.0"] <- "fw"
names(df_alc)[names(df_alc) == "4407-0.0"] <- "rwm" #mthly rw
names(df_alc)[names(df_alc) == "4418-0.0"] <- "wwm"
names(df_alc)[names(df_alc) == "4429-0.0"] <- "pintm"
names(df_alc)[names(df_alc) == "4440-0.0"] <- "spm"
names(df_alc)[names(df_alc) == "4462-0.0"] <- "otherm"
names(df_alc)[names(df_alc) == "4451-0.0"] <- "fwm"
names(df_alc)[names(df_alc) == "20117-0.0"] <- "alcohol_status"
names(df_alc)[names(df_alc) == "20406-0.0"] <- "ever_alc_addict"

# calculate units
df_alc$weekly_rw <- df_alc$rw*1.69 # convert to weekly rw
df_alc$weekly_ww <- df_alc$ww*1.69
df_alc$weekly_pint <- df_alc$pint*2.41
df_alc$weekly_sp <- df_alc$sp*1
df_alc$weekly_other <- df_alc$other*1.1
df_alc$weekly_fw <- df_alc$fw*1.19
df_alc$units_week <- rowSums(df_alc[,grep("weekly", names(df_alc))], na.rm=T) # sum weekly units

df_alc$monthly_rw <- df_alc$rwm*1.69
df_alc$monthly_ww <- df_alc$wwm*1.69
df_alc$monthly_pint <-df_alc$pintm*2.41
df_alc$monthly_sp <- df_alc$spm*1
df_alc$monthly_other <- df_alc$otherm*1.1
df_alc$monthly_units_fw <- df_alc$fwm*1.19 

df_alc$units_month <- rowSums(df_alc[,grep("monthly", names(df_alc))], na.rm=T)
df_alc$units_weekfrommonth <- df_alc$units_month/4.33 # convert the monthly units into something weekly
df_alc$units_combined <- rowSums(df_alc[,c("units_week","units_weekfrommonth")],na.rm=T) # add this to the weekly total to make a units per week grand total 

summary(df_alc$units_combined)
hist(df_alc$units_combined)
df_alc$eid <-as.character(df_alc$eid)

# merge with main df
df <- list(df, df_alc) %>% reduce(left_join, by = "eid")
#to check - do i need to merge all variables of df_alc?
summary(df$units_combined)

# define light drinking (women)
df$alcohol_status_0_0 <- ifelse(df$units_combined==0,0,
                         ifelse(df$Sex==0&df$units_combined>=0&df$units_combined<14,1,
                         ifelse(df$Sex==0&df$units_combined>=14,0,
                         ifelse(df$Sex==1&df$units_combined>=0&df$units_combined<21,1,
                         ifelse(df$Sex==1&df$units_combined>=21,0,NA)))))

View(df[ c("alcohol_status_0_0", "Sex", "units_combined")])

summary(df$alcohol_status_0_0)

# calculate points based alcohol variable
df$point_alcohol_anu_adri_score <- ifelse(df$alcohol_status_0_0==0, 0,
                                   ifelse(df$alcohol_status_0_0==1,-3,NA))
# calculate regression based alcohol variable 
df$beta_alcohol_anu_adri_score <- ifelse(df$alcohol_status_0_0==0, 0,
                         ifelse(df$alcohol_status_0_0<=1,-0.33,NA))

summary(as.factor(df$point_alcohol_anu_adri_score))
summary(as.factor(df$beta_alcohol_anu_adri_score))

# 2.12 SOCIAL ENGAGEMENT
# Create a composite SE score ranging from 0 - 3
# Here we assign a 1 if other individuals live in household, if family and friend visit and if individual participates in group activities on weekly basis 
# still need percieved quality of family/friend relationships
#df[ , 270:274]  <- lapply(df[ , 270:274] , 
#                               FUN = function(x) car::recode(x, "-3=NA;-7=NA"))
#df[ , 578:585]  <- lapply(df[ , 578:585] , 
#                          FUN = function(x) car::recode(x, "-1=NA;-7=NA"))

df$Social_engagement_0_1 <-      ifelse(df$Number_in_household_0_0>1, 1, 
                                 ifelse(df$Number_in_household_0_0<=1, 0, NA))

df$Social_engagement_0_2 <-      ifelse(df$Frequency_of_friendfam_visits_0_0>=1 & df$Frequency_of_friendfam_visits_0_0<=3, 1,
                                 ifelse(df$Frequency_of_friendfam_visits_0_0>=4 | df$Frequency_of_friendfam_visits_0_0<=7, 0, NA))

df$Social_engagement_0_3 <-      ifelse(df$Leisure_social_activities_0_0>=1 | df$Leisure_social_activities_0_1>=1 |
                                 df$Leisure_social_activities_0_2>=1 | df$Leisure_social_activities_0_3>=1 |    
                                 df$Leisure_social_activities_0_4>=1, 1,0)

df$Social_engagement_0_3[is.na(df$Leisure_social_activities_0_0)&is.na(df$Leisure_social_activities_0_1)&is.na(df$Leisure_social_activities_0_2)& is.na(df$Leisure_social_activities_0_3)&is.na(df$Leisure_social_activities_0_4)] <- 0

df$Social_engagement_0_4 <-      ifelse(df$Family_relationship_satisfaction_0_0>=1 & df$Family_relationship_satisfaction_0_0<=3|  
                                 df$Friendship_satisfaction_0_0>=1 & df$Friendship_satisfaction_0_0<=3, 1,
                                 ifelse(df$Family_relationship_satisfaction_0_0>=4 | df$Friendship_satisfaction_0_0>=4, 0, NA))


# quick check of variables
myvars <- c("Social_engagement_0_3", "Leisure_social_activities_0_0", "Leisure_social_activities_0_1", "Leisure_social_activities_0_2", "Leisure_social_activities_0_3", "Leisure_social_activities_0_4")
myvars <- c("Number_in_household_0_0", "Social_engagement_0_1", "Frequency_of_friendfam_visits_0_0", "Social_engagement_0_2", 
            "Leisure_social_activities_0_1", "Leisure_social_activities_0_2", "Leisure_social_activities_0_3", "Leisure_social_activities_0_4",
            "Social_engagement_0_3")
SE_check <- df[myvars]
rm(SE_check)

# create SE composite score
df$SE_TOTAL <- rowSums(df[,c("Social_engagement_0_1","Social_engagement_0_2", "Social_engagement_0_3", "Social_engagement_0_4")], 
                       na.rm=TRUE) 
summary(as.factor(df$SE_TOTAL))
#0      1      2      3      4 
#6183  48012 158216 214978  73872 

# current coding: 0 = lowest, 1 = low-med, 2 = med-high, >3 = highest
df$point_SE_anu_adri_recoded <- ifelse(df$SE_TOTAL==0, 6, 
                                ifelse(df$SE_TOTAL==1, 4,
                                ifelse(df$SE_TOTAL==2, 1, 
                                ifelse(df$SE_TOTAL>=3, 0, NA))))

df$beta_SE_anu_adri_recoded <- ifelse(df$SE_TOTAL==0, 0.84, 
                               ifelse(df$SE_TOTAL==1, 0.51,
                               ifelse(df$SE_TOTAL==2,0.17,
                               ifelse(df$SE_TOTAL>=3, 0, NA))))

summary(as.factor(df$point_SE_anu_adri_recoded))
summary(as.factor(df$beta_SE_anu_adri_recoded))

apply(df[,grep("Social_engagement", names(df))],2,pMiss)

# 2.13 PHYSICAL ACTIVITY
#To align with IPA questionnaire
##### high PA: vigorous activity greater than or equal to 15000 MET min/week or low/mod/vig activity greater than or equal to 3000 MET min/week
#### moderate PA: 3 or more days of vigorous intensity activity and/or walking of at least 30 min per day
#### OR 5 or more days of moderate intensity activity and/or walking for at least 30 min per day
###  OR 5 or more days of any combination of low/mod/vig to equate to 600 MET min/week
###  low PA: Not meeting the above
#df$Number_days_allPA_perweek_0_0 <-  rowSums(df[,c("Number_days_vigorousPA_perweek_0_0","Number_days_moderatePA_perweek_0_0", "Number_days_walked_perweek_0_0")], 
#                                             na.rm=TRUE)
#summary(df$Number_days_allPA_perweek_0_0)

df$point_pa_score_anu_adri_recoded <- NULL
df$point_pa_score_anu_adri_recoded <- factor(df$IPAQ_activity_group_0_0, levels=0:2, labels=c(0, -2, -3))

df$beta_pa_score_anu_adri_recoded <- NULL
df$beta_pa_score_anu_adri_recoded <- factor(df$IPAQ_activity_group_0_0, levels=0:2, labels=c(0, -0.29, -0.40))

summary(as.factor(df$point_pa_score_anu_adri_recoded))
summary(as.factor(df$beta_pa_score_anu_adri_recoded))

# 2.14 Fish intake
df$total_fish_intake_per_week_0_0 <- rowSums(df[,c("Fish_oily_intake_0_0","Fish_nonoily_intake_0_0")], 
                                             na.rm=TRUE)
# quick check
myvars <- c("total_fish_intake_per_week_0_0","Fish_oily_intake_0_0", "Fish_nonoily_intake_0_0")
fish_data <- df[myvars]
rm(fish_data)

df$point_fish_score_anu_adri_recoded <- ifelse(df$total_fish_intake_per_week_0_0==0, 0, 
                                        ifelse(df$total_fish_intake_per_week_0_0==1, -3, 
                                        ifelse(df$total_fish_intake_per_week_0_0>=2 & df$total_fish_intake_per_week_0_0<=3, -4,
                                        ifelse(df$total_fish_intake_per_week_0_0>=4, -5, NA))))

df$beta_fish_score_anu_adri_recoded <- ifelse(df$total_fish_intake_per_week_0_0==0, 0, 
                                       ifelse(df$total_fish_intake_per_week_0_0==1, -0.33, 
                                       ifelse(df$total_fish_intake_per_week_0_0>=2 & df$total_fish_intake_per_week_0_0<=3, -0.53,
                                       ifelse(df$total_fish_intake_per_week_0_0>=4, -0.62, NA))))
summary(df$total_fish_intake_per_week_0_0)
summary(df$point_fish_score_anu_adri_recoded)
summary(df$beta_fish_score_anu_adri_recoded)

# 2.15 Pesticide exposure
# recode pesticide exposure - currently 39 variables 
#first recode indivdiual variables and then make summary variable
df[ , 488:527]  <- lapply(df[ , 488:527] , FUN = function(x) car::recode(x, "NA=0;-141=2;-131=1;-212=NA"))
pesticide_vars <- df[ , 488:527] 
df$Worked_with_pesticides_total_0_0 <- apply(pesticide_vars[, -1], 1, function(x) {max(as.numeric(as.character(x)), na.rm = T)}) 

summary(as.factor(df$Worked_with_pesticides_total_0_0))
#0      1      2 
#498000   2847    414 

df$point_pesticide_anu_adri_recoded <- ifelse(df$Worked_with_pesticides_total_0_0==1, 2, 
                                       ifelse(df$Worked_with_pesticides_total_0_0==2, 2,
                                       ifelse(df$Worked_with_pesticides_total_0_0==0, 0, 0)))

df$beta_pesticide_anu_adri_recoded <- ifelse(df$Worked_with_pesticides_total_0_0==1, 0.31, 
                                      ifelse(df$Worked_with_pesticides_total_0_0==2, 0.31,
                                      ifelse(df$Worked_with_pesticides_total_0_0==0, 0, 0)))
summary(as.factor(df$point_pesticide_anu_adri_recode))
summary(as.factor(df$beta_pesticide_anu_adri_recoded))

#remove redundant dfs
rm(pesticide_vars)

#------ 3. Calculate ANU-ADRI scores ------
#point-based score

df$point_TOTAL_ANU_ADRI_SCORE <- (df$point_age_anu_adri_recoded + df$point_edu_anu_adri_recoded + df$point_BMI_anu_adri_recoded + 
                                  df$point_diabetes_anu_adri_recoded + df$point_depression_anu_adri_score + df$point_cholesterol_anu_adri_recoded + 
                                  df$point_tbi_anu_adri_score + df$point_smok_anu_adri_score + df$point_alcohol_anu_adri_score + df$point_SE_anu_adri_recoded + 
                                  df$point_fish_score_anu_adri_recoded + df$point_pesticide_anu_adri_recoded)

# excluded domain: df$point_cog_score_anu_adri_recoded 
summary(df$point_TOTAL_ANU_ADRI_SCORE)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  -8.00   -3.81   -1.00   -0.76    2.00   24.58   33919 
hist(df$point_TOTAL_ANU_ADRI_SCORE)

#beta-weighted score
df$ANU_ADRI <- df$beta_age_anu_adri_recoded + df$beta_edu_anu_adri_recoded + df$beta_BMI_anu_adri_recoded + 
               df$beta_diabetes_anu_adri_recoded + df$beta_depression_anu_adri_score + df$beta_cholesterol_anu_adri_recoded + 
               df$beta_tbi_anu_adri_score + df$beta_smok_anu_adri_score + df$beta_alcohol_anu_adri_score + df$beta_SE_anu_adri_recoded + 
               df$beta_fish_score_anu_adri_recoded + df$beta_pesticide_anu_adri_recoded

# excluded domain: df$beta_cog_score_anu_adri_recoded 
summary(df$ANU_ADRI)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  -0.95   -0.20    0.44    0.74    1.10    7.13   33919
hist(df$ANU_ADRI)

# check missingness
pMiss <- function(x){sum(is.na(x))/length(x)*100}
apply(df[,grep("ANU|anu_adri", names(df))],2,pMiss)

#save current df
save(df, file = paste0(data_pathway, "ukbdata_diagnoses_baseline_diseasestatus_baselinemedications_ANUADRI.rda"))
