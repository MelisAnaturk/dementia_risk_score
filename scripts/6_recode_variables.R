# Written by Melis Anaturk (July 2020)

# STEPS IN SCRIPT
# 1. Data set-up
# 2. Recode variables for main analysis

#------------- 1.  Data set-up ----
# 1.1 load libraries
library(data.table)

# 1.2 data pathway
data_pathway = "~/Documents/Oxford_DPhil/Biobank_Analyses/Dementia_risk_project/ukb_data/For_Rai/raw_data/"

# 1.3 set working directory 
setwd(data_pathway)

# 1.4 load csv file
load(file = paste0(data_pathway,"ukb_data_orig_merged_ANU-ADRI_CAIDE_FRS_oct22.rda"))

# 1.5 load in townsend deprivation variables
df_townsend_deprivation <- fread(paste0(data_pathway, "town_send_deprivation.csv"), header = TRUE, sep =',', stringsAsFactors = FALSE)
df_townsend_deprivation$eid <- as.character(df_townsend_deprivation$eid)

# merge dataframes 
df <- list(df,df_townsend_deprivation) %>% reduce(left_join, by = "eid")

# updated df of column names
all_col_names <- data.frame(names(df))

#------------- 2. Recode variables ----
# 2.1 Sex
summary(df$Sex)
df$Sex <- as.factor(df$Sex)

# 2.2 BMI
summary(df$BMI_0_0)
df$BMI_bin <-                ifelse(df$BMI_0_0<30, 0,
                             ifelse(df$BMI_0_0>=30, 1, NA))
df$BMI_bin <- as.factor(df$BMI_bin)
summary(df$BMI_bin)

# 2.3 Education
summary(df$education_years)

#df$education_stratified <- NULL
df$education_stratified <- ifelse(df$education_years<12,  0,
                           ifelse(df$education_years>=12, 1,NA))
df$education_stratified <- as.factor(df$education_stratified)

# 2.4 Alcohol intake
df$Alcohol_alt  <-        ifelse(df$Alcohol_status_0_0==0, 0,
                          ifelse(df$Alcohol_status_0_0==1, 1,
                          ifelse(df$Alcohol_status_0_0==2, 1, NA)))

df$Alcohol_alt<- as.factor(df$Alcohol_alt)

# based on codings of meta-analysis

# 2.5 Physical activity (IPAQ group)
summary(df$IPAQ_activity_group_0_0)
df$IPAQ_activity_group_0_0<- as.factor(df$IPAQ_activity_group_0_0)

# recode PA groups
df$IPAQ_activity_group_0_0 <- factor(df$IPAQ_activity_group_0_0, levels=0:2, labels=c("0", "1", "1"))
summary(df$IPAQ_activity_group_0_0)

# 2.6 Hearing problems
df$Hearing_prob <-      ifelse(df$Hearing_problems_0_0==1, 1, 
                        ifelse(df$Hearing_problems_0_0==0, 0,
                        ifelse(df$Hearing_problems_0_0==99,1, NA)))

df$Hearing_prob<- as.factor(df$Hearing_prob)

# 2.7 APOE e4
df$APOE_genotype_bin  <-      ifelse(df$APOE_e4==0, 0, 
                              ifelse(df$APOE_e4==1, 1,NA))
df$APOE_genotype_bin<- as.factor(df$APOE_genotype_bin)
summary(df$APOE_genotype_bin)

# 2.8 Stroke/TIA
summary(df$self_report_stroke)
summary(df$stroke_BIN_TOTAL)

df$stroke_TIA_BIN_FINAL <- ifelse(df$stroke_BIN_FINAL_0_0==1|df$TIA_BIN_FINAL_0_0==1,1,0)

# 2.9 Fish intake
df$Fish_intake_BIN <- ifelse(df$total_fish_intake_per_week_0_0==0, 0,
                      ifelse(df$total_fish_intake_per_week_0_0>=1, 1, NA))
df$Fish_intake_BIN <-as.factor(df$Fish_intake_BIN)

# 2.10 HDL Cholesterol
#https://www.heartuk.org.uk/cholesterol/getting-a-cholesterol-test
df$HDL_CHOLESTEROL_BIN <- ifelse(df$Sex==1&df$HDL_cholesterol_0_0<=1, 0, 
                          ifelse(df$Sex==1&df$HDL_cholesterol_0_0>1,1,
                          ifelse(df$Sex==0&df$HDL_cholesterol_0_0<=1.2, 0,
                          ifelse(df$Sex==0&df$HDL_cholesterol_0_0>=1.2, 0,NA))))
df$HDL_CHOLESTEROL_BIN <- as.factor(df$HDL_CHOLESTEROL_BIN)


# 2.11 Total Cholesterol
df$TOTAL_CHOLESTEROL_BIN <- ifelse(df$Total_cholesterol_0_0<=6.5, 0, 
                            ifelse(df$Total_cholesterol_0_0>6.5, 1, NA))
df$TOTAL_CHOLESTEROL_BIN <- as.factor(df$TOTAL_CHOLESTEROL_BIN)

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5240556/

# 2.12 Systolic BP
df$Systolic_BP_bin <-  ifelse(df$Systolic_BP_auto_mean<=140, 0, 
                       ifelse(df$Systolic_BP_auto_mean>140, 1 ,NA))
df$Systolic_BP_bin <- as.factor(df$Systolic_BP_bin)

# 2.13 Sleep
df$Sleep_BIN <-        ifelse(df$Sleep_duration_0_0<7, 0,
                       ifelse(df$Sleep_duration_0_0>=7 & df$Sleep_duration_0_0<=9, 1, 
                       ifelse(df$Sleep_duration_0_0>9, 2,NA)))
df$Sleep_BIN <- as.factor(df$Sleep_BIN)

#https://jcsm.aasm.org/doi/10.5664/jcsm.4758#d3e471

summary(df$Sleeplesness_insomnia_0_0)
df$Sleeplesness_insomnia_0_0_bin <-   ifelse(df$Sleeplesness_insomnia_0_0==3, 2,
                                      ifelse(df$Sleeplesness_insomnia_0_0==2, 1, 
                                      ifelse(df$Sleeplesness_insomnia_0_0==1, 0, NA)))
df$Sleeplesness_insomnia_0_0_bin<-as.factor(df$Sleeplesness_insomnia_0_0_bin)

# 2.14 Smoker status
df$Smoker_bin <-          ifelse(df$Smoking_Status_0_0==0, 0,
                          ifelse(df$Smoking_Status_0_0==1, 1,
                          ifelse(df$Smoking_Status_0_0==2, 1, NA)))
class(df$Smoker_bin)

df$Smoker_bin<-as.factor(df$Smoker_bin)

# 2.15 Social engagement
df$Social_engagement_0_1 <- as.factor(df$Social_engagement_0_1)
df$Social_engagement_0_2 <- as.factor(df$Social_engagement_0_2)
df$Social_engagement_0_3 <- as.factor(df$Social_engagement_0_3)

# 2.16 Hormone replacement therapy
#df$HRT_0_0<- ifelse(df$primary_care_prescription_for_HRTs==1|df$HRT_meds_0_0==1,1,0)
#df$HRT_0_0[is.na(df$HRT_0_0)] <- 0
df$HRT_0_0 <- df$HRT_meds_0_0
class(df$HRT_0_0)
summary(df$HRT_0_0)

# 2.17 Diastolic BP
#df$Diastolic_BP_auto_mean <- (df$Diastolic_BP_auto_0_0 + df$Diastolic_BP_auto_0_1)/2
#df$Diastolic_BP_bin <-        ifelse(df$Diastolic_BP_auto_mean<=90, 0, 
#                              ifelse(df$Diastolic_BP_auto_mean>90, 1 ,NA))
#df$Diastolic_BP_bin <- as.factor(df$Diastolic_BP_bin)

# read in LDL data
df_LDL <- read.csv(paste0(data_pathway, "LDL.csv"), header=TRUE, sep=",", stringsAsFactors = FALSE)
df_LDL[df_LDL == "NA"] <- NA
df_LDL$eid <- as.character(df_LDL$eid)

# merge LDL with main df
df <- list(df,df_LDL) %>% reduce(left_join, by = "eid")
rm(df_LDL)

df$LDL_0_0_BIN <- ifelse(df$LDL_0_0<3, 0, 
                  ifelse(df$LDL_0_0>=3, 1, NA))
df$LDL_0_0_BIN <- as.factor(df$LDL_0_0_BIN)

# 2.18 arterial stiffness
#df$Pulse_wave_art_stiff_0_0 <- ifelse(df$Pulse_wave_arterial_stiffness_0_0<=10, 0, 
#                              ifelse(df$Pulse_wave_arterial_stiffness_0_0>10, 1, NA))
#df$Pulse_wave_art_stiff_0_0<- as.factor(df$Pulse_wave_art_stiff_0_0)
#summary(df$Pulse_wave_art_stiff_0_0)

# 2.19 parental history of illness 
# for coding see: https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=1010
# 11 = parkinson's disease
# 10 = Alzheimer's disease/dementia

df$family_history_of_dementia <- apply(df[, c(596:605, 636:646)], 1, function(x) {
  if(any(x %in% c("11", "10"))) {
    return(1)
  } else {
    return(0)
  }
})
summary(as.factor(df$family_history_of_dementia))
rm(df_townsend_deprivation)
