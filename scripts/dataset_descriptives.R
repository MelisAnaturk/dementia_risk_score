#this script is used to check demographic descriptives for Table 1
library(qwraps2)

load("../../raw_data/train_data_outliers_removed.rda")
load("../../raw_data/test_data_outliers_removed.rda")

test.data$dataset <- "UKBtest"
train.data$dataset <- "UKB_train"

df <- rbind(test.data, train.data) 

#variables i need stats on
vars_demographic <- c("Age_when_attended_assesment_centre_0_0","Sex","education_years", "Townsend_deprivation_Groups_0_0")

summary_demographic <-
  list("Age" = 
         list("Age at baseline, years" = ~ qwraps2::mean_sd(Age_when_attended_assesment_centre_0_0, denote_sd = "paren")),
       "Sex" =
         list("Female" = ~ qwraps2::n_perc(Sex == 0, show_denom = "never")),
       "Education" =
         list("Education, years" = ~ qwraps2::mean_sd(education_years, denote_sd = "paren")),
       "Townsend Deprivation" = 
         list("Group 1" = ~ qwraps2::n_perc(Townsend_deprivation_Groups_0_0 == 0),
              "Group 2" = ~ qwraps2::n_perc(Townsend_deprivation_Groups_0_0 == 1),
              "Group 3" = ~ qwraps2::n_perc(Townsend_deprivation_Groups_0_0 == 2),
              "Group 4" = ~ qwraps2::n_perc(Townsend_deprivation_Groups_0_0 == 3),
              "Group 5" = ~ qwraps2::n_perc(Townsend_deprivation_Groups_0_0 == 4))
  )

df_demographic <- summary_table(df, summary_demographic, by = c("dataset"))
#~ qwraps2::mean_sd(, denote_sd = "paren")
#~ qwraps2::n_perc()

biomed_vars <- c("BMI_0_0","Systolic_BP_auto_mean","LDL_0_0","HDL_cholesterol_0_0", "Total_cholesterol_0_0")

summary_biomed <-
  list("BMI" = 
         list("BMI, kg/m2" = ~ qwraps2::mean_sd(BMI_0_0, denote_sd = "paren")),
       "Systolic BP" = 
         list("Systolic BP, mm Hg" = ~ qwraps2::mean_sd(Systolic_BP_auto_mean, denote_sd = "paren")),
       "Total Cholesterol" = 
         list("Total, mmol/L" = ~ qwraps2::mean_sd(Total_cholesterol_0_0, denote_sd = "paren")),
       "HDL" = 
         list("HDL, mmol/L" = ~ qwraps2::mean_sd(HDL_cholesterol_0_0, denote_sd = "paren")),
       "LDL" = 
         list("LDL, mmol/L" = ~ qwraps2::mean_sd(LDL_0_0, denote_sd = "paren"))
)

df_biomed <- summary_table(df, summary_biomed, by = c("dataset"))

#df_table1 <- rbind(data.frame(df_demographic), data.frame(df_biomed))

lifestyle_vars <- c("Sleep_duration_0_0", "IPAQ_activity_group_0_0", "units_combined", "Smoker_bin", "total_fish_intake_per_week_0_0","Sleeplesness_insomnia_0_0_bin")
#need to add leisure social activities...cound number of responses in Leisure_social_activities_0_*

#leisure_vars<-c("Leisure_social_activities_0_0", "Leisure_social_activities_0_1", "Leisure_social_activities_0_2", "Leisure_social_activities_0_3", "Leisure_social_activities_0_4")
df$Leisure_social_activities_0_0_BIN<-ifelse(is.na(df$Leisure_social_activities_0_0),0,
                                                    ifelse(df$Leisure_social_activities_0_0>0,1,0))
df$Leisure_social_activities_0_1_BIN<-ifelse(is.na(df$Leisure_social_activities_0_1),0,
                                                    ifelse(df$Leisure_social_activities_0_1>0,1,0))
df$Leisure_social_activities_0_2_BIN<-ifelse(is.na(df$Leisure_social_activities_0_2),0,
                                                    ifelse(df$Leisure_social_activities_0_2>0,1,0))
df$Leisure_social_activities_0_3_BIN<-ifelse(is.na(df$Leisure_social_activities_0_3),0,
                                                    ifelse(df$Leisure_social_activities_0_3>0,1,0))
df$Leisure_social_activities_0_4_BIN<-ifelse(is.na(df$Leisure_social_activities_0_4),0,
                                                    ifelse(df$Leisure_social_activities_0_4>0,1,0))

leisure_vars<-c("Leisure_social_activities_0_0_BIN", "Leisure_social_activities_0_1_BIN", "Leisure_social_activities_0_2_BIN", 
                "Leisure_social_activities_0_3_BIN", "Leisure_social_activities_0_4_BIN")
#test.data[leisure_vars][is.na(test.data$leisure_vars)] <- 0

#df$weekly_leisure_activities<-sum(df[leisure_vars])
df$weekly_leisure_activities<-rowSums(df[,leisure_vars])
View(df[,c("Leisure_social_activities_0_0", "Leisure_social_activities_0_1", 
                  "Leisure_social_activities_0_2", "Leisure_social_activities_0_3", "Leisure_social_activities_0_4","weekly_leisure_activities")])
View(df[leisure_vars])


lifestyle_vars <- c("Sleep_duration_0_0", "IPAQ_activity_group_0_0", "units_combined", "Smoker_bin", "total_fish_intake_per_week_0_0",
                    "Number_in_household_0_0", "weekly_leisure_activities","Sleeplesness_insomnia_0_0_bin")
summary_lifestyle <-
  list("Sleep duration" = 
         list("Sleep duration" = ~ qwraps2::mean_sd(Sleep_duration_0_0, denote_sd = "paren")),
       "IPAQ" = 
         list("Low" = ~ qwraps2::n_perc(IPAQ_activity_group_0_0 == 0),
              "Moderate-vigorous" = ~ qwraps2::n_perc(IPAQ_activity_group_0_0 == 1)),
       "Alcohol" = 
         list("Alcohol" = ~ qwraps2::mean_sd(units_combined, denote_sd = "paren")),
       "Smoker" =
         list("Current smoker" = ~ qwraps2::n_perc(Smoker_bin == 1, show_denom = "never")),
       "Weekly Fish Intake" = 
         list("Weekly Fish Intake" = ~ qwraps2::mean_sd(total_fish_intake_per_week_0_0, denote_sd = "paren")),
       "N. in household" = 
         list("N. in household" = ~ qwraps2::mean_sd(Number_in_household_0_0, denote_sd = "paren")),
       "N. leisure" = 
         list("N. of weekly leisure activities" = ~ qwraps2::mean_sd(weekly_leisure_activities, denote_sd = "paren")),
       "Insomnia" = 
         list("Never" = ~ qwraps2::n_perc(Sleeplesness_insomnia_0_0_bin == 0),
              "Sometimes" = ~ qwraps2::n_perc(Sleeplesness_insomnia_0_0_bin == 1),
              "Usually" = ~ qwraps2::n_perc(Sleeplesness_insomnia_0_0_bin == 2))
  )

df_lifestyle <- summary_table(df, summary_lifestyle, by = c("dataset"))

genetic_vars <- c("APOE_genotype_bin")

summary_genetic <-
  list("APOE" = 
         list("APOE4status" = ~ qwraps2::n_perc(APOE_genotype_bin == 1, show_denom = "never"))
  )
df_genetic <- summary_table(df, summary_genetic, by = c("dataset"))


med_vars <- c("current_history_depression","Diabetes_BIN_FINAL_0_0", "stroke_TIA_BIN_FINAL", "TBI_BIN_FINAL_0_0", "Atrial_Fibrillation_BIN_FINAL_0_0",
               "statins_0_0", "NSAIDs_0_0", "Aspirin_0_0","HRT_0_0","Antihypertensive_meds_0_0","family_history_of_dementia","dementia_BIN_TOTAL","years_diff_all_cause_dementia_0_0")

summary_medical <-
  list("Depression" = 
         list("Current or history of Depression" = ~ qwraps2::n_perc(current_history_depression == 1, show_denom = "never")),
       "Diabetes" = 
         list("History of Diabetes" = ~ qwraps2::n_perc(Diabetes_BIN_FINAL_0_0 == 1, show_denom = "never")),
       "Stroke/TIA" = 
         list("Stroke/TIA" = ~ qwraps2::n_perc(stroke_TIA_BIN_FINAL == 1, show_denom = "never")),
       "TBI" = 
         list("TBI" = ~ qwraps2::n_perc(TBI_BIN_FINAL_0_0 == 1, show_denom = "never")),
       "AtrialFib" = 
         list("AtrialFib" = ~ qwraps2::n_perc(Atrial_Fibrillation_BIN_FINAL_0_0 == 1, show_denom = "never")),
       "Statins" = 
         list("Statins" = ~ qwraps2::n_perc(statins_0_0 == 1, show_denom = "never")),
       "NSAIDs" = 
         list("NSAIDs" = ~ qwraps2::n_perc(NSAIDs_0_0 == 1, show_denom = "never")),
       "Aspirin" = 
         list("Aspirin" = ~ qwraps2::n_perc(Aspirin_0_0 == 1, show_denom = "never")),
       "HRT" = 
         list("HRT" = ~ qwraps2::n_perc(HRT_0_0 == 1, show_denom = "never")),
       "Antihypertensive" = 
         list("Antihypertensive" = ~ qwraps2::n_perc(Antihypertensive_meds_0_0 == 1, show_denom = "never")),
       "FamHX" = 
         list("FamHX" = ~ qwraps2::n_perc(family_history_of_dementia == 1, show_denom = "never")),
       "Incident dementia" = 
         list("Incident dementia" = ~ qwraps2::n_perc(dementia_BIN_TOTAL == 1, show_denom = "never")),
       "Years to diagnosis" = 
         list("Years to diagnosis" = ~ qwraps2::mean_sd(years_diff_all_cause_dementia_0_0, denote_sd = "paren", na_rm = TRUE))
  )

df_medical <- summary_table(df, summary_medical, by = c("dataset"))

table1 <- rbind(df_demographic, df_biomed, df_lifestyle, df_genetic, df_medical)
