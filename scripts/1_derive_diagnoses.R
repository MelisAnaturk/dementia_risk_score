# Written by Melis Anaturk

# This r script derives dates and diagnoses for dementia, stroke, ect. using: 
# gp/primary care data
# hospital inpatient data
# self-report
# death reports

#----------- 1. Data set-up -----------
# 1.1 load packages
library(tidyverse)
library(psych)

# 1.2 set working directory
setwd("./")

# 1.3 data pathway
data_pathway = "../../raw_data/"

#----------- 2. Primary care data ---------
# 2.1 identifying dementia codes based on records derived from gp practices
datalist <-list()

# 2.2 create list of diseases
disease <- c("dementia", "Depression", "Stroke", "TBI", "Diabetes","Diabetes_II", "TIA", "Atrial_fibrillation")

# 2.3 load in csv with individuals with primary care diagnoses
# code everyone with a diagnosis of interest as "1" and extract date of diagnosis
# lapply is much more efficient than a for loop in r - need to change this

for (x in disease){
  df <- read.csv(sprintf('../../raw_data/participants_with_%s.csv', x), header = TRUE, sep =',', stringsAsFactors = FALSE)
  df$eid <- as.character(df$eid)
  df[paste('primary_care_diagnosis_for', x, sep="_")] <- 1
  print(sprintf('number of individuals with %s', x))
  print(summary(as.factor(df[,sprintf('primary_care_diagnosis_for_%s', x)])))
  df <- df[c("eid", sprintf('primary_care_diagnosis_for_%s', x), "event_dt_1")]
  df[paste("primary_care_diagnosis_date_for", x, sep="_")] <- df$event_dt_1
  df <- df[c("eid", sprintf('primary_care_diagnosis_for_%s',x),sprintf('primary_care_diagnosis_date_for_%s', x))]
  datalist[[x]] <- df}

# 2.4 medication list
meds <- c("Dementia", "NSAIDs", "statins", "HRTs", "Diabetes_II")

# 2.5 code anyone with a prescription code as "1" from primary care data
for (x in meds){
  df <- read.csv(sprintf('../../raw_data/participants_%s.csv', x), header = TRUE, sep =',', stringsAsFactors = FALSE)
  df <- df[c("eid", "issue_date_1")]
  df$eid <- as.character(df$eid)
  df[paste("primary_care_prescription_for", x, sep="_")] <- 1
  print(sprintf('number of individuals with prescriptions related to (%s)', x))
  print(summary(as.factor(df[,sprintf('primary_care_prescription_for_%s', x)])))
  df[paste("primary_care_prescription_date_for", x, sep="_")]  <- df$issue_date_1
  df <- df[c("eid", sprintf('primary_care_prescription_for_%s', x),sprintf('primary_care_prescription_date_for_%s', x))]
  datalist[[sprintf("%s_1", x)]] <- df}

primary_care_df = Reduce(function(...) merge(..., all=T, by="eid"), datalist)

  #------------- 3. Secondary care data ------
# 3.1 identifying dementia codes based on HES data - hospital inpatient records
df_diagnoses <- read.csv(paste0(data_pathway,"ukb50321_diagnoses.csv"), header=TRUE, sep=",", stringsAsFactors = FALSE)
diagnoses_list <- read.csv(paste0(data_pathway,"ICD9_10_codes.csv") , header = TRUE, sep =',', stringsAsFactors = FALSE)

# 3.2 convert eid as character
df_diagnoses$eid <- as.character(df_diagnoses$eid)

# 3.4 set NAs as real NAs
diagnoses_list[diagnoses_list == ""] <- NA

# 3.5 create list of diseases
datalist <-list()
disease <- c("Dementia", "Depression", "Stroke", "TBI", "Diabetes")

# 3.6 extract ICD-9/ICD-10 codes 
for (x in disease){
  codes_to_extract <-diagnoses_list[sprintf('%s_ICD910_codes', x)]
  codes_to_extract <- codes_to_extract[complete.cases(codes_to_extract),]
  df_diagnoses[paste("secondary_care_diagnosis_of", x, sep="_")]<- apply(df_diagnoses[, -1], 1, function(x) {
    if(any(x %in% codes_to_extract)) {
      return(1)
    } else {
      return(0)
    }
  })}

# 3.7 merge primary and secondary care data using purr
merged_diagnoses_df <- list(df_diagnoses,primary_care_df) %>% reduce(left_join, by = "eid")

#------------- 4. Death reports and self report diseases -----------------------
# 4.1 identify cases based on self-report data (baseline only)
df_self_report_diagnoses <- read.csv(paste0(data_pathway,"add_vars_stroke_death.csv"), header=TRUE, sep=",", stringsAsFactors = FALSE)
cols_diagoses <- read.csv(paste0(data_pathway,"names_stroke_ICD.csv"), header=TRUE, sep=",", stringsAsFactors = FALSE)
df_self_report_diagnoses[df_self_report_diagnoses == "NA"] <- NA
df_self_report_diagnoses$eid <- as.character(df_self_report_diagnoses$eid)

# 4.2 recode column names in self report df
newcols <- as.vector(cols_diagoses$col_name)
names(df_self_report_diagnoses) <- newcols
col_names_self_report_diagnoses <- data.frame(names(df_self_report_diagnoses))

# 4.3 baseline diagnosis of stroke 
# For codes see https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=6
# 1583 = ischaemic stroke
# 1081 = stroke

# One hot encoding for individuals with self report diagnosis of stroke at baseline
df_self_report_diagnoses$self_report_stroke_0_0 <- apply(df_self_report_diagnoses[, 37:70], 1, function(x) {
  if(any(x %in% c("1583", "1081"))) {
    return(1)
  } else {
    return(0)
  }
})

# convert variable to factor
df_self_report_diagnoses$self_report_stroke_0_0 <- as.factor(df_self_report_diagnoses$self_report_stroke_0_0)
summary(df_self_report_diagnoses$self_report_stroke_0_0)
#0      1 
#495807   6699 

# 4.4 baseline diagnosis of depression
# 1286 = depression
# 1291 = mania/bipolar disorder/manic depression
# 1531 = post-natal depression

df_self_report_diagnoses$self_report_depression_0_0 <- apply(df_self_report_diagnoses[, 37:70], 1, function(x) {
  if(any(x %in% c("1286","1291", "1531"))) {
    return(1)
  } else {
    return(0)
  }
})

df_self_report_diagnoses$self_report_depression_0_0 <- as.factor(df_self_report_diagnoses$self_report_depression_0_0)
summary(df_self_report_diagnoses$self_report_depression_0_0)
#0      1 
#472738  29768

# 4.5 identify participants with self-reported history of dementia at baseline
# 1262 = parkinsons disease
# 1263 = dementia/alzheimers/cognitive impairment

df_self_report_diagnoses$self_report_dementia_0_0 <- apply(df_self_report_diagnoses[, 37:70], 1, function(x) {
  if(any(x %in% c("1262","1263"))) {
    return(1)
  } else {
    return(0)
  }
})

df_self_report_diagnoses$self_report_dementia_0_0 <- as.factor(df_self_report_diagnoses$self_report_dementia_0_0)
summary(df_self_report_diagnoses$self_report_dementia_0_0)
#0      1 
#501530    976 

# merge using purr
#df <- list(df_merge,df_diagnoses_combined) %>% reduce(left_join, by = "eid")
#why is the above line here now? 

# 4.6 One hot encoding for individuals with dementia (across all phases)
df_self_report_diagnoses$self_report_dementia <- apply(df_self_report_diagnoses[, 37:172], 1, function(x) {
  if(any(x %in% c("1262","1263"))) {
    return(1)
  } else {
    return(0)
  }
})

df_self_report_diagnoses$self_report_dementia <- as.factor(df_self_report_diagnoses$self_report_dementia)
summary(df_self_report_diagnoses$self_report_dementia)
#0      1 
#501426   1080

# 4.7 Diabetes
load(paste0(data_pathway,"df_ukb_raw.rda"))
#df_self_report_diagnoses <- list(df_self_report_diagnoses, df_merge[,c("eid","Diabetes_diagnosed_bydoctor_0_0")]) %>% reduce(left_join, by = "eid")
#rp replaced above with below
df_self_report_diagnoses <- list(df_self_report_diagnoses, df_ukb_raw[,c("eid","Diabetes_diagnosed_bydoctor_0_0")]) %>% reduce(left_join, by = "eid")

# 4.8 Diabetes II - baseline
# 1223 = 	type 2 diabetes
df_self_report_diagnoses$self_report_diabetes_II_0_0 <- apply(df_self_report_diagnoses[, 37:70], 1, function(x) {
  if(any(x %in% c("1223"))) {
    return(1)
  } else {
    return(0)
  }
})

df_self_report_diagnoses$self_report_diabetes_II_0_0 <- as.factor(df_self_report_diagnoses$self_report_diabetes_II_0_0)
summary(df_self_report_diagnoses$self_report_diabetes_II_0_0)
#0      1 
#499130   3376 

# 4.9 TIA - baseline
# 1082 = transient ischaemic attack (tia)

df_self_report_diagnoses$self_report_TIA_0_0 <- apply(df_self_report_diagnoses[, 37:70], 1, function(x) {
  if(any(x %in% c("1082"))) {
    return(1)
  } else {
    return(0)
  }
})

df_self_report_diagnoses$self_report_TIA_0_0 <- as.factor(df_self_report_diagnoses$self_report_TIA_0_0)
summary(df_self_report_diagnoses$self_report_TIA_0_0)
#0      1 
#500727   1779 

# 4.10 Atrial fibrillation - baseline
# 1471 = atrial fibrillation
# 1483 = atrial flutter
df_self_report_diagnoses$self_report_atrial_fibrillation_0_0 <- apply(df_self_report_diagnoses[, 37:70], 1, function(x) {
  if(any(x %in% c("1471", "1483"))) {
    return(1)
  } else {
    return(0)
  }
})
df_self_report_diagnoses$self_report_atrial_fibrillation_0_0 <- as.factor(df_self_report_diagnoses$self_report_atrial_fibrillation_0_0)
summary(df_self_report_diagnoses$self_report_atrial_fibrillation_0_0)  
#0      1 
#498786   3720 

#------------- 5. Death reports -------
# identify dementia cases based on death reports
#load updated death report data from data refresh
df_deathreport_refreshed <- read.csv(paste0(data_pathway,"ukb50321_deathvars.csv"), header=TRUE, sep=",", stringsAsFactors = FALSE)
df_deathreport_refreshed$eid <- as.character(df_deathreport_refreshed$eid)

#look for dementia codes in death report columns of updated HES data, is all columns except 1 (eid)
df_deathreport_refreshed$death_report_dementia <- apply(df_deathreport_refreshed[, 2:31], 1, function(x) {
  if(any(x %in% c("F00","F000","F001","F002","F009","G30","G300","G301","G308","G309","F01","F010","F011","F012","F013",
                  "F018","F019","I673","F020","G310","A810","F02","F021","F022","F023","F024","F028","F03","F051","F106","G311","G318"))) {
    return(1)
  } else {
    return(0)
  }
})
df_deathreport_refreshed$death_report_dementia <- as.factor(df_deathreport_refreshed$death_report_dementia)
summary(df_deathreport_refreshed$death_report_dementia)
#0      1 
#499746   2667 

#merge death report field back with self report df
df_self_report_diagnoses <- list(df_self_report_diagnoses, df_deathreport_refreshed[,c("eid","death_report_dementia")]) %>% reduce(left_join, by = "eid")
df_self_report_diagnoses$death_report_dementia[is.na(df_self_report_diagnoses$death_report_dementia)] <- 0
summary(df_self_report_diagnoses$death_report_dementia)
#0      1 
#499839   2667 

# merge primary care and self-report. merged_diagnoses_df contains primary and secondary care diagnosis
df_diagnoses_combined <- list(merged_diagnoses_df, df_self_report_diagnoses) %>% reduce(left_join, by = "eid")

# merge above with secondary care data. dont need this as merged_diganoses_df contins secondary care
#df_diagnoses_combined <- merge(df_diagnoses_combined, df_diagnoses, by.c="eid", all.x=TRUE)

# remove dfs
rm(df_diagnoses,df_self_report_diagnoses,primary_care_df)

#------------- 6. Create final variables -----
# 6.1 This is the original dataset downloaded from UKB, merged with additional lifestyle/genetic variables
#rp is skipping this for now. i think i already have what this was with df_ukb_raw
#load(file = paste0(data_pathway, "ukb_data_orig_merged.rda"))

# TOTAL N = 502,413


#following 3 lines adjusted by RP to get the columns needed with the column ordering I have
#names(df_diagnoses_combined[,c(1, 310:336, 508:514, 516)])
#df_diagnoses_combined <- df_diagnoses_combined[,c(1, 310:336, 508:514, 516)]
#save(df_diagnoses_combined, file = paste0(data_pathway,"merged_diagnoses.rda"))

names(df_diagnoses_combined[,c(1,310,315:340,512:520)])
#it should just be 512:520 and not redefine atrial fib self report later
df_diagnoses_combined <- df_diagnoses_combined[,c(1,310,315:340,512:520)]
#save the current diagnosis df, so that we can load and pickup from this point in future instead of redefining diagnoses
save(df_diagnoses_combined, file = paste0(data_pathway,"ukbdata_interim_diagnoses.rda"))

# 6.2 merge with purr
#df <- list(df_merge,df_diagnoses_combined) %>% reduce(left_join, by = "eid")
df <- list(df_ukb_raw,df_diagnoses_combined) %>% reduce(left_join, by = "eid")
#rp good till now
# 6.3 if a participant has a dementia code either in self-report data, secondary care data, primary care data or cause of death then code as a 1
# BUG fixed: need to code primary care variables as numeric before coding dementia cases
df$primary_care_diagnosis_for_dementia <- as.numeric(df$primary_care_diagnosis_for_dementia)
df$primary_care_prescription_for_Dementia <- as.numeric(df$primary_care_prescription_for_Dementia)

# 6.4 this variable is behaving oddly so going to recode
df$primary_care_diagnosis_for_dementia[is.na(df$primary_care_diagnosis_for_dementia)] <- 0
df$primary_care_prescription_for_Dementia[is.na(df$primary_care_prescription_for_Dementia)] <- 0

# 6.5 create final dementia variable
myvars <- c("secondary_care_diagnosis_of_Dementia", "self_report_dementia", "death_report_dementia", "primary_care_diagnosis_for_dementia", "primary_care_prescription_for_Dementia")
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
#497757   4764 
#View(df[,c("dementia_BIN_TOTAL", "primary_care_diagnosis_for_dementia", "primary_care_prescription_for_Dementia", "secondary_care_diagnosis_of_Dementia", "self_report_dementia", "death_report_dementia")])

#------ 7. Time to diagnosis and disease status for dementia -----
# 7.1 Use base R for conversion to variables with "Date" format 
df$baseline_date                        <-  as.Date(df$Date_of_assessment_0_0_new, format="%Y-%m-%d")
df$date_all_cause_dementia_0_0          <-  as.Date(df$Date_of_all_cause_dementia_report_0_0, format="%Y-%m-%d")
df$primary_care_prescription_date_for_Dementia <-  as.Date(df$primary_care_prescription_date_for_Dementia, format="%d/%m/%Y")
#07/03/2022 addition:
df$primary_care_diagnosis_date_for_dementia <- as.Date(df$primary_care_diagnosis_date_for_dementia, format="%d/%m/%Y")
df$date_of_latest_follow_up <- as.Date(df$Date_of_assessment_3_0_new, format="%Y-%m-%d")
df$date_of_death_all <- ifelse(is.na(df$death_date_1), df$death_date_2, df$death_date_1)
df$date_of_death_all <-  as.Date(df$date_of_death_all, format="%Y-%m-%d")

# 7.2 remove impossible dates
#following chunk of code up until View command was what was here initially. i dont think it works. do a simple for loop
#date_variables <- c("baseline_date", "date_all_cause_dementia_0_0", "primary_care_prescription_date_for_Dementia")

#remove_impossible_dates <- function(x){
#  print(class(df[,x]))
#  df[,x] <- as.Date(df[,x], format="%d/%m/%Y")
#  df[,x][df[,x]=="1901-01-01"] <- NA
#  df[,x][df[,x]=="1900-01-01"] <- NA
#  df[,x][df[,x]=="1902-02-02"] <- NA
#  df[,x][df[,x]=="1903-03-03"] <- NA
#  df[,x][df[,x]=="2037-07-07"] <- NA
#}

#lapply(df, remove_impossible_dates)
#summary(df[date_variables])
#View(df[date_variables])

#date_variables <- c("baseline_date", "date_all_cause_dementia_0_0", "primary_care_prescription_date_for_Dementia")
date_variables <- c("baseline_date", "date_all_cause_dementia_0_0", "primary_care_diagnosis_date_for_dementia","primary_care_prescription_date_for_Dementia")
summary(df[date_variables]) #some impossible dates (1900-01-01)
View(df[date_variables])

for (x in date_variables){
  print(class(df[,x]))
  df[,x] <- as.Date(df[,x], format="%d/%m/%Y")
  df[,x][df[,x]=="1901-01-01"] <- NA
  df[,x][df[,x]=="1900-01-01"] <- NA
  df[,x][df[,x]=="1902-02-02"] <- NA
  df[,x][df[,x]=="1903-03-03"] <- NA
  df[,x][df[,x]=="2037-07-07"] <- NA
}
summary(df[date_variables]) #no impossible dates
View(df[date_variables])

# 7.3 Dementia cases: create a variable that reflects the earliest date of diagnosis
df <- transform(df, earliest_dementia_date = pmin(date_all_cause_dementia_0_0, primary_care_diagnosis_date_for_dementia, primary_care_prescription_date_for_Dementia, na.rm=TRUE))
summary(df$earliest_dementia_date)
class(df$earliest_dementia_date)
#View(df[c("earliest_dementia_date", "date_all_cause_dementia_0_0", "primary_care_prescription_date_for_Dementia","primary_care_diagnosis_date_for_dementia")])

# 7.4 create a difference between dates variable
df$date_diff_all_cause_dementia_0_0                    <- df$earliest_dementia_date - df$baseline_date
df$date_diff_time_between_baseline_latest_follow_up    <- df$date_of_latest_follow_up - df$baseline_date
df$date_diff_time_between_baseline_death               <- df$date_of_death_all - df$baseline_date

#double check everything looks ok
View(df[c("date_diff_time_between_baseline_death", "date_of_death_all", "baseline_date")])

myvars <- c("baseline_date", "earliest_dementia_date",  "date_diff_all_cause_dementia_0_0")
summary(as.factor(df$Date_of_all_cause_dementia_report_0_0))
View(df[myvars])

# save file
#save(df, file = paste0(data_pathway,"ukb_data_orig_merged.rda"))
save(df, file = paste0(data_pathway,"ukbdata_diagnoses.rda"))

# 7.5 clear working space
rm(list=ls())

# 7.6 load in .rda file
#load(file = paste0(data_pathway,"ukb_data_orig_merged.rda"))
data_pathway = "../../raw_data/"
load(file = paste0(data_pathway,"ukbdata_diagnoses.rda"))
#load(file="../../../raw_data/ukb_data_orig_merged_final_diseases_oct22.rda")

# 7.7 create an index date for individuals without dementia (for individuals with dementia it is date of diagnosis)
df$date_diff_healthy <- ifelse(is.na(df$date_diff_time_between_baseline_latest_follow_up), df$date_diff_time_between_baseline_death,
                        df$date_diff_time_between_baseline_latest_follow_up)
summary(df$date_diff_healthy)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#      4    1263    2098    2033    2786    4885  481249 
View(df[c("date_diff_time_between_baseline_latest_follow_up", "date_diff_time_between_baseline_death", "date_diff_healthy")])

# combine into same variables
for (x in c("date_diff_all_cause_dementia_0_0", "date_diff_time_between_baseline_latest_follow_up", "date_diff_time_between_baseline_death"))
{df[,x] <- as.numeric(df[,x])
print(summary(df[,x]))
print(class(df[,x])) }

df$final_date_diff <- ifelse(is.na(df$date_diff_all_cause_dementia_0_0), df$date_diff_healthy, df$date_diff_all_cause_dementia_0_0)
View(df[c("final_date_diff", "date_diff_all_cause_dementia_0_0", "date_diff_healthy")])
summary(df$final_date_diff)

# 7.8 follow-up for non-dementia cases
df <- transform(df, latest_follow_up_date_1 = pmax(df$date_of_latest_follow_up, df$date_of_death_all, as.Date(df$Date_of_assessment_2_0_new, format="%Y-%m-%d"), as.Date(df$Date_of_assessment_1_0_new, format="%Y-%m-%d"), baseline_date, na.rm=TRUE))
df$date_diff_baseline_all <- df$latest_follow_up_date_1- df$baseline_date
summary(as.numeric(df$date_diff_baseline_all))

df$date_diff_final <- ifelse(is.na(df$date_diff_all_cause_dementia_0_0), df$date_diff_baseline_all, df$date_diff_all_cause_dementia_0_0)
summary(df$date_diff_final)

# Checks
View(df[c("date_diff_final", "date_diff_baseline_all", "date_diff_all_cause_dementia_0_0")])
#View(df[c("earliest_dementia_date", "date_all_cause_dementia", "baseline_date", "years_diff_all_time", "dementia_BIN_TOTAL", "date_diff_all_cause_dementia_0_0", "final_date_diff")])
#View(df[, c("date_diff_baseline_all","latest_follow_up_date_1", "date_of_latest_follow_up", "date_of_death_all", "Date_of_assessment_2_0_new", "Date_of_assessment_1_0_new", "baseline_date")])

# convert date diff to years
df$years_diff_all_time  <- df$date_diff_baseline_all /365.242 
summary(as.numeric(df$years_diff_all_time))

#------- 8. Begin data subsetting -----
# 8.1 exclude people without a baseline asssessment
df <- subset(df, complete.cases(baseline_date))
#  CURRENT N = 502,506 

# 8.2 exclude people with pre-existing dementia based on date_diff variable
summary(df$date_diff_all_cause_dementia)
df <- subset(df, date_diff_all_cause_dementia_0_0>0| is.na(date_diff_all_cause_dementia_0_0)) 
# CURRENT N = 501,871
summary(df$dementia_BIN_TOTAL)
#0      1 
#497742   4129 

# 8.3 identify anyone who reported dementia at baseline
summary(df$self_report_dementia_0_0)
#  0      1 
#501020    851 

# 8.4 exclude people self-reporting a diagnosis at baseline
df <- subset(df, self_report_dementia_0_0==0)

#n = 501 020

# 8.5 convert date diff to years
df$years_diff_all_cause_dementia_0_0  <- df$date_diff_all_cause_dementia_0_0 /365.242 

describe(as.numeric(df$years_diff_all_cause_dementia_0_0))
#vars    n mean   sd median trimmed  mad  min   max range  skew kurtosis   se
#X1    1 3183 5.61 2.25    5.9    5.72 2.37 0.03 10.71 10.68 -0.41    -0.56 0.04
hist(as.numeric(df$years_diff_all_cause_dementia_0_0))

# 8.6 exclude people who have developed dementia within a year of baseline (minimize reverse causation)

df <- subset(df, years_diff_all_cause_dementia_0_0>=1| is.na(date_diff_all_cause_dementia_0_0)) 
#n = 500 920

# CURRENT N = 501,261
#updated n, after incorporating pcare diagnosis date: 500 920

summary(df$dementia_BIN_TOTAL)
#0      1 
#497742   3178 

describe(as.numeric(df$years_diff_all_cause_dementia_0_0))
#vars    n mean   sd median trimmed  mad min   max range skew kurtosis   se
#X1    1 3083 5.77 2.09      6    5.85 2.25   1 10.71   9.7 -0.3    -0.69 0.04
hist(as.numeric(df$years_diff_all_cause_dementia_0_0))

#------- 9. Dates for other diseases -----
# 9.1 create vector of primary care diagnosis dates
date_variables <- c("primary_care_diagnosis_date_for_Stroke", "primary_care_diagnosis_date_for_Depression", "primary_care_diagnosis_date_for_TIA", "primary_care_diagnosis_date_for_TBI", "primary_care_diagnosis_date_for_Diabetes", "primary_care_diagnosis_date_for_Diabetes_II")
#View(df[,date_variables]) #date_variables])

#lapply(date_variables, remove_impossible_dates)
#summary(df[date_variables])
#View(df[date_variables])
sapply(df[date_variables], class) #they are characters, make them dates
for (x in date_variables){
  df[,x] <- as.Date(df[,x], format="%d/%m/%Y")
}
summary(df[date_variables]) #some impossible dates (1900-01-01, 1902-02-02, etc)
View(df[date_variables])

#remove impossible dates
for (x in date_variables){
  print(class(df[,x]))
  df[,x] <- as.Date(df[,x], format="%d/%m/%Y")
  df[,x][df[,x]=="1901-01-01"] <- NA
  df[,x][df[,x]=="1900-01-01"] <- NA
  df[,x][df[,x]=="1902-02-02"] <- NA
  df[,x][df[,x]=="1903-03-03"] <- NA
  df[,x][df[,x]=="2037-07-07"] <- NA
}
summary(df[date_variables]) #no impossible dates
View(df[date_variables])

# 9.2 Stroke
# BASELINE STROKE ONLY
df$primary_care_diagnosis_for_Stroke <- as.numeric(df$primary_care_diagnosis_for_Stroke)
df$primary_care_diagnosis_for_Stroke[is.na(df$primary_care_diagnosis_for_Stroke)] <- 0

# final stroke variable
df$date_diff_stroke_diagnosis <- df$primary_care_diagnosis_date_for_Stroke - df$baseline_date
summary(as.factor(df$date_diff_stroke_diagnosis))

# excluding cases who developed dementia during study period
df$stroke_BIN_FINAL_0_0 <- ifelse(df$date_diff_stroke_diagnosis<=0&df$primary_care_diagnosis_for_Stroke==1|df$stroke_BIN_TOTAL_0_0==1,1,0) 
df$stroke_BIN_FINAL_0_0[is.na(df$stroke_BIN_FINAL_0_0)] <- 0
df$stroke_BIN_FINAL_0_0 <- as.factor(df$stroke_BIN_FINAL_0_0)
summary(df$self_report_stroke_0_0)
#0      1 
#494298   6622 
summary(df$stroke_BIN_FINAL_0_0)
#0      1 
#493147   7773 

View(df[,c("stroke_BIN_FINAL_0_0", "primary_care_diagnosis_date_for_Stroke", "baseline_date","date_diff_stroke_diagnosis")])

# 9.3 Depression
df$primary_care_diagnosis_date_for_Depression <-  as.Date(df$primary_care_diagnosis_date_for_Depression, format="%d/%m/%Y")
df$date_diff_depression_diagnosis <- df$primary_care_diagnosis_date_for_Depression - df$baseline_date
summary(as.factor(df$date_diff_depression_diagnosis))

# primary_care depression variable
df$primary_care_diagnosis_for_Depression[is.na(df$primary_care_diagnosis_for_Depression)] <- 0
df$depression_BIN_FINAL_0_0 <- ifelse(df$date_diff_depression_diagnosis<=0&df$primary_care_diagnosis_for_Depression==1|df$self_report_depression_0_0==1,1,0) 
df$depression_BIN_FINAL_0_0 <- as.factor(df$depression_BIN_FINAL_0_0)
summary(df$depression_BIN_FINAL_0_0)
#0      1   NA's 
#443766  57146      8 

# 9.4 TBI
df$TBI_BIN_TOTAL <- ifelse(is.na(df[,c("primary_care_diagnosis_for_TBI")]),0,1)

df$TBI_BIN_TOTAL <- as.factor(df$TBI_BIN_TOTAL)
summary(df$TBI_BIN_TOTAL)
#0      1 
#484962  16228
summary(as.factor(df$primary_care_diagnosis_for_TBI))
#1   NA's 
# 16228 484692

df$primary_care_diagnosis_date_for_TBI <-  as.Date(df$primary_care_diagnosis_date_for_TBI, format="%d/%m/%Y")
df$date_diff_TBI_diagnosis <- df$primary_care_diagnosis_date_for_TBI - df$baseline_date
summary(as.factor(df$date_diff_TBI_diagnosis))

df$TBI_BIN_FINAL_0_0 <- ifelse(df$date_diff_TBI_diagnosis<=0&df$primary_care_diagnosis_for_TBI==1,1,0) 
df$TBI_BIN_FINAL_0_0[is.na(df$TBI_BIN_FINAL_0_0)] <- 0
df$TBI_BIN_FINAL_0_0 <- as.factor(df$TBI_BIN_FINAL_0_0)
summary(df$TBI_BIN_FINAL_0_0)
#0      1 
#490620  10300

# 9.5 Diabetes - BASELINE
# incident + historic diabetes
df$primary_care_diagnosis_date_for_Diabetes <-  as.Date(df$primary_care_diagnosis_date_for_Diabetes, format="%d/%m/%Y")
df$date_diff_Diabetes_diagnosis <- df$primary_care_diagnosis_date_for_Diabetes - df$baseline_date
summary(as.factor(df$date_diff_Diabetes_diagnosis))

# historic diabetes
df$Diabetes_BIN_FINAL_0_0 <- ifelse(df$date_diff_Diabetes_diagnosis<=0&df$primary_care_diagnosis_for_Diabetes==1|df$Diabetes_diagnosed_bydoctor_0_0.x==1,1,0) 
df$Diabetes_BIN_FINAL_0_0[is.na(df$Diabetes_BIN_FINAL_0_0)] <- 0
df$Diabetes_BIN_FINAL_0_0 <- as.factor(df$Diabetes_BIN_FINAL_0_0)
summary(as.factor(df$Diabetes_diagnosed_bydoctor_0_0.x))
#0      1   NA's 
#472037  26281   2602
summary(df$Diabetes_BIN_FINAL_0_0)
#0      1 
#474304  26616 

# 9.6 Type 2 Diabetes
df$primary_care_diagnosis_date_for_Diabetes_II <-  as.Date(df$primary_care_diagnosis_date_for_Diabetes_II, format="%d/%m/%Y")
df$date_diff_Diabetes_II_diagnosis <- df$primary_care_diagnosis_date_for_Diabetes_II - df$baseline_date
summary(as.factor(df$date_diff_Diabetes_II_diagnosis))

df$Diabetes_II_BIN_FINAL_0_0 <- ifelse(df$date_diff_Diabetes_II_diagnosis<=0&df$primary_care_diagnosis_for_Diabetes_II==1,1,0) 
df$Diabetes_II_BIN_FINAL_0_0[is.na(df$Diabetes_II_BIN_FINAL_0_0)] <- 0
df$Diabetes_II_BIN_FINAL_0_0 <- as.factor(df$Diabetes_II_BIN_FINAL_0_0)
summary(as.factor(df$Diabetes_diagnosed_bydoctor_0_0.x))
summary(df$Diabetes_II_BIN_FINAL_0_0)
#0      1 
#490939   9981

View(df[,c("Diabetes_II_BIN_FINAL_0_0", "primary_care_diagnosis_for_Diabetes_II", "primary_care_diagnosis_date_for_Diabetes_II", "baseline_date")])

# 9.7 TIA
df$primary_care_diagnosis_for_TIA[is.na(df$primary_care_diagnosis_for_TIA)] <-0
summary(df$primary_care_diagnosis_for_TIA)
class(df$primary_care_diagnosis_for_TIA)

df$primary_care_diagnosis_date_for_TIA <-  as.Date(df$primary_care_diagnosis_date_for_TIA, format="%d/%m/%Y")
df$date_diff_TIA_diagnosis <- df$primary_care_diagnosis_date_for_TIA - df$baseline_date
summary(as.factor(df$date_diff_TIA_diagnosis))

df$TIA_BIN_FINAL_0_0 <- ifelse(df$date_diff_TIA_diagnosis<=0&df$primary_care_diagnosis_date_for_TIA==1|df$self_report_TIA_0_0==1,1,0)
df$TIA_BIN_FINAL_0_0[is.na(df$TIA_BIN_FINAL_0_0)] <- 0
df$TIA_BIN_FINAL_0_0 <- as.factor(df$TIA_BIN_FINAL_0_0)
summary(df$TIA_BIN_FINAL_0_0)
#0      1 
#499156   1764  

# 9.8 atrial fibrillation - baseline #rp: atrial fib self report was already computed, no need to redo so 9.8, 9.9 commented out
#df_self_report_diagnoses <- read.csv('../../raw_data/add_vars_stroke_death.csv', header=TRUE, sep=",", stringsAsFactors = FALSE)
#cols_diagoses <- read.csv('../../names_stroke_ICD.csv', header=TRUE, sep=",", stringsAsFactors = FALSE)
#df_self_report_diagnoses[df_self_report_diagnoses == "NA"] <- NA
#names(df_self_report_diagnoses) <- cols_diagoses[,2]
#df_self_report_diagnoses$eid <- as.character(df_self_report_diagnoses$eid)

#df_self_report_diagnoses$self_report_atrial_fibrillation_0_0 <- apply(df_self_report_diagnoses[, 37:70], 1, function(x) {
#  if(any(x %in% c("1471", "1483"))) {
#    return(1)
#  } else {
#    return(0)
#  }
#})
#df_self_report_diagnoses$self_report_atrial_fibrillation_0_0 <- as.factor(df_self_report_diagnoses$self_report_atrial_fibrillation_0_0)
#summary(df_self_report_diagnoses$self_report_atrial_fibrillation_0_0)  

# 9.9 merge with main df
#df_self_report_diagnoses <- df_self_report_diagnoses[, c(1, 33:173)]
#df <- list(df, df_self_report_diagnoses) %>% reduce(left_join, by = "eid")
#summary(df$self_report_atrial_fibrillation_0_0)

df$primary_care_diagnosis_date_for_Atrial_fibrillation <-  as.Date(df$primary_care_diagnosis_date_for_Atrial_fibrillation, format="%d/%m/%Y")
df$date_diff_Atrial_fibrillation_diagnosis <- df$primary_care_diagnosis_date_for_Atrial_fibrillation - df$baseline_date
summary(as.factor(df$date_diff_Atrial_fibrillation_diagnosis))

df$Atrial_Fibrillation_BIN_FINAL_0_0 <- ifelse(df$date_diff_Atrial_fibrillation_diagnosis<=0&df$primary_care_diagnosis_for_Atrial_fibrillation==1|df$self_report_atrial_fibrillation_0_0==1,1,0) 
df$Atrial_Fibrillation_BIN_FINAL_0_0[is.na(df$Atrial_Fibrillation_BIN_FINAL_0_0)] <- 0
df$Atrial_Fibrillation_BIN_FINAL_0_0 <- as.factor(df$Atrial_Fibrillation_BIN_FINAL_0_0)
summary(df$Atrial_Fibrillation_BIN_FINAL_0_0)
#0      1 
#495548   5372 

# 9.10 Save file
#save(df, file = paste0(data_pathway, "ukb_data_orig_merged_final_diseases_oct22.rda"))
#load(file = paste0(data_pathway,"ukb_data_orig_merged_final_diseases_oct22.rda"))
save(df, file = paste0(data_pathway, "ukbdata_diagnoses_baseline_diseasestatus.rda"))

#--------- 10. Medications 
# read csv
df_meds <- read.csv(paste0(data_pathway, "add_vars_meds.csv"), header=TRUE, sep=",", stringsAsFactors = FALSE)
col_names_meds <- read.csv(paste0(data_pathway, "names_meds.csv"), header=TRUE, sep=",", stringsAsFactors = FALSE)
list_meds <- read.csv(paste0(data_pathway, "list_meds_new.csv"), header=TRUE, sep=",", stringsAsFactors = FALSE)

#rename df_med headings and make "NA" real NA
newcols <- as.vector(col_names_meds$col_name)
names(df_meds) <- newcols
df_meds[df_meds == "NA"] <- NA
df_meds$eid <- as.character(df_meds$eid)

df <- list(df, df_meds) %>% reduce(left_join, by = "eid")

# list of drugs
names(list_meds)
list_meds[list_meds == "NA"] <- NA
list_meds[list_meds == ""] <- NA
data.frame(names(list_meds))
View(names(df))

# 10.1 Statins (baseline only)
meds_statins <- as.vector(as.character(na.omit(list_meds[,2])))
print(meds_statins)

#**MA has columns 893:940 as baseline treatment (Treatment_medication_code_0_*)
#Rp instead has this as 753:800 because rp did not remerge cols 33:173 of df_self_report_diagnoses in 9.9
df$statins_0_0 <- apply(df[, 753:800], 1, function(x) {  # 1230:1277
  if(any(x %in% meds_statins)) {
    return(1)
  } else {
    return(0)
  }
})

df$statins_0_0 <- as.factor(df$statins_0_0)
summary(df$statins_0_0)
#0      1 
#418733  82187

# 10.2 HRT (baseline only)
meds_HRT <- as.vector(as.character(na.omit(list_meds[,4])))
print(meds_HRT)

df$HRT_meds_0_0 <- apply(df[, 753:800], 1, function(x) {
  if(any(x %in% meds_HRT)) {
    return(1)
  } else {
    return(0)
  }
})

df$HRT_meds_0_0 <- as.factor(df$HRT_meds_0_0)
summary(df$HRT_meds_0_0)
#0      1 
#482606  18314 

# 10.3 NSAIDs (Excluding Aspirin)
meds_NSAIDs <- as.vector(as.character(na.omit(list_meds[,6])))
print(meds_NSAIDs)

df$NSAIDs_0_0 <- apply(df[,  753:800], 1, function(x) {
  if(any(x %in% meds_NSAIDs)) {
    return(1)
  } else {
    return(0)
  }
})

df$NSAIDs_0_0  <- as.factor(df$NSAIDs_0_0)
summary(df$NSAIDs_0_0)
#0      1 
#404164  96756 

# 10.4 ASPIRIN
meds_Aspirin <- as.vector(as.character(na.omit(list_meds[,14])))
print(meds_Aspirin)

df$Aspirin_0_0 <- apply(df[, 753:800], 1, function(x) {
  if(any(x %in% meds_Aspirin)) {
    return(1)
  } else {
    return(0)
  }
})

df$Aspirin_0_0  <- as.factor(df$Aspirin_0_0)
summary(df$Aspirin_0_0)
#0      1 
#431456  69464

# 10.5 Antihypertensive medications
meds_Antihypertensive <- as.vector(as.character(na.omit(list_meds[,10])))
print(meds_Antihypertensive)

df$Antihypertensive_meds_0_0 <- apply(df[, 753:800], 1, function(x) {
  if(any(x %in% meds_Antihypertensive)) {
    return(1)
  } else {
    return(0)
  }
})

df$Antihypertensive_meds_0_0  <- as.factor(df$Antihypertensive_meds_0_0)
summary(df$Antihypertensive_meds_0_0)
#0      1 
#387890 113030

# 10.6 Antidepressants
meds_Antidepressants <- as.vector(as.character(na.omit(list_meds[,12])))
print(meds_Antidepressants)

df$Antidepressant_meds_0_0 <- apply(df[, 753:800], 1, function(x) {
  if(any(x %in% meds_Antidepressants)) {
    return(1)
  } else {
    return(0)
  }
})

df$Antidepressant_meds_0_0  <- as.factor(df$Antidepressant_meds_0_0)
summary(df$Antidepressant_meds_0_0)
#0      1 
#463523  37397 

rm(col_names_meds,df_meds,list_meds)

length(unique(df$eid))

# Save file
#save(df, file = paste0(data_pathway, "ukb_data_orig_merged_final_diseases_oct22.rda"))
#load(file = paste0(data_pathway, "ukb_data_orig_merged_final_diseases_oct22.rda"))
save(df, file = paste0(data_pathway, "ukbdata_diagnoses_baseline_diseasestatus_baselinemedications.rda"))
