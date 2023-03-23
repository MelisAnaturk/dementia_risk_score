#this script derives diagnoses of hypertension and high cholesterol

library(ggplot2)
library(tidyverse)
library(broom)
library(data.table)
library(dplyr)
library(tidyr)
library(cowplot)
library(knitr)
library(kableExtra)
library(MASS)
library(glmnet)
library(lm.beta)
library(extrafont)
library(pROC)
library(caret)
library(psych)
library(corrplot)

data_pathway = "../../raw_data/"
model_pathway = "../models/"
save_pathway = "../results/"

# load .rda file
load(file = paste0(data_pathway,"ukbdata_diagnoses_baseline_diseasestatus_baselinemedications_ANUADRI_CAIDE_FRS_recoded_DRS_fiftyplusnoapoe.rda"))
names(df)

#we'll need the self report fields
df_self_report_diagnoses <- read.csv(paste0(data_pathway,"add_vars_stroke_death.csv"), header=TRUE, sep=",", stringsAsFactors = FALSE)
cols_diagoses <- read.csv(paste0(data_pathway,"names_stroke_ICD.csv"), header=TRUE, sep=",", stringsAsFactors = FALSE)
df_self_report_diagnoses[df_self_report_diagnoses == "NA"] <- NA
df_self_report_diagnoses$eid <- as.character(df_self_report_diagnoses$eid)

#recode column names in self report df
newcols <- as.vector(cols_diagoses$col_name)
names(df_self_report_diagnoses) <- newcols
names(df_self_report_diagnoses)
#we only need the 224k at study
df_self_report_diagnoses <- list(df[,c("eid","Sex")], df_self_report_diagnoses) %>% reduce(left_join, by="eid")
names(df_self_report_diagnoses)
#self report at baseline is columns 38 - 71
names(df_self_report_diagnoses)[37]
names(df_self_report_diagnoses)[38]
names(df_self_report_diagnoses)[70]
names(df_self_report_diagnoses)[71]
names(df_self_report_diagnoses)[72]

#### self report ####
#make a self reported hypertension var
df_self_report_diagnoses$self_report_hypertension_0_0 <- apply(df_self_report_diagnoses[, 38:71], 1, function(x) {
  if(any(x %in% c("1065","1072", "1073"))) {
    return(1)
  } else {
    return(0)
  }
})
summary(as.factor(df_self_report_diagnoses$self_report_hypertension_0_0))
#0      1 
#156728  68034 
df_self_report_diagnoses$self_report_hypertension_0_0<-as.factor(df_self_report_diagnoses$self_report_hypertension_0_0)

#and cholesterol
df_self_report_diagnoses$self_report_cholesterol_0_0 <- apply(df_self_report_diagnoses[, 38:71], 1, function(x) {
  if(any(x %in% c("1473"))) {
    return(1)
  } else {
    return(0)
  }
})
summary(as.factor(df_self_report_diagnoses$self_report_cholesterol_0_0))
#0      1 
#191550  33212  
df_self_report_diagnoses$self_report_cholesterol_0_0<-as.factor(df_self_report_diagnoses$self_report_cholesterol_0_0)

#how many of these people are on meds

df <- list(df, df_self_report_diagnoses[,c("eid","self_report_hypertension_0_0","self_report_cholesterol_0_0")]) %>% reduce(left_join, by="eid")
summary(subset(df, df$self_report_hypertension_0_0==1)$Antihypertensive_meds_0_0)
#0     1 
#17276 50758  25% who sr hypertension dont sr hypertnsive meds

summary(subset(df, df$Antihypertensive_meds_0_0==1)$self_report_hypertension_0_0)
#0     1 
#8382 50758  #14% who sr no meds sr yes hypertension

summary(subset(df, df$self_report_cholesterol_0_0==1)$statins_0_0)
#0     1 
#4062 29150   12% who sr cholesterol dont sr statins

summary(subset(df, df$statins_0_0==1)$self_report_cholesterol_0_0)
#0     1 
#15468 29150 #34% who sr no statins sr yes cholesterol

#### primary care ####
datalist <-list()
for (ltc in c("hypertension","hyperlipidemia","hyperchol")){
  
  #load in the subjects with an identified pcare diagnosis. this is in two splits

  df_pcare<-read.csv(paste('../../raw_data/participants_with_',ltc,'.csv',sep=""))
  df_ltc <- df_pcare[,c("eid","event_dt_1")]
  
  df_ltc[paste("pcare_diagnosis_of",ltc,sep="_")] <- 1 #everyone here has a diagnosis
  df_ltc[paste("pcare_diagnosis_date_of",ltc,sep="_")] <- as.Date(df_ltc$event_dt_1, format="%d/%m/%Y") #format date
  df_ltc<-df_ltc[which(df_ltc[,paste("pcare_diagnosis_date_of",ltc,sep="_")]>as.Date("1950-01-01")),] #drop invalid dates
  df_ltc<-df_ltc[c("eid", sprintf('pcare_diagnosis_of_%s',ltc),sprintf('pcare_diagnosis_date_of_%s', ltc))]
  df_ltc$eid<-as.character(df_ltc$eid)
  datalist[[ltc]] <- df_ltc #add to list, merge later
}
#merge list into one df
ukb_ltcdiagnoses_pcare = Reduce(function(...) merge(..., all=T, by="eid"), datalist)

#change NAs to 0
for(ltc in c("hypertension","hyperlipidemia","hyperchol")){
  ukb_ltcdiagnoses_pcare[is.na(ukb_ltcdiagnoses_pcare[paste("pcare_diagnosis_of",ltc,sep="_")]),paste("pcare_diagnosis_of",ltc,sep="_")]<-0
}

names(ukb_ltcdiagnoses_pcare)

#now which are at baseline
ukb_ltcdiagnoses_pcare <- list(df[,c("eid","Date_of_assessment_0_0_new")], ukb_ltcdiagnoses_pcare) %>% reduce(inner_join, by="eid")

ukb_ltcdiagnoses_pcare$Date_of_assessment_0_0_new <- as.Date(ukb_ltcdiagnoses_pcare$Date_of_assessment_0_0_new, format="%Y-%m-%d")
ukb_ltcdiagnoses_pcare$pcare_diagnosis_date_of_hypertension <- as.Date(ukb_ltcdiagnoses_pcare$pcare_diagnosis_date_of_hypertension, format="%Y-%m-%d")
ukb_ltcdiagnoses_pcare$pcare_diagnosis_date_of_hyperlipidemia <- as.Date(ukb_ltcdiagnoses_pcare$pcare_diagnosis_date_of_hyperlipidemia, format="%Y-%m-%d")
ukb_ltcdiagnoses_pcare$pcare_diagnosis_date_of_hyperchol <- as.Date(ukb_ltcdiagnoses_pcare$pcare_diagnosis_date_of_hyperchol, format="%Y-%m-%d")

for (ltc in c("hypertension","hyperlipidemia","hyperchol")){
  ukb_ltcdiagnoses_pcare[paste("pcare_diagnosis_of", ltc, "0_0", sep="_")]<-as.factor(ifelse(is.na(
    (ukb_ltcdiagnoses_pcare[paste("pcare_diagnosis_date_of", ltc, sep="_")] - ukb_ltcdiagnoses_pcare["Date_of_assessment_0_0_new"])),0,
    ifelse((ukb_ltcdiagnoses_pcare[paste("pcare_diagnosis_date_of", ltc, sep="_")] - ukb_ltcdiagnoses_pcare["Date_of_assessment_0_0_new"]) <=0,1,0)))
}
summary(as.factor(ukb_ltcdiagnoses_pcare$pcare_diagnosis_of_hypertension))
#0     1 
#4056 15109 
summary(as.factor(ukb_ltcdiagnoses_pcare$pcare_diagnosis_of_hypertension_0_0))
#0     1 
#8237 10928 

df <- list(df, ukb_ltcdiagnoses_pcare) %>% reduce(left_join, by="eid")
#lets merge hyperchol and hyperlipi

df$pcare_diagnosis_of_cholesterol_0_0 <- as.factor(ifelse(df$pcare_diagnosis_of_hyperchol_0_0==1 |
                                                            df$pcare_diagnosis_of_hyperlipidemia_0_0==1,1,0))

summary(df$pcare_diagnosis_of_cholesterol_0_0)
#0      1   NA's 
# 14732   4433 205597 
summary(df$pcare_diagnosis_of_hyperlipidemia_0_0)
#     0      1   NA's 
#16111   3054 205597 
summary(df$pcare_diagnosis_of_hyperchol_0_0)
#0      1   NA's 
# 17599   1566 205597 

df$pcare_diagnosis_of_cholesterol_0_0[is.na(df$pcare_diagnosis_of_cholesterol_0_0)] <- 0
summary(df$pcare_diagnosis_of_cholesterol_0_0)
#0      1 
#220329   4433 

df$pcare_diagnosis_of_hypertension_0_0[is.na(df$pcare_diagnosis_of_hypertension_0_0)] <- 0


#### derive variable ####
#ok now lets make variable based on
#self report + med use + pcare
#### hypertension ####
summary(df$self_report_hypertension_0_0)
#0      1 
#156728  68034 
summary(df$Antihypertensive_meds_0_0)
#0      1 
#165622  59140 
summary(df$pcare_diagnosis_of_hypertension_0_0)
#0      1 
#213834  10928 

df$hypertensive <- as.factor(ifelse( df$self_report_hypertension_0_0==1 |
                                       df$Antihypertensive_meds_0_0==1 |
                                       df$pcare_diagnosis_of_hypertension_0_0==1,1,0))
summary(df$hypertensive)
#0      1 
#147849  76913 
summary(subset(df, df$hypertensive==0)$Systolic_BP_auto_mean)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#72.0   124.5   135.5   137.1   148.0   253.0 
summary(subset(df, df$hypertensive==1)$Systolic_BP_auto_mean)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#79.0   133.5   145.5   146.2   158.0   252.5 

#thats a big difference
summary(subset(df, df$hypertensive==0)$dementia_BIN_TOTAL)
#0      1 
#145870   1979 #11 %
summary(subset(df, df$hypertensive==1)$dementia_BIN_TOTAL)
#0     1 
#74937  1976 #20%

summary(subset(df, df$Antihypertensive_meds_0_0==0)$dementia_BIN_TOTAL)
#0      1 
#163328   2294 #12%
summary(subset(df, df$Antihypertensive_meds_0_0==1)$dementia_BIN_TOTAL)
#0     1 
#57479  1661 #22%
#similar, but stronger spread in hypertensive


#### cholesterol ####
summary(df$self_report_cholesterol_0_0)
#0      1 
#191550  33212 
summary(df$statins_0_0)
#0      1 
#180144  44618 
summary(df$pcare_diagnosis_of_cholesterol_0_0)
#0      1 
#220329   4433 

df$cholesterol <- as.factor(ifelse( df$self_report_cholesterol_0_0==1 |
                                      df$statins_0_0==1 |
                                      df$pcare_diagnosis_of_cholesterol_0_0==1,1,0))
summary(df$cholesterol)
#0      1 
#175013  49749 

summary(subset(df, df$cholesterol==0)$dementia_BIN_TOTAL)
#0      1 
#172558   2455 #1.4%
summary(subset(df, df$cholesterol==1)$dementia_BIN_TOTAL)
#0     1 
#48249  1500  #3%

#save
save(df, file=paste0(data_pathway,"ukbdata_diagnoses_baseline_diseasestatus_baselinemedications_ANUADRI_CAIDE_FRS_recoded_DRS_derivedhythchol_fiftyplusnoapoe.rda"))
