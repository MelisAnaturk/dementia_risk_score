#compute discrimination and calibration of models in specific age ranges
library(tidyverse)
library(broom)
library(rcompanion)
library(data.table)

library(survival)
library(cmprsk)
library(riskRegression)

library(ggplot2)
library(tidyverse)
library(broom)
library(data.table)
library(dplyr)
library(tidyr)


load(file="../../raw_data/modelvar/12_train_data_outliers_removed_cox_crr_fitted.rda")
load(file="../../raw_data/modelvar/12_test_data_outliers_removed_cox_crr_fitted.rda")
set.seed(9)


#define status
#dementia
df_dementia <- train.data[which(train.data$dementia_BIN_TOTAL==1),]
df_dementia$crr_status <- 1

#deaths
df_deaths <- train.data[which( (train.data$dementia_BIN_TOTAL==0) & (train.data$death_record_refresh==1) ),]
df_deaths$crr_status <- 2

#healthy
df_healthy <- train.data[which( (train.data$dementia_BIN_TOTAL==0) & (train.data$death_record_refresh==0) ),]
df_healthy$crr_status<-0

train.data <- rbind(df_dementia, df_healthy)
train.data<-rbind(train.data, df_deaths)
rm(df_dementia, df_deaths, df_healthy)


df_dementia <- test.data[which(test.data$dementia_BIN_TOTAL==1),]
df_dementia$crr_status <- 1

#deaths
df_deaths <- test.data[which( (test.data$dementia_BIN_TOTAL==0) & (test.data$death_record_refresh==1) ),]
df_deaths$crr_status <- 2

#healthy
df_healthy <- test.data[which( (test.data$dementia_BIN_TOTAL==0) & (test.data$death_record_refresh==0) ),]
df_healthy$crr_status<-0

test.data <- rbind(df_dementia, df_healthy)
test.data<-rbind(test.data, df_deaths)
rm(df_dementia, df_deaths, df_healthy)

train.data$y <- ifelse(train.data$dementia_BIN_TOTAL==1,1,0)
test.data$y <- ifelse(test.data$dementia_BIN_TOTAL==1,1,0)

train.data$dataset <- "train"
test.data$dataset <- "test"

df_test<-rbind(train.data, test.data)
rm(train.data, test.data)
#### CAIDE ####

#compute the auc of CAIDE and UKB-DRS in age specific subset of train and test data
#age range for CAIDE is 39-64

models <- c("age_only", "UKBDRS_LASSO", "CAIDE")

#subset data into caide age range
caide_data <- subset(df_test, Age_at_recruitment_0_0<65)

data <- subset(caide_data, dataset=="train") #n=133965
auc_test<-Score(list('UKBDRS'=data$UKBDRS_LASSO_crr_predicted_prob,
                     'DRS'=data$DRS_predicted_prob,
                     'Age_only'=data$age_only_crr_predicted_prob,
                     'CAIDE'=data$CAIDE_predicted_prob),
                formula=Hist(time_at_risk,crr_status)~1,
                data = data,
                null.model = FALSE,
                conf.int = TRUE,
                times = c(365.25*14.2),
                plots="ROC",
                metrics="AUC",
                cens.model = "cox",
                conservative=FALSE,
                censoring.save.memory=FALSE,
                contrasts = list(c(1,2,3,4)))
auc_test$AUC$score
# model   times       AUC           se     lower     upper
# 1:   UKBDRS 5186.55 0.7641147 0.0004465031 0.7632396 0.7649898
# 2:      DRS 5186.55 0.7214535 0.0004520721 0.7205674 0.7223395
# 3: Age_only 5186.55 0.7219440 0.0004285201 0.7211041 0.7227839
# 4:    CAIDE 5186.55 0.6279616 0.0005195334 0.6269433 0.6289799

data <- subset(caide_data, dataset=="test") #n=33476
auc_test<-Score(list('UKBDRS'=data$UKBDRS_LASSO_crr_predicted_prob,
                     'DRS'=data$DRS_predicted_prob,
                     'Age_only'=data$age_only_crr_predicted_prob,
                     'CAIDE'=data$CAIDE_predicted_prob),
                formula=Hist(time_at_risk,crr_status)~1,
                data = data,
                null.model = FALSE,
                conf.int = TRUE,
                times = c(365.25*14.2),
                plots="ROC",
                metrics="AUC",
                cens.model = "cox",
                conservative=FALSE,
                censoring.save.memory=FALSE,
                contrasts = list(c(1,2,3,4)))
auc_test$AUC$score
# model   times       AUC         se     lower     upper
# 1:   UKBDRS 5186.55 0.7704222 0.01578605 0.7394821 0.8013623
# 2:      DRS 5186.55 0.7142099 0.02194836 0.6711919 0.7572279
# 3: Age_only 5186.55 0.7182549 0.01987973 0.6792913 0.7572184
# 4:    CAIDE 5186.55 0.5965240 0.02737856 0.5428630 0.6501850

rm(caide_data)

#### DRS ####

#compute the auc of DRS and UKB-DRS in age specific subset of train and test data
#age range for DRS is 60-79

#subset data into drs age range
drs_data <- subset(df_test, (Age_at_recruitment_0_0>=60 & Age_at_recruitment_0_0<=79))
summary(as.factor(drs_data$dataset))
# test train 
# 24794 98697 

data <- subset(drs_data, dataset=="train") 
auc_test<-Score(list('UKBDRS'=data$UKBDRS_LASSO_crr_predicted_prob,
                     'DRS'=data$DRS_predicted_prob,
                     'Age_only'=data$age_only_crr_predicted_prob,
                     'CAIDE'=data$CAIDE_predicted_prob),
                formula=Hist(time_at_risk,crr_status)~1,
                data = data,
                null.model = FALSE,
                conf.int = TRUE,
                times = c(365.25*14.2),
                plots="ROC",
                metrics="AUC",
                cens.model = "cox",
                conservative=FALSE,
                censoring.save.memory=FALSE,
                contrasts = list(c(1,2,3,4)))
auc_test$AUC$score
# model   times       AUC          se     lower     upper
# 1:   UKBDRS 5186.55 0.6953488 0.001084213 0.6932238 0.6974738
# 2:      DRS 5186.55 0.6478113 0.001106159 0.6456433 0.6499793
# 3: Age_only 5186.55 0.6469872 0.001102189 0.6448270 0.6491475
# 4:    CAIDE 5186.55 0.5547889 0.001143404 0.5525479 0.5570299

data <- subset(drs_data, dataset=="test") #n=33476
auc_test<-Score(list('UKBDRS'=data$UKBDRS_LASSO_crr_predicted_prob,
                     'DRS'=data$DRS_predicted_prob,
                     'Age_only'=data$age_only_crr_predicted_prob,
                     'CAIDE'=data$CAIDE_predicted_prob),
                formula=Hist(time_at_risk,crr_status)~1,
                data = data,
                null.model = FALSE,
                conf.int = TRUE,
                times = c(365.25*14.2),
                plots="ROC",
                metrics="AUC",
                cens.model = "cox",
                conservative=FALSE,
                censoring.save.memory=FALSE,
                contrasts = list(c(1,2,3,4)))
auc_test$AUC$score
# model   times       AUC         se     lower     upper
# 1:   UKBDRS 5186.55 0.7098931 0.01670061 0.6771605 0.7426257
# 2:      DRS 5186.55 0.6590201 0.01743919 0.6248399 0.6932003
# 3: Age_only 5186.55 0.6767813 0.01619058 0.6450483 0.7085142
# 4:    CAIDE 5186.55 0.5466845 0.01726683 0.5128421 0.5805269

rm(drs_data)
gc()


#### ANU ADRI (MAP) ####

#compute the auc of DRS and UKB-DRS in age specific subset of train and test data
#age range for DRS is 60-79

#subset data into anu map age range
anu_data <- subset(df_test, Age_at_recruitment_0_0>=54 & !is.na(df_test$ANU_ADRI))
summary(as.factor(anu_data$dataset))
# test  train 
# 36912 147277 


data <- subset(anu_data, dataset=="train") 
auc_test<-Score(list('UKBDRS'=data$UKBDRS_LASSO_crr_predicted_prob,
                     'anu'=data$ANU_ADRI),
                formula=Hist(time_at_risk,crr_status)~1,
                data = data,
                null.model = FALSE,
                conf.int = TRUE,
                times = c(365.25*14.2),
                plots="ROC",
                metrics="AUC",
                cens.model = "cox",
                conservative=FALSE,
                censoring.save.memory=FALSE)
auc_test$AUC$score
# model   times       AUC           se     lower     upper
# 1: UKBDRS 5186.55 0.7585621 0.0002230804 0.7581248 0.7589993
# 2:    anu 5186.55 0.5666151 0.0002570924 0.5661112 0.5671190

data <- subset(anu_data, dataset=="test") 
auc_test<-Score(list('UKBDRS'=data$UKBDRS_LASSO_crr_predicted_prob,
                     'anu'=data$ANU_ADRI),
                formula=Hist(time_at_risk,crr_status)~1,
                data = data,
                null.model = FALSE,
                conf.int = TRUE,
                times = c(365.25*14.2),
                plots="ROC",
                metrics="AUC",
                cens.model = "cox",
                conservative=FALSE,
                censoring.save.memory=FALSE)
auc_test$AUC$score
# model   times       AUC         se     lower     upper
# 1: UKBDRS 5186.55 0.7685242 0.01383943 0.7413994 0.7956490
# 2:    anu 5186.55 0.5516272 0.01759415 0.5171433 0.5861111


rm(anu_data)



#### ANU ADRI (CVHS) ####

#subset data into age range
anu_data <- subset(df_test, Age_at_recruitment_0_0>=62 & !is.na(df_test$ANU_ADRI))
summary(as.factor(anu_data$dataset))
# test train 
# 18948 75401 

data <- subset(anu_data, dataset=="train") 
auc_test<-Score(list('UKBDRS'=data$UKBDRS_LASSO_crr_predicted_prob,
                     'anu'=data$ANU_ADRI),
                formula=Hist(time_at_risk,crr_status)~1,
                data = data,
                null.model = FALSE,
                conf.int = TRUE,
                times = c(365.25*14.2),
                plots="ROC",
                metrics="AUC",
                cens.model = "cox",
                conservative=FALSE,
                censoring.save.memory=FALSE)
auc_test$AUC$score
# model   times       AUC          se     lower     upper
# 1: UKBDRS 5186.55 0.6660904 0.002588970 0.6610161 0.6711647
# 2:    anu 5186.55 0.5783159 0.002390164 0.5736313 0.5830005

data <- subset(anu_data, dataset=="test") 
auc_test<-Score(list('UKBDRS'=data$UKBDRS_LASSO_crr_predicted_prob,
                     'anu'=data$ANU_ADRI),
                formula=Hist(time_at_risk,crr_status)~1,
                data = data,
                null.model = FALSE,
                conf.int = TRUE,
                times = c(365.25*14.2),
                plots="ROC",
                metrics="AUC",
                cens.model = "cox",
                conservative=FALSE,
                censoring.save.memory=FALSE)
auc_test$AUC$score
# model   times       AUC         se     lower     upper
# 1: UKBDRS 5186.55 0.6707152 0.02120909 0.6291461 0.7122842
# 2:    anu 5186.55 0.5704525 0.02124239 0.5288182 0.6120868


rm(anu_data)
