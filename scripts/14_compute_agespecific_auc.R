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

#subset data into caide age range
caide_data <- subset(df_test, Age_at_recruitment_0_0<65)
summary(as.factor(caide_data$dataset))
# test  train 
# 33476 133965 

data <- subset(caide_data, dataset=="train") 
auc_test<-Score(list('UKBDRS'=data$UKBDRS_LASSO_crr_predicted_prob,
                     'CAIDE'=data$CAIDE_predicted_prob),
                formula=Hist(time_at_risk,crr_status)~1,
                data = data,
                null.model = FALSE,
                conf.int = TRUE,
                times = c(365.25*14),
                plots="ROC",
                metrics="AUC",
                cens.model = "cox",
                conservative=FALSE,
                censoring.save.memory=FALSE,
                contrasts = list(c(1,2)))
auc_test$AUC$score
# model  times       AUC           se     lower     upper
# 1: UKBDRS 5113.5 0.7659735 0.0003962542 0.7651969 0.7667502
# 2:  CAIDE 5113.5 0.6084587 0.0004486838 0.6075793 0.6093381

data <- subset(caide_data, dataset=="test") 
auc_test<-Score(list('UKBDRS'=data$UKBDRS_LASSO_crr_predicted_prob,
                     'CAIDE'=data$CAIDE_predicted_prob),
                formula=Hist(time_at_risk,crr_status)~1,
                data = data,
                null.model = FALSE,
                conf.int = TRUE,
                times = c(365.25*14),
                plots="ROC",
                metrics="AUC",
                cens.model = "cox",
                conservative=FALSE,
                censoring.save.memory=FALSE,
                contrasts = list(c(1,2)))
auc_test$AUC$score
# model  times       AUC         se     lower     upper
# 1: UKBDRS 5113.5 0.7726198 0.01539741 0.7424415 0.8027982
# 2:  CAIDE 5113.5 0.6233200 0.01876640 0.5865385 0.6601014

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
                     'DRS'=data$DRS_predicted_prob),
                formula=Hist(time_at_risk,crr_status)~1,
                data = data,
                null.model = FALSE,
                conf.int = TRUE,
                times = c(365.25*14),
                plots="ROC",
                metrics="AUC",
                cens.model = "cox",
                conservative=FALSE,
                censoring.save.memory=FALSE,
                contrasts = list(c(1,2)))
auc_test$AUC$score
# model  times       AUC           se     lower     upper
# 1: UKBDRS 5113.5 0.6934775 0.0008575114 0.6917968 0.6951582
# 2:    DRS 5113.5 0.6520244 0.0008892761 0.6502814 0.6537673

data <- subset(drs_data, dataset=="test") 
auc_test<-Score(list('UKBDRS'=data$UKBDRS_LASSO_crr_predicted_prob,
                     'DRS'=data$DRS_predicted_prob),
                formula=Hist(time_at_risk,crr_status)~1,
                data = data,
                null.model = FALSE,
                conf.int = TRUE,
                times = c(365.25*14),
                plots="ROC",
                metrics="AUC",
                cens.model = "cox",
                conservative=FALSE,
                censoring.save.memory=FALSE,
                contrasts = list(c(1,2)))
auc_test$AUC$score
# model  times       AUC         se     lower     upper
# 1: UKBDRS 5113.5 0.7132067 0.01262436 0.6884634 0.7379500
# 2:    DRS 5113.5 0.6673618 0.01409002 0.6397459 0.6949778

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
                times = c(365.25*14),
                plots="ROC",
                metrics="AUC",
                cens.model = "cox",
                conservative=FALSE,
                censoring.save.memory=FALSE)
auc_test$AUC$score
# model  times      AUC           se     lower     upper
# 1: UKBDRS 5113.5 0.758261 0.0001405946 0.7579855 0.7585366
# 2:    anu 5113.5 0.562036 0.0001555521 0.5617312 0.5623409

data <- subset(anu_data, dataset=="test") 
auc_test<-Score(list('UKBDRS'=data$UKBDRS_LASSO_crr_predicted_prob,
                     'anu'=data$ANU_ADRI),
                formula=Hist(time_at_risk,crr_status)~1,
                data = data,
                null.model = FALSE,
                conf.int = TRUE,
                times = c(365.25*14),
                plots="ROC",
                metrics="AUC",
                cens.model = "cox",
                conservative=FALSE,
                censoring.save.memory=FALSE)
auc_test$AUC$score
# model  times       AUC          se     lower     upper
# 1: UKBDRS 5113.5 0.7711868 0.009912036 0.7517596 0.7906140
# 2:    anu 5113.5 0.5683208 0.013890332 0.5410962 0.5955453


rm(anu_data)



#### ANU ADRI (CVHS) ####

#subset data into age range
anu_data <- subset(df_test, Age_at_recruitment_0_0>=60 & !is.na(df_test$ANU_ADRI))
summary(as.factor(anu_data$dataset))
# test train 
# 24794 98697 

data <- subset(anu_data, dataset=="train") 
auc_test<-Score(list('UKBDRS'=data$UKBDRS_LASSO_crr_predicted_prob,
                     'anu'=data$ANU_ADRI),
                formula=Hist(time_at_risk,crr_status)~1,
                data = data,
                null.model = FALSE,
                conf.int = TRUE,
                times = c(365.25*14),
                plots="ROC",
                metrics="AUC",
                cens.model = "cox",
                conservative=FALSE,
                censoring.save.memory=FALSE)
auc_test$AUC$score
# model  times       AUC           se     lower     upper
# 1: UKBDRS 5113.5 0.6934775 0.0008575114 0.6917968 0.6951582
# 2:    anu 5113.5 0.5817251 0.0008861294 0.5799883 0.5834619

data <- subset(anu_data, dataset=="test") 
auc_test<-Score(list('UKBDRS'=data$UKBDRS_LASSO_crr_predicted_prob,
                     'anu'=data$ANU_ADRI),
                formula=Hist(time_at_risk,crr_status)~1,
                data = data,
                null.model = FALSE,
                conf.int = TRUE,
                times = c(365.25*14),
                plots="ROC",
                metrics="AUC",
                cens.model = "cox",
                conservative=FALSE,
                censoring.save.memory=FALSE)
auc_test$AUC$score
# model  times       AUC         se     lower     upper
# 1: UKBDRS 5113.5 0.6846833 0.01441278 0.6564348 0.7129319
# 2:    anu 5113.5 0.5817918 0.01677286 0.5489176 0.6146660


rm(anu_data)
