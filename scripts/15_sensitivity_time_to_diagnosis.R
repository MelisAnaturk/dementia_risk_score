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


#lets look at performance in predicting
#dementia in 1-5 years
#dementia 5-10 years
#dementia 10-14

#early
nrow(test.data[which( (test.data$dementia_BIN_TOTAL==1) & (test.data$time_at_risk<=(362.25*5))   ),])
#99
auc_test<-Score(list('UKBDRS'=test.data$UKBDRS_LASSO_crr_predicted_prob,
                     'DRS'=test.data$DRS_predicted_prob,
                     'Age_only'=test.data$age_only_crr_predicted_prob,
                     'CAIDE'=test.data$CAIDE_predicted_prob),
                formula=Hist(time_at_risk,crr_status)~1,
                data = test.data,
                null.model = FALSE,
                conf.int = TRUE,
                times = c(365.25*5),
                plots="ROC",
                metrics="AUC",
                cens.model = "cox",
                conservative=FALSE,
                censoring.save.memory=FALSE,
                contrasts = list(c(1,2,3,4)))
auc_test$AUC$score
# model   times       AUC         se     lower     upper
# 1:   UKBDRS 1826.25 0.7539711 0.02331617 0.7082722 0.7996699
# 2:      DRS 1826.25 0.7323078 0.02296532 0.6872966 0.7773190
# 3: Age_only 1826.25 0.7204349 0.02365019 0.6740814 0.7667884
# 4:    CAIDE 1826.25 0.5592098 0.02640792 0.5074513 0.6109684


#anu
anu_test <- test.data[which(!is.na(test.data$ANU_ADRI)),]
auc_test<-Score(list('ANU'=anu_test$ANU_ADRI),
                formula=Hist(time_at_risk,crr_status)~1,
                data = anu_test,
                null.model = FALSE,
                conf.int = TRUE,
                times = c(365.25*5),
                plots="ROC",
                metrics="AUC",
                cens.model = "cox",
                conservative=FALSE,
                censoring.save.memory=FALSE)
auc_test$AUC$score
# model   times       AUC         se     lower     upper
# 1:   ANU 1826.25 0.5355414 0.02792355 0.4808123 0.5902706


#mid
nrow(test.data[which( (test.data$dementia_BIN_TOTAL==1) & (test.data$time_at_risk<=(362.25*10))   ),])
#441
auc_test<-Score(list('UKBDRS'=test.data$UKBDRS_LASSO_crr_predicted_prob,
                     'DRS'=test.data$DRS_predicted_prob,
                     'Age_only'=test.data$age_only_crr_predicted_prob,
                     'CAIDE'=test.data$CAIDE_predicted_prob),
                formula=Hist(time_at_risk,crr_status)~1,
                data = test.data,
                null.model = FALSE,
                conf.int = TRUE,
                times = c(365.25*10),
                plots="ROC",
                metrics="AUC",
                cens.model = "cox",
                conservative=FALSE,
                censoring.save.memory=FALSE,
                contrasts = list(c(1,2,3,4)))
auc_test$AUC$score
# model  times       AUC         se     lower     upper
# 1:   UKBDRS 3652.5 0.7578182 0.01062555 0.7369925 0.7786439
# 2:      DRS 3652.5 0.7313246 0.01097558 0.7098128 0.7528363
# 3: Age_only 3652.5 0.7327285 0.01096800 0.7112316 0.7542254
# 4:    CAIDE 3652.5 0.5792560 0.01251398 0.5547290 0.6037829


#anu
auc_test<-Score(list('ANU'=anu_test$ANU_ADRI),
                formula=Hist(time_at_risk,crr_status)~1,
                data = anu_test,
                null.model = FALSE,
                conf.int = TRUE,
                times = c(365.25*10),
                plots="ROC",
                metrics="AUC",
                cens.model = "cox",
                conservative=FALSE,
                censoring.save.memory=FALSE)
auc_test$AUC$score
# model  times       AUC         se     lower     upper
# 1:   ANU 3652.5 0.5776875 0.01318589 0.5518436 0.6035314

