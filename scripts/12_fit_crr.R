#load required packages
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
library(survival)
library(cmprsk)

# specify data, model and results pathway
data_pathway = "../../raw_data/"
model_pathway = "../models/"
save_pathway = "../results/"

load(file= paste0(data_pathway, "11_train_data_outliers_removed_cox_fitted.rda"))
load(file= paste0(data_pathway, "11_test_data_outliers_removed_cox_fitted.rda"))

#update death records
df_death<-read.csv("../../raw_data/ukb50321_deathvars_refresh.csv", stringsAsFactors = FALSE)
names(df_death) #date is X40000.0.0
df_death <- df_death[,c("eid","X40000.0.0")]
names(df_death) <- c("eid","death_date_refresh")
df_death$eid <- as.character(df_death$eid)
df_death$death_date_refresh <- as.Date(df_death$death_date_refresh, format="%Y-%m-%d")
summary(df_death$death_date_refresh)
# Min.      1st Qu.       Median         Mean      3rd Qu.         Max.         NA's 
# "2006-05-10" "2014-04-18" "2017-06-22" "2016-11-21" "2019-12-21" "2021-11-12"     "464516" 
df_death$death_record_refresh <- ifelse( !is.na(df_death$death_date_refresh) &
                                           df_death$death_date_refresh<as.Date("2021-10-31",format="%Y-%m-%d"),1,0)
summary(as.factor(df_death$death_record_refresh))
# 0      1 
# 464532  37881 

#merge
train.data <- list(train.data, df_death) %>% reduce(left_join, by="eid")
test.data <- list(test.data, df_death) %>% reduce(left_join, by="eid")
rm(df_death)



#### competing risk ####
#apply competing risk adjustment to both models to correct coefficients

#we need a new failcode status to separate censor vs death vs dementia

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

crr_status <- train.data$crr_status
time_at_risk <- train.data$time_at_risk

#### UKBDRS_LASSO ####
UKBDRS_LASSO  <-  paste("Surv(time_at_risk, dementia_BIN_surv) ~  Age_when_attended_assesment_centre_0_0 +  family_history_of_dementia +
                            education_years + Diabetes_BIN_FINAL_0_0  + Townsend_deprivation_Groups_0_0 +
                            current_history_depression + stroke_TIA_BIN_FINAL +  
                            hypertensive + cholesterol + household_occupancy + Sex")
modelvars <- c("Age_when_attended_assesment_centre_0_0","family_history_of_dementia","education_years","Diabetes_BIN_FINAL_0_0",
            "Townsend_deprivation_Groups_0_0","current_history_depression","stroke_TIA_BIN_FINAL","hypertensive","cholesterol",
            "household_occupancy","Sex","dementia_BIN_surv")

ukbdrs.cr.train <- crr(ftime = time_at_risk, fstatus = crr_status, failcode=1, cencode = 0,
                       cov1=model.matrix(as.formula(UKBDRS_LASSO), train.data[modelvars])[,-1], variance=TRUE)
save(ukbdrs.cr.train, file=paste0(save_pathway,"ukbdrs.cr.train.rda"))
summary(ukbdrs.cr.train)
# coef exp(coef) se(coef)      z p-value
# Age_when_attended_assesment_centre_0_0  0.1760     1.192  0.00480 36.696 0.0e+00
# family_history_of_dementia1             0.4389     1.551  0.04294 10.221 0.0e+00
# education_years                        -0.0415     0.959  0.00648 -6.407 1.5e-10
# Diabetes_BIN_FINAL_0_01                 0.5177     1.678  0.05775  8.964 0.0e+00
# Townsend_deprivation_Groups_0_01       -0.0359     0.965  0.05983 -0.600 5.5e-01
# Townsend_deprivation_Groups_0_02        0.0188     1.019  0.05876  0.321 7.5e-01
# Townsend_deprivation_Groups_0_03        0.0438     1.045  0.05918  0.740 4.6e-01
# Townsend_deprivation_Groups_0_04        0.2402     1.272  0.05707  4.209 2.6e-05
# current_history_depression1             0.5517     1.736  0.04623 11.933 0.0e+00
# stroke_TIA_BIN_FINAL1                   0.6394     1.895  0.07784  8.215 2.2e-16
# hypertensive1                           0.1581     1.171  0.04093  3.862 1.1e-04
# cholesterol1                            0.1046     1.110  0.04509  2.319 2.0e-02
# household_occupancy1                    0.1204     1.128  0.04429  2.718 6.6e-03
# household_occupancy2                   -0.0553     0.946  0.05579 -0.990 3.2e-01
# Sex1                                    0.1656     1.180  0.03799  4.359 1.3e-05
# 
# exp(coef) exp(-coef)  2.5% 97.5%
#   Age_when_attended_assesment_centre_0_0     1.192      0.839 1.181 1.204
# family_history_of_dementia1                1.551      0.645 1.426 1.687
# education_years                            0.959      1.042 0.947 0.972
# Diabetes_BIN_FINAL_0_01                    1.678      0.596 1.499 1.879
# Townsend_deprivation_Groups_0_01           0.965      1.037 0.858 1.085
# Townsend_deprivation_Groups_0_02           1.019      0.981 0.908 1.143
# Townsend_deprivation_Groups_0_03           1.045      0.957 0.930 1.173
# Townsend_deprivation_Groups_0_04           1.272      0.786 1.137 1.422
# current_history_depression1                1.736      0.576 1.586 1.901
# stroke_TIA_BIN_FINAL1                      1.895      0.528 1.627 2.208
# hypertensive1                              1.171      0.854 1.081 1.269
# cholesterol1                               1.110      0.901 1.016 1.213
# household_occupancy1                       1.128      0.887 1.034 1.230
# household_occupancy2                       0.946      1.057 0.848 1.056
# Sex1                                       1.180      0.847 1.095 1.271
# 
# Num. cases = 176611
# Pseudo Log-likelihood = -34944 
# Pseudo likelihood ratio test = 3001  on 15 df,

#save the coefficients
ukbdrs_coefs <- cbind(summary(ukbdrs.cr.train)$coef,summary(ukbdrs.cr.train)$conf.int[,c(3,4)])
df_ukbdrs_coefficients <- as.data.frame(ukbdrs_coefs) #5 columns
write.csv(df_ukbdrs_coefficients, paste0(save_pathway,"crr_ukbdrs_coefficients.csv"))



#### UKBDRS_APOE_LASSO ####
UKBDRS_APOE_LASSO  <-  paste("Surv(time_at_risk, dementia_BIN_surv) ~  Age_when_attended_assesment_centre_0_0 +  family_history_of_dementia +
                            education_years + Diabetes_BIN_FINAL_0_0  + Townsend_deprivation_Groups_0_0 +
                            current_history_depression + stroke_TIA_BIN_FINAL +  
                            hypertensive + cholesterol + household_occupancy + Sex + APOE_genotype_bin")

apoe.train.data <- train.data[which(!is.na(train.data$APOE_genotype_bin)),]

crr_status <- apoe.train.data$crr_status
time_at_risk <- apoe.train.data$time_at_risk

ukbdrs.apoe.cr.train <- crr(ftime = time_at_risk, fstatus = crr_status,
                            failcode=1, cencode = 0,
                            cov1=model.matrix(as.formula(UKBDRS_APOE_LASSO), apoe.train.data[modelvars])[,-1], variance=TRUE)
save(ukbdrs.apoe.cr.train, file=paste0(save_pathway,"ukbdrs.apoe.cr.train.rda"))
summary(ukbdrs.apoe.cr.train)
# coef exp(coef) se(coef)       z p-value
# Age_when_attended_assesment_centre_0_0  0.1840     1.202  0.00558 32.9450 0.0e+00
# family_history_of_dementia1             0.3248     1.384  0.04913  6.6101 3.8e-11
# education_years                        -0.0383     0.962  0.00749 -5.1176 3.1e-07
# Diabetes_BIN_FINAL_0_01                 0.5295     1.698  0.06886  7.6899 1.5e-14
# Townsend_deprivation_Groups_0_01       -0.0767     0.926  0.06832 -1.1220 2.6e-01
# Townsend_deprivation_Groups_0_02       -0.0377     0.963  0.06776 -0.5564 5.8e-01
# Townsend_deprivation_Groups_0_03       -0.0040     0.996  0.06809 -0.0587 9.5e-01
# Townsend_deprivation_Groups_0_04        0.2359     1.266  0.06530  3.6125 3.0e-04
# current_history_depression1             0.5622     1.755  0.05385 10.4403 0.0e+00
# stroke_TIA_BIN_FINAL1                   0.6126     1.845  0.09303  6.5849 4.6e-11
# hypertensive1                           0.1900     1.209  0.04737  4.0108 6.1e-05
# cholesterol1                            0.0221     1.022  0.05195  0.4258 6.7e-01
# household_occupancy1                    0.1275     1.136  0.05143  2.4796 1.3e-02
# household_occupancy2                   -0.0108     0.989  0.06370 -0.1698 8.7e-01
# Sex1                                    0.1598     1.173  0.04368  3.6581 2.5e-04
# APOE_genotype_bin1                      1.1285     3.091  0.04242 26.6019 0.0e+00
# 
# exp(coef) exp(-coef)  2.5% 97.5%
#   Age_when_attended_assesment_centre_0_0     1.202      0.832 1.189 1.215
# family_history_of_dementia1                1.384      0.723 1.257 1.524
# education_years                            0.962      1.039 0.948 0.977
# Diabetes_BIN_FINAL_0_01                    1.698      0.589 1.484 1.944
# Townsend_deprivation_Groups_0_01           0.926      1.080 0.810 1.059
# Townsend_deprivation_Groups_0_02           0.963      1.038 0.843 1.100
# Townsend_deprivation_Groups_0_03           0.996      1.004 0.872 1.138
# Townsend_deprivation_Groups_0_04           1.266      0.790 1.114 1.439
# current_history_depression1                1.755      0.570 1.579 1.950
# stroke_TIA_BIN_FINAL1                      1.845      0.542 1.538 2.214
# hypertensive1                              1.209      0.827 1.102 1.327
# cholesterol1                               1.022      0.978 0.923 1.132
# household_occupancy1                       1.136      0.880 1.027 1.256
# household_occupancy2                       0.989      1.011 0.873 1.121
# Sex1                                       1.173      0.852 1.077 1.278
# APOE_genotype_bin1                         3.091      0.324 2.844 3.359
# 
# Num. cases = 125701
# Pseudo Log-likelihood = -24970 
# Pseudo likelihood ratio test = 2994  on 16 df,

#save the coefficients
ukbdrs_apoecoefs <- cbind(summary(ukbdrs.apoe.cr.train)$coef,summary(ukbdrs.apoe.cr.train)$conf.int[,c(3,4)])
df_ukbdrs_apoecoefficients <- as.data.frame(ukbdrs_apoecoefs) 
write.csv(df_ukbdrs_apoecoefficients, paste0(save_pathway,"crr_ukbdrs_apoe_coefficients.csv"))

#formatted coefs
x <- rbind(df_ukbdrs_coefficients, df_ukbdrs_apoecoefficients)
x$FDR_BH <- p.adjust(x$`p-value`, method="BH") 

df_table1_ukbdrs <- x %>% mutate_if(is.numeric, ~round(., 3))
write.csv(df_table1_ukbdrs, paste0(save_pathway,"crr_coefficients_formatted.csv"))
rm(apoe.train.data)




#### AGE ONLY ####
#compute cr coefficients for the baseline age only model
age_only <-      paste("Surv(time_at_risk, dementia_BIN_surv) ~  Age_when_attended_assesment_centre_0_0")

crr_status <- train.data$crr_status
time_at_risk <- train.data$time_at_risk

modelvars <- c("Age_when_attended_assesment_centre_0_0","dementia_BIN_surv")
covs<-model.matrix(as.formula(age_only), train.data[modelvars])

ageonly.cr.train <- crr(ftime = time_at_risk, fstatus = crr_status,
                            failcode=1, cencode = 0,
                            cov1=covs[,-1], variance=TRUE)
save(ageonly.cr.train, file=paste0(save_pathway,"ageonly.cr.train.rda"))
summary(ageonly.cr.train)
# Call:
#   crr(ftime = time_at_risk, fstatus = crr_status, cov1 = covs[, 
#                                                               -1], failcode = 1, cencode = 0, variance = TRUE)
# 
# coef exp(coef) se(coef)    z p-value
# covs[, -1]1 0.188      1.21  0.00446 42.2       0
# 
# exp(coef) exp(-coef) 2.5% 97.5%
#   covs[, -1]1      1.21      0.828  1.2  1.22
# 
# Num. cases = 176611
# Pseudo Log-likelihood = -35262 
# Pseudo likelihood ratio test = 2365  on 1 df,
ageonly.cr.train$coef
# covs[, -1]1 
# 0.1883235 


#### FIT LP ####
#using the cr adjusted coefficients, build the linear predictor and predicted probablity for each model

#convert categoricals to numerics to multiply by coefficients
train.data$UKBDRS_familyhistory <- ifelse(train.data$family_history_of_dementia==0,0,1)
train.data$UKBDRS_diabetes <- ifelse(train.data$Diabetes_BIN_FINAL_0_0==0,0,1)
train.data$UKBDRS_depression <- ifelse(train.data$current_history_depression==0,0,1)
train.data$UKBDRS_stroke <- ifelse(train.data$stroke_TIA_BIN_FINAL==0,0,1)
train.data$UKBDRS_hypertensive <- ifelse(train.data$hypertensive==0,0,1)
train.data$UKBDRS_cholesterol <- ifelse(train.data$cholesterol==0,0,1)
train.data$UKBDRS_townsend_group1 <- ifelse(train.data$Townsend_deprivation_Groups_0_0==0,1,0)
train.data$UKBDRS_townsend_group2 <- ifelse(train.data$Townsend_deprivation_Groups_0_0==1,1,0)
train.data$UKBDRS_townsend_group3 <- ifelse(train.data$Townsend_deprivation_Groups_0_0==2,1,0)
train.data$UKBDRS_townsend_group4 <- ifelse(train.data$Townsend_deprivation_Groups_0_0==3,1,0)
train.data$UKBDRS_townsend_group5 <- ifelse(train.data$Townsend_deprivation_Groups_0_0==4,1,0)
train.data$UKBDRS_livesalone <- ifelse(train.data$household_occupancy==1,1,0)
train.data$UKBDRS_liveswithmultiple <- ifelse(train.data$household_occupancy==2,1,0)
train.data$UKBDRS_APOE <- ifelse(train.data$APOE_genotype_bin==1,1,0)
train.data$UKBDRS_sex_score <- ifelse(train.data$Sex==0,0,1)


test.data$UKBDRS_familyhistory <- ifelse(test.data$family_history_of_dementia==0,0,1)
test.data$UKBDRS_diabetes <- ifelse(test.data$Diabetes_BIN_FINAL_0_0==0,0,1)
test.data$UKBDRS_depression <- ifelse(test.data$current_history_depression==0,0,1)
test.data$UKBDRS_stroke <- ifelse(test.data$stroke_TIA_BIN_FINAL==0,0,1)
test.data$UKBDRS_hypertensive <- ifelse(test.data$hypertensive==0,0,1)
test.data$UKBDRS_cholesterol <- ifelse(test.data$cholesterol==0,0,1)
test.data$UKBDRS_townsend_group1 <- ifelse(test.data$Townsend_deprivation_Groups_0_0==0,1,0)
test.data$UKBDRS_townsend_group2 <- ifelse(test.data$Townsend_deprivation_Groups_0_0==1,1,0)
test.data$UKBDRS_townsend_group3 <- ifelse(test.data$Townsend_deprivation_Groups_0_0==2,1,0)
test.data$UKBDRS_townsend_group4 <- ifelse(test.data$Townsend_deprivation_Groups_0_0==3,1,0)
test.data$UKBDRS_townsend_group5 <- ifelse(test.data$Townsend_deprivation_Groups_0_0==4,1,0)
test.data$UKBDRS_livesalone <- ifelse(test.data$household_occupancy==1,1,0)
test.data$UKBDRS_liveswithmultiple <- ifelse(test.data$household_occupancy==2,1,0)
test.data$UKBDRS_APOE <- ifelse(test.data$APOE_genotype_bin==1,1,0)
test.data$UKBDRS_sex_score <- ifelse(test.data$Sex==0,0,1)

#build linear predictor for each score using cr adjusted coeffients
#for continuous predictors, subtract the mean of the train data vals

#### UKBDRS LASSO ####
train.data$UKBDRS_LASSO_crr_linear_predictor <- 0.176031632032254*(train.data$Age_when_attended_assesment_centre_0_0 - mean(train.data$Age_when_attended_assesment_centre_0_0)) +
  0.438900573093536*train.data$UKBDRS_familyhistory + -0.0415438411765326*(train.data$education_years - mean(train.data$education_years)) +
  0.517727093577745*train.data$UKBDRS_diabetes + 
  -0.0358789897498946*train.data$UKBDRS_townsend_group2 + 0.0188491412350755*train.data$UKBDRS_townsend_group3 +
  0.0438160634033426*train.data$UKBDRS_townsend_group4 + 0.240197314505592*train.data$UKBDRS_townsend_group5 +
  0.551685657152675*train.data$UKBDRS_depression + 0.639444965849634*train.data$UKBDRS_stroke +
  0.158067938257923*train.data$UKBDRS_hypertensive + 0.104581315791995*train.data$UKBDRS_cholesterol +
  0.120392245428737*train.data$UKBDRS_livesalone + -0.0552531234441752*train.data$UKBDRS_liveswithmultiple +
  0.165615246946691*train.data$UKBDRS_sex_score
summary(train.data$UKBDRS_LASSO_crr_linear_predictor)  
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -2.3179 -0.4498  0.4574  0.4063  1.2480  4.2345 

train.data$UKBDRS_LASSO_crr_predicted_prob <- 1 - 0.9911458^exp(train.data$UKBDRS_LASSO_crr_linear_predictor)  
summary(train.data$UKBDRS_LASSO_crr_predicted_prob)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0008755 0.0056559 0.0139532 0.0226691 0.0305044 0.4587754 


test.data$UKBDRS_LASSO_crr_linear_predictor <- 0.176031632032254*(test.data$Age_when_attended_assesment_centre_0_0 - mean(train.data$Age_when_attended_assesment_centre_0_0)) +
  0.438900573093536*test.data$UKBDRS_familyhistory + -0.0415438411765326*(test.data$education_years - mean(train.data$education_years)) +
  0.517727093577745*test.data$UKBDRS_diabetes + 
  -0.0358789897498946*test.data$UKBDRS_townsend_group2 + 0.0188491412350755*test.data$UKBDRS_townsend_group3 +
  0.0438160634033426*test.data$UKBDRS_townsend_group4 + 0.240197314505592*test.data$UKBDRS_townsend_group5 +
  0.551685657152675*test.data$UKBDRS_depression + 0.639444965849634*test.data$UKBDRS_stroke +
  0.158067938257923*test.data$UKBDRS_hypertensive + 0.104581315791995*test.data$UKBDRS_cholesterol +
  0.120392245428737*test.data$UKBDRS_livesalone + -0.0552531234441752*test.data$UKBDRS_liveswithmultiple +
  0.165615246946691*test.data$UKBDRS_sex_score
summary(test.data$UKBDRS_LASSO_crr_linear_predictor)  
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -2.3642 -0.4465  0.4548  0.4093  1.2478  4.1141  

test.data$UKBDRS_LASSO_crr_predicted_prob <- 1 - 0.9911458^exp(test.data$UKBDRS_LASSO_crr_linear_predictor)  
summary(test.data$UKBDRS_LASSO_crr_predicted_prob)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0008359 0.0056744 0.0139174 0.0227658 0.0304991 0.4197438 


#### UKBDRS APOE ####
train.data$UKBDRS_APOE_LASSO_crr_linear_predictor <- 0.183970644979205*(train.data$Age_when_attended_assesment_centre_0_0 - mean(train.data$Age_when_attended_assesment_centre_0_0)) +
  0.32476842583806*train.data$UKBDRS_familyhistory + -0.0383496766100552*(train.data$education_years - mean(train.data$education_years)) +
  0.529530544064499*train.data$UKBDRS_diabetes + 
  -0.0766572689764079*train.data$UKBDRS_townsend_group2 + -0.0377024242266824*train.data$UKBDRS_townsend_group3 +
  -0.00399517369427498*train.data$UKBDRS_townsend_group4 + 0.235901509809486*train.data$UKBDRS_townsend_group5 +
  0.562230274804497*train.data$UKBDRS_depression + 0.612567091982896*train.data$UKBDRS_stroke +
  0.190005466407646*train.data$UKBDRS_hypertensive + 0.022123262331353*train.data$UKBDRS_cholesterol +
  0.127522987028888*train.data$UKBDRS_livesalone + -0.0108170960528888*train.data$UKBDRS_liveswithmultiple +
  0.159769970963938*train.data$UKBDRS_sex_score + 1.1284610152188*train.data$UKBDRS_APOE
summary(train.data$UKBDRS_APOE_LASSO_crr_linear_predictor)  
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   -2.29   -0.21    0.72    0.70    1.59    5.33   50910 

train.data$UKBDRS_APOE_LASSO_crr_predicted_prob <- 1 - 0.9941616^exp(train.data$UKBDRS_APOE_LASSO_crr_linear_predictor)  
summary(train.data$UKBDRS_APOE_LASSO_crr_predicted_prob)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.00    0.00    0.01    0.02    0.03    0.70   50910 


test.data$UKBDRS_APOE_LASSO_crr_linear_predictor <- 0.183970644979205*(test.data$Age_when_attended_assesment_centre_0_0 - mean(train.data$Age_when_attended_assesment_centre_0_0)) +
  0.32476842583806*test.data$UKBDRS_familyhistory + -0.0383496766100552*(test.data$education_years - mean(train.data$education_years)) +
  0.529530544064499*test.data$UKBDRS_diabetes + 
  -0.0766572689764079*test.data$UKBDRS_townsend_group2 + -0.0377024242266824*test.data$UKBDRS_townsend_group3 +
  -0.00399517369427498*test.data$UKBDRS_townsend_group4 + 0.235901509809486*test.data$UKBDRS_townsend_group5 +
  0.562230274804497*test.data$UKBDRS_depression + 0.612567091982896*test.data$UKBDRS_stroke +
  0.190005466407646*test.data$UKBDRS_hypertensive + 0.022123262331353*test.data$UKBDRS_cholesterol +
  0.127522987028888*test.data$UKBDRS_livesalone + -0.0108170960528888*test.data$UKBDRS_liveswithmultiple +
  0.159769970963938*test.data$UKBDRS_sex_score + 1.1284610152188*test.data$UKBDRS_APOE
summary(test.data$UKBDRS_APOE_LASSO_crr_linear_predictor)  
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  -2.393  -0.204   0.716   0.698   1.596   5.238   12762 

test.data$UKBDRS_APOE_LASSO_crr_predicted_prob <- 1 - 0.9941616^exp(test.data$UKBDRS_APOE_LASSO_crr_linear_predictor)  
summary(test.data$UKBDRS_APOE_LASSO_crr_predicted_prob)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.001   0.005   0.012   0.024   0.028   0.668   12762 

#### age only ####
train.data$age_only_crr_linear_predictor <- 0.1883235*(train.data$Age_when_attended_assesment_centre_0_0 - mean(train.data$Age_when_attended_assesment_centre_0_0))
train.data$age_only_crr_predicted_prob <- 1 - 0.9874905^exp(train.data$age_only_crr_linear_predictor)  
test.data$age_only_crr_linear_predictor <- 0.1883235*(test.data$Age_when_attended_assesment_centre_0_0 - mean(train.data$Age_when_attended_assesment_centre_0_0))
test.data$age_only_crr_predicted_prob <- 1 - 0.9874905^exp(test.data$age_only_crr_linear_predictor)  
summary(train.data$age_only_crr_predicted_prob)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.001924 0.005942 0.012580 0.019659 0.026530 0.136222 
summary(test.data$age_only_crr_predicted_prob)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.001924 0.005942 0.012580 0.019730 0.026530 0.095597 

#save train and test, containing lps and predicted probs, for future auc tests
save(train.data, file=paste0(data_pathway,"12_train_data_outliers_removed_cox_crr_fitted.rda"))
save(test.data, file=paste0(data_pathway,"12_test_data_outliers_removed_cox_crr_fitted.rda"))

