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
data_pathway = "../../raw_data/modelvar/"
model_pathway = "../models/"
save_pathway = "../results/"

load(file= paste0(data_pathway, "11_train_data_outliers_removed_cox_fitted.rda"))
load(file= paste0(data_pathway, "11_test_data_outliers_removed_cox_fitted.rda"))



#### competing risk ####
#apply competing risk adjustment to both models to correct coefficients

#we need a new failcode status to separate censor vs death vs dementia

#dementia
df_dementia <- train.data[which(train.data$dementia_BIN_TOTAL==1),]
df_dementia$crr_status <- 1

#deaths
df_deaths <- train.data[which( (train.data$dementia_BIN_TOTAL==0) & (train.data$has_death_record==1) ),]
df_deaths$crr_status <- 2

#healthy
df_healthy <- train.data[which( (train.data$dementia_BIN_TOTAL==0) & (train.data$has_death_record==0) ),]
df_healthy$crr_status<-0

train.data <- rbind(df_dementia, df_healthy)
train.data<-rbind(train.data, df_deaths)
rm(df_dementia, df_deaths, df_healthy)

crr_status <- train.data$crr_status
time_at_risk <- train.data$time_at_risk

#### UKBDRS_LASSO ####
UKBDRS_LASSO  <-  paste("Surv(time_at_risk, dementia_BIN_surv) ~  Age_when_attended_assesment_centre_0_0 +  family_history_of_dementia +
                            education_years + Diabetes_BIN_FINAL_0_0  + Townsend_deprivation_modelvar +
                            current_history_depression + stroke_TIA_BIN_FINAL +  
                            hypertensive + cholesterol + livesalone + Sex")
#*crr is diff formula but we only need right side to build model matrix, so can repeat cox form here

modelvars <- c("Age_when_attended_assesment_centre_0_0","family_history_of_dementia","education_years","Diabetes_BIN_FINAL_0_0",
            "Townsend_deprivation_modelvar","current_history_depression","stroke_TIA_BIN_FINAL","hypertensive","cholesterol",
            "livesalone","Sex","dementia_BIN_surv")

ukbdrs.cr.train <- crr(ftime = time_at_risk, fstatus = crr_status, failcode=1, cencode = 0,
                       cov1=model.matrix(as.formula(UKBDRS_LASSO), train.data[modelvars])[,-1], variance=TRUE)
save(ukbdrs.cr.train, file=paste0(save_pathway,"ukbdrs.cr.train.rda"))
summary(ukbdrs.cr.train)
# coef exp(coef) se(coef)     z p-value
# Age_when_attended_assesment_centre_0_0  0.178      1.19  0.00470 37.80 0.0e+00
# family_history_of_dementia1             0.431      1.54  0.04293 10.04 0.0e+00
# education_years                        -0.041      0.96  0.00645 -6.35 2.1e-10
# Diabetes_BIN_FINAL_0_01                 0.536      1.71  0.05756  9.32 0.0e+00
# Townsend_deprivation_modelvar1          0.228      1.26  0.04324  5.26 1.4e-07
# current_history_depression1             0.556      1.74  0.04618 12.04 0.0e+00
# stroke_TIA_BIN_FINAL1                   0.655      1.92  0.07784  8.41 0.0e+00
# hypertensive1                           0.159      1.17  0.04113  3.87 1.1e-04
# cholesterol1                            0.104      1.11  0.04528  2.30 2.2e-02
# livesalone1                             0.141      1.15  0.04313  3.27 1.1e-03
# Sex1                                    0.169      1.18  0.03788  4.45 8.4e-06
# 
# exp(coef) exp(-coef)  2.5% 97.5%
#   Age_when_attended_assesment_centre_0_0      1.19      0.837 1.184 1.206
# family_history_of_dementia1                 1.54      0.650 1.415 1.674
# education_years                             0.96      1.042 0.948 0.972
# Diabetes_BIN_FINAL_0_01                     1.71      0.585 1.528 1.914
# Townsend_deprivation_modelvar1              1.26      0.796 1.153 1.367
# current_history_depression1                 1.74      0.573 1.593 1.909
# stroke_TIA_BIN_FINAL1                       1.92      0.520 1.652 2.242
# hypertensive1                               1.17      0.853 1.082 1.271
# cholesterol1                                1.11      0.901 1.015 1.213
# livesalone1                                 1.15      0.869 1.058 1.253
# Sex1                                        1.18      0.845 1.099 1.275
# 
# Num. cases = 176611
# Pseudo Log-likelihood = -34924 
# Pseudo likelihood ratio test = 3020  on 11 df,
#save the coefficients
ukbdrs_coefs <- cbind(summary(ukbdrs.cr.train)$coef,summary(ukbdrs.cr.train)$conf.int[,c(3,4)])
df_ukbdrs_coefficients <- as.data.frame(ukbdrs_coefs) #5 columns
write.csv(df_ukbdrs_coefficients, paste0(save_pathway,"crr_ukbdrs_coefficients.csv"))


#get baseline cif
modelvars <- c("Age_when_attended_assesment_centre_0_0","family_history_of_dementia","education_years","Diabetes_BIN_FINAL_0_0",
               "Townsend_deprivation_modelvar","current_history_depression","stroke_TIA_BIN_FINAL","hypertensive","cholesterol",
               "livesalone","Sex")

df_base <- train.data[1,modelvars]
df_base$Age_when_attended_assesment_centre_0_0=mean(train.data$Age_when_attended_assesment_centre_0_0)
df_base$family_history_of_dementia <-0
df_base$Diabetes_BIN_FINAL_0_0<-0
df_base$current_history_depression<-0
df_base$stroke_TIA_BIN_FINAL<-0
df_base$hypertensive<-0
df_base$cholesterol<-0
df_base$Townsend_deprivation_modelvar<-0
df_base$education_years<-mean(train.data$education_years)
df_base$Sex<-0
df_base$livesalone<-0

cif_baseline <- predict(ukbdrs.cr.train, df_base)
#what time do we want?
cif_baseline[2000:2060,1] #2053 is 13.99 years
cif_baseline[2053,]
#5.117000e+03 8.3804896e-03
#baseline cif for 14 years is 0.008304896
baseline_surv <- 1 - cif_baseline[2053,2]
baseline_surv #0.9916195



#### UKBDRS_APOE_LASSO ####
UKBDRS_APOE_LASSO  <-  paste("Surv(time_at_risk, dementia_BIN_surv) ~  Age_when_attended_assesment_centre_0_0 +  family_history_of_dementia +
                            education_years + Diabetes_BIN_FINAL_0_0  + Townsend_deprivation_modelvar +
                            current_history_depression + stroke_TIA_BIN_FINAL +  
                            hypertensive + cholesterol + livesalone + Sex + APOE_genotype_bin")

modelvars <- c("Age_when_attended_assesment_centre_0_0","family_history_of_dementia","education_years","Diabetes_BIN_FINAL_0_0",
               "Townsend_deprivation_modelvar","current_history_depression","stroke_TIA_BIN_FINAL","hypertensive","cholesterol",
               "livesalone","Sex","APOE_genotype_bin","dementia_BIN_surv")

apoe.train.data <- train.data[which(!is.na(train.data$APOE_genotype_bin)),]

crr_status <- apoe.train.data$crr_status
time_at_risk <- apoe.train.data$time_at_risk

ukbdrs.apoe.cr.train <- crr(ftime = time_at_risk, fstatus = crr_status,
                            failcode=1, cencode = 0,
                            cov1=model.matrix(as.formula(UKBDRS_APOE_LASSO), apoe.train.data[modelvars])[,-1], variance=TRUE)
save(ukbdrs.apoe.cr.train, file=paste0(save_pathway,"ukbdrs.apoe.cr.train.rda"))
summary(ukbdrs.apoe.cr.train)
# coef exp(coef) se(coef)      z p-value
# Age_when_attended_assesment_centre_0_0  0.1854     1.204  0.00548 33.833 0.0e+00
# family_history_of_dementia1             0.3113     1.365  0.04920  6.327 2.5e-10
# education_years                        -0.0376     0.963  0.00745 -5.048 4.5e-07
# Diabetes_BIN_FINAL_0_01                 0.5260     1.692  0.06858  7.670 1.7e-14
# Townsend_deprivation_modelvar1          0.2470     1.280  0.05012  4.928 8.3e-07
# current_history_depression1             0.5669     1.763  0.05376 10.545 0.0e+00
# stroke_TIA_BIN_FINAL1                   0.6427     1.902  0.09271  6.932 4.1e-12
# hypertensive1                           0.1895     1.209  0.04738  4.000 6.3e-05
# cholesterol1                            0.0271     1.028  0.05192  0.523 6.0e-01
# livesalone1                             0.1366     1.146  0.05004  2.729 6.4e-03
# Sex1                                    0.1644     1.179  0.04355  3.775 1.6e-04
# APOE_genotype_bin1                      1.1287     3.091  0.04242 26.609 0.0e+00
# 
# exp(coef) exp(-coef)  2.5% 97.5%
#   Age_when_attended_assesment_centre_0_0     1.204      0.831 1.191 1.217
# family_history_of_dementia1                1.365      0.733 1.240 1.503
# education_years                            0.963      1.038 0.949 0.977
# Diabetes_BIN_FINAL_0_01                    1.692      0.591 1.479 1.936
# Townsend_deprivation_modelvar1             1.280      0.781 1.160 1.412
# current_history_depression1                1.763      0.567 1.587 1.959
# stroke_TIA_BIN_FINAL1                      1.902      0.526 1.586 2.280
# hypertensive1                              1.209      0.827 1.101 1.326
# cholesterol1                               1.028      0.973 0.928 1.138
# livesalone1                                1.146      0.872 1.039 1.264
# Sex1                                       1.179      0.848 1.082 1.284
# APOE_genotype_bin1                         3.091      0.323 2.845 3.359
# 
# Num. cases = 125844
# Pseudo Log-likelihood = -24958 
# Pseudo likelihood ratio test = 3009  on 12 df,
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


#get baseline cif
modelvars <- c("Age_when_attended_assesment_centre_0_0","family_history_of_dementia","education_years","Diabetes_BIN_FINAL_0_0",
               "Townsend_deprivation_modelvar","current_history_depression","stroke_TIA_BIN_FINAL","hypertensive","cholesterol",
               "livesalone","Sex", "APOE_genotype_bin")

df_base <- apoe.train.data[1,modelvars]
df_base$Age_when_attended_assesment_centre_0_0=mean(train.data$Age_when_attended_assesment_centre_0_0)
df_base$family_history_of_dementia <-0
df_base$Diabetes_BIN_FINAL_0_0<-0
df_base$current_history_depression<-0
df_base$stroke_TIA_BIN_FINAL<-0
df_base$hypertensive<-0
df_base$cholesterol<-0
df_base$Townsend_deprivation_modelvar<-0
df_base$education_years<-mean(train.data$education_years)
df_base$Sex<-0
df_base$livesalone<-0
df_base$APOE_genotype_bin<-0

cif_baseline <- predict(ukbdrs.apoe.cr.train, df_base)
#what time do we want?
cif_baseline[1660:1677,1] #1673 is 13.99 years
cif_baseline[1673,]
#5.112000e+03 5.447631e-03
#baseline cif for 14 years is 0.005447631
baseline_surv <- 1 - cif_baseline[1673,2]
baseline_surv #0.9945524




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
# covs[, -1]1 0.189      1.21  0.00447 42.3       0
# 
# exp(coef) exp(-coef) 2.5% 97.5%
#   covs[, -1]1      1.21      0.828  1.2  1.22
# 
# Num. cases = 176611
# Pseudo Log-likelihood = -35244 
# Pseudo likelihood ratio test = 2380  on 1 df,
ageonly.cr.train$coef
# covs[, -1]1 
# 0.1891088 

#get baseline cif
modelvars <- c("Age_when_attended_assesment_centre_0_0", "dementia_BIN_surv")
covs<-model.matrix(as.formula(age_only), train.data[modelvars])

covs[1,2]<-mean(train.data$Age_when_attended_assesment_centre_0_0)
covs[1,]

cif_baseline <- predict(ageonly.cr.train, covs[1,][-1])
#what time do we want?
cif_baseline[2000:2060,1] #2053 is 13.99 years
cif_baseline[2053,]
#5.117000e+03 1.381633e-02
#baseline cif for 14 years is 0.01381633
baseline_surv <- 1 - cif_baseline[2053,2]
baseline_surv #0.9861837


#### FIT LP ####
#using the cr adjusted coefficients, build the linear predictor and predicted probablity for each model

#convert categoricals to numerics to multiply by coefficients
train.data$UKBDRS_APOE <- ifelse(train.data$APOE_genotype_bin==1,1,0)
test.data$UKBDRS_APOE <- ifelse(test.data$APOE_genotype_bin==1,1,0)


#build linear predictor for each score using cr adjusted coeffients
#for continuous predictors, subtract the mean of the train data vals

#### UKBDRS LASSO ####
train.data$UKBDRS_LASSO_crr_linear_predictor <- 0.177692959251294*(train.data$Age_when_attended_assesment_centre_0_0 - mean(train.data$Age_when_attended_assesment_centre_0_0)) +
  0.430976053237681*train.data$UKBDRS_familyhistory + -0.0409925728007542*(train.data$education_years - mean(train.data$education_years)) +
  0.536470125918681*train.data$UKBDRS_diabetes + 
  0.227535362193283*train.data$UKBDRS_townsend_group5 +
  0.556016474679774*train.data$UKBDRS_depression + 0.654852204288696*train.data$UKBDRS_stroke +
  0.159181761673448*train.data$UKBDRS_hypertensive + 0.103934423846803*train.data$UKBDRS_cholesterol +
  0.140912692178426*train.data$UKBDRS_livesalone +
  0.168726275413058*train.data$UKBDRS_sex_score
summary(train.data$UKBDRS_LASSO_crr_linear_predictor)  
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -2.3051 -0.4333  0.4651  0.4189  1.2613  4.2976 

train.data$UKBDRS_LASSO_crr_predicted_prob <- 1 - 0.9916195^exp(train.data$UKBDRS_LASSO_crr_linear_predictor)  
summary(train.data$UKBDRS_LASSO_crr_predicted_prob)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0008391 0.0054418 0.0133106 0.0217947 0.0292720 0.4613910 


test.data$UKBDRS_LASSO_crr_linear_predictor <- 0.177692959251294*(test.data$Age_when_attended_assesment_centre_0_0 - mean(train.data$Age_when_attended_assesment_centre_0_0)) +
  0.430976053237681*test.data$UKBDRS_familyhistory + -0.0409925728007542*(test.data$education_years - mean(train.data$education_years)) +
  0.536470125918681*test.data$UKBDRS_diabetes + 
  0.227535362193283*test.data$UKBDRS_townsend_group5 +
  0.556016474679774*test.data$UKBDRS_depression + 0.654852204288696*test.data$UKBDRS_stroke +
  0.159181761673448*test.data$UKBDRS_hypertensive + 0.103934423846803*test.data$UKBDRS_cholesterol +
  0.140912692178426*test.data$UKBDRS_livesalone +
  0.168726275413058*test.data$UKBDRS_sex_score
summary(test.data$UKBDRS_LASSO_crr_linear_predictor)  
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -2.1543 -0.4334  0.4651  0.4245  1.2613  4.2566 

test.data$UKBDRS_LASSO_crr_predicted_prob <- 1 - 0.9916195^exp(test.data$UKBDRS_LASSO_crr_linear_predictor)  
summary(test.data$UKBDRS_LASSO_crr_predicted_prob)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0009756 0.0054411 0.0133106 0.0219884 0.0292720 0.4478378 



#### UKBDRS APOE ####
train.data$UKBDRS_APOE_LASSO_crr_linear_predictor <- 0.185448583498348*(train.data$Age_when_attended_assesment_centre_0_0 - mean(train.data$Age_when_attended_assesment_centre_0_0)) +
  0.311285747731805*train.data$UKBDRS_familyhistory + -0.037579578137673*(train.data$education_years - mean(train.data$education_years)) +
  0.526047157958713*train.data$UKBDRS_diabetes +
  0.246972315824607*train.data$UKBDRS_townsend_group5 +
  0.566892403366833*train.data$UKBDRS_depression + 0.64267917239979*train.data$UKBDRS_stroke +
  0.189525433227791*train.data$UKBDRS_hypertensive + 0.0271391344013689*train.data$UKBDRS_cholesterol +
  0.136570007308125*train.data$UKBDRS_livesalone +
  0.164374628844355*train.data$UKBDRS_sex_score + 1.12865113064495*train.data$UKBDRS_APOE
summary(train.data$UKBDRS_APOE_LASSO_crr_linear_predictor)  


train.data$UKBDRS_APOE_LASSO_crr_predicted_prob <- 1 - 0.9945524^exp(train.data$UKBDRS_APOE_LASSO_crr_linear_predictor)  
summary(train.data$UKBDRS_APOE_LASSO_crr_predicted_prob)
 

test.data$UKBDRS_APOE_LASSO_crr_linear_predictor <- 0.185448583498348*(test.data$Age_when_attended_assesment_centre_0_0 - mean(train.data$Age_when_attended_assesment_centre_0_0)) +
  0.311285747731805*test.data$UKBDRS_familyhistory + -0.037579578137673*(test.data$education_years - mean(train.data$education_years)) +
  0.526047157958713*test.data$UKBDRS_diabetes +
  0.246972315824607*test.data$UKBDRS_townsend_group5 +
  0.566892403366833*test.data$UKBDRS_depression + 0.64267917239979*test.data$UKBDRS_stroke +
  0.189525433227791*test.data$UKBDRS_hypertensive + 0.0271391344013689*test.data$UKBDRS_cholesterol +
  0.136570007308125*test.data$UKBDRS_livesalone +
  0.164374628844355*test.data$UKBDRS_sex_score + 1.12865113064495*test.data$UKBDRS_APOE
summary(test.data$UKBDRS_APOE_LASSO_crr_linear_predictor)  


test.data$UKBDRS_APOE_LASSO_crr_predicted_prob <- 1 - 0.9945524^exp(test.data$UKBDRS_APOE_LASSO_crr_linear_predictor)  
summary(test.data$UKBDRS_APOE_LASSO_crr_predicted_prob)


#### age only ####
train.data$age_only_crr_linear_predictor <- 0.1891088*(train.data$Age_when_attended_assesment_centre_0_0 - mean(train.data$Age_when_attended_assesment_centre_0_0))
train.data$age_only_crr_predicted_prob <- 1 - 0.9861837^exp(train.data$age_only_crr_linear_predictor)  
test.data$age_only_crr_linear_predictor <- 0.1891088*(test.data$Age_when_attended_assesment_centre_0_0 - mean(train.data$Age_when_attended_assesment_centre_0_0))
test.data$age_only_crr_predicted_prob <- 1 - 0.9861837^exp(test.data$age_only_crr_linear_predictor)  
summary(train.data$age_only_crr_predicted_prob)
summary(test.data$age_only_crr_predicted_prob)

#save train and test, containing lps and predicted probs, for future auc tests
save(train.data, file=paste0(data_pathway,"12_train_data_outliers_removed_cox_crr_fitted.rda"))
save(test.data, file=paste0(data_pathway,"12_test_data_outliers_removed_cox_crr_fitted.rda"))

