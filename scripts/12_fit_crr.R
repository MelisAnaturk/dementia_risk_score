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
                            education_years + Diabetes_BIN_FINAL_0_0  + Townsend_deprivation_Groups_0_0 +
                            current_history_depression + stroke_TIA_BIN_FINAL +  
                            hypertensive + cholesterol + household_occupancy + Sex")
#*crr is diff formula but we only need right side to build model matrix, so can repeat cox form here

modelvars <- c("Age_when_attended_assesment_centre_0_0","family_history_of_dementia","education_years","Diabetes_BIN_FINAL_0_0",
            "Townsend_deprivation_Groups_0_0","current_history_depression","stroke_TIA_BIN_FINAL","hypertensive","cholesterol",
            "household_occupancy","Sex","dementia_BIN_surv")

ukbdrs.cr.train <- crr(ftime = time_at_risk, fstatus = crr_status, failcode=1, cencode = 0,
                       cov1=model.matrix(as.formula(UKBDRS_LASSO), train.data[modelvars])[,-1], variance=TRUE)
save(ukbdrs.cr.train, file=paste0(save_pathway,"ukbdrs.cr.train.rda"))
summary(ukbdrs.cr.train)
# coef exp(coef) se(coef)      z p-value
# Age_when_attended_assesment_centre_0_0  0.1767     1.193  0.00481 36.733 0.0e+00
# family_history_of_dementia1             0.4307     1.538  0.04293 10.032 0.0e+00
# education_years                        -0.0405     0.960  0.00646 -6.269 3.6e-10
# Diabetes_BIN_FINAL_0_01                 0.5376     1.712  0.05762  9.329 0.0e+00
# Townsend_deprivation_Groups_0_01       -0.0330     0.968  0.05981 -0.551 5.8e-01
# Townsend_deprivation_Groups_0_02        0.0104     1.010  0.05874  0.176 8.6e-01
# Townsend_deprivation_Groups_0_03        0.0491     1.050  0.05917  0.830 4.1e-01
# Townsend_deprivation_Groups_0_04        0.2382     1.269  0.05706  4.174 3.0e-05
# current_history_depression1             0.5543     1.741  0.04618 12.002 0.0e+00
# stroke_TIA_BIN_FINAL1                   0.6546     1.924  0.07781  8.413 0.0e+00
# hypertensive1                           0.1585     1.172  0.04111  3.855 1.2e-04
# cholesterol1                            0.1031     1.109  0.04529  2.276 2.3e-02
# household_occupancy1                    0.1251     1.133  0.04433  2.822 4.8e-03
# household_occupancy2                   -0.0615     0.940  0.05575 -1.103 2.7e-01
# Sex1                                    0.1712     1.187  0.03790  4.517 6.3e-06
# 
# exp(coef) exp(-coef)  2.5% 97.5%
#   Age_when_attended_assesment_centre_0_0     1.193      0.838 1.182 1.205
# family_history_of_dementia1                1.538      0.650 1.414 1.673
# education_years                            0.960      1.041 0.948 0.973
# Diabetes_BIN_FINAL_0_01                    1.712      0.584 1.529 1.916
# Townsend_deprivation_Groups_0_01           0.968      1.034 0.861 1.088
# Townsend_deprivation_Groups_0_02           1.010      0.990 0.901 1.134
# Townsend_deprivation_Groups_0_03           1.050      0.952 0.935 1.179
# Townsend_deprivation_Groups_0_04           1.269      0.788 1.135 1.419
# current_history_depression1                1.741      0.574 1.590 1.906
# stroke_TIA_BIN_FINAL1                      1.924      0.520 1.652 2.241
# hypertensive1                              1.172      0.853 1.081 1.270
# cholesterol1                               1.109      0.902 1.014 1.211
# household_occupancy1                       1.133      0.882 1.039 1.236
# household_occupancy2                       0.940      1.063 0.843 1.049
# Sex1                                       1.187      0.843 1.102 1.278
# 
# Num. cases = 176611
# Pseudo Log-likelihood = -34922 
# Pseudo likelihood ratio test = 3023  on 15 df,
#save the coefficients
ukbdrs_coefs <- cbind(summary(ukbdrs.cr.train)$coef,summary(ukbdrs.cr.train)$conf.int[,c(3,4)])
df_ukbdrs_coefficients <- as.data.frame(ukbdrs_coefs) #5 columns
write.csv(df_ukbdrs_coefficients, paste0(save_pathway,"crr_ukbdrs_coefficients.csv"))



#### UKBDRS_APOE_LASSO ####
UKBDRS_APOE_LASSO  <-  paste("Surv(time_at_risk, dementia_BIN_surv) ~  Age_when_attended_assesment_centre_0_0 +  family_history_of_dementia +
                            education_years + Diabetes_BIN_FINAL_0_0  + Townsend_deprivation_Groups_0_0 +
                            current_history_depression + stroke_TIA_BIN_FINAL +  
                            hypertensive + cholesterol + household_occupancy + Sex + APOE_genotype_bin")

modelvars <- c("Age_when_attended_assesment_centre_0_0","family_history_of_dementia","education_years","Diabetes_BIN_FINAL_0_0",
               "Townsend_deprivation_Groups_0_0","current_history_depression","stroke_TIA_BIN_FINAL","hypertensive","cholesterol",
               "household_occupancy","Sex","APOE_genotype_bin","dementia_BIN_surv")

apoe.train.data <- train.data[which(!is.na(train.data$APOE_genotype_bin)),]

crr_status <- apoe.train.data$crr_status
time_at_risk <- apoe.train.data$time_at_risk

ukbdrs.apoe.cr.train <- crr(ftime = time_at_risk, fstatus = crr_status,
                            failcode=1, cencode = 0,
                            cov1=model.matrix(as.formula(UKBDRS_APOE_LASSO), apoe.train.data[modelvars])[,-1], variance=TRUE)
save(ukbdrs.apoe.cr.train, file=paste0(save_pathway,"ukbdrs.apoe.cr.train.rda"))
summary(ukbdrs.apoe.cr.train)
# coef exp(coef) se(coef)      z p-value
# Age_when_attended_assesment_centre_0_0  0.18532     1.204  0.00561 33.011 0.0e+00
# family_history_of_dementia1             0.31118     1.365  0.04920  6.324 2.5e-10
# education_years                        -0.03769     0.963  0.00746 -5.054 4.3e-07
# Diabetes_BIN_FINAL_0_01                 0.52603     1.692  0.06873  7.654 2.0e-14
# Townsend_deprivation_Groups_0_01       -0.07897     0.924  0.06831 -1.156 2.5e-01
# Townsend_deprivation_Groups_0_02       -0.05246     0.949  0.06772 -0.775 4.4e-01
# Townsend_deprivation_Groups_0_03       -0.00443     0.996  0.06806 -0.065 9.5e-01
# Townsend_deprivation_Groups_0_04        0.21368     1.238  0.06541  3.267 1.1e-03
# current_history_depression1             0.56648     1.762  0.05373 10.543 0.0e+00
# stroke_TIA_BIN_FINAL1                   0.64312     1.902  0.09271  6.937 4.0e-12
# hypertensive1                           0.18874     1.208  0.04738  3.983 6.8e-05
# cholesterol1                            0.02750     1.028  0.05195  0.529 6.0e-01
# household_occupancy1                    0.13371     1.143  0.05151  2.596 9.4e-03
# household_occupancy2                   -0.00922     0.991  0.06369 -0.145 8.8e-01
# Sex1                                    0.16451     1.179  0.04358  3.775 1.6e-04
# APOE_genotype_bin1                      1.12876     3.092  0.04244 26.599 0.0e+00
# 
# exp(coef) exp(-coef)  2.5% 97.5%
#   Age_when_attended_assesment_centre_0_0     1.204      0.831 1.190 1.217
# family_history_of_dementia1                1.365      0.733 1.240 1.503
# education_years                            0.963      1.038 0.949 0.977
# Diabetes_BIN_FINAL_0_01                    1.692      0.591 1.479 1.936
# Townsend_deprivation_Groups_0_01           0.924      1.082 0.808 1.056
# Townsend_deprivation_Groups_0_02           0.949      1.054 0.831 1.084
# Townsend_deprivation_Groups_0_03           0.996      1.004 0.871 1.138
# Townsend_deprivation_Groups_0_04           1.238      0.808 1.089 1.408
# current_history_depression1                1.762      0.568 1.586 1.958
# stroke_TIA_BIN_FINAL1                      1.902      0.526 1.586 2.281
# hypertensive1                              1.208      0.828 1.101 1.325
# cholesterol1                               1.028      0.973 0.928 1.138
# household_occupancy1                       1.143      0.875 1.033 1.264
# household_occupancy2                       0.991      1.009 0.875 1.123
# Sex1                                       1.179      0.848 1.082 1.284
# APOE_genotype_bin1                         3.092      0.323 2.845 3.360
# 
# Num. cases = 125844
# Pseudo Log-likelihood = -24957 
# Pseudo likelihood ratio test = 3011  on 16 df,

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
train.data$UKBDRS_LASSO_crr_linear_predictor <- 0.176669418489078*(train.data$Age_when_attended_assesment_centre_0_0 - mean(train.data$Age_when_attended_assesment_centre_0_0)) +
  0.430676926035003*train.data$UKBDRS_familyhistory + -0.0405038840167207*(train.data$education_years - mean(train.data$education_years)) +
  0.537562983032388*train.data$UKBDRS_diabetes + 
  -0.0329769531560162*train.data$UKBDRS_townsend_group2 + 0.0103579985413302*train.data$UKBDRS_townsend_group3 +
  0.0491156700850136*train.data$UKBDRS_townsend_group4 + 0.238156169012864*train.data$UKBDRS_townsend_group5 +
  0.554267780232139*train.data$UKBDRS_depression + 0.654638552618771*train.data$UKBDRS_stroke +
  0.158474158854365*train.data$UKBDRS_hypertensive + 0.10307003740093*train.data$UKBDRS_cholesterol +
  0.125104378688588*train.data$UKBDRS_livesalone + -0.0614962253749352*train.data$UKBDRS_liveswithmultiple +
  0.171202344829797*train.data$UKBDRS_sex_score
summary(train.data$UKBDRS_LASSO_crr_linear_predictor)  
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -2.3511 -0.4497  0.4570  0.4076  1.2534  4.2815

train.data$UKBDRS_LASSO_crr_predicted_prob <- 1 - 0.9910584^exp(train.data$UKBDRS_LASSO_crr_linear_predictor)  
summary(train.data$UKBDRS_LASSO_crr_predicted_prob)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0008553 0.0057122 0.0140844 0.0230130 0.0309656 0.4778577 


test.data$UKBDRS_LASSO_crr_linear_predictor <- 0.176669418489078*(test.data$Age_when_attended_assesment_centre_0_0 - mean(train.data$Age_when_attended_assesment_centre_0_0)) +
  0.430676926035003*test.data$UKBDRS_familyhistory + -0.0405038840167207*(test.data$education_years - mean(train.data$education_years)) +
  0.537562983032388*test.data$UKBDRS_diabetes + 
  -0.0329769531560162*test.data$UKBDRS_townsend_group2 + 0.0103579985413302*test.data$UKBDRS_townsend_group3 +
  0.0491156700850136*test.data$UKBDRS_townsend_group4 + 0.238156169012864*test.data$UKBDRS_townsend_group5 +
  0.554267780232139*test.data$UKBDRS_depression + 0.654638552618771*test.data$UKBDRS_stroke +
  0.158474158854365*test.data$UKBDRS_hypertensive + 0.10307003740093*test.data$UKBDRS_cholesterol +
  0.125104378688588*test.data$UKBDRS_livesalone + -0.0614962253749352*test.data$UKBDRS_liveswithmultiple +
  0.171202344829797*test.data$UKBDRS_sex_score
summary(test.data$UKBDRS_LASSO_crr_linear_predictor)  
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -2.1414 -0.4490  0.4625  0.4132  1.2544  4.2410 

test.data$UKBDRS_LASSO_crr_predicted_prob <- 1 - 0.9910584^exp(test.data$UKBDRS_LASSO_crr_linear_predictor)  
summary(test.data$UKBDRS_LASSO_crr_predicted_prob)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.001055 0.005717 0.014162 0.023220 0.030996 0.464214 



#### UKBDRS APOE ####
train.data$UKBDRS_APOE_LASSO_crr_linear_predictor <- 0.185318600297277*(train.data$Age_when_attended_assesment_centre_0_0 - mean(train.data$Age_when_attended_assesment_centre_0_0)) +
  0.311180078767657*train.data$UKBDRS_familyhistory + -0.0376858334633324*(train.data$education_years - mean(train.data$education_years)) +
  0.526030705230127*train.data$UKBDRS_diabetes + 
  -0.0789680050990133*train.data$UKBDRS_townsend_group2 + -0.0524563711914478*train.data$UKBDRS_townsend_group3 +
  -0.00442555339500229*train.data$UKBDRS_townsend_group4 + 0.213675130036749*train.data$UKBDRS_townsend_group5 +
  0.566480952398541*train.data$UKBDRS_depression + 0.64312025662005*train.data$UKBDRS_stroke +
  0.188735504751471*train.data$UKBDRS_hypertensive + 0.0274992857607639*train.data$UKBDRS_cholesterol +
  0.133706416811094*train.data$UKBDRS_livesalone + -0.00922351693031058*train.data$UKBDRS_liveswithmultiple +
  0.164505689768259*train.data$UKBDRS_sex_score + 1.12876295214559*train.data$UKBDRS_APOE
summary(train.data$UKBDRS_APOE_LASSO_crr_linear_predictor)  
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   -2.39   -0.22    0.71    0.69    1.59    4.99   50767

train.data$UKBDRS_APOE_LASSO_crr_predicted_prob <- 1 - 0.9941051^exp(train.data$UKBDRS_APOE_LASSO_crr_linear_predictor)  
summary(train.data$UKBDRS_APOE_LASSO_crr_predicted_prob)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#    0.00    0.00    0.01    0.02    0.03    0.58   50767 


test.data$UKBDRS_APOE_LASSO_crr_linear_predictor <- 0.185318600297277*(test.data$Age_when_attended_assesment_centre_0_0 - mean(train.data$Age_when_attended_assesment_centre_0_0)) +
  0.311180078767657*test.data$UKBDRS_familyhistory + -0.0376858334633324*(test.data$education_years - mean(train.data$education_years)) +
  0.526030705230127*test.data$UKBDRS_diabetes + 
  -0.0789680050990133*test.data$UKBDRS_townsend_group2 + -0.0524563711914478*test.data$UKBDRS_townsend_group3 +
  -0.00442555339500229*test.data$UKBDRS_townsend_group4 + 0.213675130036749*test.data$UKBDRS_townsend_group5 +
  0.566480952398541*test.data$UKBDRS_depression + 0.64312025662005*test.data$UKBDRS_stroke +
  0.188735504751471*test.data$UKBDRS_hypertensive + 0.0274992857607639*test.data$UKBDRS_cholesterol +
  0.133706416811094*test.data$UKBDRS_livesalone + -0.00922351693031058*test.data$UKBDRS_liveswithmultiple +
  0.164505689768259*test.data$UKBDRS_sex_score + 1.12876295214559*test.data$UKBDRS_APOE
summary(test.data$UKBDRS_APOE_LASSO_crr_linear_predictor)  
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  -2.211  -0.207   0.723   0.700   1.589   5.362   12905 

test.data$UKBDRS_APOE_LASSO_crr_predicted_prob <- 1 - 0.9941051^exp(test.data$UKBDRS_APOE_LASSO_crr_linear_predictor)  
summary(test.data$UKBDRS_APOE_LASSO_crr_predicted_prob)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   0.001   0.005   0.012   0.024   0.029   0.716   12905

#### age only ####
train.data$age_only_crr_linear_predictor <- 0.1891088*(train.data$Age_when_attended_assesment_centre_0_0 - mean(train.data$Age_when_attended_assesment_centre_0_0))
train.data$age_only_crr_predicted_prob <- 1 - 0.9846693^exp(train.data$age_only_crr_linear_predictor)  
test.data$age_only_crr_linear_predictor <- 0.1891088*(test.data$Age_when_attended_assesment_centre_0_0 - mean(train.data$Age_when_attended_assesment_centre_0_0))
test.data$age_only_crr_predicted_prob <- 1 - 0.9846693^exp(test.data$age_only_crr_linear_predictor)  
summary(train.data$age_only_crr_predicted_prob)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.002342 0.007267 0.015420 0.024098 0.032569 0.139566 
summary(test.data$age_only_crr_predicted_prob)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.002342 0.007267 0.015420 0.024267 0.032569 0.166074 

#save train and test, containing lps and predicted probs, for future auc tests
save(train.data, file=paste0(data_pathway,"12_train_data_outliers_removed_cox_crr_fitted.rda"))
save(test.data, file=paste0(data_pathway,"12_test_data_outliers_removed_cox_crr_fitted.rda"))

