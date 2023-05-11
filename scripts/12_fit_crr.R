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
                            cov1=covs[,-1], variance=FALSE)
save(ageonly.cr.train, file=paste0(save_pathway,"ageonly.cr.train.rda"))
summary(ageonly.cr.train)
# Call:
#   crr(ftime = time_at_risk, fstatus = crr_status, cov1 = covs[, 
#                                                               -1], failcode = 1, cencode = 0, variance = FALSE)
# 
# coef exp(coef) se(coef)  z p-value
# covs[, -1]1 0.189      1.21       NA NA      NA
# 
# exp(coef) exp(-coef) 2.5% 97.5%
#   covs[, -1]1      1.21      0.828   NA    NA
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
train.data$age_only_crr_linear_predictor <- 0.1891088*(train.data$Age_when_attended_assesment_centre_0_0 - mean(train.data$Age_when_attended_assesment_centre_0_0))
train.data$age_only_crr_predicted_prob <- 1 - 0.9846693^exp(train.data$age_only_crr_linear_predictor)  
test.data$age_only_crr_linear_predictor <- 0.1891088*(test.data$Age_when_attended_assesment_centre_0_0 - mean(train.data$Age_when_attended_assesment_centre_0_0))
test.data$age_only_crr_predicted_prob <- 1 - 0.9846693^exp(test.data$age_only_crr_linear_predictor)  
summary(train.data$age_only_crr_predicted_prob)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.001924 0.005942 0.012580 0.019659 0.026530 0.136222 
summary(test.data$age_only_crr_predicted_prob)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.001924 0.005942 0.012580 0.019730 0.026530 0.095597 

#save train and test, containing lps and predicted probs, for future auc tests
save(train.data, file=paste0(data_pathway,"12_train_data_outliers_removed_cox_crr_fitted.rda"))
save(test.data, file=paste0(data_pathway,"12_test_data_outliers_removed_cox_crr_fitted.rda"))

