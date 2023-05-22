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

# specify data, model and results pathway
data_pathway = "../../raw_data/modelvar/"
model_pathway = "../models/"
save_pathway = "../results/"

load(file= paste0(data_pathway, "train_data_outliers_removed_postlasso.rda"))
load(file= paste0(data_pathway, "test_data_outliers_removed_postlasso.rda"))


#----- 1. COX REGRESSION ------------------------------------------
# cox regression is now run to calculate the beta-weights for each of the 
#components in our risk score (note, we also add in age, sex and education)

#### 1.1 Test beta coefficients ####
#based on the lasso selected vars, compute model coefficients in train data
#use 2 models - one with apoe, one without. 

#### UKBDRS_LASSO ####
UKBDRS_LASSO  <-  paste("Surv(time_at_risk, dementia_BIN_surv) ~  Age_when_attended_assesment_centre_0_0 +  family_history_of_dementia +
                            education_years + Diabetes_BIN_FINAL_0_0  + Townsend_deprivation_Groups_0_0 +
                            current_history_depression + stroke_TIA_BIN_FINAL +  
                            hypertensive + cholesterol + household_occupancy + Sex")

ukbdrs.cox <- coxph(as.formula(UKBDRS_LASSO), data = train.data)
summary(ukbdrs.cox)

#only most deprived sig
#only living alone sig

train.data$Townsend_deprivation_modelvar<-ifelse(train.data$Townsend_deprivation_Groups_0_0==4,1,0)
train.data$Townsend_deprivation_modelvar<-as.factor(train.data$Townsend_deprivation_modelvar)
test.data$Townsend_deprivation_modelvar<-ifelse(test.data$Townsend_deprivation_Groups_0_0==4,1,0)
test.data$Townsend_deprivation_modelvar<-as.factor(test.data$Townsend_deprivation_modelvar)

train.data$livesalone <- ifelse(train.data$household_occupancy==1,1,0)
train.data$livesalone <- as.factor(train.data$livesalone)
test.data$livesalone <- ifelse(test.data$household_occupancy==1,1,0)
test.data$livesalone <- as.factor(test.data$livesalone)

UKBDRS_LASSO  <-  paste("Surv(time_at_risk, dementia_BIN_surv) ~  Age_when_attended_assesment_centre_0_0 +  family_history_of_dementia +
                            education_years + Diabetes_BIN_FINAL_0_0  + Townsend_deprivation_modelvar +
                            current_history_depression + stroke_TIA_BIN_FINAL +  
                            hypertensive + cholesterol + livesalone + Sex")

ukbdrs.cox <- coxph(as.formula(UKBDRS_LASSO), data = train.data)
summary(ukbdrs.cox)

# coxph(formula = as.formula(UKBDRS_LASSO), data = train.data)
# 
# n= 176611, number of events= 3051 
# 
# coef exp(coef)  se(coef)      z Pr(>|z|)    
# Age_when_attended_assesment_centre_0_0  0.183109  1.200946  0.004532 40.400  < 2e-16 ***
#   family_history_of_dementia1             0.424164  1.528312  0.042978  9.869  < 2e-16 ***
#   education_years                        -0.043349  0.957577  0.006183 -7.011 2.36e-12 ***
#   Diabetes_BIN_FINAL_0_01                 0.596603  1.815939  0.057535 10.369  < 2e-16 ***
#   Townsend_deprivation_modelvar1          0.259727  1.296576  0.043365  5.989 2.11e-09 ***
#   current_history_depression1             0.571791  1.771437  0.046180 12.382  < 2e-16 ***
#   stroke_TIA_BIN_FINAL1                   0.717527  2.049358  0.076593  9.368  < 2e-16 ***
#   hypertensive1                           0.176537  1.193079  0.040813  4.326 1.52e-05 ***
#   cholesterol1                            0.103220  1.108736  0.044515  2.319 0.020408 *  
#   livesalone1                             0.167627  1.182495  0.043227  3.878 0.000105 ***
#   Sex1                                    0.215354  1.240300  0.037756  5.704 1.17e-08 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# exp(coef) exp(-coef) lower .95 upper .95
# Age_when_attended_assesment_centre_0_0    1.2009     0.8327     1.190    1.2117
# family_history_of_dementia1               1.5283     0.6543     1.405    1.6626
# education_years                           0.9576     1.0443     0.946    0.9693
# Diabetes_BIN_FINAL_0_01                   1.8159     0.5507     1.622    2.0327
# Townsend_deprivation_modelvar1            1.2966     0.7713     1.191    1.4116
# current_history_depression1               1.7714     0.5645     1.618    1.9393
# stroke_TIA_BIN_FINAL1                     2.0494     0.4880     1.764    2.3813
# hypertensive1                             1.1931     0.8382     1.101    1.2924
# cholesterol1                              1.1087     0.9019     1.016    1.2098
# livesalone1                               1.1825     0.8457     1.086    1.2870
# Sex1                                      1.2403     0.8063     1.152    1.3356
# 
# Concordance= 0.777  (se = 0.004 )
# Likelihood ratio test= 3241  on 11 df,   p=<2e-16
# Wald test            = 2795  on 11 df,   p=<2e-16
# Score (logrank) test = 3230  on 11 df,   p=<2e-16
#save the coefficients
ukbdrs_coefs <- cbind(summary(ukbdrs.cox)$coef,summary(ukbdrs.cox)$conf.int[,c(3,4)])
df_ukbdrs_coefficients <- as.data.frame(ukbdrs_coefs) #5 columns
write.csv(df_ukbdrs_coefficients, paste0(save_pathway,"cox_ukbdrs_coefficients.csv"))

#baseline survival
#The baseline survival is the distribution of the predicted survival for the patient whose predictor values are either the average or 0 (or the reference group for categorical predictors) across the complete follow-up time under study.
#https://www.acpjournals.org/doi/full/10.7326/M22-0844
df_base <- data.frame(Age_when_attended_assesment_centre_0_0 = mean(train.data$Age_when_attended_assesment_centre_0_0),
                      family_history_of_dementia = as.factor(0), Diabetes_BIN_FINAL_0_0 = as.factor(0),
                      current_history_depression = as.factor(0), stroke_TIA_BIN_FINAL = as.factor(0),
                      hypertensive = as.factor(0), cholesterol = as.factor(0), Townsend_deprivation_modelvar = as.factor(0),
                      education_years = mean(train.data$education_years), Sex = as.factor(0), livesalone = as.factor(0))

ukbdrs.surv.baseline <- survfit(ukbdrs.cox, newdata = df_base)
length(ukbdrs.surv.baseline$time) #there are 4880 timepoints
#survival over entire time window:
ukbdrs.surv.baseline$surv[4880]
#0.9910855




#### UKBDRS_APOE_LASSO ####
UKBDRS_APOE_LASSO  <-  paste("Surv(time_at_risk, dementia_BIN_surv) ~  Age_when_attended_assesment_centre_0_0 +  family_history_of_dementia +
                            education_years + Diabetes_BIN_FINAL_0_0  + Townsend_deprivation_modelvar +
                            current_history_depression + stroke_TIA_BIN_FINAL +  
                            hypertensive + cholesterol + livesalone + Sex + APOE_genotype_bin")

ukbdrs.apoe.cox <- coxph(as.formula(UKBDRS_APOE_LASSO), data = train.data)
summary(ukbdrs.apoe.cox)
# coef exp(coef)  se(coef)      z Pr(>|z|)    
# Age_when_attended_assesment_centre_0_0  0.191047  1.210516  0.005282 36.170  < 2e-16 ***
#   family_history_of_dementia1             0.304597  1.356078  0.049486  6.155 7.50e-10 ***
#   education_years                        -0.039661  0.961115  0.007109 -5.579 2.42e-08 ***
#   Diabetes_BIN_FINAL_0_01                 0.589732  1.803504  0.068541  8.604  < 2e-16 ***
#   Townsend_deprivation_modelvar1          0.279947  1.323059  0.050158  5.581 2.39e-08 ***
#   current_history_depression1             0.579893  1.785847  0.053606 10.818  < 2e-16 ***
#   stroke_TIA_BIN_FINAL1                   0.704237  2.022304  0.090659  7.768 7.97e-15 ***
#   hypertensive1                           0.205654  1.228328  0.047098  4.367 1.26e-05 ***
#   cholesterol1                            0.022284  1.022534  0.051211  0.435  0.66346    
# livesalone1                             0.163070  1.177120  0.050204  3.248  0.00116 ** 
#   Sex1                                    0.209128  1.232603  0.043614  4.795 1.63e-06 ***
#   APOE_genotype_bin1                      1.132189  3.102440  0.042597 26.579  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# exp(coef) exp(-coef) lower .95 upper .95
# Age_when_attended_assesment_centre_0_0    1.2105     0.8261    1.1980    1.2231
# family_history_of_dementia1               1.3561     0.7374    1.2307    1.4942
# education_years                           0.9611     1.0405    0.9478    0.9746
# Diabetes_BIN_FINAL_0_01                   1.8035     0.5545    1.5768    2.0628
# Townsend_deprivation_modelvar1            1.3231     0.7558    1.1992    1.4597
# current_history_depression1               1.7858     0.5600    1.6077    1.9837
# stroke_TIA_BIN_FINAL1                     2.0223     0.4945    1.6931    2.4155
# hypertensive1                             1.2283     0.8141    1.1200    1.3471
# cholesterol1                              1.0225     0.9780    0.9249    1.1305
# livesalone1                               1.1771     0.8495    1.0668    1.2988
# Sex1                                      1.2326     0.8113    1.1316    1.3426
# APOE_genotype_bin1                        3.1024     0.3223    2.8539    3.3726
# 
# Concordance= 0.808  (se = 0.005 )
# Likelihood ratio test= 3173  on 12 df,   p=<2e-16
# Wald test            = 2787  on 12 df,   p=<2e-16
# Score (logrank) test = 3195  on 12 df,   p=<2e-16

#save the coefficients
ukbdrs_apoecoefs <- cbind(summary(ukbdrs.apoe.cox)$coef,summary(ukbdrs.apoe.cox)$conf.int[,c(3,4)])
df_ukbdrs_apoecoefficients <- as.data.frame(ukbdrs_apoecoefs) 
write.csv(df_ukbdrs_apoecoefficients, paste0(save_pathway,"cox_ukbdrs_apoe_coefficients.csv"))

#formatted coefs
x <- rbind(df_ukbdrs_coefficients, df_ukbdrs_apoecoefficients)
x$FDR_BH <- p.adjust(x[,5], method="BH") 

df_table1_ukbdrs <- x %>% mutate_if(is.numeric, ~round(., 3))
write.csv(df_table1_ukbdrs, paste0(save_pathway,"cox_coefficients_formatted.csv"))

#baseline survival
df_base <- data.frame(Age_when_attended_assesment_centre_0_0 = mean(train.data$Age_when_attended_assesment_centre_0_0),
                      family_history_of_dementia = as.factor(0), Diabetes_BIN_FINAL_0_0 = as.factor(0),
                      current_history_depression = as.factor(0), stroke_TIA_BIN_FINAL = as.factor(0),
                      hypertensive = as.factor(0), cholesterol = as.factor(0), Townsend_deprivation_modelvar = as.factor(0),
                      education_years = mean(train.data$education_years), Sex = as.factor(0), livesalone = as.factor(0),
                      APOE_genotype_bin = as.factor(0))

ukbdrs.apoe.surv.baseline <- survfit(ukbdrs.apoe.cox, newdata = df_base)
length(ukbdrs.apoe.surv.baseline$time) #there are 4210 timepoints
#survival over entire time window:
ukbdrs.apoe.surv.baseline$surv[4556]
#0.9942871


#### AGE ONLY ####
#compute cox and cr coefficients for the baseline age only model
age_only <-      paste("Surv(time_at_risk, dementia_BIN_surv) ~  Age_when_attended_assesment_centre_0_0")

ageonly.cox <- coxph(as.formula(age_only), data = train.data)
summary(ageonly.cox)
# coef exp(coef) se(coef)     z Pr(>|z|)    
# Age_when_attended_assesment_centre_0_0 0.194889  1.215176 0.004404 44.25   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# exp(coef) exp(-coef) lower .95 upper .95
# Age_when_attended_assesment_centre_0_0     1.215     0.8229     1.205     1.226
# 
# Concordance= 0.75  (se = 0.004 )
# Likelihood ratio test= 2517  on 1 df,   p=<2e-16
# Wald test            = 1958  on 1 df,   p=<2e-16
# Score (logrank) test = 2263  on 1 df,   p=<2e-16

#baseline survival
df_base <- data.frame(Age_when_attended_assesment_centre_0_0 = mean(train.data$Age_when_attended_assesment_centre_0_0))

ageonly.surv.baseline <- survfit(ageonly.cox, newdata = df_base)
length(ageonly.surv.baseline$time) #there are 4880 timepoints
#survival over entire time window:
ageonly.surv.baseline$surv[4880]
#0.9846693

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
test.data$UKBDRS_APOE <- ifelse(test.data$APOE_genotype_bin==1,1,0)
test.data$UKBDRS_sex_score <- ifelse(test.data$Sex==0,0,1)

#build linear predictor for each score using cr adjusted coeffients
#for continuous predictors, subtract the mean of the train data vals

#### UKBDRS LASSO ####
train.data$UKBDRS_LASSO_cox_linear_predictor <- 0.183109361720272*(train.data$Age_when_attended_assesment_centre_0_0 - mean(train.data$Age_when_attended_assesment_centre_0_0)) +
  0.424163934579772*train.data$UKBDRS_familyhistory + -0.0433488369350041*(train.data$education_years - mean(train.data$education_years)) +
  0.596602945300926*train.data$UKBDRS_diabetes + 
  0.259726771758075*train.data$UKBDRS_townsend_group5 +
  0.571791052595081*train.data$UKBDRS_depression + 0.717526717536598*train.data$UKBDRS_stroke +
  0.176537055345449*train.data$UKBDRS_hypertensive + 0.103220277285387*train.data$UKBDRS_cholesterol +
  0.167626941229058*train.data$UKBDRS_livesalone + 0.215353609282428*train.data$UKBDRS_sex_score
summary(train.data$UKBDRS_LASSO_cox_linear_predictor)  
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -2.3712 -0.4286  0.5078  0.4643  1.3388  4.6156 

train.data$UKBDRS_LASSO_cox_predicted_prob <- 1 - 0.9910855^exp(train.data$UKBDRS_LASSO_cox_linear_predictor)  
summary(train.data$UKBDRS_LASSO_cox_predicted_prob)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0008357 0.0058160 0.0147690 0.0252238 0.0335814 0.5954036 


test.data$UKBDRS_LASSO_cox_linear_predictor <- 0.183109361720272*(test.data$Age_when_attended_assesment_centre_0_0 - mean(train.data$Age_when_attended_assesment_centre_0_0)) +
  0.424163934579772*test.data$UKBDRS_familyhistory + -0.0433488369350041*(test.data$education_years - mean(train.data$education_years)) +
  0.596602945300926*test.data$UKBDRS_diabetes + 
  0.259726771758075*test.data$UKBDRS_townsend_group5 +
  0.571791052595081*test.data$UKBDRS_depression + 0.717526717536598*test.data$UKBDRS_stroke +
  0.176537055345449*test.data$UKBDRS_hypertensive + 0.103220277285387*test.data$UKBDRS_cholesterol +
  0.167626941229058*test.data$UKBDRS_livesalone + 0.215353609282428*test.data$UKBDRS_sex_score
summary(test.data$UKBDRS_LASSO_cox_linear_predictor)  
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -2.1934 -0.4299  0.5144  0.4702  1.3375  4.5723 

test.data$UKBDRS_LASSO_cox_predicted_prob <- 1 - 0.9910855^exp(test.data$UKBDRS_LASSO_cox_linear_predictor)  
summary(test.data$UKBDRS_LASSO_cox_predicted_prob)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0009982 0.0058085 0.0148656 0.0254664 0.0335355 0.5795705 


#### UKBDRS APOE ####
train.data$UKBDRS_APOE_LASSO_cox_linear_predictor <- 0.191046740820465*(train.data$Age_when_attended_assesment_centre_0_0 - mean(train.data$Age_when_attended_assesment_centre_0_0)) +
  0.304596870777095*train.data$UKBDRS_familyhistory + -0.0396609431935823*(train.data$education_years - mean(train.data$education_years)) +
  0.589731633009964*train.data$UKBDRS_diabetes + 
  0.279946576466443*train.data$UKBDRS_townsend_group5 +
  0.579892973305337*train.data$UKBDRS_depression + 0.704237268964376*train.data$UKBDRS_stroke +
  0.205653692009378*train.data$UKBDRS_hypertensive + 0.0222843040547898*train.data$UKBDRS_cholesterol +
  0.163070371206403*train.data$UKBDRS_livesalone +
  0.209128446128077*train.data$UKBDRS_sex_score + 1.13218887704916*train.data$UKBDRS_APOE
summary(train.data$UKBDRS_APOE_LASSO_cox_linear_predictor)  
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   -2.43   -0.20    0.75    0.74    1.67    5.27   50767 

train.data$UKBDRS_APOE_LASSO_cox_predicted_prob <- 1 - 0.9942871^exp(train.data$UKBDRS_APOE_LASSO_cox_linear_predictor)  
summary(train.data$UKBDRS_APOE_LASSO_cox_predicted_prob)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#    0.00    0.00    0.01    0.03    0.03    0.68   50767 


test.data$UKBDRS_APOE_LASSO_cox_linear_predictor <- 0.191046740820465*(test.data$Age_when_attended_assesment_centre_0_0 - mean(train.data$Age_when_attended_assesment_centre_0_0)) +
  0.304596870777095*test.data$UKBDRS_familyhistory + -0.0396609431935823*(test.data$education_years - mean(train.data$education_years)) +
  0.589731633009964*test.data$UKBDRS_diabetes + 
  0.279946576466443*test.data$UKBDRS_townsend_group5 +
  0.579892973305337*test.data$UKBDRS_depression + 0.704237268964376*test.data$UKBDRS_stroke +
  0.205653692009378*test.data$UKBDRS_hypertensive + 0.0222843040547898*test.data$UKBDRS_cholesterol +
  0.163070371206403*test.data$UKBDRS_livesalone +
  0.209128446128077*test.data$UKBDRS_sex_score + 1.13218887704916*test.data$UKBDRS_APOE
summary(test.data$UKBDRS_APOE_LASSO_cox_linear_predictor)  
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  -2.208  -0.146   0.810   0.780   1.690   5.712   12905 

test.data$UKBDRS_APOE_LASSO_cox_predicted_prob <- 1 - 0.9942871^exp(test.data$UKBDRS_APOE_LASSO_cox_linear_predictor)  
summary(test.data$UKBDRS_APOE_LASSO_cox_predicted_prob)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.001   0.005   0.013   0.026   0.031   0.823   12905 




#### age only ####
train.data$age_only_cox_linear_predictor <- 0.194889*(train.data$Age_when_attended_assesment_centre_0_0 - mean(train.data$Age_when_attended_assesment_centre_0_0))
train.data$age_only_cox_predicted_prob <- 1 - 0.9846693^exp(train.data$age_only_cox_linear_predictor)  
test.data$age_only_cox_linear_predictor <- 0.194889*(test.data$Age_when_attended_assesment_centre_0_0 - mean(train.data$Age_when_attended_assesment_centre_0_0))
test.data$age_only_cox_predicted_prob <- 1 - 0.9846693^exp(test.data$age_only_cox_linear_predictor)  
summary(train.data$age_only_cox_predicted_prob)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.002211 0.007103 0.015423 0.024713 0.033324 0.148831 
summary(test.data$age_only_cox_predicted_prob)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.002211 0.007103 0.015423 0.024893 0.033324 0.177838 

#save train and test, containing lps and predicted probs, for future auc tests
save(train.data, file=paste0(data_pathway,"11_train_data_outliers_removed_cox_fitted.rda"))
save(test.data, file=paste0(data_pathway,"11_test_data_outliers_removed_cox_fitted.rda"))




