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
data_pathway = "../../raw_data/"
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
# coxph(formula = as.formula(UKBDRS_LASSO), data = train.data)
# 
# n= 176611, number of events= 3051 
# 
# coef exp(coef)  se(coef)      z Pr(>|z|)    
# Age_when_attended_assesment_centre_0_0  0.180696  1.198050  0.004643 38.914  < 2e-16 ***
#   family_history_of_dementia1             0.433780  1.543079  0.042985 10.091  < 2e-16 ***
#   education_years                        -0.043925  0.957026  0.006211 -7.072 1.53e-12 ***
#   Diabetes_BIN_FINAL_0_01                 0.568808  1.766161  0.057648  9.867  < 2e-16 ***
#   Townsend_deprivation_Groups_0_01       -0.036469  0.964188  0.059765 -0.610  0.54173    
# Townsend_deprivation_Groups_0_02        0.024363  1.024663  0.058884  0.414  0.67906    
# Townsend_deprivation_Groups_0_03        0.054756  1.056283  0.059248  0.924  0.35538    
# Townsend_deprivation_Groups_0_04        0.271078  1.311378  0.057395  4.723 2.32e-06 ***
#   current_history_depression1             0.560606  1.751735  0.046207 12.132  < 2e-16 ***
#   stroke_TIA_BIN_FINAL1                   0.695670  2.005053  0.076570  9.085  < 2e-16 ***
#   hypertensive1                           0.170602  1.186019  0.040738  4.188 2.82e-05 ***
#   cholesterol1                            0.105227  1.110962  0.044406  2.370  0.01780 *  
#   household_occupancy1                    0.144535  1.155502  0.044497  3.248  0.00116 ** 
#   household_occupancy2                   -0.053093  0.948292  0.056009 -0.948  0.34316    
# Sex1                                    0.204492  1.226901  0.037846  5.403 6.54e-08 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# exp(coef) exp(-coef) lower .95 upper .95
# Age_when_attended_assesment_centre_0_0    1.1981     0.8347    1.1872    1.2090
# family_history_of_dementia1               1.5431     0.6481    1.4184    1.6787
# education_years                           0.9570     1.0449    0.9454    0.9687
# Diabetes_BIN_FINAL_0_01                   1.7662     0.5662    1.5775    1.9774
# Townsend_deprivation_Groups_0_01          0.9642     1.0371    0.8576    1.0840
# Townsend_deprivation_Groups_0_02          1.0247     0.9759    0.9130    1.1500
# Townsend_deprivation_Groups_0_03          1.0563     0.9467    0.9405    1.1863
# Townsend_deprivation_Groups_0_04          1.3114     0.7626    1.1719    1.4675
# current_history_depression1               1.7517     0.5709    1.6001    1.9178
# stroke_TIA_BIN_FINAL1                     2.0051     0.4987    1.7256    2.3297
# hypertensive1                             1.1860     0.8432    1.0950    1.2846
# cholesterol1                              1.1110     0.9001    1.0184    1.2120
# household_occupancy1                      1.1555     0.8654    1.0590    1.2608
# household_occupancy2                      0.9483     1.0545    0.8497    1.0583
# Sex1                                      1.2269     0.8151    1.1392    1.3214
# 
# Concordance= 0.776  (se = 0.004 )
# Likelihood ratio test= 3189  on 15 df,   p=<2e-16
# Wald test            = 2746  on 15 df,   p=<2e-16
# Score (logrank) test = 3176  on 15 df,   p=<2e-16

#save the coefficients
summary(ukbdrs.cox)

ukbdrs_coefs <- cbind(summary(ukbdrs.cox)$coef,summary(ukbdrs.cox)$conf.int[,c(3,4)])
df_ukbdrs_coefficients <- as.data.frame(ukbdrs_coefs) #5 columns
write.csv(df_ukbdrs_coefficients, paste0(save_pathway,"cox_ukbdrs_coefficients.csv"))

#baseline survival
#The baseline survival is the distribution of the predicted survival for the patient whose predictor values are either the average or 0 (or the reference group for categorical predictors) across the complete follow-up time under study.
#https://www.acpjournals.org/doi/full/10.7326/M22-0844
df_base <- data.frame(Age_when_attended_assesment_centre_0_0 = mean(train.data$Age_when_attended_assesment_centre_0_0),
                      family_history_of_dementia = as.factor(0), Diabetes_BIN_FINAL_0_0 = as.factor(0),
                      current_history_depression = as.factor(0), stroke_TIA_BIN_FINAL = as.factor(0),
                      hypertensive = as.factor(0), cholesterol = as.factor(0), Townsend_deprivation_Groups_0_0 = as.factor(0),
                      education_years = mean(train.data$education_years), Sex = as.factor(0), household_occupancy = as.factor(0))

ukbdrs.surv.baseline <- survfit(ukbdrs.cox, newdata = df_base)
length(ukbdrs.surv.baseline$time) #there are 4614 timepoints
#survival over entire time window:
ukbdrs.surv.baseline$surv[4614]
#0.9911458




#### UKBDRS_APOE_LASSO ####
UKBDRS_APOE_LASSO  <-  paste("Surv(time_at_risk, dementia_BIN_surv) ~  Age_when_attended_assesment_centre_0_0 +  family_history_of_dementia +
                            education_years + Diabetes_BIN_FINAL_0_0  + Townsend_deprivation_Groups_0_0 +
                            current_history_depression + stroke_TIA_BIN_FINAL +  
                            hypertensive + cholesterol + household_occupancy + Sex + APOE_genotype_bin")

ukbdrs.apoe.cox <- coxph(as.formula(UKBDRS_APOE_LASSO), data = train.data)
summary(ukbdrs.apoe.cox)
# coef exp(coef)  se(coef)      z Pr(>|z|)    
# Age_when_attended_assesment_centre_0_0  0.188811  1.207813  0.005397 34.987  < 2e-16 ***
#   family_history_of_dementia1             0.320784  1.378208  0.049463  6.485 8.85e-11 ***
#   education_years                        -0.040396  0.960409  0.007146 -5.653 1.58e-08 ***
#   Diabetes_BIN_FINAL_0_01                 0.578210  1.782845  0.068723  8.414  < 2e-16 ***
#   Townsend_deprivation_Groups_0_01       -0.076188  0.926642  0.068207 -1.117  0.26399    
# Townsend_deprivation_Groups_0_02       -0.034438  0.966149  0.067806 -0.508  0.61153    
# Townsend_deprivation_Groups_0_03        0.004195  1.004204  0.068055  0.062  0.95084    
# Townsend_deprivation_Groups_0_04        0.266213  1.305013  0.065484  4.065 4.80e-05 ***
#   current_history_depression1             0.568130  1.764963  0.053668 10.586  < 2e-16 ***
#   stroke_TIA_BIN_FINAL1                   0.670637  1.955482  0.090726  7.392 1.45e-13 ***
#   hypertensive1                           0.204445  1.226844  0.047102  4.340 1.42e-05 ***
#   cholesterol1                            0.017352  1.017504  0.051247  0.339  0.73491    
# household_occupancy1                    0.149911  1.161731  0.051722  2.898  0.00375 ** 
#   household_occupancy2                   -0.011820  0.988250  0.063830 -0.185  0.85309    
# Sex1                                    0.197152  1.217929  0.043690  4.512 6.41e-06 ***
#   APOE_genotype_bin1                      1.132391  3.103066  0.042591 26.587  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# exp(coef) exp(-coef) lower .95 upper .95
# Age_when_attended_assesment_centre_0_0    1.2078     0.8279    1.1951     1.221
# family_history_of_dementia1               1.3782     0.7256    1.2509     1.519
# education_years                           0.9604     1.0412    0.9471     0.974
# Diabetes_BIN_FINAL_0_01                   1.7828     0.5609    1.5582     2.040
# Townsend_deprivation_Groups_0_01          0.9266     1.0792    0.8107     1.059
# Townsend_deprivation_Groups_0_02          0.9661     1.0350    0.8459     1.103
# Townsend_deprivation_Groups_0_03          1.0042     0.9958    0.8788     1.147
# Townsend_deprivation_Groups_0_04          1.3050     0.7663    1.1478     1.484
# current_history_depression1               1.7650     0.5666    1.5887     1.961
# stroke_TIA_BIN_FINAL1                     1.9555     0.5114    1.6369     2.336
# hypertensive1                             1.2268     0.8151    1.1187     1.345
# cholesterol1                              1.0175     0.9828    0.9203     1.125
# household_occupancy1                      1.1617     0.8608    1.0497     1.286
# household_occupancy2                      0.9882     1.0119    0.8720     1.120
# Sex1                                      1.2179     0.8211    1.1180     1.327
# APOE_genotype_bin1                        3.1031     0.3223    2.8545     3.373
# 
# Concordance= 0.807  (se = 0.005 )
# Likelihood ratio test= 3134  on 16 df,   p=<2e-16
# Wald test            = 2761  on 16 df,   p=<2e-16
# Score (logrank) test = 3168  on 16 df,   p=<2e-16

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
                      hypertensive = as.factor(0), cholesterol = as.factor(0), Townsend_deprivation_Groups_0_0 = as.factor(0),
                      education_years = mean(train.data$education_years), Sex = as.factor(0), household_occupancy = as.factor(0),
                      APOE_genotype_bin = as.factor(0))

ukbdrs.apoe.surv.baseline <- survfit(ukbdrs.apoe.cox, newdata = df_base)
length(ukbdrs.apoe.surv.baseline$time) #there are 4210 timepoints
#survival over entire time window:
ukbdrs.apoe.surv.baseline$surv[4210]
#0.9941616


#### AGE ONLY ####
#compute cox and cr coefficients for the baseline age only model
age_only <-      paste("Surv(time_at_risk, dementia_BIN_surv) ~  Age_when_attended_assesment_centre_0_0")

ageonly.cox <- coxph(as.formula(age_only), data = train.data)
summary(ageonly.cox)
# coef exp(coef) se(coef)     z Pr(>|z|)    
# Age_when_attended_assesment_centre_0_0 0.193307  1.213255 0.004396 43.98   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# exp(coef) exp(-coef) lower .95 upper .95
# Age_when_attended_assesment_centre_0_0     1.213     0.8242     1.203     1.224
# 
# Concordance= 0.749  (se = 0.004 )
# Likelihood ratio test= 2482  on 1 df,   p=<2e-16
# Wald test            = 1934  on 1 df,   p=<2e-16
# Score (logrank) test = 2232  on 1 df,   p=<2e-16

#baseline survival
df_base <- data.frame(Age_when_attended_assesment_centre_0_0 = mean(train.data$Age_when_attended_assesment_centre_0_0))

ageonly.surv.baseline <- survfit(ageonly.cox, newdata = df_base)
length(ageonly.surv.baseline$time) #there are 4614 timepoints
#survival over entire time window:
ageonly.surv.baseline$surv[4614]
#0.9849508

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
train.data$UKBDRS_LASSO_cox_linear_predictor <- 0.180695590861614*(train.data$Age_when_attended_assesment_centre_0_0 - mean(train.data$Age_when_attended_assesment_centre_0_0)) +
  0.433779744348196*train.data$UKBDRS_familyhistory + -0.0439250682154443*(train.data$education_years - mean(train.data$education_years)) +
  0.568808164096848*train.data$UKBDRS_diabetes + 
  -0.0364688189274691*train.data$UKBDRS_townsend_group2 + 0.0243634555633536*train.data$UKBDRS_townsend_group3 +
  0.0547564068743358*train.data$UKBDRS_townsend_group4 + 0.271078208031246*train.data$UKBDRS_townsend_group5 +
  0.560606497261319*train.data$UKBDRS_depression + 0.695670361087043*train.data$UKBDRS_stroke +
  0.170602226203275*train.data$UKBDRS_hypertensive + 0.105226745067214*train.data$UKBDRS_cholesterol +
  0.144535045635076*train.data$UKBDRS_livesalone + -0.053093199368648*train.data$UKBDRS_liveswithmultiple +
  0.204491820563287*train.data$UKBDRS_sex_score
summary(train.data$UKBDRS_LASSO_cox_linear_predictor)  
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -2.3967 -0.4345  0.4972  0.4484  1.3137  4.5084 

train.data$UKBDRS_LASSO_cox_predicted_prob <- 1 - 0.9911458^exp(train.data$UKBDRS_LASSO_cox_linear_predictor)  
summary(train.data$UKBDRS_LASSO_cox_predicted_prob)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0008091 0.0057430 0.0145154 0.0244633 0.0325438 0.5539517 


test.data$UKBDRS_LASSO_cox_linear_predictor <- 0.180695590861614*(test.data$Age_when_attended_assesment_centre_0_0 - mean(train.data$Age_when_attended_assesment_centre_0_0)) +
  0.433779744348196*test.data$UKBDRS_familyhistory + -0.0439250682154443*(test.data$education_years - mean(train.data$education_years)) +
  0.568808164096848*test.data$UKBDRS_diabetes + 
  -0.0364688189274691*test.data$UKBDRS_townsend_group2 + 0.0243634555633536*test.data$UKBDRS_townsend_group3 +
  0.0547564068743358*test.data$UKBDRS_townsend_group4 + 0.271078208031246*test.data$UKBDRS_townsend_group5 +
  0.560606497261319*test.data$UKBDRS_depression + 0.695670361087043*test.data$UKBDRS_stroke +
  0.170602226203275*test.data$UKBDRS_hypertensive + 0.105226745067214*test.data$UKBDRS_cholesterol +
  0.144535045635076*test.data$UKBDRS_livesalone + -0.053093199368648*test.data$UKBDRS_liveswithmultiple +
  0.204491820563287*test.data$UKBDRS_sex_score
summary(test.data$UKBDRS_LASSO_cox_linear_predictor)  
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -2.4332 -0.4298  0.4945  0.4513  1.3147  4.3639 

test.data$UKBDRS_LASSO_cox_predicted_prob <- 1 - 0.9911458^exp(test.data$UKBDRS_LASSO_cox_linear_predictor)  
summary(test.data$UKBDRS_LASSO_cox_predicted_prob)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0007802 0.0057697 0.0144764 0.0245681 0.0325748 0.5027595




#### UKBDRS APOE ####
train.data$UKBDRS_APOE_LASSO_cox_linear_predictor <- 0.188811098734031*(train.data$Age_when_attended_assesment_centre_0_0 - mean(train.data$Age_when_attended_assesment_centre_0_0)) +
  0.320784372273051*train.data$UKBDRS_familyhistory + -0.0403960218653043*(train.data$education_years - mean(train.data$education_years)) +
  0.578210406879922*train.data$UKBDRS_diabetes + 
  -0.0761880900198712*train.data$UKBDRS_townsend_group2 + -0.034437726599586*train.data$UKBDRS_townsend_group3 +
  0.00419541159982712*train.data$UKBDRS_townsend_group4 + 0.266212713779133*train.data$UKBDRS_townsend_group5 +
  0.568129949510332*train.data$UKBDRS_depression + 0.670636737314582*train.data$UKBDRS_stroke +
  0.204445092809709*train.data$UKBDRS_hypertensive + 0.0173523036725743*train.data$UKBDRS_cholesterol +
  0.149911167963347*train.data$UKBDRS_livesalone + -0.011819789317642*train.data$UKBDRS_liveswithmultiple +
  0.197151795600108*train.data$UKBDRS_sex_score + 1.13239067599009*train.data$UKBDRS_APOE
summary(train.data$UKBDRS_APOE_LASSO_cox_linear_predictor)  
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   -2.37   -0.20    0.75    0.73    1.66    5.59   50910 

train.data$UKBDRS_APOE_LASSO_cox_predicted_prob <- 1 - 0.9941616^exp(train.data$UKBDRS_APOE_LASSO_cox_linear_predictor)  
summary(train.data$UKBDRS_APOE_LASSO_cox_predicted_prob)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#    0.00    0.00    0.01    0.02    0.03    0.71   50910


test.data$UKBDRS_APOE_LASSO_cox_linear_predictor <- 0.188811098734031*(test.data$Age_when_attended_assesment_centre_0_0 - mean(train.data$Age_when_attended_assesment_centre_0_0)) +
  0.320784372273051*test.data$UKBDRS_familyhistory + -0.0403960218653043*(test.data$education_years - mean(train.data$education_years)) +
  0.578210406879922*test.data$UKBDRS_diabetes + 
  -0.0761880900198712*test.data$UKBDRS_townsend_group2 + -0.034437726599586*test.data$UKBDRS_townsend_group3 +
  0.00419541159982712*test.data$UKBDRS_townsend_group4 + 0.266212713779133*test.data$UKBDRS_townsend_group5 +
  0.568129949510332*test.data$UKBDRS_depression + 0.670636737314582*test.data$UKBDRS_stroke +
  0.204445092809709*test.data$UKBDRS_hypertensive + 0.0173523036725743*test.data$UKBDRS_cholesterol +
  0.149911167963347*test.data$UKBDRS_livesalone + -0.011819789317642*test.data$UKBDRS_liveswithmultiple +
  0.197151795600108*test.data$UKBDRS_sex_score + 1.13239067599009*test.data$UKBDRS_APOE
summary(test.data$UKBDRS_APOE_LASSO_cox_linear_predictor)  
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  -2.438  -0.193   0.752   0.737   1.659   5.483   12762 

test.data$UKBDRS_APOE_LASSO_cox_predicted_prob <- 1 - 0.9941616^exp(test.data$UKBDRS_APOE_LASSO_cox_linear_predictor)  
summary(test.data$UKBDRS_APOE_LASSO_cox_predicted_prob)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#    0.00    0.00    0.01    0.02    0.03    0.71   50910




#### age only ####
train.data$age_only_cox_linear_predictor <- 0.193307*(train.data$Age_when_attended_assesment_centre_0_0 - mean(train.data$Age_when_attended_assesment_centre_0_0))
train.data$age_only_cox_predicted_prob <- 1 - 0.9874905^exp(train.data$age_only_cox_linear_predictor)  
test.data$age_only_cox_linear_predictor <- 0.193307*(test.data$Age_when_attended_assesment_centre_0_0 - mean(train.data$Age_when_attended_assesment_centre_0_0))
test.data$age_only_cox_predicted_prob <- 1 - 0.9874905^exp(test.data$age_only_cox_linear_predictor)  
summary(train.data$age_only_cox_predicted_prob)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.001830 0.005826 0.012582 0.020094 0.027061 0.144666 
summary(test.data$age_only_cox_predicted_prob)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.001830 0.005826 0.012582 0.020168 0.027061 0.100717 

#save train and test, containing lps and predicted probs, for future auc tests
save(train.data, file=paste0(data_pathway,"11_train_data_outliers_removed_cox_fitted.rda"))
save(test.data, file=paste0(data_pathway,"11_test_data_outliers_removed_cox_fitted.rda"))




