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


#----- 2. COX REGRESSION ------------------------------------------
# cox regression is now run to calculate the beta-weights for each of the 
#components in our risk score (note, we also add in age, sex and education)

#coxph will want 1,2 as outcome
df_filtered$dementia_BIN_surv <- as.numeric(df_filtered$dementia_BIN_TOTAL)
summary(as.factor(df_filtered$dementia_BIN_surv))
#1      2 
#216949   3813 

# conduct train/test split (same partion as above conducted but more columns are retained)
train.data  <- df_filtered[training.samples, ]
test.data <- df_filtered[-training.samples, ]

#membership in highest deprived group (4, on scale of 0-4) was selected, not others
#create binary version of deprivation to use in model
train.data$Townsend_deprivation_modelvar<-ifelse(train.data$Townsend_deprivation_Groups_0_0==4,1,0)
train.data$Townsend_deprivation_modelvar<-as.factor(train.data$Townsend_deprivation_modelvar)
test.data$Townsend_deprivation_modelvar<-ifelse(test.data$Townsend_deprivation_Groups_0_0==4,1,0)
test.data$Townsend_deprivation_modelvar<-as.factor(test.data$Townsend_deprivation_modelvar)

save(train.data, file = paste0(save_pathway, "train_data_outliers_removed_postlasso.rda"))
save(test.data, file = paste0(save_pathway, "test_data_outliers_removed_postlasso.rda"))

load(file = paste0(save_pathway, "train_data_outliers_removed_postlasso.rda"))
load(file = paste0(save_pathway, "test_data_outliers_removed_postlasso.rda"))

#### 2.1 Test beta coefficients ####
#based on the lasso selected vars, compute model coefficients in train data
#use 2 models - one with apoe, one without. 

#### Test models ####

#### UKBDRS Linear ####
#TESTED THIS IN WH AND IT WAS OK
UKBDRS_LASSO  <-  paste("Surv(time_at_risk, dementia_BIN_surv) ~  Age_when_attended_assesment_centre_0_0 +  family_history_of_dementia +
                            education_years + Diabetes_BIN_FINAL_0_0  + Townsend_deprivation_Groups_0_0 +
                            current_history_depression + stroke_TIA_BIN_FINAL +  
                            hypertensive + cholesterol + household_occupancy + Sex")

ukbdrs.cox <- coxph(as.formula(UKBDRS_LASSO), data = train.data)
summary(ukbdrs.cox)
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

cox.zph(ukbdrs.cox)
# Age_when_attended_assesment_centre_0_0  9.0199  1 0.00267
# family_history_of_dementia              0.3350  1 0.56273
# education_years                         0.0316  1 0.85885
# Diabetes_BIN_FINAL_0_0                  0.8196  1 0.36529
# Townsend_deprivation_Groups_0_0         8.0088  4 0.09126
# current_history_depression             17.3982  1   3e-05
# stroke_TIA_BIN_FINAL                    1.5327  1 0.21570
# hypertensive                            1.1205  1 0.28981
# cholesterol                             0.6421  1 0.42297
# household_occupancy                     0.3183  2 0.85288
# Sex                                     0.7551  1 0.38488
# GLOBAL                                 40.4159 15 0.00039

#age, depr and townsend violate
ph<-cox.zph(ukbdrs.cox)


#baseline survival
#The baseline survival is the distribution of the predicted survival for the patient whose predictor values are either the average or 0 (or the reference group for categorical predictors) across the complete follow-up time under study.
#https://www.acpjournals.org/doi/full/10.7326/M22-0844
survv <- survfit(ukbdrs.cox)
summary(survfit(ukbdrs.cox))

df_base <- data.frame(Age_when_attended_assesment_centre_0_0 = mean(train.data$Age_when_attended_assesment_centre_0_0),
                      family_history_of_dementia = as.factor(0), Diabetes_BIN_FINAL_0_0 = as.factor(0),
                      current_history_depression = as.factor(0), stroke_TIA_BIN_FINAL = as.factor(0),
                      hypertensive = as.factor(0), cholesterol = as.factor(0), Townsend_deprivation_Groups_0_0 = as.factor(0),
                      education_years = mean(train.data$education_years), Sex = as.factor(0), household_occupancy = as.factor(0))

survv_baseline <- survfit(ukbdrs.cox, newdata = df_base)
length(survv_baseline$time) #there are 4614 timepoints
#survival over entire time window:
survv_baseline$surv[4614]
#0.9911458

