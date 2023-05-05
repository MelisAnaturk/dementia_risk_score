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
# Call:
#   coxph(formula = as.formula(UKBDRS_LASSO), data = train.data)
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
# Call:
#   coxph(formula = as.formula(UKBDRS_APOE_LASSO), data = train.data)
# 
# n= 125701, number of events= 2280 
# (50910 observations deleted due to missingness)
# 
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
modelvars <- c("Age_when_attended_assesment_centre_0_0","family_history_of_dementia","education_years","Diabetes_BIN_FINAL_0_0",
            "Townsend_deprivation_Groups_0_0","current_history_depression","stroke_TIA_BIN_FINAL","hypertensive","cholesterol",
            "household_occupancy","Sex","dementia_BIN_surv")

#covs <- model.matrix(as.formula(UKBDRS_LASSO), train.data[myvars])[,-1]

ukbdrs.cr.train <- crr(ftime = time_at_risk, fstatus = crr_status, failcode=1, cencode = 0,
                       cov1=model.matrix(as.formula(UKBDRS_LASSO), train.data[modelvars])[,-1], variance=FALSE)
summary(ukbdrs.cr.train)
# coef exp(coef) se(coef)  z p-value
# Age_when_attended_assesment_centre_0_0  0.1766     1.193       NA NA      NA
# family_history_of_dementia1             0.4401     1.553       NA NA      NA
# education_years                        -0.0416     0.959       NA NA      NA
# Diabetes_BIN_FINAL_0_01                 0.5212     1.684       NA NA      NA
# Townsend_deprivation_Groups_0_01       -0.0345     0.966       NA NA      NA
# Townsend_deprivation_Groups_0_02        0.0203     1.021       NA NA      NA
# Townsend_deprivation_Groups_0_03        0.0459     1.047       NA NA      NA
# Townsend_deprivation_Groups_0_04        0.2443     1.277       NA NA      NA
# current_history_depression1             0.5528     1.738       NA NA      NA
# stroke_TIA_BIN_FINAL1                   0.6447     1.905       NA NA      NA
# hypertensive1                           0.1593     1.173       NA NA      NA
# cholesterol1                            0.1051     1.111       NA NA      NA
# household_occupancy1                    0.1222     1.130       NA NA      NA
# household_occupancy2                   -0.0545     0.947       NA NA      NA
# Sex1                                    0.1691     1.184       NA NA      NA
# 
# exp(coef) exp(-coef) 2.5% 97.5%
#   Age_when_attended_assesment_centre_0_0     1.193      0.838   NA    NA
# family_history_of_dementia1                1.553      0.644   NA    NA
# education_years                            0.959      1.042   NA    NA
# Diabetes_BIN_FINAL_0_01                    1.684      0.594   NA    NA
# Townsend_deprivation_Groups_0_01           0.966      1.035   NA    NA
# Townsend_deprivation_Groups_0_02           1.021      0.980   NA    NA
# Townsend_deprivation_Groups_0_03           1.047      0.955   NA    NA
# Townsend_deprivation_Groups_0_04           1.277      0.783   NA    NA
# current_history_depression1                1.738      0.575   NA    NA
# stroke_TIA_BIN_FINAL1                      1.905      0.525   NA    NA
# hypertensive1                              1.173      0.853   NA    NA
# cholesterol1                               1.111      0.900   NA    NA
# household_occupancy1                       1.130      0.885   NA    NA
# household_occupancy2                       0.947      1.056   NA    NA
# Sex1                                       1.184      0.844   NA    NA
# 
# Num. cases = 176611
# Pseudo Log-likelihood = -34921 
# Pseudo likelihood ratio test = 3021  on 15 df,

#save the coefficients
ukbdrs_coefs <- summary(ukbdrs.cr.train)$coef
df_ukbdrs_coefficients <- as.data.frame(ukbdrs_coefs)

#formatted coefs
df_table1 <- as.data.frame(cbind(rownames(ukbdrs_coefs),
                                 round(ukbdrs_coefs[,1],3),
                                 round(ukbdrs_coefs[,2],3)))
names(df_table1) <- c("Predictor","coef","exp(coef)")


#### UKBDRS_APOE_LASSO ####
modelvars <- c("Age_when_attended_assesment_centre_0_0","family_history_of_dementia","education_years","Diabetes_BIN_FINAL_0_0",
               "Townsend_deprivation_Groups_0_0","current_history_depression","stroke_TIA_BIN_FINAL","hypertensive","cholesterol",
               "household_occupancy","Sex","APOE_genotype_bin","dementia_BIN_surv")

apoe.train.data <- train.data[which(!is.na(train.data$APOE_genotype_bin)),]

crr_status <- apoe.train.data$crr_status
time_at_risk <- apoe.train.data$time_at_risk

ukbdrs.apoe.cr.train <- crr(ftime = time_at_risk, fstatus = crr_status,
                            failcode=1, cencode = 0,
                            cov1=model.matrix(as.formula(UKBDRS_APOE_LASSO), apoe.train.data[modelvars])[,-1], variance=FALSE)
summary(ukbdrs.apoe.cr.train)
# coef exp(coef) se(coef)  z p-value
# Age_when_attended_assesment_centre_0_0  0.18447     1.203       NA NA      NA
# family_history_of_dementia1             0.32582     1.385       NA NA      NA
# education_years                        -0.03832     0.962       NA NA      NA
# Diabetes_BIN_FINAL_0_01                 0.53316     1.704       NA NA      NA
# Townsend_deprivation_Groups_0_01       -0.07518     0.928       NA NA      NA
# Townsend_deprivation_Groups_0_02       -0.03608     0.965       NA NA      NA
# Townsend_deprivation_Groups_0_03       -0.00177     0.998       NA NA      NA
# Townsend_deprivation_Groups_0_04        0.24007     1.271       NA NA      NA
# current_history_depression1             0.56333     1.757       NA NA      NA
# stroke_TIA_BIN_FINAL1                   0.61739     1.854       NA NA      NA
# hypertensive1                           0.19161     1.211       NA NA      NA
# cholesterol1                            0.02223     1.022       NA NA      NA
# household_occupancy1                    0.12935     1.138       NA NA      NA
# household_occupancy2                   -0.00981     0.990       NA NA      NA
# Sex1                                    0.16291     1.177       NA NA      NA
# APOE_genotype_bin1                      1.12849     3.091       NA NA      NA
# 
# exp(coef) exp(-coef) 2.5% 97.5%
#   Age_when_attended_assesment_centre_0_0     1.203      0.832   NA    NA
# family_history_of_dementia1                1.385      0.722   NA    NA
# education_years                            0.962      1.039   NA    NA
# Diabetes_BIN_FINAL_0_01                    1.704      0.587   NA    NA
# Townsend_deprivation_Groups_0_01           0.928      1.078   NA    NA
# Townsend_deprivation_Groups_0_02           0.965      1.037   NA    NA
# Townsend_deprivation_Groups_0_03           0.998      1.002   NA    NA
# Townsend_deprivation_Groups_0_04           1.271      0.787   NA    NA
# current_history_depression1                1.757      0.569   NA    NA
# stroke_TIA_BIN_FINAL1                      1.854      0.539   NA    NA
# hypertensive1                              1.211      0.826   NA    NA
# cholesterol1                               1.022      0.978   NA    NA
# household_occupancy1                       1.138      0.879   NA    NA
# household_occupancy2                       0.990      1.010   NA    NA
# Sex1                                       1.177      0.850   NA    NA
# APOE_genotype_bin1                         3.091      0.324   NA    NA
# 
# Num. cases = 125701
# Pseudo Log-likelihood = -24954 
# Pseudo likelihood ratio test = 3008  on 16 df,

#save the coefficients
ukbdrs_apoe_coefs <- summary(ukbdrs.apoe.cr.train)$coef
df_ukbdrs_apoe_coefficients <- as.data.frame(ukbdrs_coefs)

#formatted coefs
x<-as.data.frame(cbind(rownames(ukbdrs_apoe_coefs),
                      round(ukbdrs_apoe_coefs[,1],3),
                      round(ukbdrs_apoe_coefs[,2],3)))
names(x) <- c("Predictor","coef","exp(coef)")
df_table1 <- rbind(df_table1,x)

#save
write.csv(df_ukbdrs_coefficients, paste0(save_pathway,"cox_ukbdrs_coefficients.csv"))
write.csv(df_ukbdrs_apoe_coefficients, paste0(save_pathway,"cox_ukbdrs_apoe_coefficients.csv"))

write.csv(df_table1, paste0(save_pathway,"cox_coefficients_formatted.csv"), row.names = FALSE)
