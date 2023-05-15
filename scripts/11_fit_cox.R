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
# Age_when_attended_assesment_centre_0_0  0.182118  1.199755  0.004649 39.171  < 2e-16 ***
#   family_history_of_dementia1             0.423763  1.527699  0.042980  9.860  < 2e-16 ***
#   education_years                        -0.042787  0.958115  0.006197 -6.905 5.02e-12 ***
#   Diabetes_BIN_FINAL_0_01                 0.597528  1.817620  0.057579 10.377  < 2e-16 ***
#   Townsend_deprivation_Groups_0_01       -0.031091  0.969387  0.059760 -0.520 0.602875    
# Townsend_deprivation_Groups_0_02        0.019710  1.019906  0.058878  0.335 0.737801    
# Townsend_deprivation_Groups_0_03        0.063439  1.065495  0.059236  1.071 0.284186    
# Townsend_deprivation_Groups_0_04        0.277068  1.319256  0.057384  4.828 1.38e-06 ***
#   current_history_depression1             0.570034  1.768327  0.046192 12.340  < 2e-16 ***
#   stroke_TIA_BIN_FINAL1                   0.717111  2.048507  0.076594  9.363  < 2e-16 ***
#   hypertensive1                           0.175820  1.192223  0.040815  4.308 1.65e-05 ***
#   cholesterol1                            0.102253  1.107664  0.044524  2.297 0.021642 *  
#   household_occupancy1                    0.150660  1.162601  0.044529  3.383 0.000716 ***
#   household_occupancy2                   -0.062270  0.939629  0.055974 -1.112 0.265930    
# Sex1                                    0.217989  1.243573  0.037796  5.768 8.04e-09 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# exp(coef) exp(-coef) lower .95 upper .95
# Age_when_attended_assesment_centre_0_0    1.1998     0.8335    1.1889    1.2107
# family_history_of_dementia1               1.5277     0.6546    1.4043    1.6620
# education_years                           0.9581     1.0437    0.9465    0.9698
# Diabetes_BIN_FINAL_0_01                   1.8176     0.5502    1.6236    2.0348
# Townsend_deprivation_Groups_0_01          0.9694     1.0316    0.8622    1.0898
# Townsend_deprivation_Groups_0_02          1.0199     0.9805    0.9087    1.1447
# Townsend_deprivation_Groups_0_03          1.0655     0.9385    0.9487    1.1967
# Townsend_deprivation_Groups_0_04          1.3193     0.7580    1.1789    1.4763
# current_history_depression1               1.7683     0.5655    1.6153    1.9359
# stroke_TIA_BIN_FINAL1                     2.0485     0.4882    1.7630    2.3803
# hypertensive1                             1.1922     0.8388    1.1006    1.2915
# cholesterol1                              1.1077     0.9028    1.0151    1.2087
# household_occupancy1                      1.1626     0.8601    1.0654    1.2686
# household_occupancy2                      0.9396     1.0642    0.8420    1.0486
# Sex1                                      1.2436     0.8041    1.1548    1.3392
# 
# Concordance= 0.777  (se = 0.004 )
# Likelihood ratio test= 3245  on 15 df,   p=<2e-16
# Wald test            = 2792  on 15 df,   p=<2e-16
# Score (logrank) test = 3240  on 15 df,   p=<2e-16

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
ukbdrs.surv.baseline$surv[4880]
#0.9910584




#### UKBDRS_APOE_LASSO ####
UKBDRS_APOE_LASSO  <-  paste("Surv(time_at_risk, dementia_BIN_surv) ~  Age_when_attended_assesment_centre_0_0 +  family_history_of_dementia +
                            education_years + Diabetes_BIN_FINAL_0_0  + Townsend_deprivation_Groups_0_0 +
                            current_history_depression + stroke_TIA_BIN_FINAL +  
                            hypertensive + cholesterol + household_occupancy + Sex + APOE_genotype_bin")

ukbdrs.apoe.cox <- coxph(as.formula(UKBDRS_APOE_LASSO), data = train.data)
summary(ukbdrs.apoe.cox)
# coef exp(coef)  se(coef)      z Pr(>|z|)    
# Age_when_attended_assesment_centre_0_0  0.190923  1.210366  0.005414 35.266  < 2e-16 ***
#   family_history_of_dementia1             0.304407  1.355821  0.049488  6.151 7.69e-10 ***
#   education_years                        -0.039707  0.961071  0.007126 -5.572 2.51e-08 ***
#   Diabetes_BIN_FINAL_0_01                 0.589402  1.802910  0.068621  8.589  < 2e-16 ***
#   Townsend_deprivation_Groups_0_01       -0.075792  0.927009  0.068217 -1.111 0.266552    
# Townsend_deprivation_Groups_0_02       -0.045304  0.955707  0.067796 -0.668 0.503974    
# Townsend_deprivation_Groups_0_03        0.006939  1.006963  0.068052  0.102 0.918783    
# Townsend_deprivation_Groups_0_04        0.252315  1.287001  0.065548  3.849 0.000118 ***
#   current_history_depression1             0.579418  1.785000  0.053622 10.806  < 2e-16 ***
#   stroke_TIA_BIN_FINAL1                   0.704846  2.023534  0.090668  7.774 7.61e-15 ***
#   hypertensive1                           0.204735  1.227200  0.047107  4.346 1.39e-05 ***
#   cholesterol1                            0.022505  1.022760  0.051233  0.439 0.660476    
# household_occupancy1                    0.159235  1.172614  0.051768  3.076 0.002099 ** 
#   household_occupancy2                   -0.010453  0.989601  0.063845 -0.164 0.869947    
# Sex1                                    0.209449  1.232999  0.043659  4.797 1.61e-06 ***
#   APOE_genotype_bin1                      1.132261  3.102665  0.042599 26.579  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# exp(coef) exp(-coef) lower .95 upper .95
# Age_when_attended_assesment_centre_0_0    1.2104     0.8262    1.1976    1.2233
# family_history_of_dementia1               1.3558     0.7376    1.2305    1.4939
# education_years                           0.9611     1.0405    0.9477    0.9746
# Diabetes_BIN_FINAL_0_01                   1.8029     0.5547    1.5760    2.0625
# Townsend_deprivation_Groups_0_01          0.9270     1.0787    0.8110    1.0596
# Townsend_deprivation_Groups_0_02          0.9557     1.0463    0.8368    1.0915
# Townsend_deprivation_Groups_0_03          1.0070     0.9931    0.8812    1.1506
# Townsend_deprivation_Groups_0_04          1.2870     0.7770    1.1318    1.4634
# current_history_depression1               1.7850     0.5602    1.6069    1.9828
# stroke_TIA_BIN_FINAL1                     2.0235     0.4942    1.6941    2.4171
# hypertensive1                             1.2272     0.8149    1.1190    1.3459
# cholesterol1                              1.0228     0.9777    0.9250    1.1308
# household_occupancy1                      1.1726     0.8528    1.0595    1.2978
# household_occupancy2                      0.9896     1.0105    0.8732    1.1215
# Sex1                                      1.2330     0.8110    1.1319    1.3432
# APOE_genotype_bin1                        3.1027     0.3223    2.8541    3.3728
# 
# Concordance= 0.808  (se = 0.005 )
# Likelihood ratio test= 3175  on 16 df,   p=<2e-16
# Wald test            = 2789  on 16 df,   p=<2e-16
# Score (logrank) test = 3207  on 16 df,   p=<2e-16
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
ukbdrs.apoe.surv.baseline$surv[4556]
#0.9941051


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
train.data$UKBDRS_LASSO_cox_linear_predictor <- 0.182117705330166*(train.data$Age_when_attended_assesment_centre_0_0 - mean(train.data$Age_when_attended_assesment_centre_0_0)) +
  0.423762823699024*train.data$UKBDRS_familyhistory + -0.0427871312029261*(train.data$education_years - mean(train.data$education_years)) +
  0.597528076668602*train.data$UKBDRS_diabetes + 
  -0.0310910753295164*train.data$UKBDRS_townsend_group2 + 0.0197102410695051*train.data$UKBDRS_townsend_group3 +
  0.0634390217533084*train.data$UKBDRS_townsend_group4 + 0.277067904026735*train.data$UKBDRS_townsend_group5 +
  0.570033639803766*train.data$UKBDRS_depression + 0.71711132183806*train.data$UKBDRS_stroke +
  0.175819933388762*train.data$UKBDRS_hypertensive + 0.102253387457651*train.data$UKBDRS_cholesterol +
  0.150660021680297*train.data$UKBDRS_livesalone + -0.0622701288907812*train.data$UKBDRS_liveswithmultiple +
  0.217988544195546*train.data$UKBDRS_sex_score
summary(train.data$UKBDRS_LASSO_cox_linear_predictor)  
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -2.4309 -0.4307  0.5048  0.4591  1.3344  4.6048 

train.data$UKBDRS_LASSO_cox_predicted_prob <- 1 - 0.9910584^exp(train.data$UKBDRS_LASSO_cox_linear_predictor)  
summary(train.data$UKBDRS_LASSO_cox_predicted_prob)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0007897 0.0058217 0.0147699 0.0252265 0.0335347 0.5925421 


test.data$UKBDRS_LASSO_cox_linear_predictor <- 0.182117705330166*(test.data$Age_when_attended_assesment_centre_0_0 - mean(train.data$Age_when_attended_assesment_centre_0_0)) +
  0.423762823699024*test.data$UKBDRS_familyhistory + -0.0427871312029261*(test.data$education_years - mean(train.data$education_years)) +
  0.597528076668602*test.data$UKBDRS_diabetes + 
  -0.0310910753295164*test.data$UKBDRS_townsend_group2 + 0.0197102410695051*test.data$UKBDRS_townsend_group3 +
  0.0634390217533084*test.data$UKBDRS_townsend_group4 + 0.277067904026735*test.data$UKBDRS_townsend_group5 +
  0.570033639803766*test.data$UKBDRS_depression + 0.71711132183806*test.data$UKBDRS_stroke +
  0.175819933388762*test.data$UKBDRS_hypertensive + 0.102253387457651*test.data$UKBDRS_cholesterol +
  0.150660021680297*test.data$UKBDRS_livesalone + -0.0622701288907812*test.data$UKBDRS_liveswithmultiple +
  0.217988544195546*test.data$UKBDRS_sex_score
summary(test.data$UKBDRS_LASSO_cox_linear_predictor)  
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -2.1975 -0.4307  0.5111  0.4649  1.3302  4.5620 

test.data$UKBDRS_LASSO_cox_predicted_prob <- 1 - 0.9910584^exp(test.data$UKBDRS_LASSO_cox_linear_predictor)  
summary(test.data$UKBDRS_LASSO_cox_predicted_prob)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0009972 0.0058215 0.0148625 0.0254736 0.0333959 0.5769280 


#### UKBDRS APOE ####
train.data$UKBDRS_APOE_LASSO_cox_linear_predictor <- 0.190922647325301*(train.data$Age_when_attended_assesment_centre_0_0 - mean(train.data$Age_when_attended_assesment_centre_0_0)) +
  0.304407328114425*train.data$UKBDRS_familyhistory + -0.0397071055598187*(train.data$education_years - mean(train.data$education_years)) +
  0.58940191568914*train.data$UKBDRS_diabetes + 
  -0.0757918760267292*train.data$UKBDRS_townsend_group2 + -0.0453043718914188*train.data$UKBDRS_townsend_group3 +
  0.00693901461196727*train.data$UKBDRS_townsend_group4 + 0.252314521192837*train.data$UKBDRS_townsend_group5 +
  0.579418204294867*train.data$UKBDRS_depression + 0.704845544426178*train.data$UKBDRS_stroke +
  0.204735145667873*train.data$UKBDRS_hypertensive + 0.0225045388404213*train.data$UKBDRS_cholesterol +
  0.159235467633398*train.data$UKBDRS_livesalone + -0.0104530896743045*train.data$UKBDRS_liveswithmultiple +
  0.209449140689626*train.data$UKBDRS_sex_score + 1.13226135309444*train.data$UKBDRS_APOE
summary(train.data$UKBDRS_APOE_LASSO_cox_linear_predictor)  
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   -2.43   -0.20    0.75    0.74    1.67    5.27   50767 

train.data$UKBDRS_APOE_LASSO_cox_predicted_prob <- 1 - 0.9941051^exp(train.data$UKBDRS_APOE_LASSO_cox_linear_predictor)  
summary(train.data$UKBDRS_APOE_LASSO_cox_predicted_prob)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#    0.00    0.00    0.01    0.03    0.03    0.68   50767 


test.data$UKBDRS_APOE_LASSO_cox_linear_predictor <- 0.190922647325301*(test.data$Age_when_attended_assesment_centre_0_0 - mean(train.data$Age_when_attended_assesment_centre_0_0)) +
  0.304407328114425*test.data$UKBDRS_familyhistory + -0.0397071055598187*(test.data$education_years - mean(train.data$education_years)) +
  0.58940191568914*test.data$UKBDRS_diabetes + 
  -0.0757918760267292*test.data$UKBDRS_townsend_group2 + -0.0453043718914188*test.data$UKBDRS_townsend_group3 +
  0.00693901461196727*test.data$UKBDRS_townsend_group4 + 0.252314521192837*test.data$UKBDRS_townsend_group5 +
  0.579418204294867*test.data$UKBDRS_depression + 0.704845544426178*test.data$UKBDRS_stroke +
  0.204735145667873*test.data$UKBDRS_hypertensive + 0.0225045388404213*test.data$UKBDRS_cholesterol +
  0.159235467633398*test.data$UKBDRS_livesalone + -0.0104530896743045*test.data$UKBDRS_liveswithmultiple +
  0.209449140689626*test.data$UKBDRS_sex_score + 1.13226135309444*test.data$UKBDRS_APOE
summary(test.data$UKBDRS_APOE_LASSO_cox_linear_predictor)  
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  -2.283  -0.189   0.770   0.748   1.662   5.679   12905 

test.data$UKBDRS_APOE_LASSO_cox_predicted_prob <- 1 - 0.9941051^exp(test.data$UKBDRS_APOE_LASSO_cox_linear_predictor)  
summary(test.data$UKBDRS_APOE_LASSO_cox_predicted_prob)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   0.001   0.005   0.013   0.026   0.031   0.823   12905




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




