#assess diagnostics of the fitted coxph model

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

load(file=paste0(data_pathway,"11_train_data_outliers_removed_fitted.rda"))
load(file=paste0(data_pathway,"11_test_data_outliers_removed_fitted.rda"))


#### PH ####
UKBDRS_LASSO  <-  paste("Surv(time_at_risk, dementia_BIN_surv) ~  Age_when_attended_assesment_centre_0_0 +  family_history_of_dementia +
                            education_years + Diabetes_BIN_FINAL_0_0  + Townsend_deprivation_Groups_0_0 +
                            current_history_depression + stroke_TIA_BIN_FINAL +  
                            hypertensive + cholesterol + household_occupancy + Sex")

ukbdrs.cox <- coxph(as.formula(UKBDRS_LASSO), data = train.data)
summary(ukbdrs.cox)

ph <- cox.zph(ukbdrs.cox)
ph
# chisq df       p
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
#age, depression appear to violate

plot(ph) #plots schoenfield residuals vs time for each variable
#non zero slope of time vs residual indicates violation of ph

#save plots
vars<-names(as.data.frame(ph$y))
for (v in seq(1,11)){
  png(file=paste0("../results/diagnostics/ph/",vars[v],"_schoenfeld.png",sep=""),
      height = 400, width = 600)
  plot(ph[v])
  dev.off()
}

#depression slope is minimally non zero..
#same with age

#depression - lets look at surv fcn vs surv time
#parallel curves indicate satisfaction of ph
depression_fit <- survfit(Surv(time_at_risk, dementia_BIN_surv) ~ current_history_depression,
                          data = train.data)

library(ggfortify)
#curves are pretty similar
plot(log(depression_fit$time), log(-log(depression_fit$surv)))
png(file="../results/diagnostics/ph/depression_loglog.png",
    height = 400, width = 600)
plot(log(depression_fit$time), log(-log(depression_fit$surv)))
dev.off()

depression.plot <- autoplot(depression_fit)
depression.plot <- depression.plot + 
  ggtitle("Depression Survival") +
  labs(x = "Time", y = "Survival Probability") +
  guides(fill="none") +
  labs(colour = "Depression") +
  scale_color_manual(labels = c("No", "Yes"), values = c(1,2))
print(depression.plot)
ggsave("../results/diagnostics/ph/depression_survplot.png", 
       plot = depression.plot, 
       device = "png", 
       dpi = 320)

#lets compare with a categorical that was not sig, family history
famhx_fit <- survfit(Surv(time_at_risk, dementia_BIN_surv) ~ family_history_of_dementia,
                          data = train.data)

famhx.plot <- autoplot(famhx_fit)
famhx.plot <- famhx.plot + 
  ggtitle("Family History Survival") +
  labs(x = "Time", y = "Survival Probability") +
  guides(fill="none") +
  labs(colour = "Family History") +
  scale_color_manual(labels = c("No", "Yes"), values = c(1,2))
print(famhx.plot)
ggsave("../results/diagnostics/ph/famhx_survplot.png", 
       plot = famhx.plot, 
       device = "png", 
       dpi = 320)
#curves are pretty similar
plot(log(famhx_fit$time), log(-log(famhx_fit$surv)))
png(file="../results/diagnostics/ph/famhistloglog.png",
    height = 400, width = 600)
plot(log(famhx_fit$time), log(-log(famhx_fit$surv)))
dev.off()


diabetes_fit <- survfit(Surv(time_at_risk, dementia_BIN_surv) ~ Diabetes_BIN_FINAL_0_0,
                     data = train.data)

diabetes.plot <- autoplot(diabetes_fit)
diabetes.plot <- diabetes.plot + 
  ggtitle("Diabetes Survival") +
  labs(x = "Time", y = "Survival Probability") +
  guides(fill="none") +
  labs(colour = "Diabetes") +
  scale_color_manual(labels = c("No", "Yes"), values = c(1,2))
print(diabetes.plot)
ggsave("../results/diagnostics/ph/diabetes_survplot.png", 
       plot = diabetes.plot, 
       device = "png", 
       dpi = 320)
#curves are pretty similar
plot(log(diabetes_fit$time), log(-log(diabetes_fit$surv)))
png(file="../results/diagnostics/ph/diabetes_loglog.png",
    height = 400, width = 600)
plot(log(diabetes_fit$time), log(-log(diabetes_fit$surv)))
dev.off()

#curve of depression is not graphically much diff than those for famhx, diabetes, which satisfied hazards
#move on

#age plot is similarly close to non zero
#move on






#### Non linear ####
#look at the martingale residuals vs predictor of interest to assess linearity
#base/simple cox
UKBDRS_LASSO  <-  paste("Surv(time_at_risk, dementia_BIN_surv) ~  Age_when_attended_assesment_centre_0_0 +  family_history_of_dementia +
                            education_years + Diabetes_BIN_FINAL_0_0  + Townsend_deprivation_Groups_0_0 +
                            current_history_depression + stroke_TIA_BIN_FINAL +  
                            hypertensive + cholesterol + household_occupancy + Sex")

ukbdrs.cox <- coxph(as.formula(UKBDRS_LASSO), data = train.data)
summary(ukbdrs.cox)
martin.resids <- as.data.frame(residuals(ukbdrs.cox, type="martingale"))

#now plot residuals vs age
age.loess <- loess(martin.resids[,] ~ train.data$Age_when_attended_assesment_centre_0_0)
plot(train.data$Age_when_attended_assesment_centre_0_0, martin.resids[,])
j <- order(train.data$Age_when_attended_assesment_centre_0_0)
lines(train.data$Age_when_attended_assesment_centre_0_0[j],age.loess$fitted[j])
#FLAT
#save it
png(file="../results/diagnostics/lin/age.vs.martin.png",
    height = 400, width = 600)
plot(train.data$Age_when_attended_assesment_centre_0_0, martin.resids[,], xlab="Age",ylab="Martingale Residual")
j <- order(train.data$Age_when_attended_assesment_centre_0_0)
lines(train.data$Age_when_attended_assesment_centre_0_0[j],age.loess$fitted[j])
dev.off()

#now plot residuals vs education
education.loess <- loess(martin.resids[,] ~ train.data$education_years)
plot(train.data$education_years, martin.resids[,])
j <- order(train.data$education_years)
lines(train.data$education_years[j],education.loess$fitted[j])

png(file="../results/diagnostics/lin/edu.vs.martin.png",
    height = 400, width = 600)
plot(train.data$education_years, martin.resids[,], xlab="Education (years)",ylab="Martingale Residual")
j <- order(train.data$education_years)
lines(train.data$education_years[j],education.loess$fitted[j])
dev.off()


#nontheless, lets test nonlinear impact on auc
age2 = (train.data$Age_when_attended_assesment_centre_0_0 - mean(train.data$Age_when_attended_assesment_centre_0_0))^2
age3 = (train.data$Age_when_attended_assesment_centre_0_0 - mean(train.data$Age_when_attended_assesment_centre_0_0))^3

train.data$age2o = resid(lm(age2 ~ train.data$Age_when_attended_assesment_centre_0_0))
train.data$age3o = resid(lm(age3 ~ train.data$Age_when_attended_assesment_centre_0_0))

education2 = (train.data$education_years - mean(train.data$education_years))^2
education3 = (train.data$education_years - mean(train.data$education_years))^3

train.data$education2o = resid(lm(education2 ~ train.data$education_years))
train.data$education3o = resid(lm(education3 ~ train.data$education_years))

UKBDRS_LASSO_nonlinear  <-  paste("Surv(time_at_risk, dementia_BIN_surv) ~  Age_when_attended_assesment_centre_0_0 + family_history_of_dementia +
                            age2o + age3o + education2o + education3o +
                            education_years + Diabetes_BIN_FINAL_0_0  + Townsend_deprivation_Groups_0_0 +
                            current_history_depression + stroke_TIA_BIN_FINAL +  
                            hypertensive + cholesterol + household_occupancy + Sex")

ukbdrs.nlin.cox <- coxph(as.formula(UKBDRS_LASSO_nonlinear), data = train.data)
summary(ukbdrs.nlin.cox)
# coef  exp(coef)   se(coef)      z Pr(>|z|)    
# Age_when_attended_assesment_centre_0_0  1.652e-01  1.180e+00  5.963e-03 27.701  < 2e-16 ***
#   family_history_of_dementia1             4.397e-01  1.552e+00  4.302e-02 10.220  < 2e-16 ***
#   age2o                                   3.852e-03  1.004e+00  9.040e-04  4.260 2.04e-05 ***
#   age3o                                  -6.082e-04  9.994e-01  1.496e-04 -4.067 4.76e-05 ***
#   education2o                             3.695e-03  1.004e+00  1.117e-03  3.308 0.000941 ***
#   education3o                             1.626e-04  1.000e+00  8.753e-05  1.858 0.063131 .  
# education_years                        -4.276e-02  9.581e-01  6.142e-03 -6.962 3.35e-12 ***
#   Diabetes_BIN_FINAL_0_01                 5.675e-01  1.764e+00  5.764e-02  9.845  < 2e-16 ***
#   Townsend_deprivation_Groups_0_01       -3.875e-02  9.620e-01  5.977e-02 -0.648 0.516756    
# Townsend_deprivation_Groups_0_02        1.988e-02  1.020e+00  5.890e-02  0.338 0.735714    
# Townsend_deprivation_Groups_0_03        4.964e-02  1.051e+00  5.926e-02  0.838 0.402242    
# Townsend_deprivation_Groups_0_04        2.611e-01  1.298e+00  5.749e-02  4.542 5.56e-06 ***
#   current_history_depression1             5.610e-01  1.752e+00  4.621e-02 12.139  < 2e-16 ***
#   stroke_TIA_BIN_FINAL1                   6.917e-01  1.997e+00  7.657e-02  9.034  < 2e-16 ***
#   hypertensive1                           1.676e-01  1.183e+00  4.076e-02  4.113 3.91e-05 ***
#   cholesterol1                            1.027e-01  1.108e+00  4.440e-02  2.312 0.020782 *  
#   household_occupancy1                    1.427e-01  1.153e+00  4.451e-02  3.205 0.001351 ** 
#   household_occupancy2                   -6.256e-02  9.394e-01  5.646e-02 -1.108 0.267862    
# Sex1                                    2.014e-01  1.223e+00  3.786e-02  5.319 1.04e-07 ***
#nlin terms are sig but not as strong

nlinvars <- c("Age_when_attended_assesment_centre_0_0","age2o","age3o","family_history_of_dementia","education_years",
              "education2o","education3o","Diabetes_BIN_FINAL_0_0",
              "Townsend_deprivation_Groups_0_0","current_history_depression","stroke_TIA_BIN_FINAL","hypertensive","cholesterol",
              "household_occupancy","Sex","dementia_BIN_surv")
crr_status <- train.data$crr_status
time_at_risk <- train.data$time_at_risk
nlincovs <- model.matrix(as.formula(UKBDRS_LASSO_nonlinear), train.data[nlinvars])[,-1]

cr.train.nlin <- crr(ftime = time_at_risk, fstatus = crr_status, failcode=1, cencode = 0, cov1=nlincovs, variance=TRUE)
save(cr.train.nlin, file=paste0(save_pathway,"ukbdrs.nonlin.cr.train.rda"))
summary(cr.train.nlin)
# coef exp(coef) se(coef)      z p-value
# Age_when_attended_assesment_centre_0_0  0.162334     1.176 5.99e-03 27.089 0.0e+00
# family_history_of_dementia1             0.445328     1.561 4.30e-02 10.357 0.0e+00
# age2o                                   0.003643     1.004 8.92e-04  4.085 4.4e-05
# age3o                                  -0.000621     0.999 1.50e-04 -4.144 3.4e-05
# education2o                             0.003635     1.004 1.04e-03  3.494 4.8e-04
# education3o                             0.000144     1.000 8.14e-05  1.767 7.7e-02
# education_years                        -0.040395     0.960 6.20e-03 -6.515 7.3e-11
# Diabetes_BIN_FINAL_0_01                 0.519960     1.682 5.77e-02  9.010 0.0e+00
# Townsend_deprivation_Groups_0_01       -0.036296     0.964 5.98e-02 -0.607 5.4e-01
# Townsend_deprivation_Groups_0_02        0.016565     1.017 5.88e-02  0.282 7.8e-01
# Townsend_deprivation_Groups_0_03        0.041226     1.042 5.92e-02  0.697 4.9e-01
# Townsend_deprivation_Groups_0_04        0.235028     1.265 5.71e-02  4.117 3.8e-05
# current_history_depression1             0.552648     1.738 4.62e-02 11.950 0.0e+00
# stroke_TIA_BIN_FINAL1                   0.641151     1.899 7.78e-02  8.246 2.2e-16
# hypertensive1                           0.156354     1.169 4.10e-02  3.815 1.4e-04
# cholesterol1                            0.102545     1.108 4.51e-02  2.275 2.3e-02
# household_occupancy1                    0.120786     1.128 4.43e-02  2.725 6.4e-03
# household_occupancy2                   -0.062682     0.939 5.65e-02 -1.109 2.7e-01
# Sex1                                    0.165576     1.180 3.80e-02  4.358 1.3e-05
# 
# exp(coef) exp(-coef)  2.5% 97.5%
#   Age_when_attended_assesment_centre_0_0     1.176      0.850 1.163 1.190
# family_history_of_dementia1                1.561      0.641 1.435 1.698
# age2o                                      1.004      0.996 1.002 1.005
# age3o                                      0.999      1.001 0.999 1.000
# education2o                                1.004      0.996 1.002 1.006
# education3o                                1.000      1.000 1.000 1.000
# education_years                            0.960      1.041 0.949 0.972
# Diabetes_BIN_FINAL_0_01                    1.682      0.595 1.502 1.883
# Townsend_deprivation_Groups_0_01           0.964      1.037 0.858 1.084
# Townsend_deprivation_Groups_0_02           1.017      0.984 0.906 1.141
# Townsend_deprivation_Groups_0_03           1.042      0.960 0.928 1.170
# Townsend_deprivation_Groups_0_04           1.265      0.791 1.131 1.415
# current_history_depression1                1.738      0.575 1.587 1.903
# stroke_TIA_BIN_FINAL1                      1.899      0.527 1.630 2.211
# hypertensive1                              1.169      0.855 1.079 1.267
# cholesterol1                               1.108      0.903 1.014 1.210
# household_occupancy1                       1.128      0.886 1.035 1.231
# household_occupancy2                       0.939      1.065 0.841 1.049
# Sex1                                       1.180      0.847 1.095 1.271
# 
# Num. cases = 176611
# Pseudo Log-likelihood = -34904 
# Pseudo likelihood ratio test = 3054  on 19 df,

#mape age2o..etc into test
load(file=paste0(data_pathway,"11_test_data_outliers_removed_fitted.rda"))
train.data$age2 = (train.data$Age_when_attended_assesment_centre_0_0 - mean(train.data$Age_when_attended_assesment_centre_0_0))^2
train.data$age3 = (train.data$Age_when_attended_assesment_centre_0_0 - mean(train.data$Age_when_attended_assesment_centre_0_0))^3

testage2 = (test.data$Age_when_attended_assesment_centre_0_0 - mean(train.data$Age_when_attended_assesment_centre_0_0))^2
testage3 = (test.data$Age_when_attended_assesment_centre_0_0 - mean(train.data$Age_when_attended_assesment_centre_0_0))^3

age2mod <-  lm(age2 ~ Age_when_attended_assesment_centre_0_0, data = train.data)
age3mod <-  lm(age3 ~ Age_when_attended_assesment_centre_0_0, data = train.data)

test.data$age2o = testage2 - predict(age2mod, test.data)
test.data$age3o = testage3 - predict(age3mod, test.data)

train.data$education2 = (train.data$education_years - mean(train.data$education_years))^2
train.data$education3 = (train.data$education_years - mean(train.data$education_years))^3

testeducation2 = (test.data$education_years - mean(train.data$education_years))^2
testeducation3 = (test.data$education_years - mean(train.data$education_years))^3

education2mod <-  lm(education2 ~ education_years, data = train.data)
education3mod <-  lm(education3 ~ education_years, data = train.data)

test.data$education2o = testeducation2 - predict(education2mod, test.data)
test.data$education3o = testeducation3 - predict(education3mod, test.data)


#baseline survival of nlin model
df_base <- data.frame(Age_when_attended_assesment_centre_0_0 = mean(train.data$Age_when_attended_assesment_centre_0_0),
                      family_history_of_dementia = as.factor(0), Diabetes_BIN_FINAL_0_0 = as.factor(0),
                      current_history_depression = as.factor(0), stroke_TIA_BIN_FINAL = as.factor(0),
                      hypertensive = as.factor(0), cholesterol = as.factor(0), Townsend_deprivation_Groups_0_0 = as.factor(0),
                      education_years = mean(train.data$education_years), Sex = as.factor(0), household_occupancy = as.factor(0),
                      age2o = 0, age3o = 0, education2o = 0, education3o = 0)

nlinsurvv_baseline <- survfit(ukbdrs.nlin.cox, newdata = df_base)
length(nlinsurvv_baseline$time) #there are 4614 timepoints
#survival over entire time window:
nlinsurvv_baseline$surv[4614]
# 0.9908608

#build linear predictor, predicted probability for train and test

train.data$NLIN_UKBDRS_Tsend_score <- ifelse(train.data$Townsend_deprivation_Groups_0_0==0,0,
                                             ifelse(train.data$Townsend_deprivation_Groups_0_0==1,-0.0362962405,
                                                    ifelse(train.data$Townsend_deprivation_Groups_0_0==2,0.0165652972,
                                                           ifelse(train.data$Townsend_deprivation_Groups_0_0==3,0.0412255997,0.2350284))))
train.data$NLIN_UKBDRS_household_score <- ifelse(train.data$household_occupancy==0,0,
                                                 ifelse(train.data$household_occupancy==1,0.1207861952,-0.062682))

train.data$NLIN_UKBDRS_LASSO_linear_predictor <- 0.1623336309*(train.data$Age_when_attended_assesment_centre_0_0 - mean(train.data$Age_when_attended_assesment_centre_0_0)) +
  0.4453283967*train.data$UKBDRS_familyhistory - 0.0403948791*(train.data$education_years - mean(train.data$education_years)) +
  0.5199598637*train.data$UKBDRS_diabetes + train.data$NLIN_UKBDRS_Tsend_score +
  0.5526476070*train.data$UKBDRS_depression+ 0.6411513996*train.data$UKBDRS_stroke +
  0.1563540387*train.data$UKBDRS_hypertensive + 0.1025447161*train.data$UKBDRS_cholesterol +
  train.data$NLIN_UKBDRS_household_score + 0.1655764383*train.data$UKBDRS_sex_score +
  0.0036433708*train.data$age2o - 0.0006210342*train.data$age3o +
  0.0036349084*train.data$education2o + 0.0001438530*train.data$education3o 
summary(train.data$NLIN_UKBDRS_LASSO_linear_predictor)  
#note age2o, age3o, education2o, education3o have mean 0 so no need to subtract training data mean

train.data$NLIN_UKBDRS_LASSO_predicted_prob <- 1 - 0.9908608^exp(train.data$NLIN_UKBDRS_LASSO_linear_predictor)  
summary(train.data$NLIN_UKBDRS_LASSO_predicted_prob)
plot(train.data$UKBDRS_LASSO_predicted_prob, train.data$NLIN_UKBDRS_LASSO_predicted_prob)


#test
test.data$NLIN_UKBDRS_Tsend_score <- ifelse(test.data$Townsend_deprivation_Groups_0_0==0,0,
                                            ifelse(test.data$Townsend_deprivation_Groups_0_0==1,-0.0362962405,
                                                   ifelse(test.data$Townsend_deprivation_Groups_0_0==2,0.0165652972,
                                                          ifelse(test.data$Townsend_deprivation_Groups_0_0==3,0.0412255997,0.2350284))))
test.data$NLIN_UKBDRS_household_score <- ifelse(test.data$household_occupancy==0,0,
                                                ifelse(test.data$household_occupancy==1,0.1207861952,-0.062682))
## LEFT OFF HERE
test.data$NLIN_UKBDRS_LASSO_linear_predictor <- 0.1623336309*(test.data$Age_when_attended_assesment_centre_0_0 - mean(train.data$Age_when_attended_assesment_centre_0_0)) +
  0.4453283967*test.data$UKBDRS_familyhistory - 0.0403948791*(test.data$education_years - mean(train.data$education_years)) +
  0.5199598637*test.data$UKBDRS_diabetes + test.data$NLIN_UKBDRS_Tsend_score +
  0.5526476070*test.data$UKBDRS_depression+ 0.6411513996*test.data$UKBDRS_stroke +
  0.1563540387*test.data$UKBDRS_hypertensive + 0.1025447161*test.data$UKBDRS_cholesterol +
  test.data$NLIN_UKBDRS_household_score + 0.1655764383*test.data$UKBDRS_sex_score +
  0.0036433708*test.data$age2o - 0.0006210342*test.data$age3o +
  0.0036349084*test.data$education2o + 0.0001438530*test.data$education3o 
summary(test.data$NLIN_UKBDRS_LASSO_linear_predictor)  

test.data$NLIN_UKBDRS_LASSO_predicted_prob <- 1 - 0.9908608^exp(test.data$NLIN_UKBDRS_LASSO_linear_predictor)  
summary(test.data$NLIN_UKBDRS_LASSO_predicted_prob)


#save train and test, containing lps and predicted probs, for future auc tests
save(train.data, file=paste0(data_pathway,"12_train_data_outliers_removed_fitted.rda"))
save(test.data, file=paste0(data_pathway,"12_test_data_outliers_removed_fitted.rda"))
