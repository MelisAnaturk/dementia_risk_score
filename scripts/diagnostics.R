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
data_pathway = "../../raw_data/modelvar/"
model_pathway = "../models/"
save_pathway = "../results/"

load(file=paste0(data_pathway,"12_train_data_outliers_removed_cox_crr_fitted.rda"))
load(file=paste0(data_pathway,"12_test_data_outliers_removed_cox_crr_fitted.rda"))

#define status
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


df_dementia <- test.data[which(test.data$dementia_BIN_TOTAL==1),]
df_dementia$crr_status <- 1

#deaths
df_deaths <- test.data[which( (test.data$dementia_BIN_TOTAL==0) & (test.data$death_record_refresh==1) ),]
df_deaths$crr_status <- 2

#healthy
df_healthy <- test.data[which( (test.data$dementia_BIN_TOTAL==0) & (test.data$death_record_refresh==0) ),]
df_healthy$crr_status<-0

test.data <- rbind(df_dementia, df_healthy)
test.data<-rbind(test.data, df_deaths)
rm(df_dementia, df_deaths, df_healthy)

train.data$y <- ifelse(train.data$dementia_BIN_TOTAL==1,1,0)
test.data$y <- ifelse(test.data$dementia_BIN_TOTAL==1,1,0)


#### PH ####
UKBDRS_LASSO  <-  paste("Surv(time_at_risk, dementia_BIN_surv) ~  Age_when_attended_assesment_centre_0_0 +  family_history_of_dementia +
                            education_years + Diabetes_BIN_FINAL_0_0  + Townsend_deprivation_modelvar +
                            current_history_depression + stroke_TIA_BIN_FINAL +  
                            hypertensive + cholesterol + livesalone + Sex")

ukbdrs.cox <- coxph(as.formula(UKBDRS_LASSO), data = train.data)
summary(ukbdrs.cox)

ph <- cox.zph(ukbdrs.cox)
ph
# chisq df       p
# Age_when_attended_assesment_centre_0_0 11.275  1 0.00079
# family_history_of_dementia              0.522  1 0.47012
# education_years                         0.110  1 0.74003
# Diabetes_BIN_FINAL_0_0                  1.987  1 0.15864
# Townsend_deprivation_modelvar           3.318  1 0.06852
# current_history_depression             16.852  1 4.0e-05
# stroke_TIA_BIN_FINAL                    0.887  1 0.34628
# hypertensive                            2.131  1 0.14432
# cholesterol                             1.305  1 0.25324
# livesalone                              0.469  1 0.49351
# Sex                                     0.283  1 0.59449
# GLOBAL                                 37.559 11 9.3e-05
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
                            education_years + Diabetes_BIN_FINAL_0_0  + Townsend_deprivation_modelvar +
                            current_history_depression + stroke_TIA_BIN_FINAL +  
                            hypertensive + cholesterol + livesalone + Sex")


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
                            education_years + Diabetes_BIN_FINAL_0_0  + Townsend_deprivation_modelvar +
                            current_history_depression + stroke_TIA_BIN_FINAL +  
                            hypertensive + cholesterol + livesalone + Sex")

ukbdrs.nlin.cox <- coxph(as.formula(UKBDRS_LASSO_nonlinear), data = train.data)
summary(ukbdrs.nlin.cox)
# Age_when_attended_assesment_centre_0_0  1.675e-01  1.182e+00  5.806e-03 28.852  < 2e-16 ***
#   family_history_of_dementia1             4.305e-01  1.538e+00  4.301e-02 10.009  < 2e-16 ***
#   age2o                                   3.855e-03  1.004e+00  9.034e-04  4.267 1.98e-05 ***
#   age3o                                  -5.899e-04  9.994e-01  1.496e-04 -3.943 8.04e-05 ***
#   education2o                             3.527e-03  1.004e+00  1.123e-03  3.140 0.001689 ** 
#   education3o                             1.932e-04  1.000e+00  8.797e-05  2.196 0.028068 *  
#   education_years                        -4.259e-02  9.583e-01  6.122e-03 -6.957 3.48e-12 ***
#   Diabetes_BIN_FINAL_0_01                 5.944e-01  1.812e+00  5.754e-02 10.330  < 2e-16 ***
#   Townsend_deprivation_modelvar1          2.527e-01  1.287e+00  4.344e-02  5.817 5.99e-09 ***
#   current_history_depression1             5.714e-01  1.771e+00  4.619e-02 12.368  < 2e-16 ***
#   stroke_TIA_BIN_FINAL1                   7.141e-01  2.042e+00  7.659e-02  9.324  < 2e-16 ***
#   hypertensive1                           1.737e-01  1.190e+00  4.084e-02  4.253 2.11e-05 ***
#   cholesterol1                            1.011e-01  1.106e+00  4.452e-02  2.271 0.023162 *  
#   livesalone1                             1.670e-01  1.182e+00  4.322e-02  3.865 0.000111 ***
#   Sex1                                    2.122e-01  1.236e+00  3.778e-02  5.618 1.94e-08 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# exp(coef) exp(-coef) lower .95 upper .95
# Age_when_attended_assesment_centre_0_0    1.1824     0.8458    1.1690    1.1959
# family_history_of_dementia1               1.5380     0.6502    1.4137    1.6733
# age2o                                     1.0039     0.9962    1.0021    1.0056
# age3o                                     0.9994     1.0006    0.9991    0.9997
# education2o                               1.0035     0.9965    1.0013    1.0057
# education3o                               1.0002     0.9998    1.0000    1.0004
# education_years                           0.9583     1.0435    0.9469    0.9699
# Diabetes_BIN_FINAL_0_01                   1.8119     0.5519    1.6187    2.0282
# Townsend_deprivation_modelvar1            1.2875     0.7767    1.1824    1.4019
# current_history_depression1               1.7707     0.5648    1.6174    1.9385
# stroke_TIA_BIN_FINAL1                     2.0424     0.4896    1.7577    2.3733
# hypertensive1                             1.1897     0.8405    1.0982    1.2889
# cholesterol1                              1.1064     0.9038    1.0139    1.2073
# livesalone1                               1.1818     0.8462    1.0858    1.2863
# Sex1                                      1.2364     0.8088    1.1482    1.3314
# 
# Concordance= 0.778  (se = 0.004 )
# Likelihood ratio test= 3276  on 15 df,   p=<2e-16
# Wald test            = 2904  on 15 df,   p=<2e-16
# Score (logrank) test = 3883  on 15 df,   p=<2e-16

nlinvars <- c("Age_when_attended_assesment_centre_0_0","age2o","age3o","family_history_of_dementia","education_years",
              "education2o","education3o","Diabetes_BIN_FINAL_0_0",
              "Townsend_deprivation_modelvar","current_history_depression","stroke_TIA_BIN_FINAL","hypertensive","cholesterol",
              "livesalone","Sex","dementia_BIN_surv")
crr_status <- train.data$crr_status
time_at_risk <- train.data$time_at_risk
nlincovs <- model.matrix(as.formula(UKBDRS_LASSO_nonlinear), train.data[nlinvars])[,-1]

cr.train.nlin <- crr(ftime = time_at_risk, fstatus = crr_status, failcode=1, cencode = 0, cov1=nlincovs, variance=TRUE)
save(cr.train.nlin, file=paste0(save_pathway,"ukbdrs.nonlin.cr.train.rda"))
summary(cr.train.nlin)
# coef exp(coef) se(coef)      z p-value
# Age_when_attended_assesment_centre_0_0  0.163572     1.178 5.82e-03 28.11 0.0e+00
# family_history_of_dementia1             0.436469     1.547 4.30e-02 10.15 0.0e+00
# age2o                                   0.003606     1.004 8.91e-04  4.05 5.2e-05
# age3o                                  -0.000600     0.999 1.50e-04 -4.00 6.2e-05
# education2o                             0.003490     1.003 1.05e-03  3.31 9.4e-04
# education3o                             0.000172     1.000 8.24e-05  2.09 3.7e-02
# education_years                        -0.040186     0.961 6.18e-03 -6.50 8.0e-11
# Diabetes_BIN_FINAL_0_01                 0.534182     1.706 5.75e-02  9.29 0.0e+00
# Townsend_deprivation_modelvar1          0.220785     1.247 4.33e-02  5.10 3.4e-07
# current_history_depression1             0.554889     1.742 4.62e-02 12.01 0.0e+00
# stroke_TIA_BIN_FINAL1                   0.651614     1.919 7.78e-02  8.37 0.0e+00
# hypertensive1                           0.156246     1.169 4.12e-02  3.79 1.5e-04
# cholesterol1                            0.101889     1.107 4.53e-02  2.25 2.4e-02
# livesalone1                             0.140500     1.151 4.31e-02  3.26 1.1e-03
# Sex1                                    0.165356     1.180 3.79e-02  4.36 1.3e-05

#mape age2o..etc into test
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


#build linear predictor, predicted probability for train and test
train.data$NLIN_UKBDRS_LASSO_linear_predictor <- 0.163572*(train.data$Age_when_attended_assesment_centre_0_0 - mean(train.data$Age_when_attended_assesment_centre_0_0)) +
  0.436469*train.data$UKBDRS_familyhistory + -0.040186*(train.data$education_years - mean(train.data$education_years)) +
  0.534182*train.data$UKBDRS_diabetes + 0.220785*train.data$UKBDRS_townsend_group5 +
  0.554889*train.data$UKBDRS_depression+ 0.651614*train.data$UKBDRS_stroke +
  0.156246*train.data$UKBDRS_hypertensive + 0.101889*train.data$UKBDRS_cholesterol +
  0.140500*train.data$UKBDRS_livesalone + 0.165356*train.data$UKBDRS_sex_score +
  0.003606*train.data$age2o + -0.000600*train.data$age3o +
  0.003490*train.data$education2o + 0.000172*train.data$education3o 
summary(train.data$NLIN_UKBDRS_LASSO_linear_predictor)  
#note age2o, age3o, education2o, education3o have mean 0 so no need to subtract training data mean

train.data$NLIN_UKBDRS_LASSO_predicted_prob <- 1 - 0.9916195^exp(train.data$NLIN_UKBDRS_LASSO_linear_predictor)  
summary(train.data$NLIN_UKBDRS_LASSO_predicted_prob)
plot(train.data$UKBDRS_LASSO_crr_predicted_prob, train.data$NLIN_UKBDRS_LASSO_predicted_prob)


#test
test.data$NLIN_UKBDRS_LASSO_linear_predictor <- 0.163572*(test.data$Age_when_attended_assesment_centre_0_0 - mean(train.data$Age_when_attended_assesment_centre_0_0)) +
  0.436469*test.data$UKBDRS_familyhistory + -0.040186*(test.data$education_years - mean(train.data$education_years)) +
  0.534182*test.data$UKBDRS_diabetes + 0.220785*test.data$UKBDRS_townsend_group5 +
  0.554889*test.data$UKBDRS_depression+ 0.651614*test.data$UKBDRS_stroke +
  0.156246*test.data$UKBDRS_hypertensive + 0.101889*test.data$UKBDRS_cholesterol +
  0.140500*test.data$UKBDRS_livesalone + 0.165356*test.data$UKBDRS_sex_score +
  0.003606*test.data$age2o + -0.000600*test.data$age3o +
  0.003490*test.data$education2o + 0.000172*test.data$education3o 
summary(test.data$NLIN_UKBDRS_LASSO_linear_predictor)  
#note age2o, age3o, education2o, education3o have mean 0 so no need to subtract training data mean

test.data$NLIN_UKBDRS_LASSO_predicted_prob <- 1 - 0.9916195^exp(test.data$NLIN_UKBDRS_LASSO_linear_predictor)  
summary(test.data$NLIN_UKBDRS_LASSO_predicted_prob)

library(riskRegression)
auc_test<-Score(list('UKBDRS'=test.data$UKBDRS_LASSO_crr_predicted_prob,
                      'nlin'=test.data$NLIN_UKBDRS_LASSO_predicted_prob),
                 formula=Hist(time_at_risk,crr_status)~1,
                 data = test.data,
                 null.model = FALSE,
                 conf.int = TRUE,
                 times = c(365.25*14),
                 plots="ROC",
                 metrics="AUC",
                 cens.model = "cox",
                 conservative=FALSE,
                 censoring.save.memory=FALSE,
                 contrasts = list(c(1,2)))
auc_test$AUC$score
# model  times       AUC          se     lower     upper
# 1: UKBDRS 5113.5 0.8010068 0.008737752 0.7838811 0.8181325
# 2:   nlin 5113.5 0.8003450 0.008779824 0.7831368 0.8175531
auc_test$AUC$contrasts
# times model reference     delta.AUC          se        lower       upper         p
# 1: 5113.5  nlin    UKBDRS -0.0006618058 0.000882729 -0.002391923 0.001068311 0.4534191
