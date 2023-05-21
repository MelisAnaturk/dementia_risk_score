library(tidyverse)
library(broom)
library(rcompanion)
library(data.table)

library(survival)
library(cmprsk)
library(riskRegression)

library(ggplot2)
library(tidyverse)
library(broom)
library(data.table)
library(dplyr)
library(tidyr)


load(file="../../raw_data/modelvar/12_train_data_outliers_removed_cox_crr_fitted.rda")
load(file="../../raw_data/modelvar/12_test_data_outliers_removed_cox_crr_fitted.rda")
set.seed(9)


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

#### Train auc ####
auc_train<-Score(list('UKBDRS'=train.data$UKBDRS_LASSO_crr_predicted_prob,
                'DRS'=train.data$DRS_predicted_prob,
                'Age_only'=train.data$age_only_crr_predicted_prob,
               'CAIDE'=train.data$CAIDE_predicted_prob),
           formula=Hist(time_at_risk,crr_status)~1,
           data = train.data,
           null.model = FALSE,
           conf.int = TRUE,
           times = c(365.25*14),
           plots="ROC",
           metrics="AUC",
           cens.model = "cox",
           conservative=FALSE,
           censoring.save.memory=FALSE,
           contrasts = list(c(1,2,3,4)))
auc_train$AUC$score
# model  times       AUC           se     lower     upper
# 1:   UKBDRS 5113.5 0.7853066 0.0002246764 0.7848663 0.7857470
# 2:      DRS 5113.5 0.7595194 0.0002254273 0.7590775 0.7599612
# 3: Age_only 5113.5 0.7539931 0.0002227548 0.7535565 0.7544297
# 4:    CAIDE 5113.5 0.5998509 0.0002563156 0.5993485 0.6003532

auc_train$AUC$contrasts
# times    model reference   delta.AUC           se       lower       upper p
# 1: 5113.5      DRS    UKBDRS -0.02578724 1.114554e-04 -0.02600569 -0.02556879 0
# 2: 5113.5 Age_only    UKBDRS -0.03131353 9.962497e-05 -0.03150879 -0.03111827 0
# 3: 5113.5    CAIDE    UKBDRS -0.18545574 2.862850e-04 -0.18601685 -0.18489463 0

df_auc <- as.data.frame(auc_train$AUC$score)[,c("model","AUC","lower","upper")]

#anu subset
anu_train <- train.data[which(!is.na(train.data$ANU_ADRI)),] #n=174982
auc_train_anu<-Score(list('ukbdrs_anu'=anu_train$UKBDRS_LASSO_crr_predicted_prob,
                      'ANUADRI'=anu_train$ANU_ADRI),
                 formula=Hist(time_at_risk,crr_status)~1,
                 data = anu_train,
                 null.model = FALSE,
                 conf.int = TRUE,
                 times = c(365.25*14),
                 plots="ROC",
                 metrics="AUC",
                 cens.model = "cox",
                 conservative=FALSE,
                 censoring.save.memory=FALSE,
                 contrasts = list(c(1,2)))
auc_train_anu$AUC$score
auc_train_anu$AUC$contrasts

df_auc <- rbind(df_auc, as.data.frame(auc_train_anu$AUC$score)[2,c("model","AUC","lower","upper")])
write.csv(df_auc, file="../results/crr_auc_ukb_train.csv")

df_train_auc_compare <- as.data.frame(auc_train$AUC$contrasts[,c("model","reference","delta.AUC",
                                                                 "lower","upper","p")])

df_train_auc_compare <- rbind(df_train_auc_compare, as.data.frame(auc_train_anu$AUC$contrasts[1,c("model","reference","delta.AUC",
                                                                                             "lower","upper","p")]))
df_train_auc_compare$p_fdr <- p.adjust(df_train_auc_compare$p, method="BH")
write.csv(df_train_auc_compare, file="../results/auc_compare_training.csv")


#apoe
train.apoe <- train.data[which(!is.na(train.data$UKBDRS_APOE_LASSO_crr_predicted_prob)),]
auc_train_apoe<-Score(list('ukbdrs_apoe'=train.apoe$UKBDRS_APOE_LASSO_crr_predicted_prob),
                     formula=Hist(time_at_risk,crr_status)~1,
                     data = train.apoe,
                     null.model = FALSE,
                     conf.int = TRUE,
                     times = c(365.25*14),
                     plots="ROC",
                     metrics="AUC",
                     cens.model = "cox",
                     conservative=FALSE,
                     censoring.save.memory=FALSE)
auc_train_apoe$AUC$score
# model  times       AUC           se     lower     upper
# 1: ukbdrs_anu 5113.5 0.8115939 0.0005506572 0.8105146 0.8126731

#### test data ####
auc_test<-Score(list('UKBDRS'=test.data$UKBDRS_LASSO_crr_predicted_prob,
                      'DRS'=test.data$DRS_predicted_prob,
                      'Age_only'=test.data$age_only_crr_predicted_prob,
                      'CAIDE'=test.data$CAIDE_predicted_prob),
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
                 contrasts = list(c(1,2,3,4)))
auc_test$AUC$score
auc_test$AUC$contrasts
# times    model reference   delta.AUC          se       lower       upper            p
# 1: 5113.5      DRS    UKBDRS -0.02717938 0.005432618 -0.03782711 -0.01653164 5.644562e-07
# 2: 5113.5 Age_only    UKBDRS -0.02772356 0.006262003 -0.03999686 -0.01545026 9.543466e-06
# 3: 5113.5    CAIDE    UKBDRS -0.20020435 0.014674655 -0.22896614 -0.17144255 2.226380e-42

df_auc <- as.data.frame(auc_test$AUC$score)[,c("model","AUC","lower","upper")]

#anu subset
anu_test <- test.data[which(!is.na(test.data$ANU_ADRI)),] #n = 43749
auc_test_anu<-Score(list('ukbdrs_anu'=anu_test$UKBDRS_LASSO_crr_predicted_prob,
                          'ANUADRI'=anu_test$ANU_ADRI),
                     formula=Hist(time_at_risk,crr_status)~1,
                     data = anu_test,
                     null.model = FALSE,
                     conf.int = TRUE,
                     times = c(365.25*14),
                     plots="ROC",
                     metrics="AUC",
                     cens.model = "cox",
                     conservative=FALSE,
                     censoring.save.memory=FALSE,
                     contrasts = list(c(1,2)))
auc_test_anu$AUC$score
auc_test_anu$AUC$contrasts

df_auc <- rbind(df_auc, as.data.frame(auc_test_anu$AUC$score)[2,c("model","AUC","lower","upper")])
write.csv(df_auc, file="../results/crr_auc_ukb_test.csv")

df_test_auc_compare <- as.data.frame(auc_test$AUC$contrasts[,c("model","reference","delta.AUC",
                                                                 "lower","upper","p")])

df_test_auc_compare <- rbind(df_test_auc_compare, as.data.frame(auc_test_anu$AUC$contrasts[1,c("model","reference","delta.AUC",
                                                                                                  "lower","upper","p")]))
df_test_auc_compare$p_fdr <- p.adjust(df_test_auc_compare$p, method="BH")
write.csv(df_test_auc_compare, file="../results/auc_compare_testing.csv")


plot_auc_obj <- auc_test
plot_auc_obj$ROC$plotframe <- rbind(plot_auc_obj$ROC$plotframe, auc_test_anu$ROC$plotframe)
plot_auc_obj$AUC$score <- rbind(plot_auc_obj$AUC$score , auc_test_anu$AUC$score )
plotROC(plot_auc_obj, models = c("UKBDRS","DRS","Age_only","CAIDE","ANUADRI"))

plotROC(plot_auc_obj, models = c("UKBDRS","DRS","Age_only","CAIDE","ANUADRI"),auc.in.legend = FALSE)

png(file="../results/roc_plotted_all.png",
    height = 600, width = 600, res=120)
plotROC(plot_auc_obj, models = c("UKBDRS","DRS","Age_only","CAIDE","ANUADRI"),auc.in.legend = FALSE)
dev.off()



test.apoe <- test.data[which(!is.na(test.data$UKBDRS_APOE_LASSO_crr_predicted_prob)),]
auc_test_apoe<-Score(list('ukbdrs_apoe'=test.apoe$UKBDRS_APOE_LASSO_crr_predicted_prob),
                      formula=Hist(time_at_risk,crr_status)~1,
                      data = test.apoe,
                      null.model = FALSE,
                      conf.int = TRUE,
                      times = c(365.25*14),
                      plots="ROC",
                      metrics="AUC",
                      cens.model = "cox",
                      conservative=FALSE,
                      censoring.save.memory=FALSE)
auc_test_apoe$AUC$score
# model  times       AUC          se     lower     upper
# 1: ukbdrs_anu 5113.5 0.8260386 0.009626646 0.8071707 0.8449065



#### calibration ####
#use Score and plotCalibration
calib_test<-Score(list('UKBDRS'=test.data$UKBDRS_LASSO_crr_predicted_prob),
                formula=Hist(time_at_risk,crr_status)~1,
                data = test.data,
                null.model = FALSE,
                conf.int = TRUE,
                times = c(365.25*14),
                plots="calibration",
                metrics="brier",
                cens.model = "cox",
                conservative=FALSE,
                censoring.save.memory=FALSE,)
plotCalibration(calib_test,cens.method = "local", method = "quantile", q=10,xlim = (c(0,0.15)), ylim=c(0,0.15))
#line matches diagonal

png(file="../results/ukbdrs_calib.png",
    height = 800, width = 800, res=120)
plotCalibration(calib_test,cens.method = "local", method = "quantile", q=10,xlim = (c(0,0.15)), ylim=c(0,0.15))
dev.off()

#plotCalibration can return plottable object
#get the points plotted, model an int and lope
x<-plotCalibration(calib_test,cens.method = "local", method = "quantile", q=10,xlim = (c(0,0.5)), ylim=c(0,0.5), plot=FALSE)

mod <- lm(x$plotFrames$UKBDRS$Obs ~ x$plotFrames$UKBDRS$Pred)
summary(mod)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)              -0.0008998  0.0013042   -0.69     0.51    
# x$plotFrames$UKBDRS$Pred  1.0250295  0.0411657   24.90 7.24e-09 ***

confint(mod)
# 2.5 %      97.5 %
#   (Intercept)              -0.003907166 0.002107643
# x$plotFrames$UKBDRS$Pred  0.930101201 1.119957834

#now do apoe
calib_test_apoe<-Score(list('UKBDRS'=test.apoe$UKBDRS_APOE_LASSO_crr_predicted_prob),
                  formula=Hist(time_at_risk,crr_status)~1,
                  data = test.apoe,
                  null.model = FALSE,
                  conf.int = TRUE,
                  times = c(365.25*14),
                  plots="calibration",
                  metrics="brier",
                  cens.model = "cox",
                  conservative=FALSE,
                  censoring.save.memory=FALSE,)
plotCalibration(calib_test_apoe,cens.method = "local", method = "quantile", q=10,xlim = (c(0,0.15)), ylim=c(0,0.15))
#line matches diagonal

png(file="../results/ukbdrs.apoe_calib.png",
    height = 800, width = 800, res=120)
plotCalibration(calib_test_apoe,cens.method = "local", method = "quantile", q=10,xlim = (c(0,0.15)), ylim=c(0,0.15))
dev.off()


#plotCalibration can return plottable object
#get the points plotted, model an int and lope
x<-plotCalibration(calib_test_apoe,cens.method = "local", method = "quantile", q=10,xlim = (c(0,0.5)), ylim=c(0,0.5), plot=FALSE)

mod <- lm(x$plotFrames$UKBDRS$Obs ~ x$plotFrames$UKBDRS$Pred)
summary(mod)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)              0.001312   0.001876   0.699    0.504    
# x$plotFrames$UKBDRS$Pred 0.968970   0.051149  18.944 6.24e-08 ***

confint(mod)
# 2.5 %      97.5 %
#   (Intercept)              -0.00301432 0.005637377
# x$plotFrames$UKBDRS$Pred  0.85101978 1.086920921


rm(anu_train, anu_test)
#sens spec thresholds
#these can be found in the ROC plotframe of the auc output
df_sens_spec <- rbind(auc_train$ROC$plotframe, auc_train_apoe$ROC$plotframe)

#only care about ukbdrs
df_sens_spec <- df_sens_spec[which(df_sens_spec$model=="UKBDRS" |
                                     df_sens_spec$model=="ukbdrs_apoe")]
#we have sensitivity already as TPR
#we have FPR which can give specificity
df_sens_spec$sensitivity <- df_sens_spec$TPR
df_sens_spec$specificity <- 1 - df_sens_spec$FPR

#how about ppv npv..
summary(train.data$dementia_BIN_TOTAL)
prev <- 3051/176611

df_sens_spec$ppv <- (df_sens_spec$sensitivity*prev) /
  (df_sens_spec$sensitivity*prev + (1 - df_sens_spec$specificity)*(1-prev))

df_sens_spec$npv <- (df_sens_spec$specificity*(1-prev)) /
  (df_sens_spec$specificity*(1-prev) + (1 - df_sens_spec$sensitivity)*prev)


#Sensitivity = 80
threshold <-  df_sens_spec[df_sens_spec$sensitivity >= .799 & df_sens_spec$sensitivity <= 0.801, ] %>% mutate_if(is.numeric, ~round(., 6))
print(threshold[which(threshold$model=="UKBDRS"),c("model","risk","sensitivity","specificity","ppv","npv")])
#UKBDRS 0.018432    0.800150    0.643462 0.037954 0.994570
print(threshold[which(threshold$model=="ukbdrs_apoe"),c("model","risk","sensitivity","specificity","ppv","npv")])
#ukbdrs_apoe ukbdrs_apoe 0.018620    0.800089    0.676096 0.041615 0.994829

#Sensitivity = 85
threshold <-  df_sens_spec[df_sens_spec$sensitivity >= .849 & df_sens_spec$sensitivity <= 0.8501, ] %>% mutate_if(is.numeric, ~round(., 6))
print(threshold[which(threshold$model=="UKBDRS"),c("model","risk","sensitivity","specificity","ppv","npv")])
#UKBDRS 0.014506    0.849899    0.569271 0.033523 0.995386
print(threshold[which(threshold$model=="ukbdrs_apoe"),c("model","risk","sensitivity","specificity","ppv","npv")])
#ukbdrs_apoe  0.014577    0.849780    0.606564 0.036580 0.995665

#Sensitivity = 90
threshold <-  df_sens_spec[df_sens_spec$sensitivity >= .8995 & df_sens_spec$sensitivity <= 0.901, ] %>% mutate_if(is.numeric, ~round(.,6))
print(threshold[which(threshold$model=="UKBDRS"),c("model","risk","sensitivity","specificity","ppv","npv")])
#UKBDRS 0.011402    0.900092    0.499999 0.030675 0.996500
print(threshold[which(threshold$model=="ukbdrs_apoe"),c("model","risk","sensitivity","specificity","ppv","npv")])
#ukbdrs_apoe0.010014    0.900023    0.504218 0.030925 0.996527

#Sensitivity = 95
threshold <-  df_sens_spec[df_sens_spec$sensitivity >= .949 & df_sens_spec$sensitivity <= 0.9501, ] %>% mutate_if(is.numeric, ~round(., 6))
print(threshold[which(threshold$model=="UKBDRS"),c("model","risk","sensitivity","specificity","ppv","npv")])
#UKBDRS 0.006714    0.949919    0.344727 0.024850 0.997453
print(threshold[which(threshold$model=="ukbdrs_apoe"),c("model","risk","sensitivity","specificity","ppv","npv")])
#ukbdrs_apoe 0.005672    0.949932    0.342535 0.024770 0.997437







#Specificity = 80
threshold <-  df_sens_spec[df_sens_spec$specificity >= .799 & df_sens_spec$specificity <= 0.801, ] %>% mutate_if(is.numeric, ~round(., 5))
print(threshold[which(threshold$model=="UKBDRS"),c("model","risk","sensitivity","specificity","ppv","npv")])
#UKBDRS 0.03201     0.58450     0.80015 0.04890 0.99095
print(threshold[which(threshold$model=="ukbdrs_apoe"),c("model","risk","sensitivity","specificity","ppv","npv")])
#ukbdrs_apoe 0.03038     0.66103     0.80004 0.05492 0.99261

#Specificity = 85
threshold <-  df_sens_spec[df_sens_spec$specificity >= .849 & df_sens_spec$specificity <= 0.8505, ] %>% mutate_if(is.numeric, ~round(., 5))
print(threshold[which(threshold$model=="UKBDRS"),c("model","risk","sensitivity","specificity","ppv","npv")])
#UKBDRS 0.03848     0.49284     0.85033 0.05472 0.98962
print(threshold[which(threshold$model=="ukbdrs_apoe"),c("model","risk","sensitivity","specificity","ppv","npv")])
#kbdrs_apoe  0.03746     0.59061     0.85001 0.06474 0.99160

#Specificity = 90
threshold <-  df_sens_spec[df_sens_spec$specificity >= .8995 & df_sens_spec$specificity <= 0.901, ] %>% mutate_if(is.numeric, ~round(., 5))
print(threshold[which(threshold$model=="UKBDRS"),c("model","risk","sensitivity","specificity","ppv","npv")])
#UKBDRS 0..04804     0.36777     0.90009 0.06077 0.98780
print(threshold[which(threshold$model=="ukbdrs_apoe"),c("model","risk","sensitivity","specificity","ppv","npv")])
#ukbdrs_apoe 0.05052     0.47431     0.90005 0.07700 0.98984

#Specificity = 95
threshold <-  df_sens_spec[df_sens_spec$specificity >= .949 & df_sens_spec$specificity <= 0.9501, ] %>% mutate_if(is.numeric, ~round(., 5))
print(threshold[which(threshold$model=="UKBDRS"),c("model","risk","sensitivity","specificity","ppv","npv")])
#UKBDRS 0.06373     0.23446     0.95000 0.07616 0.98603
print(threshold[which(threshold$model=="ukbdrs_apoe"),c("model","risk","sensitivity","specificity","ppv","npv")])
#ukbdrs_apoe 0.07563     0.31919     0.94998 0.10086 0.98756
