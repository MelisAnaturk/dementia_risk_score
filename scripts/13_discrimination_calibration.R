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
auc_train_apoe<-Score(list('ukbdrs_anu'=train.apoe$UKBDRS_APOE_LASSO_crr_predicted_prob),
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
anu_test <- test.data[which(!is.na(test.data$ANU_ADRI)),]
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
auc_test_apoe<-Score(list('ukbdrs_anu'=test.apoe$UKBDRS_APOE_LASSO_crr_predicted_prob),
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

#plotCalibration can return plottable object
#get the points plotted, model an int and lope
x<-plotCalibration(calib_test,cens.method = "local", method = "quantile", q=10,xlim = (c(0,0.5)), ylim=c(0,0.5), plot=FALSE)

mod <- lm(x$plotFrames$UKBDRS$Obs ~ x$plotFrames$UKBDRS$Pred)
summary(mod)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)              -0.000932   0.001298  -0.718    0.493      #intercept
# x$plotFrames$UKBDRS$Pred  0.966367   0.038600  25.035 6.93e-09 ***  #slope

confint(mod)
# 2.5 %     97.5 %
#   (Intercept)              -0.003925378 0.00206128
# x$plotFrames$UKBDRS$Pred  0.877355155 1.05537918

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
plotCalibration(calib_test_apoe,cens.method = "jackknife", method = "quantile", q=10,xlim = (c(0,0.5)), ylim=c(0,0.5))
#line matches diagonal

#plotCalibration can return plottable object
#get the points plotted, model an int and lope
x<-plotCalibration(calib_test_apoe,cens.method = "jackknife", method = "quantile", q=10,xlim = (c(0,0.5)), ylim=c(0,0.5), plot=FALSE)

mod <- lm(x$plotFrames$UKBDRS$Obs ~ x$plotFrames$UKBDRS$Pred)
summary(mod)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)              0.001276   0.001865   0.684    0.513    
# x$plotFrames$UKBDRS$Pred 0.926982   0.048621  19.065 5.93e-08 ***