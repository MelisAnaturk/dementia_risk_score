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


load(file="../../raw_data/12_train_data_outliers_removed_cox_crr_fitted.rda")
load(file="../../raw_data/12_test_data_outliers_removed_cox_crr_fitted.rda")
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
           times = c(365.25*14.2),
           plots="ROC",
           metrics="AUC",
           cens.model = "cox",
           conservative=FALSE,
           censoring.save.memory=FALSE,
           contrasts = list(c(1,2,3,4)))
auc_train$AUC$score
# model   times       AUC           se     lower     upper
# 1: ukbdrs 5186.55 0.7878599 0.0003003121 0.7872713 0.7884485
# 2:    drs 5186.55 0.7606576 0.0002831628 0.7601027 0.7612126
# 3:    age 5186.55 0.7599956 0.0002882494 0.7594306 0.7605605
# 4:  caide 5186.55 0.6178982 0.0003268553 0.6172576 0.6185388

auc_train$AUC$contrasts
# times model reference   delta.AUC           se       lower       upper p
# 1: 5186.55   drs    ukbdrs -0.02720225 0.0001339829 -0.02746486 -0.02693965 0
# 2: 5186.55   age    ukbdrs -0.02786432 0.0001130753 -0.02808595 -0.02764270 0
# 3: 5186.55 caide    ukbdrs -0.16996169 0.0003427209 -0.17063341 -0.16928996 0

df_auc <- as.data.frame(auc_train$AUC$score)[,c("model","AUC","lower","upper")]

#anu subset
anu_train <- train.data[which(!is.na(train.data$ANU_ADRI)),]
auc_train_anu<-Score(list('ukbdrs_anu'=anu_train$UKBDRS_LASSO_crr_predicted_prob,
                      'ANUADRI'=anu_train$ANU_ADRI),
                 formula=Hist(time_at_risk,crr_status)~1,
                 data = anu_train,
                 null.model = FALSE,
                 conf.int = TRUE,
                 times = c(365.25*14.2),
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



#### test data ####
auc_test<-Score(list('UKBDRS'=test.data$UKBDRS_LASSO_crr_predicted_prob,
                      'DRS'=test.data$DRS_predicted_prob,
                      'Age_only'=test.data$age_only_crr_predicted_prob,
                      'CAIDE'=test.data$CAIDE_predicted_prob),
                 formula=Hist(time_at_risk,crr_status)~1,
                 data = test.data,
                 null.model = FALSE,
                 conf.int = TRUE,
                 times = c(365.25*14.2),
                 plots="ROC",
                 metrics="AUC",
                 cens.model = "cox",
                 conservative=FALSE,
                 censoring.save.memory=FALSE,
                 contrasts = list(c(1,2,3,4)))
auc_test$AUC$score
auc_test$AUC$contrasts
# times model reference   delta.AUC           se       lower       upper p
# 1: 5186.55   drs    ukbdrs -0.02720225 0.0001339829 -0.02746486 -0.02693965 0
# 2: 5186.55   age    ukbdrs -0.02786432 0.0001130753 -0.02808595 -0.02764270 0
# 3: 5186.55 caide    ukbdrs -0.16996169 0.0003427209 -0.17063341 -0.16928996 0

df_auc <- as.data.frame(auc_test$AUC$score)[,c("model","AUC","lower","upper")]

#anu subset
anu_test <- test.data[which(!is.na(test.data$ANU_ADRI)),]
auc_test_anu<-Score(list('ukbdrs_anu'=anu_test$UKBDRS_LASSO_crr_predicted_prob,
                          'ANUADRI'=anu_test$ANU_ADRI),
                     formula=Hist(time_at_risk,crr_status)~1,
                     data = anu_test,
                     null.model = FALSE,
                     conf.int = TRUE,
                     times = c(365.25*14.2),
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

