library(tidyverse)
library(pROC)
library(CalibrationCurves)


load(file="../../raw_data/train_data_outliers_removed_fiftyplusnoapoe_sexstratify.rda")
load(file="../../raw_data/test_data_outliers_removed_fiftyplusnoapoe_sexstratify.rda")


#lets look at performance in predicting
#dementia in 1-5 years
#dementia 5-10 years
#dementia 10-14

summary(as.factor(test.data$dementia_BIN_TOTAL))
#0     1 
#43389   762 

summary(test.data$years_diff_all_cause_dementia_0_0)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   1.10    6.79    9.24    8.81   11.33   14.19   43389 

hist(train.data$years_diff_all_cause_dementia_0_0)
hist(test.data$years_diff_all_cause_dementia_0_0)

#### definitions ####
#### early ####
test.data$early_dementia<-ifelse(test.data$dementia_BIN_TOTAL==1,
                                 ifelse(test.data$years_diff_all_cause_dementia_0_0<=5,1,0),0)
summary(as.factor(test.data$early_dementia))                                  
#0     1 
#44050   101 
summary(test.data[which(test.data$early_dementia==1),"years_diff_all_cause_dementia_0_0"])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.103   2.642   3.406   3.374   4.205   4.997 
hist(test.data[which(test.data$early_dementia==1),"years_diff_all_cause_dementia_0_0"])

#### mid ####
test.data$mid_dementia<-ifelse(test.data$dementia_BIN_TOTAL==1,
                                 ifelse(test.data$years_diff_all_cause_dementia_0_0>5 &
                                          test.data$years_diff_all_cause_dementia_0_0<=10,1,0),0)
summary(as.factor(test.data$mid_dementia))                                  
#0     1 
#43802   349 
summary(test.data[which(test.data$mid_dementia==1),"years_diff_all_cause_dementia_0_0"])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#5.060   6.763   8.096   7.863   9.052   9.991 
hist(test.data[which(test.data$mid_dementia==1),"years_diff_all_cause_dementia_0_0"])


#### late ####
test.data$late_dementia<-ifelse(test.data$dementia_BIN_TOTAL==1,
                               ifelse(test.data$years_diff_all_cause_dementia_0_0>10,1,0),0)
summary(as.factor(test.data$late_dementia))                                  
#0     1 
#43839   312 
summary(test.data[which(test.data$late_dementia==1),"years_diff_all_cause_dementia_0_0"])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#10.01   10.89   11.58   11.64   12.25   14.19 
hist(test.data[which(test.data$late_dementia==1),"years_diff_all_cause_dementia_0_0"])



#### assess performance ####
models <- c("age_only", "UKBDRS_LASSO", "UKBDRS_APOE_LASSO", "CAIDE", "DRS")

#dataframe for storing/organizing auc results for each model/time frame
df_auc_timetodx<-data.frame(matrix(ncol=4))
names(df_auc_timetodx)<-c("Model","early_dementia","mid_dementia","late_dementia")

#cycle through train/test for each model, store auc
for (m in models){
  df_auc<-data.frame(matrix(ncol=4))
  names(df_auc)<-c("Model","early_dementia","mid_dementia","late_dementia")
  df_auc$Model<-m
  for (t in c("early_dementia","mid_dementia","late_dementia")){
    data <- test.data
    print(paste0('computing AUC for ', m, ' in ', t))
    # AUC part
    roc <- pROC::roc(data[, c(t)], data[, paste(m, "predicted_prob", sep="_")], plot=TRUE, smooth = FALSE, ci=TRUE)
    print(roc$auc)
    print(roc$ci)
    df_auc[t]<-paste(round(roc$auc[1],5), " [", round(roc$ci[1],5), ",", round(roc$ci[3],5),"]",sep="")
    print('========================================================')
    print('========================================================')
  }
  df_auc_timetodx<-rbind(df_auc_timetodx, df_auc)
}

# for anu_adri only
df_auc<-data.frame(matrix(ncol=4))
names(df_auc)<-c("Model","early_dementia","mid_dementia","late_dementia")
df_auc$Model<-"ANU_ADRI"
for (t in c("early_dementia","mid_dementia","late_dementia")){
  data <- test.data
  print(paste0('AUC for anu-adri in ', t))  
  roc <- pROC::roc(data[, c(t)], data[, c("ANU_ADRI")], plot=TRUE, smooth = FALSE, ci=TRUE)
  print(roc$auc)
  print(roc$ci)
  df_auc[t]<-paste(round(roc$auc[1],5), " [", round(roc$ci[1],5), ",", round(roc$ci[3],5),"]",sep="")
}
df_auc_timetodx<-rbind(df_auc_timetodx, df_auc)

