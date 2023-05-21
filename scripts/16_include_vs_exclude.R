library(broom)
library(data.table)
library(tidyverse)

#compare demographics in excluded vs included data

#get raw data
load("../../raw_data/df_ukb_raw.rda")
rawdata_demo<-df_ukb_raw[,c("eid","Age_when_attended_assesment_centre_0_0","Sex","education_years")]
rm(df_ukb_raw)

#subset to 50+
rawdata_demo<-subset(rawdata_demo, rawdata_demo$Age_when_attended_assesment_centre_0_0>=50)
#n=384642

#get analysis sample
load("../../raw_data/modelvar/12_test_data_outliers_removed_cox_crr_fitted.rda")
load("../../raw_data/modelvar/12_train_data_outliers_removed_cox_crr_fitted.rda")
df_analysis<-rbind(train.data[,names(test.data)],test.data)
rm(train.data, test.data)

analysis_sample<-data.frame(df_analysis[,c("eid")])
rm(df_analysis)
analysis_sample$dataset<-"included"
names(analysis_sample)<-c("eid","dataset")
analysis_sample$eid<-as.character(analysis_sample$eid)

#join full sample and analysis on id
df_merged<- list(rawdata_demo, analysis_sample) %>% reduce(left_join, by="eid")
#check
View(df_merged)
View(analysis_sample)

df_merged$dataset[is.na(df_merged$dataset)]<-"excluded"
df_merged$dataset<-as.factor(df_merged$dataset)
summary(df_merged$dataset)
#excluded included 
#163880   220762 
#cool

excluded_sample <- subset(df_merged, df_merged$dataset=="excluded")
included_sample <- subset(df_merged, df_merged$dataset=="included")
age_test<-t.test(Age_when_attended_assesment_centre_0_0 ~ dataset, data = df_merged, var.equal = FALSE)
#data:  Age_when_attended_assesment_centre_0_0 by dataset
#t = 13.31, df = 352022, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  0.2017481 0.2714272
#sample estimates:
#  mean in group excluded mean in group included 
#60.21101               59.97442 
sd(excluded_sample$Age_when_attended_assesment_centre_0_0) #5.46
sd(included_sample$Age_when_attended_assesment_centre_0_0) #5.43


edu_test<-t.test(education_years ~ dataset, data = df_merged, var.equal = FALSE)
edu_test
#data:  education_years by dataset
#t = -55.771, df = 339144, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.5825166 -0.5429639
#sample estimates:
#  mean in group excluded mean in group included 
#12.98370               13.54644  
sd(excluded_sample$education_years, na.rm=TRUE) #3.07
sd(included_sample$education_years) #3.05

df_merged$Sex<-as.factor(df_merged$Sex)
sex_test<-chisq.test(df_merged$dataset, df_merged$Sex, correct=FALSE)
sex_test
#Pearson's Chi-squared test
#
#data:  df_merged$dataset and df_merged$Sex
#X-squared = 1819.4, df = 1, p-value < 2.2e-16
summary(as.factor(excluded_sample$Sex))
#0      1 
#95444 68436 
summary(as.factor(included_sample$Sex))
#0      1 
#113276 107486 