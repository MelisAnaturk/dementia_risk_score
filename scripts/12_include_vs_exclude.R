library(broom)
library(data.table)
library(tidyverse)

#compare demographics in excluded vs included data

#get raw data
load("../../raw_data/df_ukb_raw.rda")
rawdata_demo<-df_ukb_raw[,c("eid","Age_when_attended_assesment_centre_0_0","Sex","education_years")]
rm(df_ukb_raw)

#get analysis sample
load("../../raw_data/train_data_outliers_removed.rda")
load("../../raw_data/test_data_outliers_removed.rda")
df_analysis<-rbind(train.data,test.data)
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
#294550   207971 
#cool

excluded_sample <- subset(df_merged, df_merged$dataset=="excluded")
included_sample <- subset(df_merged, df_merged$dataset=="included")
age_test<-t.test(Age_when_attended_assesment_centre_0_0 ~ dataset, data = df_merged, var.equal = FALSE)
#data:  Age_when_attended_assesment_centre_0_0 by dataset
#t = 19.779, df = 445169, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  0.4136474 0.5046456
#sample estimates:
#  mean in group excluded mean in group included 
#56.71851               56.25936 
sd(excluded_sample$Age_when_attended_assesment_centre_0_0) #8.06
sd(included_sample$Age_when_attended_assesment_centre_0_0) #8.14


edu_test<-t.test(education_years ~ dataset, data = df_merged, var.equal = FALSE)
edu_test
#data:  education_years by dataset
#t = -56.42, df = 451183, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.5127060 -0.4782804
#sample estimates:
#  mean in group excluded mean in group included 
#13.30331               13.79880 
sd(excluded_sample$education_years, na.rm=TRUE) #3.07
sd(included_sample$education_years) #3.03

df_merged$Sex<-as.factor(df_merged$Sex)
sex_test<-chisq.test(df_merged$dataset, df_merged$Sex, correct=FALSE)
sex_test
#Pearson's Chi-squared test

#data:  df_merged$dataset and df_merged$Sex
#X-squared = 653.57, df = 1, p-value < 2.2e-16
summary(as.factor(excluded_sample$Sex))
#0      1 
#164694 129856 
summary(as.factor(included_sample$Sex))
#0      1 
#108700  99271 