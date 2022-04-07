# Created by Raihaan Patel
# Modified by Melis Anaturk (Feb 2020)

# This r script calculates the BDSI score

# STEPS IN SCRIPT
# 1. Data set up
# 2. Data recoding
# 3. BDSI caculation

#------ 1. DATA SET UP -----
# 1.1 load packages
#library(car)
#library(sjstats)
#library(pROC)
#library(ggplot2)
#library(janitor)
#library(dplyr)

# 1.2 data pathway
data_pathway = "../../raw_data/"

# 1.3 read in csv files
#df <- read.csv('ukb_demo-phys-cog-dementia-depression-tbi-stroke-edu_clean.csv', header=TRUE, sep=",", stringsAsFactors = FALSE)
#cols <- read.csv('column_names.csv')

#additional variables
#df_add <- read.csv('../subsets_new/add_vars2.csv', header=TRUE, sep=",", stringsAsFactors = FALSE)
#cols_add <- read.csv('../subsets_new/add_vars_head.csv', header=TRUE, sep=",", stringsAsFactors = FALSE)

#add meaningful headings
#original dataset
#newcols <- as.vector(cols$UID_text)
#names(df) <- newcols
#names(df)

#additional dataset
#newcols <- as.vector(cols_add$col_name)
#names(df_add) <- newcols
#names(df_add)

# make NAs real NAs
#df[df == "NA"] <- NA
#df_add[df_add == "NA"] <- NA

#remove empty cols
#df <- remove_empty(df, which = c("cols"))
#493 vars ---> 487

#merge datasets
#df <- merge(df, df_add, by.c="eid", all.x=TRUE)

#------ 2. Data recoding ----
# 2.1 Age
#subset df into ppts who are between 65-79
df_bdsi = data.frame(subset(df,
                       df$Age_when_attended_assesment_centre_0_0>=65 & 
                       df$Age_when_attended_assesment_centre_0_0<= 79))

#recode age
df_bdsi$bdsi_age_score<-ifelse(df_bdsi$Age_when_attended_assesment_centre_0_0==65, 0, 
                        ifelse(df_bdsi$Age_when_attended_assesment_centre_0_0>65, 
                        round(df_bdsi$Age_when_attended_assesment_centre_0_0 - 65,0 ), NA))

# 2.2 Education
df_bdsi$bdsi_edu_score<-ifelse(df_bdsi$education_years<12, 9, ifelse(df_bdsi$education_years>=12, 0, NA))
summary(as.factor(df_bdsi$bdsi_edu_score))

# 2.3 BMI
df_bdsi$bdsi_bmi_score<-ifelse(df_bdsi$BMI_0_0<18.5, 8, ifelse(df_bdsi$BMI_0_0>=18.5, 0, NA))
summary(as.factor(df_bdsi$bdsi_bmi_score))

# 2.4 DIABETES
df_bdsi$bdsi_diab_score<-ifelse(df_bdsi$Diabetes_BIN_TOTAL==1, 3, ifelse(df_bdsi$Diabetes_BIN_TOTAL==0, 0, NA))
summary(as.factor(df_bdsi$bdsi_diab_score))

# 2.5 STROKE
df_bdsi$bdsi_stroke_score <- ifelse(df_bdsi$stroke_BIN_TOTAL==1, 6, df_bdsi$stroke_BIN_TOTAL)
summary(as.factor(df_bdsi$bdsi_stroke_score))

# 2.6 DEPRESSION
#df_bdsi$bdsi_depr_score<-ifelse(df_bdsi$depression_BIN==1, 6, 0)
#summary(as.factor(df_bdsi$bdsi_depr_score))

df_bdsi$bdsi_depr_score<-ifelse(df_bdsi$depression_BIN_TOTAL==1, 6, 0)
summary(as.factor(df_bdsi$bdsi_depr_score))

#------ 3. Calculate BDSI score ------
#point-based score
df_bdsi$point_TOTAL_BDSI_SCORE <- (df_bdsi$bdsi_age_score + df_bdsi$bdsi_edu_score + df_bdsi$bdsi_bmi_score + df_bdsi$bdsi_diab_score + 
                                   df_bdsi$bdsi_stroke_score + df_bdsi$bdsi_depr_score)

# excluded domain: df_bdsi$point_cog_score_anu_adri_recoded 
summary(as.factor(df_bdsi$point_TOTAL_BDSI_SCORE))
hist(df_bdsi$point_TOTAL_BDSI_SCORE)

#combine with main dataset
myvars <- c("eid", "point_TOTAL_BDSI_SCORE")
df_bdsi <- df_bdsi[myvars]

df <- merge(df, df_bdsi, by.c="eid", all.x=TRUE)
View(df[myvars])
rm(df_bdsi)
