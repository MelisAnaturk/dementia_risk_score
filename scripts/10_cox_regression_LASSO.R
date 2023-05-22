# Created by Melis Anaturk (Feb 2020)

# This r script calculates runs LASSO regression

# STEPS IN SCRIPT
# 1. Data set up
# 2. LASSO regression - feature selection
# 3. cox regression - determine beta-weights

#LASSO/Ridge regression
# http://www.sthda.com/english/articles/36-classification-methods-essentials/149-penalized-logistic-regression-essentials-in-r-ridge-lasso-and-elastic-net/
# https://james-cole.github.io/UKBiobank-Brain-Age/
# https://eight2late.wordpress.com/2017/07/11/a-gentle-introduction-to-logistic-regression-and-lasso-regularisation-using-r/
# https://drsimonj.svbtle.com/ridge-regression-with-glmnet
# https://stats.stackexchange.com/questions/244729/lasso-with-interaction-terms-is-it-okay-if-main-effects-are-shrunk-to-zero
# https://www.surveypractice.org/article/2716-using-lasso-to-model-interactions-and-nonlinearities-in-survey-data 


#------ 1. DATA SET UP -----------------------------------------------------------------------------------
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

# load .rda file
load(file = paste0(data_pathway,"ukbdata_diagnoses_baseline_diseasestatus_baselinemedications_ANUADRI_CAIDE_FRS_recoded_DRS_derivedhythchol_fiftyplusnoapoe.rda"))

# Recode variables
df$current_history_depression <- ifelse(df$depression_BIN_FINAL_0_0==1|df$Antidepressant_meds_0_0==1, 1,0)
summary(as.factor(df$current_history_depression))

df$household_occupancy<-as.factor(ifelse(df$Number_in_household_0_0==2,0,
                                         ifelse(df$Number_in_household_0_0==1,1,
                                                ifelse(df$Number_in_household_0_0>2,2,NA))))

#CATEGORICAL & CONTINOUS PREDICTORS 
myvars <- c("Age_when_attended_assesment_centre_0_0","education_years", "Townsend_deprivation_Groups_0_0", "BMI_0_0",
            "Sex", "Sleeplesness_insomnia_0_0_bin", "family_history_of_dementia", 
            "Diabetes_BIN_FINAL_0_0", "LDL_0_0","HDL_cholesterol_0_0", "Total_cholesterol_0_0",
            "current_history_depression","TBI_BIN_FINAL_0_0", "stroke_TIA_BIN_FINAL", "Smoker_bin", "units_combined",
            "Systolic_BP_auto_mean", "IPAQ_activity_group_0_0", "Hearing_prob", "Sleep_duration_0_0", 
            "total_fish_intake_per_week_0_0", "Social_engagement_0_2", "Atrial_Fibrillation_BIN_FINAL_0_0",
            "household_occupancy","dementia_BIN_TOTAL", "APOE_genotype_bin", "NSAIDs_0_0", "HRT_0_0",
            "hypertensive","cholesterol")


# number of candidate risk factors
length(myvars)-1 # to exclude dementia variable
# 29 variables

# examine missing data
pMiss <- function(x){sum(is.na(x))/length(x)*100}
apply(df[myvars],2,pMiss)

# examine class of each input feature
lapply(df[myvars], class) 
lapply(df[myvars], summary) 

# a few variables that should be categorical are encoded as continuous - convert these
numeric_vars <- c("Townsend_deprivation_Groups_0_0", "family_history_of_dementia", "current_history_depression",
                  "stroke_TIA_BIN_FINAL")

df[,numeric_vars] <- lapply(df[,numeric_vars], as.factor)

# how many individuals with dementia remaining?
summary(df$dementia_BIN_TOTAL)
#   0      1 
#220807   3955

# look at correlations between predictors
continuous_vars <- c("Age_when_attended_assesment_centre_0_0","education_years","BMI_0_0",
                     "LDL_0_0","HDL_cholesterol_0_0", "Total_cholesterol_0_0",
                     "units_combined",
                     "Systolic_BP_auto_mean", "Sleep_duration_0_0", 
                     "total_fish_intake_per_week_0_0")

df_continuous <- df[,continuous_vars]
corr_continuous<-cor(df_continuous)
mypalette = colorRampPalette(c("#4477AA", "#77AADD", "#FFFFFF", "#EE9988", "#BB4444"))
continuous_vars_names <- c("Age","Education","BMI",
                           "LDL","HDL", "Total Cholesterol",
                           "Alcohol", "SBP", "Sleep", "Fish")
rownames(corr_continuous)<-continuous_vars_names
colnames(corr_continuous)<-continuous_vars_names

corrplot(corr_continuous, method="color", addCoef.col = "black", order = "hclust",
         col=mypalette(200), type="lower", cl.cex= 1.2, number.cex = 0.75)

jpeg(paste0(save_pathway,"continuousvars_corrplot.jpg"),
     quality = 100, res=300, width = 6, height = 6, units="in")
corrplot(corr_continuous, method="color", addCoef.col = "black", order = "hclust",
         col=mypalette(200), type="lower", cl.cex= 1.2, number.cex = 0.75)
dev.off()
rm(df_continuous, corr_continuous,mypalette,continuous_vars_names,continuous_vars)
# we're now excluding total cholesterol from the set of predictors fed into LASSO due to high correlations

# 2.1 Remove outliers from dataframe before train/test split

#' Detect outliers using IQR method
#' 
#' @param x A numeric vector
#' @param na.rm Whether to exclude NAs when computing quantiles
#' 
is_outlier <- function(x, na.rm = FALSE) {
  qs = quantile(x, probs = c(0.25, 0.75), na.rm = na.rm)
  
  lowerq <- qs[1]
  upperq <- qs[2]
  iqr = upperq - lowerq 
  
  extreme.threshold.upper = (iqr * 3) + upperq
  extreme.threshold.lower = lowerq - (iqr * 3)
  
  # Return logical vector
  x > extreme.threshold.upper | x < extreme.threshold.lower
}

#' Remove rows with outliers in given columns
#' 
#' Any row with at least 1 outlier will be removed
#' 
#' @param df A data.frame
#' @param cols Names of the columns of interest. Defaults to all columns.
#' 
#' 
remove_outliers <- function(df, cols = names(df)) {
  for (col in cols) {
    cat("Removing outliers in column: ", col, " \n")
    df <- df[!is_outlier(df[[col]]),]
  }
  df
}

# create vector of continuous variables (note total cholesterol has been removed)
cont_vars <- c("Age_when_attended_assesment_centre_0_0", "education_years", "LDL_0_0",  "HDL_cholesterol_0_0", "Systolic_BP_auto_mean", "Sleep_duration_0_0", "BMI_0_0", "total_fish_intake_per_week_0_0", "units_combined")

# apply remove_outliers function to df
df_filtered <- remove_outliers(df, cont_vars)
rm(df)

#n dementia?
summary(df_filtered$dementia_BIN_TOTAL)
#216949   3813  

#we need time to death/censor/dementia
#update death records
df_death<-read.csv("../../raw_data/ukb50321_deathvars_refresh.csv", stringsAsFactors = FALSE)
names(df_death) #date is X40000.0.0
df_death <- df_death[,c("eid","X40000.0.0")]
names(df_death) <- c("eid","death_date_refresh")
df_death$eid <- as.character(df_death$eid)
df_death$death_date_refresh <- as.Date(df_death$death_date_refresh, format="%Y-%m-%d")
summary(df_death$death_date_refresh)
# Min.      1st Qu.       Median         Mean      3rd Qu.         Max.         NA's 
# "2006-05-10" "2014-04-18" "2017-06-22" "2016-11-21" "2019-12-21" "2021-11-12"     "464516" 
df_death$death_record_refresh <- ifelse( !is.na(df_death$death_date_refresh) &
                                           df_death$death_date_refresh<as.Date("2021-10-31",format="%Y-%m-%d"),1,0)
summary(as.factor(df_death$death_record_refresh))
# 0      1 
# 464532  37881 

#merge
df_filtered <- list(df_filtered, df_death) %>% reduce(left_join, by="eid")
rm(df_death)






df_filtered$Date_of_assessment_0_0_new.x <- as.Date(df_filtered$Date_of_assessment_0_0_new.x, format="%Y-%m-%d")
df_filtered$has_death_record <- ifelse(df_filtered$death_record_refresh==1,1,0)

#dementia
df_dementia <- df_filtered[which(df_filtered$dementia_BIN_TOTAL==1),]
df_dementia$time_at_risk <- df_dementia$earliest_dementia_date - df_dementia$Date_of_assessment_0_0_new.x
summary(as.numeric(df_dementia$time_at_risk))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#373    2503    3426    3239    4106    5196 

#deaths
df_deaths <- df_filtered[which( (df_filtered$dementia_BIN_TOTAL==0) & (df_filtered$has_death_record==1) ),]
df_deaths$time_at_risk <- df_deaths$death_date_refresh - df_deaths$Date_of_assessment_0_0_new.x
summary(as.numeric(df_deaths$time_at_risk))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 4    1859    2989    2836    3912    5254 

#healthy
df_healthy <- df_filtered[which( (df_filtered$dementia_BIN_TOTAL==0) & (df_filtered$has_death_record==0) ),]
df_healthy$time_at_risk <- as.Date("2021-10-31",format="%Y-%m-%d") - df_healthy$Date_of_assessment_0_0_new.x
summary(as.numeric(df_healthy$time_at_risk))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#4048    4386    4626    4630    4850    5430  

df_filtered <- rbind(df_dementia, df_healthy)
df_filtered<-rbind(df_filtered, df_deaths)
df_filtered$time_at_risk <- as.numeric(df_filtered$time_at_risk)
summary(df_filtered$time_at_risk)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#4    4346    4604    4494    4836    5430 

rm(df_dementia, df_deaths, df_healthy)

#save
save(df_filtered, file="../../raw_data/df_filtered_prelasso.rda")

#----- 2. LASSO REGRESSION ------------------------------------------
# LASSO cox regression is first run to identify the subset of predictors to be used in our risk score
# https://glmnet.stanford.edu/articles/Coxnet.html 

#load if needed
#load("../../raw_data/df_filtered_prelasso.rda")

#define variables of interest
myvars <- c("time_at_risk","Age_when_attended_assesment_centre_0_0","education_years", "Townsend_deprivation_Groups_0_0", "BMI_0_0",
            "Sex", "Sleeplesness_insomnia_0_0_bin", "family_history_of_dementia", 
            "Diabetes_BIN_FINAL_0_0", "LDL_0_0","HDL_cholesterol_0_0",
            "current_history_depression","TBI_BIN_FINAL_0_0", "stroke_TIA_BIN_FINAL", "Smoker_bin", "units_combined",
            "Systolic_BP_auto_mean", "IPAQ_activity_group_0_0", "Hearing_prob", "Sleep_duration_0_0", 
            "total_fish_intake_per_week_0_0", "Social_engagement_0_2", "Atrial_Fibrillation_BIN_FINAL_0_0",
            "household_occupancy","dementia_BIN_TOTAL", "NSAIDs_0_0", "HRT_0_0",
            "hypertensive","cholesterol")



# split df into train and test sets (we're using 80/20 split here)
set.seed(2003)
training.samples <- df_filtered$dementia_BIN_TOTAL %>% createDataPartition(p = 0.8, list = FALSE)

# conduct train/test split
train.data  <- df_filtered[training.samples, ][myvars]
test.data <- df_filtered[-training.samples, ][myvars]

save(train.data, file = "../../raw_data/modelvar/train_data_outliers_removed_prelasso.rda")
save(test.data, file = "../../raw_data/modelvar/test_data_outliers_removed_prelasso.rda")


# OR ALTERNATIVELY YOU CAN load train/test data
#load(file = "../../raw_data/train_data_outliers_removed_prelasso.rda")
#load(file = "../../raw_data/test_data_outliers_removed_prelasso.rda")
#train.data <- train.data[myvars]
#test.data <- test.data[myvars]

# Create vector of only continous variables
cont_vars <- c("Age_when_attended_assesment_centre_0_0",  "education_years", "LDL_0_0","HDL_cholesterol_0_0","Systolic_BP_auto_mean","Sleep_duration_0_0", "BMI_0_0", "total_fish_intake_per_week_0_0", "units_combined")

# scale the data
scaled.train.data <- scale(train.data[, cont_vars], scale = TRUE, center = TRUE)
scaling.parameters.center <- attr(scaled.train.data, "scaled:center")
scaling.parameters.scale <- attr(scaled.train.data, "scaled:scale")
scaled.train.data <- as.data.frame(scaled.train.data)
scaled.test_data <- as.data.frame(scale(test.data[, cont_vars], scaling.parameters.center, scaling.parameters.scale))
var_to_drop <- cont_vars
train.data <- train.data[, !(names(train.data) %in% var_to_drop)]
test.data <- test.data[, !(names(test.data) %in% var_to_drop)]
train.data <- cbind(train.data, scaled.train.data)
test.data <- cbind(test.data, scaled.test_data)
names(train.data)

# specify predictors, excluding intercept (-1) and time at risk (-2)
x <- model.matrix(dementia_BIN_TOTAL~., train.data)[,-seq(1,2)]
#check
names(as.data.frame(x))

# Convert the outcome (dementia) to survival obj
y <- Surv(train.data$time_at_risk, train.data$dementia_BIN_TOTAL, type="right")
#set type to right cens
attr(y,"type") <- "right"

# use cross-validation to identify optimal lambda
set.seed(2003)
lasso.fit.cv <-  cv.glmnet(x, y, family = "cox", alpha = 1, lambda = NULL, nfolds = 10)
plot(lasso.fit.cv)

#save fit
save(lasso.fit.cv, file = paste0(save_pathway, "lasso_fit_cox_cv.rda"))

# Or otherwise load model into R
load(file = paste0(save_pathway, "lasso_fit_cox_cv.rda"))

# Fit the final model on the training data    
# NOTE: dont use lambda.min - use lambda.1se
set.seed(2003)
lasso.final.1 <- glmnet(x, y, alpha = 1, family = "cox", lambda = lasso.fit.cv$lambda.1se) #lambda.1se
coef(lasso.final.1, s = lasso.fit.cv$lambda.1se) #lambda.1se

#Output 
# Townsend_deprivation_Groups_0_04        0.123164107
# Sex1                                    0.025075906
# family_history_of_dementia1             0.222155411
# Diabetes_BIN_FINAL_0_01                 0.530047863
# current_history_depression1             0.374105120
# stroke_TIA_BIN_FINAL1                   0.610589869
# household_occupancy1                    0.009732028
# hypertensive1                           0.105376272
# cholesterol1                            0.094987826
# Age_when_attended_assesment_centre_0_0  0.880455593
# education_years                        -0.063917765


# tidy output
View(tidy(lasso.final.1))

# save the list of selected predictors (to carry forward to cox regression analyses)
write.csv(tidy(lasso.final.1), paste0(save_pathway,"cox_LASSO_results.csv"))

# save the model
save(lasso.final.1, file = paste0(model_pathway, "lasso_fit_cox_final.rda"))


#put original data back into train/test
#coxph will want 1,2 as outcome
df_filtered$dementia_BIN_surv <- as.numeric(df_filtered$dementia_BIN_TOTAL)
summary(as.factor(df_filtered$dementia_BIN_surv))
#1      2 
#216949   3813 

train.data  <- df_filtered[training.samples, ]
test.data <- df_filtered[-training.samples, ]

save(train.data, file = paste0("../../raw_data/modelvar/train_data_outliers_removed_postlasso.rda"))
save(test.data, file = paste0("../../raw_data/modelvar/test_data_outliers_removed_postlasso.rda"))

