# Created by Melis Anaturk (Feb 2020)

# This r script calculates runs LASSO regression

# STEPS IN SCRIPT
# 1. Data set up
# 2. LASSO regression - feature selection
# 3. Logistic regression - determine beta-weights

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

# specify data, model and results pathway
data_pathway = "../../raw_data/"
model_pathway = "../models/"
save_pathway = "../results/"

# load .rda file
load(file = paste0(data_pathway,"ukbdata_diagnoses_baseline_diseasestatus_baselinemedications_ANUADRI_CAIDE_FRS_recoded_DRS.rda"))

# Recode variables
df$current_history_depression <- ifelse(df$depression_BIN_FINAL_0_0==1|df$Antidepressant_meds_0_0==1, 1,0)
summary(as.factor(df$current_history_depression))

#CATEGORICAL & CONTINOUS PREDICTORS 
myvars <- c("Age_when_attended_assesment_centre_0_0","education_years", "Townsend_deprivation_Groups_0_0", "BMI_0_0",
            "Sex", "Sleeplesness_insomnia_0_0_bin", "family_history_of_dementia", 
            "Diabetes_BIN_FINAL_0_0", "LDL_0_0","HDL_cholesterol_0_0", "Total_cholesterol_0_0",
            "current_history_depression","TBI_BIN_FINAL_0_0", "stroke_TIA_BIN_FINAL", "Smoker_bin", "units_combined",
            "Systolic_BP_auto_mean", "IPAQ_activity_group_0_0", "Hearing_prob", "Sleep_duration_0_0", "Antihypertensive_meds_0_0",
            "total_fish_intake_per_week_0_0", "Social_engagement_0_2", "Atrial_Fibrillation_BIN_FINAL_0_0",
            "Number_in_household_0_0","dementia_BIN_TOTAL", "APOE_genotype_bin", "NSAIDs_0_0", "HRT_0_0", "statins_0_0", "Aspirin_0_0")

# number of candidate risk factors
length(myvars)-1 # to exclude dementia variable
# 30 variables

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
#211040   1325

# look at correlations between predictors


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
cont_vars <- c("Age_when_attended_assesment_centre_0_0", "education_years", "LDL_0_0",  "HDL_cholesterol_0_0", "Systolic_BP_auto_mean", "Sleep_duration_0_0", "BMI_0_0", "total_fish_intake_per_week_0_0", "Number_in_household_0_0", "units_combined")

# apply remove_outliers function to df
df_filtered <- remove_outliers(df, cont_vars)

# exclude people who only have one assessment date

#----- 2. LASSO REGRESSION ------------------------------------------
# LASSO logistic regression is first run to identify the subset of predictors to be used in our risk score
# IF you would like to run a cox proportional hazard equivalent, check out:
# https://glmnet.stanford.edu/articles/Coxnet.html 
# or the r package coxphMIC

# split df into train and test sets (we're using 80/20 split here)
set.seed(2003)
training.samples <- df_filtered$dementia_BIN_TOTAL %>% createDataPartition(p = 0.8, list = FALSE)

# conduct train/test split
train.data  <- df_filtered[training.samples, ][myvars]
test.data <- df_filtered[-training.samples, ][myvars]

# OR ALTERNATIVELY YOU CAN load train/test data
#load(file = paste0(data_pathway,"train_data_outliers_removed.rda"))
#load(file = paste0(data_pathway,"test_data_outliers_removed.rda"))
#train.data <- train.data[myvars]
#test.data <- test.data[myvars]

# Create vector of only continous variables
cont_vars <- c("Age_when_attended_assesment_centre_0_0",  "education_years", "LDL_0_0",  "HDL_cholesterol_0_0", "Systolic_BP_auto_mean", "Sleep_duration_0_0", "BMI_0_0", "total_fish_intake_per_week_0_0", "Number_in_household_0_0", "units_combined")

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

# specify predictors, excluding intercept
x <- model.matrix(dementia_BIN_TOTAL~., train.data)[,-1]

# Convert the outcome (class) to a numerical variable
y <- ifelse(train.data$dementia_BIN_TOTAL==1, 1, 0)

new.y <- ifelse(test.data$dementia_BIN_TOTAL==1, 1, 0)
new.x <- model.matrix(dementia_BIN_TOTAL~., test.data)[,-1]

# use cross-validation to identify optimal lambda
set.seed(2003)
lasso.fit.cv <-  cv.glmnet(x, y, family = "binomial", alpha = 1, lambda = NULL, nfolds = 10)
plot(lasso.fit.cv)

#save fit
save(lasso.fit.cv, file = paste0(save_pathway, "lasso_fit_cv.rda"))

# Or otherwise load model into R
load(file = paste0(save_pathway, "lasso_fit_cv.rda"))

# Fit the final model on the training data    
# NOTE: dont use lambda.min - use lambda.1se
set.seed(2003)
lasso.final.1 <- glmnet(x, y, alpha = 1, family = "binomial", lambda = lasso.fit.cv$lambda.1se) #lambda.1se
coef(lasso.final.1, s = lasso.fit.cv$lambda.1se) #lambda.1se

#Output 
# 34 x 1 sparse Matrix of class "dgCMatrix"
# 1
# (Intercept)                            -5.6537125
# Townsend_deprivation_Groups_0_01        .        
# Townsend_deprivation_Groups_0_02        .        
# Townsend_deprivation_Groups_0_03        .        
# Townsend_deprivation_Groups_0_04        .        
# Sex1                                    .        
# Sleeplesness_insomnia_0_0_bin1          .        
# Sleeplesness_insomnia_0_0_bin2          .        
# family_history_of_dementia1             .        
# Diabetes_BIN_FINAL_0_01                 0.1138954
# current_history_depression1             0.2844328
# TBI_BIN_FINAL_0_01                      .        
# stroke_TIA_BIN_FINAL1                   0.6099854
# Smoker_bin1                             .        
# IPAQ_activity_group_0_01                .        
# Hearing_prob1                           .        
# Antihypertensive_meds_0_01              .        
# Social_engagement_0_21                  .        
# Atrial_Fibrillation_BIN_FINAL_0_01      .        
# APOE_genotype_bin1                      0.5757621
# NSAIDs_0_01                             .        
# HRT_0_01                                .        
# statins_0_01                            .        
# Aspirin_0_01                            .        
# Age_when_attended_assesment_centre_0_0  0.8049626
# education_years                         .        
# LDL_0_0                                 .        
# HDL_cholesterol_0_0                     .        
# Systolic_BP_auto_mean                   .        
# Sleep_duration_0_0                      .        
# BMI_0_0                                 .        
# total_fish_intake_per_week_0_0          .        
# Number_in_household_0_0                 .        
# units_combined                          .  


# tidy output
View(tidy(lasso.final.1))

# save the list of selected predictors (to carry forward to logistic regression analyses)
#write.csv(tidy(lasso.final.1), paste0(model_pathway,"logistic_LASSO_results.csv"))
#ma was writing out to the model folder, think it should be in results folder
write.csv(tidy(lasso.final.1), paste0(save_pathway,"logistic_LASSO_results.csv"))

# save the model
save(lasso.final.1, file = paste0(model_pathway, "lasso_fit_final.rda"))

# Or otherwise load model into R
#load(file = paste0(model_pathway, "lasso_fit_final.rda""))

#----- 2. LOGISTIC REGRESSION ------------------------------------------
# Simple logistic regression is now run to calculate the beta-weights for each of the 
#components in our risk score (note, we also add in age, sex and education)

# conduct train/test split (same partion as above conducted but more columns are retained)
train.data  <- df_filtered[training.samples, ]
test.data <- df_filtered[-training.samples, ]

# calculate years of follow-up
train.data$years_diff_baseline_all <- as.numeric(train.data$date_diff_baseline_all/365.242)
test.data$years_diff_baseline_all <- as.numeric(test.data$date_diff_baseline_all/365.242)

# describe
describe(as.numeric(train.data$years_diff_all_time))
describe(as.numeric(test.data$years_diff_all_time))

# specify age only and various versions of UKB-DRS (see manuscript for details)
age_only <-      paste("dementia_BIN_TOTAL~Age_when_attended_assesment_centre_0_0")

UKBDRS_LASSO  <-            paste("dementia_BIN_TOTAL ~  Age_when_attended_assesment_centre_0_0 +
                                  Diabetes_BIN_FINAL_0_0  +  current_history_depression + stroke_TIA_BIN_FINAL")

UKBDRS_LASSO_MAN  <-  paste("dementia_BIN_TOTAL ~  Age_when_attended_assesment_centre_0_0 +  Sex + education_years +
                            Diabetes_BIN_FINAL_0_0  +  current_history_depression + stroke_TIA_BIN_FINAL + 
                            family_history_of_dementia")

UKBDRS_APOE_LASSO <-      paste("dementia_BIN_TOTAL ~  Age_when_attended_assesment_centre_0_0 +
                                Diabetes_BIN_FINAL_0_0  +  current_history_depression + stroke_TIA_BIN_FINAL + APOE_genotype_bin")

UKBDRS_APOE_LASSO_MAN <-   paste("dementia_BIN_TOTAL ~  Age_when_attended_assesment_centre_0_0 +  Sex + education_years +
                                 Diabetes_BIN_FINAL_0_0  +  current_history_depression + stroke_TIA_BIN_FINAL + family_history_of_dementia + APOE_genotype_bin")


#models <- c(age_only, UKBDRS_APOE_LASSO_models,UKBDRS_LASSO_models, UKBDRS_APOE_LASSO_MAN_models, UKBDRS_APOE_LASSO_MAN_models)


# We've previously computed the linear predictor and predicted probabilties for CAIDE, FRS, ANU-ADRI and DRS
# Here we derive the linear predictor and predicted probabilites for age-only, UKB-DRS+APOE, UKB-DRS 
# using a for loop

models <- c("age_only", "UKBDRS_LASSO", "UKBDRS_LASSO_MAN", "UKBDRS_APOE_LASSO", "UKBDRS_APOE_LASSO_MAN")

for (m in models){
  print(paste0('applying logistic regression model for ', m))
  model <- glm(as.formula(m), data=train.data, family="binomial")
 
  print(paste0('UKB training set : calculating linear predictor and predicted probabilities for ', m))
  train.data[paste(m, "linear_predictor", sep="_")] <- predict(model, train.data)
  train.data[paste(m, "predicted_prob", sep="_")] <-   1/(1+exp(-train.data[paste(m, "linear_predictor", sep="_")])) # can also be computed with predict(model, type='response')
  
  print(paste0('UKB test set : calculating linear predictor and predicted probabilities for ', m))
  test.data[paste(m, "linear_predictor", sep="_")] <- predict(model, test.data)
  test.data[paste(m, "predicted_prob", sep="_")] <-   1/(1+exp(-test.data[paste(m, "linear_predictor", sep="_")])) # can also be computed with predict(model, newdata= test.data, type='response') 
}

#save train and test data
save(train.data, file = paste0(data_pathway, "train_data_outliers_removed.rda"))
save(test.data, file = paste0(data_pathway, "test_data_outliers_removed.rda"))
# Above loop rewritten as function
# model_discrimination <- function(m){
#   print(paste0('applying logistic regression model for ', m))
#   model <- glm(as.formula(m), data=train.data, family="binomial")
#   
#   print(paste0('UKB training set : calculating linear predictor and predicted probabilities for ', m))
#   train.data[paste(m, "linear_predictor", sep="_")] <- predict(model, train.data)
#   train.data[paste(m, "predicted_prob", sep="_")] <-   1/(1+exp(-train.data[paste(m, "linear_predictor", sep="_")])) # can also be computed with predict(model, type='response')
#   
#   print(paste0('UKB test set : calculating linear predictor and predicted probabilities for ', m))
#   test.data[paste(m, "linear_predictor", sep="_")] <- predict(model, test.data)
#   test.data[paste(m, "predicted_prob", sep="_")] <-   1/(1+exp(-test.data[paste(m, "linear_predictor", sep="_")])) # can also be computed with predict(model, newdata= test.data, type='response') 
# return(list(train.data, test.data))
#        }

# apply function to calculate linear predictor and predicted probabilities (still working on this..)
#final_df <- lapply(models, model_discrimination)
