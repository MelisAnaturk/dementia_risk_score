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
library(corrplot)

# specify data, model and results pathway
data_pathway = "../../raw_data/"
model_pathway = "../models/"
save_pathway = "../results/"

# load .rda file
load(file = paste0(data_pathway,"ukbdata_diagnoses_baseline_diseasestatus_baselinemedications_ANUADRI_CAIDE_FRS_recoded_DRS_fiftyplusnoapoe.rda"))

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
            "Systolic_BP_auto_mean", "IPAQ_activity_group_0_0", "Hearing_prob", "Sleep_duration_0_0", "Antihypertensive_meds_0_0",
            "total_fish_intake_per_week_0_0", "Social_engagement_0_2", "Atrial_Fibrillation_BIN_FINAL_0_0",
            "household_occupancy","dementia_BIN_TOTAL", "APOE_genotype_bin", "NSAIDs_0_0", "HRT_0_0", "statins_0_0", "Aspirin_0_0")

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
#220807   3955

# look at correlations between predictors
continuous_vars <- c("Age_when_attended_assesment_centre_0_0","education_years","BMI_0_0",
                     "LDL_0_0","HDL_cholesterol_0_0", "Total_cholesterol_0_0",
                     "units_combined",
                     "Systolic_BP_auto_mean", "Sleep_duration_0_0", 
                     "total_fish_intake_per_week_0_0", 
                     "Number_in_household_0_0")

df_continuous <- df[,continuous_vars]
corr_continuous<-cor(df_continuous)
mypalette = colorRampPalette(c("#4477AA", "#77AADD", "#FFFFFF", "#EE9988", "#BB4444"))
continuous_vars_names <- c("Age","Education","BMI",
                           "LDL","HDL", "Total Cholesterol",
                           "Alcohol", "SBP", "Sleep", "Fish", "N. Household")
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

#n dementia?
summary(df_filtered$dementia_BIN_TOTAL)
#216949   3813  

# exclude people who only have one assessment date
#redefine myvars to not include total cholesterol
myvars <- c("Age_when_attended_assesment_centre_0_0","education_years", "Townsend_deprivation_Groups_0_0", "BMI_0_0",
            "Sex", "Sleeplesness_insomnia_0_0_bin", "family_history_of_dementia", 
            "Diabetes_BIN_FINAL_0_0", "LDL_0_0","HDL_cholesterol_0_0",
            "current_history_depression","TBI_BIN_FINAL_0_0", "stroke_TIA_BIN_FINAL", "Smoker_bin", "units_combined",
            "Systolic_BP_auto_mean", "IPAQ_activity_group_0_0", "Hearing_prob", "Sleep_duration_0_0", "Antihypertensive_meds_0_0",
            "total_fish_intake_per_week_0_0", "Social_engagement_0_2", "Atrial_Fibrillation_BIN_FINAL_0_0",
            "household_occupancy","dementia_BIN_TOTAL", "NSAIDs_0_0", "HRT_0_0", "statins_0_0", "Aspirin_0_0")
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
cont_vars <- c("Age_when_attended_assesment_centre_0_0",  "education_years", "LDL_0_0",  "HDL_cholesterol_0_0", "Systolic_BP_auto_mean", "Sleep_duration_0_0", "BMI_0_0", "total_fish_intake_per_week_0_0", "units_combined")

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
#(Intercept)                            -4.58422068
#Townsend_deprivation_Groups_0_04        0.06618836
#family_history_of_dementia1             0.19324274
#Diabetes_BIN_FINAL_0_01                 0.44319669
#current_history_depression1             0.33999776
#stroke_TIA_BIN_FINAL1                   0.52360243
#Antihypertensive_meds_0_01              0.05892152
#statins_0_01                            0.06841070
#Aspirin_0_01                            0.09320814
#Age_when_attended_assesment_centre_0_0  0.84523138
#education_years                        -0.06562031  

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
#vars      n mean  sd median trimmed mad min   max range skew kurtosis   se
#X1    1 176611 1.13 2.78      0    0.32   0   0 12.27 12.27 2.33     4.01 0.01
describe(as.numeric(test.data$years_diff_all_time))
#vars     n mean   sd median trimmed mad min   max range skew kurtosis   se
#X1    1 44151 1.15 2.79      0    0.34   0   0 12.29 12.29 2.29     3.82 0.01

#membership in highest deprived group (4, on scale of 0-4) was selected, not others
#create binary version of deprivation to use in model
train.data$Townsend_deprivation_modelvar<-ifelse(train.data$Townsend_deprivation_Groups_0_0==4,1,0)
train.data$Townsend_deprivation_modelvar<-as.factor(train.data$Townsend_deprivation_modelvar)
test.data$Townsend_deprivation_modelvar<-ifelse(test.data$Townsend_deprivation_Groups_0_0==4,1,0)
test.data$Townsend_deprivation_modelvar<-as.factor(test.data$Townsend_deprivation_modelvar)


#### 2.1 Test beta coefficients ####
#based on the lasso selected vars, compute model coefficients in train data
#use 2 models - one with apoe, one without. include sex as a manually selected variable as well

# specify various versions of UKB-DRS (see manuscript for details)
UKBDRS_LASSO  <-  paste("dementia_BIN_TOTAL ~  Age_when_attended_assesment_centre_0_0 +  family_history_of_dementia +
                            education_years + Townsend_deprivation_modelvar +Diabetes_BIN_FINAL_0_0  +
                            current_history_depression + stroke_TIA_BIN_FINAL +  
                            Antihypertensive_meds_0_0 +  
                            statins_0_0 + Aspirin_0_0")

UKBDRS_APOE_LASSO <-    paste("dementia_BIN_TOTAL ~  Age_when_attended_assesment_centre_0_0 +  family_history_of_dementia +
                            education_years + Townsend_deprivation_modelvar +Diabetes_BIN_FINAL_0_0  +
                            current_history_depression + stroke_TIA_BIN_FINAL +  
                            Antihypertensive_meds_0_0 +  
                            statins_0_0 + Aspirin_0_0 + APOE_genotype_bin")

#check model coefficeints when applied to training data
format_modelcoefs <-function(lr_out){
  betas<-coef(lr_out)
  odds<-exp(coef(lr_out))
  ci<-confint(lr_out)
  odds_ci<-exp(confint(lr_out))
  df_model<-data.frame(cbind(names(betas),
                             as.vector(round(betas,5)),
                             as.vector(round(ci,3)[,1]),
                             as.vector(round(ci,3)[,2]),
                             paste(as.vector(round(odds,2)), " [",as.vector(round(odds_ci,2)[,1]), ", ",as.vector(round(odds_ci,2)[,2]), "]",sep=""),
                             as.vector(summary(model)$coefficients[,4])))
  
  names(df_model)<-c("Predictor","beta","lower","upper","OR","p")
  return(df_model)
}

table2_models <- c("UKBDRS_LASSO","UKBDRS_APOE_LASSO")

df_table2<-data.frame(matrix(ncol=6))
names(df_table2)<-c("Predictor","beta","lower","upper","OR","p")

for (m in table2_models){
  print(paste0('applying logistic regression model for ', m))
  model <- glm(as.formula(m), data=train.data, family="binomial")
  print(sprintf("storing results from %s",m))
  df_table2<-rbind(df_table2, format_modelcoefs(model))
}

df_table2$FDR_BH = p.adjust(df_table2$p, method = "BH")
write.csv(df_table2, file=paste0(save_pathway, "lr_betas_initial.csv"))
# both statins and aspirin have non sig betas, so remove them from model
rm(table2_models, df_table2, UKBDRS_LASSO, UKBDRS_APOE_LASSO)


#### 2.2 compute predicted probability/lin predictor of ukb drs models ####
# We've previously computed the linear predictor and predicted probabilties for CAIDE, FRS, ANU-ADRI and DRS
# Here we derive the linear predictor and predicted probabilites for age-only, UKB-DRS+APOE, UKB-DRS 

# specify age only and various versions of UKB-DRS (see manuscript for details)
age_only <-      paste("dementia_BIN_TOTAL~Age_when_attended_assesment_centre_0_0")

UKBDRS_LASSO  <-  paste("dementia_BIN_TOTAL ~  Age_when_attended_assesment_centre_0_0 +  family_history_of_dementia +
                            education_years + Townsend_deprivation_modelvar +Diabetes_BIN_FINAL_0_0  +
                            current_history_depression + stroke_TIA_BIN_FINAL +  
                            Antihypertensive_meds_0_0 +  
                            Aspirin_0_0")

UKBDRS_APOE_LASSO <-    paste("dementia_BIN_TOTAL ~  Age_when_attended_assesment_centre_0_0 +  family_history_of_dementia +
                            education_years + Townsend_deprivation_modelvar +Diabetes_BIN_FINAL_0_0  +
                            current_history_depression + stroke_TIA_BIN_FINAL +  
                            Antihypertensive_meds_0_0 +  
                            Aspirin_0_0 + APOE_genotype_bin")

models <- c("age_only", "UKBDRS_LASSO", "UKBDRS_APOE_LASSO")

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
save(train.data, file = paste0(data_pathway, "train_data_outliers_removed_fiftyplusnoapoe.rda"))
save(test.data, file = paste0(data_pathway, "test_data_outliers_removed_fiftyplusnoapoe.rda"))

#load train and test data, if necessary
load(file = paste0(data_pathway,"train_data_outliers_removed_fiftyplusnoapoe.rda"))
load(file = paste0(data_pathway,"test_data_outliers_removed_fiftyplusnoapoe.rda"))

table2_models <- c("age_only", "UKBDRS_LASSO", "UKBDRS_APOE_LASSO")

df_table2<-data.frame(matrix(ncol=6))
names(df_table2)<-c("Predictor","beta","lower","upper","OR","p")

for (m in table2_models){
  print(paste0('applying logistic regression model for ', m))
  model <- glm(as.formula(m), data=train.data, family="binomial")
  print(sprintf("storing results from %s",m))
  df_table2<-rbind(df_table2, format_modelcoefs(model))
}

df_table2$FDR_BH = p.adjust(df_table2$p, method = "BH")
write.csv(df_table2, file=paste0(save_pathway, "lr_betas_final.csv"))

#checking for correlations among med vars
df$statins_num<-as.numeric(df$statins_0_0)
df$antihypertensive_num<-as.numeric(df$Antihypertensive_meds_0_0)
df$aspirin_num<-as.numeric(df$Aspirin_0_0)
df$nsaid_num<-as.numeric(df$NSAIDs_0_0)
df$hrt_num<-as.numeric(df$HRT_0_0)

corr_res<-tetrachoric(df[,c("statins_num","antihypertensive_num","aspirin_num","nsaid_num","hrt_num")])
#statins_num antihypertensive_num  aspirin_num    nsaid_num     hrt_num
#statins_num           1.00000000           0.72534925  0.709486575 -0.095383210 -0.12891023
#antihypertensive_num  0.72534925           1.00000000  0.622530919 -0.073335583 -0.02766018
#aspirin_num           0.70948658           0.62253092  1.000000000 -0.009758151 -0.07794382
#nsaid_num            -0.09538321          -0.07333558 -0.009758151  1.000000000  0.13481686
#hrt_num              -0.12891023          -0.02766018 -0.077943819  0.134816862  1.00000000

#statins, aspirin, and antihypertensive use are correlated


#### Sex stratify ####
#fit LR in males and females separately to check if betas differ greatly
#load train and test data, if necessary
load(file = paste0(data_pathway,"train_data_outliers_removed_fiftyplusnoapoe.rda"))
load(file = paste0(data_pathway,"test_data_outliers_removed_fiftyplusnoapoe.rda"))

UKBDRS_LASSO_initial  <-  paste("dementia_BIN_TOTAL ~  Age_when_attended_assesment_centre_0_0 +  family_history_of_dementia +
                            education_years + Townsend_deprivation_modelvar +Diabetes_BIN_FINAL_0_0  +
                            current_history_depression + stroke_TIA_BIN_FINAL +  
                            Antihypertensive_meds_0_0 + statins_0_0 + 
                            Aspirin_0_0")

df_lr_sexstratify<-data.frame(matrix(ncol=6))
names(df_lr_sexstratify)<-c("Predictor","beta","lower","upper","OR","p")

#females
model <- glm(as.formula(UKBDRS_LASSO_initial), data=subset(train.data, train.data$Sex==0), family="binomial")
df_lr_sexstratify<-rbind(df_lr_sexstratify, format_modelcoefs(model))
#males
model <- glm(as.formula(UKBDRS_LASSO_initial), data=subset(train.data, train.data$Sex==1), family="binomial")
df_lr_sexstratify<-rbind(df_lr_sexstratify, format_modelcoefs(model))
df_lr_sexstratify$FDR_BH = p.adjust(df_lr_sexstratify$p, method = "BH")
                         
write.csv(df_lr_sexstratify, file="../results/lr_betas_sexstratify_initial.csv")

#females: age + famhx + education + townsend + diabetes + depression + stroke + aspirin
#males: age + famhx + education + townsend + diabetes + depression + stroke + antihyper

#create sexstratified predicted prob
UKBDRS_LASSO_female  <-  paste("dementia_BIN_TOTAL ~  Age_when_attended_assesment_centre_0_0 +  family_history_of_dementia +
                            education_years + Townsend_deprivation_modelvar +Diabetes_BIN_FINAL_0_0  +
                            current_history_depression + stroke_TIA_BIN_FINAL +  
                            Aspirin_0_0")

UKBDRS_LASSO_male  <-  paste("dementia_BIN_TOTAL ~  Age_when_attended_assesment_centre_0_0 +  family_history_of_dementia +
                            education_years + Townsend_deprivation_modelvar +Diabetes_BIN_FINAL_0_0  +
                            current_history_depression + stroke_TIA_BIN_FINAL +  
                            Antihypertensive_meds_0_0")

for (m in c("UKBDRS_LASSO_female","UKBDRS_LASSO_male")){
  print(paste0('applying logistic regression model for ', m))
  model <- glm(as.formula(m), data=train.data, family="binomial")
  
  print(paste0('UKB training set : calculating linear predictor and predicted probabilities for ', m))
  train.data[paste(m, "linear_predictor", sep="_")] <- predict(model, train.data)
  train.data[paste(m, "predicted_prob", sep="_")] <-   1/(1+exp(-train.data[paste(m, "linear_predictor", sep="_")])) # can also be computed with predict(model, type='response')
  
  print(paste0('UKB test set : calculating linear predictor and predicted probabilities for ', m))
  test.data[paste(m, "linear_predictor", sep="_")] <- predict(model, test.data)
  test.data[paste(m, "predicted_prob", sep="_")] <-   1/(1+exp(-test.data[paste(m, "linear_predictor", sep="_")])) # can also be computed with predict(model, newdata= test.data, type='response') 
}

#now create final UKBDRS_LASSO_sexstratify linear predictor
train.data$UKBDRS_LASSO_sexstratify_linear_predictor<-ifelse(train.data$Sex==0, train.data$UKBDRS_LASSO_female_linear_predictor,
                                                             ifelse(train.data$Sex==1, train.data$UKBDRS_LASSO_male_linear_predictor,NA))
test.data$UKBDRS_LASSO_sexstratify_linear_predictor<-ifelse(test.data$Sex==0, test.data$UKBDRS_LASSO_female_linear_predictor,
                                                             ifelse(test.data$Sex==1, test.data$UKBDRS_LASSO_male_linear_predictor,NA))

#now create final UKBDRS_LASSO_sexstratify predicted prob
train.data$UKBDRS_LASSO_sexstratify_predicted_prob<-ifelse(train.data$Sex==0, train.data$UKBDRS_LASSO_female_predicted_prob,
                                                             ifelse(train.data$Sex==1, train.data$UKBDRS_LASSO_male_predicted_prob,NA))
test.data$UKBDRS_LASSO_sexstratify_predicted_prob<-ifelse(test.data$Sex==0, test.data$UKBDRS_LASSO_female_predicted_prob,
                                                            ifelse(test.data$Sex==1, test.data$UKBDRS_LASSO_male_predicted_prob,NA))
summary(test.data$UKBDRS_LASSO_predicted_prob)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#0.0006103 0.0043762 0.0107550 0.0172858 0.0235672 0.2792533 
summary(test.data$UKBDRS_LASSO_sexstratify_predicted_prob)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#0.0006096 0.0043077 0.0108167 0.0172328 0.0233764 0.3156873 
plot(test.data$UKBDRS_LASSO_sexstratify_predicted_prob, test.data$UKBDRS_LASSO_predicted_prob)

#save train and test data now with sex stratified model as well, to compare in next script
save(train.data, file = paste0(data_pathway, "train_data_outliers_removed_fiftyplusnoapoe.rda"))
save(test.data, file = paste0(data_pathway, "test_data_outliers_removed_fiftyplusnoapoe.rda"))

#get final betas and write out, to test in wh
df_lr_sexstratify<-data.frame(matrix(ncol=6))
names(df_lr_sexstratify)<-c("Predictor","beta","lower","upper","OR","p")

#females
model <- glm(as.formula(UKBDRS_LASSO_female), data=subset(train.data, train.data$Sex==0), family="binomial")
df_lr_sexstratify<-rbind(df_lr_sexstratify, format_modelcoefs(model))
#males
model <- glm(as.formula(UKBDRS_LASSO_male), data=subset(train.data, train.data$Sex==1), family="binomial")
df_lr_sexstratify<-rbind(df_lr_sexstratify, format_modelcoefs(model))
df_lr_sexstratify$FDR_BH = p.adjust(df_lr_sexstratify$p, method = "BH")

write.csv(df_lr_sexstratify, file="../results/lr_betas_sexstratify_final.csv")
