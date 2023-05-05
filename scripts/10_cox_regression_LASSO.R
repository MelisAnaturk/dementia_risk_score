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

#----- 2. LASSO REGRESSION ------------------------------------------
# LASSO cox regression is first run to identify the subset of predictors to be used in our risk score
# https://glmnet.stanford.edu/articles/Coxnet.html 

#we need time to death/censor/dementia
df_filtered$Date_of_assessment_0_0_new.x <- as.Date(df_filtered$Date_of_assessment_0_0_new.x, format="%Y-%m-%d")
df_filtered$has_death_record <- ifelse(is.na(df_filtered$date_of_death_all),0,1)

#dementia
df_dementia <- df_filtered[which(df_filtered$dementia_BIN_TOTAL==1),]
df_dementia$time_at_risk <- df_dementia$earliest_dementia_date - df_dementia$Date_of_assessment_0_0_new.x
summary(as.numeric(df_dementia$time_at_risk))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#373    2503    3426    3239    4106    5196 

#deaths
df_deaths <- df_filtered[which( (df_filtered$dementia_BIN_TOTAL==0) & (df_filtered$has_death_record==1) ),]
df_deaths$time_at_risk <- df_deaths$date_of_death_all - df_deaths$Date_of_assessment_0_0_new.x
summary(as.numeric(df_deaths$time_at_risk))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#4    1212    2005    1931    2663    3927 

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

#load if needed
load("../../raw_data/df_filtered_prelasso.rda")

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

save(train.data, file = paste0(save_pathway, "train_data_outliers_removed_prelasso.rda"))
save(test.data, file = paste0(save_pathway, "test_data_outliers_removed_prelasso.rda"))




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
# 33 x 1 sparse Matrix of class "dgCMatrix"
#(Intercept)                            -4.59166742

#OLD, LOGISTIC      
#Townsend_deprivation_Groups_0_04        0.06505642
#family_history_of_dementia1             0.19243398
#Diabetes_BIN_FINAL_0_01                 0.45423371
#current_history_depression1             0.33689819
#stroke_TIA_BIN_FINAL1                   0.54779739
#hypertensive1                           0.09190174
#cholesterol1                            0.08886551
#Age_when_attended_assesment_centre_0_0  0.84722457
#education_years                        -0.06501424

#NEW, COX
# Townsend_deprivation_Groups_0_04        0.13382243
# Sex1                                    0.03082154
# family_history_of_dementia1             0.25174957
# Diabetes_BIN_FINAL_0_01                 0.51051696
# current_history_depression1             0.38285509
# stroke_TIA_BIN_FINAL1                   0.60957378
# household_occupancy1                    0.01811039
# hypertensive1                           0.10868495
# cholesterol1                            0.09704808
# Age_when_attended_assesment_centre_0_0  0.88294603
# education_years                        -0.07321108
#new vars are sex (male) and household occupancy (1 = living alone)


# tidy output
View(tidy(lasso.final.1))

# save the list of selected predictors (to carry forward to cox regression analyses)
write.csv(tidy(lasso.final.1), paste0(save_pathway,"cox_LASSO_results.csv"))

# save the model
save(lasso.final.1, file = paste0(model_pathway, "lasso_fit_cox_final.rda"))

# Or otherwise load model into R
#load(file = paste0(model_pathway, "lasso_fit_final.rda""))

#----- 2. COX REGRESSION ------------------------------------------
# cox regression is now run to calculate the beta-weights for each of the 
#components in our risk score (note, we also add in age, sex and education)

#coxph will want 1,2 as outcome
df_filtered$dementia_BIN_surv <- as.numeric(df_filtered$dementia_BIN_TOTAL)
summary(as.factor(df_filtered$dementia_BIN_surv))
#1      2 
#216949   3813 

# conduct train/test split (same partion as above conducted but more columns are retained)
train.data  <- df_filtered[training.samples, ]
test.data <- df_filtered[-training.samples, ]

#membership in highest deprived group (4, on scale of 0-4) was selected, not others
#create binary version of deprivation to use in model
train.data$Townsend_deprivation_modelvar<-ifelse(train.data$Townsend_deprivation_Groups_0_0==4,1,0)
train.data$Townsend_deprivation_modelvar<-as.factor(train.data$Townsend_deprivation_modelvar)
test.data$Townsend_deprivation_modelvar<-ifelse(test.data$Townsend_deprivation_Groups_0_0==4,1,0)
test.data$Townsend_deprivation_modelvar<-as.factor(test.data$Townsend_deprivation_modelvar)

save(train.data, file = paste0(save_pathway, "train_data_outliers_removed_postlasso.rda"))
save(test.data, file = paste0(save_pathway, "test_data_outliers_removed_postlasso.rda"))

load(file = paste0(save_pathway, "train_data_outliers_removed_postlasso.rda"))
load(file = paste0(save_pathway, "test_data_outliers_removed_postlasso.rda"))

#### 2.1 Test beta coefficients ####
#based on the lasso selected vars, compute model coefficients in train data
#use 2 models - one with apoe, one without. 

#### Test models ####

#### UKBDRS Linear ####
#TESTED THIS IN WH AND IT WAS OK
UKBDRS_LASSO  <-  paste("Surv(time_at_risk, dementia_BIN_surv) ~  Age_when_attended_assesment_centre_0_0 +  family_history_of_dementia +
                            education_years + Diabetes_BIN_FINAL_0_0  + Townsend_deprivation_Groups_0_0 +
                            current_history_depression + stroke_TIA_BIN_FINAL +  
                            hypertensive + cholesterol + household_occupancy + Sex")

ukbdrs.cox <- coxph(as.formula(UKBDRS_LASSO), data = train.data)
summary(ukbdrs.cox)
# coef exp(coef)  se(coef)      z Pr(>|z|)    
# Age_when_attended_assesment_centre_0_0  0.180696  1.198050  0.004643 38.914  < 2e-16 ***
#   family_history_of_dementia1             0.433780  1.543079  0.042985 10.091  < 2e-16 ***
#   education_years                        -0.043925  0.957026  0.006211 -7.072 1.53e-12 ***
#   Diabetes_BIN_FINAL_0_01                 0.568808  1.766161  0.057648  9.867  < 2e-16 ***
#   Townsend_deprivation_Groups_0_01       -0.036469  0.964188  0.059765 -0.610  0.54173    
# Townsend_deprivation_Groups_0_02        0.024363  1.024663  0.058884  0.414  0.67906    
# Townsend_deprivation_Groups_0_03        0.054756  1.056283  0.059248  0.924  0.35538    
# Townsend_deprivation_Groups_0_04        0.271078  1.311378  0.057395  4.723 2.32e-06 ***
#   current_history_depression1             0.560606  1.751735  0.046207 12.132  < 2e-16 ***
#   stroke_TIA_BIN_FINAL1                   0.695670  2.005053  0.076570  9.085  < 2e-16 ***
#   hypertensive1                           0.170602  1.186019  0.040738  4.188 2.82e-05 ***
#   cholesterol1                            0.105227  1.110962  0.044406  2.370  0.01780 *  
#   household_occupancy1                    0.144535  1.155502  0.044497  3.248  0.00116 ** 
#   household_occupancy2                   -0.053093  0.948292  0.056009 -0.948  0.34316    
# Sex1                                    0.204492  1.226901  0.037846  5.403 6.54e-08 ***

cox.zph(ukbdrs.cox)
# Age_when_attended_assesment_centre_0_0  9.0199  1 0.00267
# family_history_of_dementia              0.3350  1 0.56273
# education_years                         0.0316  1 0.85885
# Diabetes_BIN_FINAL_0_0                  0.8196  1 0.36529
# Townsend_deprivation_Groups_0_0         8.0088  4 0.09126
# current_history_depression             17.3982  1   3e-05
# stroke_TIA_BIN_FINAL                    1.5327  1 0.21570
# hypertensive                            1.1205  1 0.28981
# cholesterol                             0.6421  1 0.42297
# household_occupancy                     0.3183  2 0.85288
# Sex                                     0.7551  1 0.38488
# GLOBAL                                 40.4159 15 0.00039

#age, depr and townsend violate
ph<-cox.zph(ukbdrs.cox)


#baseline survival
#The baseline survival is the distribution of the predicted survival for the patient whose predictor values are either the average or 0 (or the reference group for categorical predictors) across the complete follow-up time under study.
#https://www.acpjournals.org/doi/full/10.7326/M22-0844
survv <- survfit(ukbdrs.cox)
summary(survfit(ukbdrs.cox))

df_base <- data.frame(Age_when_attended_assesment_centre_0_0 = mean(train.data$Age_when_attended_assesment_centre_0_0),
                      family_history_of_dementia = as.factor(0), Diabetes_BIN_FINAL_0_0 = as.factor(0),
                      current_history_depression = as.factor(0), stroke_TIA_BIN_FINAL = as.factor(0),
                      hypertensive = as.factor(0), cholesterol = as.factor(0), Townsend_deprivation_Groups_0_0 = as.factor(0),
                      education_years = mean(train.data$education_years), Sex = as.factor(0), household_occupancy = as.factor(0))

survv_baseline <- survfit(ukbdrs.cox, newdata = df_base)
length(survv_baseline$time) #there are 4614 timepoints
#survival over entire time window:
survv_baseline$surv[4614]
#0.9911458

