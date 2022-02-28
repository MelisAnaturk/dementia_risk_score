#compute discrimination and calibration of models in specific age ranges

#load required packages
library(pROC)

load(file="../../raw_data/train_data_outliers_removed.rda")
load(file="../../raw_data/test_data_outliers_removed.rda")

#### compute UKB-DRS ####
#calculate linear predictor and predicted probabilities for the age-only and UKB-DRS models
#using 9_logistic_regression_LASSO.R as sample, from line 258 and on

#define models
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
models <- c("age_only", "UKBDRS_LASSO", "UKBDRS_LASSO_MAN", "UKBDRS_APOE_LASSO", "UKBDRS_APOE_LASSO_MAN")

#apply logistic regression for each model
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

test.data$dataset <- "test"
train.data$dataset <- "train"

df_test <- rbind(test.data, train.data) 
datasets<-c("train","test")
rm(train.data, test.data)

#### CAIDE ####

#compute the auc of CAIDE and UKB-DRS in age specific subset of train and test data
#age range for CAIDE is 39-64

models <- c("age_only", "UKBDRS_APOE_LASSO", "UKBDRS_LASSO", 
            "UKBDRS_APOE_LASSO_MAN", "UKBDRS_LASSO_MAN", "CAIDE_APOE")

#define dataframe to store results
caide_auc=data.frame(matrix(ncol = 2, nrow = 6))
rownames(caide_auc)<-models
colnames(caide_auc)<-datasets

#subset data into caide age range
caide_data <- subset(df_test, Age_at_recruitment_0_0<65)

#compute auc for each model in both train and test data
#store auc and conf interval in caide_auc table
for (d in datasets){
  for (m in models){
    data <- subset(caide_data, dataset==d)
    print(paste0('computing AUC for ', m, ' in ', d))
    # AUC part
    roc <- pROC::roc(data[, c("dementia_BIN_TOTAL")], data[, paste(m, "predicted_prob", sep="_")], plot=TRUE, smooth = FALSE, ci=TRUE)
    print(roc$auc)
    print(roc$ci)
    roc_results<-paste(round(roc$auc,2), "[", round(roc$ci[1],2), round(roc$ci[3],2), "]",sep=" ")
    caide_auc[m,d]<-roc_results
  }
}
rm(caide_data)
gc()
caide_auc
#train               test
#age_only               0.7 [ 0.68 0.72 ]  0.7 [ 0.65 0.74 ]
#UKBDRS_APOE_LASSO     0.75 [ 0.72 0.77 ] 0.73 [ 0.68 0.78 ]
#UKBDRS_LASSO          0.73 [ 0.71 0.75 ] 0.73 [ 0.69 0.78 ]
#UKBDRS_APOE_LASSO_MAN 0.76 [ 0.73 0.78 ] 0.73 [ 0.69 0.78 ]
#UKBDRS_LASSO_MAN      0.75 [ 0.73 0.77 ] 0.74 [ 0.69 0.78 ]
#CAIDE_APOE             0.68 [ 0.66 0.7 ]  0.64 [ 0.6 0.69 ]


#### DRS ####

#compute the auc of DRS and UKB-DRS in age specific subset of train and test data
#age range for DRS is 60-79

models <- c("age_only", "UKBDRS_APOE_LASSO", "UKBDRS_LASSO", 
            "UKBDRS_APOE_LASSO_MAN", "UKBDRS_LASSO_MAN", "DRS")

#define dataframe to store results
drs_auc=data.frame(matrix(ncol = 2, nrow = 6))
rownames(drs_auc)<-models
colnames(drs_auc)<-datasets

#subset data into drs age range
drs_data <- subset(df_test, (Age_at_recruitment_0_0>=60 & Age_at_recruitment_0_0<=79))

#compute auc for each model in both train and test data
#store auc and conf interval in drs_auc table
for (d in datasets){
  for (m in models){
    data <- subset(drs_data, dataset==d)
    print(paste0('computing AUC for ', m, ' in ', d))
    # AUC part
    roc <- pROC::roc(data[, c("dementia_BIN_TOTAL")], data[, paste(m, "predicted_prob", sep="_")], plot=TRUE, smooth = FALSE, ci=TRUE)
    print(roc$auc)
    print(roc$ci)
    roc_results<-paste(round(roc$auc,2), "[", round(roc$ci[1],2), round(roc$ci[3],2), "]",sep=" ")
    drs_auc[m,d]<-roc_results
  }
}
rm(drs_data)
gc()
drs_auc
#train               test
#age_only              0.67 [ 0.65 0.69 ] 0.66 [ 0.62 0.69 ]
#UKBDRS_APOE_LASSO     0.74 [ 0.72 0.75 ] 0.75 [ 0.71 0.78 ]
#UKBDRS_LASSO           0.68 [ 0.66 0.7 ]  0.7 [ 0.67 0.74 ]
#UKBDRS_APOE_LASSO_MAN 0.74 [ 0.72 0.76 ] 0.74 [ 0.71 0.78 ]
#UKBDRS_LASSO_MAN       0.7 [ 0.68 0.71 ] 0.71 [ 0.67 0.74 ]
#DRS                   0.66 [ 0.65 0.68 ] 0.68 [ 0.64 0.71 ]
