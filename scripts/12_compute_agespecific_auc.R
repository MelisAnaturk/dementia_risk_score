#compute discrimination and calibration of models in specific age ranges

#load required packages
library(pROC)

load(file="../../raw_data/train_data_outliers_removed_fiftyplusnoapoe.rda")
load(file="../../raw_data/test_data_outliers_removed_fiftyplusnoapoe.rda")

#### compute UKB-DRS ####
#calculate linear predictor and predicted probabilities for the age-only and UKB-DRS models
#using 9_logistic_regression_LASSO.R as sample, from line 258 and on

#define models
#age_only <-      paste("dementia_BIN_TOTAL~Age_when_attended_assesment_centre_0_0")

#UKBDRS_LASSO  <-            paste("dementia_BIN_TOTAL ~  Age_when_attended_assesment_centre_0_0 +
#                                  Diabetes_BIN_FINAL_0_0  +  current_history_depression + stroke_TIA_BIN_FINAL")

#UKBDRS_LASSO_MAN  <-  paste("dementia_BIN_TOTAL ~  Age_when_attended_assesment_centre_0_0 +  Sex + education_years +
#                            Diabetes_BIN_FINAL_0_0  +  current_history_depression + stroke_TIA_BIN_FINAL + 
#                            family_history_of_dementia")

#UKBDRS_APOE_LASSO <-      paste("dementia_BIN_TOTAL ~  Age_when_attended_assesment_centre_0_0 +
#                                Diabetes_BIN_FINAL_0_0  +  current_history_depression + stroke_TIA_BIN_FINAL + APOE_genotype_bin")

#UKBDRS_APOE_LASSO_MAN <-   paste("dementia_BIN_TOTAL ~  Age_when_attended_assesment_centre_0_0 +  Sex + education_years +
#                                 Diabetes_BIN_FINAL_0_0  +  current_history_depression + stroke_TIA_BIN_FINAL + family_history_of_dementia + APOE_genotype_bin")
#models <- c("age_only", "UKBDRS_LASSO", "UKBDRS_LASSO_MAN", "UKBDRS_APOE_LASSO", "UKBDRS_APOE_LASSO_MAN")

test.data$dataset <- "test"
train.data$dataset <- "train"

df_test <- rbind(test.data, train.data) 
datasets<-c("train","test")
rm(train.data, test.data)


#apply logistic regression for each model
#for (m in models){
#  print(paste0('applying logistic regression model for ', m))
#  model <- glm(as.formula(m), data=train.data, family="binomial")
  
#  print(paste0('UKB training set : calculating linear predictor and predicted probabilities for ', m))
#  train.data[paste(m, "linear_predictor", sep="_")] <- predict(model, train.data)
#  train.data[paste(m, "predicted_prob", sep="_")] <-   1/(1+exp(-train.data[paste(m, "linear_predictor", sep="_")])) # can also be computed with predict(model, type='response')
  
#  print(paste0('UKB test set : calculating linear predictor and predicted probabilities for ', m))
#  test.data[paste(m, "linear_predictor", sep="_")] <- predict(model, test.data)
#  test.data[paste(m, "predicted_prob", sep="_")] <-   1/(1+exp(-test.data[paste(m, "linear_predictor", sep="_")])) # can also be computed with predict(model, newdata= test.data, type='response') 
#}

#test.data$dataset <- "test"
#train.data$dataset <- "train"

#df_test <- rbind(test.data, train.data) 
#datasets<-c("train","test")
#rm(train.data, test.data)

#### CAIDE ####

#compute the auc of CAIDE and UKB-DRS in age specific subset of train and test data
#age range for CAIDE is 39-64

models <- c("age_only", "UKBDRS_APOE_LASSO", "UKBDRS_LASSO", "CAIDE")

#define dataframe to store results
caide_auc=data.frame(matrix(ncol = 2, nrow = 4))
rownames(caide_auc)<-models
colnames(caide_auc)<-datasets

#subset data into caide age range
caide_data <- subset(df_test, Age_at_recruitment_0_0<65)
summary(as.factor(caide_data$dataset))
#test  train 
#33318 134123 

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
    roc_results<-paste(round(roc$auc,2), " [", round(roc$ci[1],2), ", ",round(roc$ci[3],2), "]",sep="")
    caide_auc[m,d]<-roc_results
  }
}
rm(caide_data)
caide_auc
#train               test
#age_only           0.69 [0.68, 0.7]  0.7 [0.67, 0.72]
#UKBDRS_APOE_LASSO 0.75 [0.74, 0.77]  0.77 [0.74, 0.8]
#UKBDRS_LASSO      0.73 [0.72, 0.74]  0.73 [0.7, 0.76]
#CAIDE             0.65 [0.64, 0.67] 0.66 [0.63, 0.69]


#### DRS ####

#compute the auc of DRS and UKB-DRS in age specific subset of train and test data
#age range for DRS is 60-79

models <- c("age_only", "UKBDRS_APOE_LASSO", "UKBDRS_LASSO", "DRS")

#define dataframe to store results
drs_auc=data.frame(matrix(ncol = 2, nrow = 4))
rownames(drs_auc)<-models
colnames(drs_auc)<-datasets

#subset data into drs age range
drs_data <- subset(df_test, (Age_at_recruitment_0_0>=60 & Age_at_recruitment_0_0<=79))
summary(as.factor(drs_data$dataset))
#test  train 
#24811 98680  

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
    roc_results<-paste(round(roc$auc,2), " [", round(roc$ci[1],2), ", ",round(roc$ci[3],2), "]",sep="")
    drs_auc[m,d]<-roc_results
  }
}
rm(drs_data)
gc()
drs_auc
#train               test
#age_only          0.66 [0.65, 0.67] 0.67 [0.65, 0.69]
#UKBDRS_APOE_LASSO 0.75 [0.74, 0.76] 0.74 [0.72, 0.76]
#UKBDRS_LASSO       0.69 [0.68, 0.7] 0.69 [0.67, 0.71]
#DRS               0.65 [0.64, 0.67] 0.65 [0.63, 0.67]


#### ANU ADRI (MAP) ####

#compute the auc of DRS and UKB-DRS in age specific subset of train and test data
#age range for DRS is 60-79

models <- c("age_only", "UKBDRS_APOE_LASSO", "UKBDRS_LASSO","ANU_ADRI")

#define dataframe to store results
anu_map_auc=data.frame(matrix(ncol = 2, nrow = 4))
rownames(anu_map_auc)<-models
colnames(anu_map_auc)<-datasets

#subset data into anu map age range
anu_data <- subset(df_test, Age_at_recruitment_0_0>=54)
summary(as.factor(anu_data$dataset))
#test  train 
#37062 148397 
#compute auc for each model in both train and test data
#store auc and conf interval in anu_map_auc table
for (d in datasets){
  for (m in models[1:3]){
    data <- subset(anu_data, dataset==d)
    print(paste0('computing AUC for ', m, ' in ', d))
    # AUC part
    roc <- pROC::roc(data[, c("dementia_BIN_TOTAL")], data[, paste(m, "predicted_prob", sep="_")], plot=TRUE, smooth = FALSE, ci=TRUE)
    print(roc$auc)
    print(roc$ci)
    roc_results<-paste(round(roc$auc,2), " [", round(roc$ci[1],2), ", ",round(roc$ci[3],2), "]",sep="")
    anu_map_auc[m,d]<-roc_results
  }
}

# for anu_adri only
for (d in datasets){
  for (m in models[4]){
    data <- subset(anu_data, dataset==d)
    print(paste0('AUC for anu-adri in ', d))  
    roc <- pROC::roc(data[, c("dementia_BIN_TOTAL")], data[, c("ANU_ADRI")], plot=TRUE, smooth = FALSE, ci=TRUE)
    print(roc$auc)
    print(roc$ci)
    roc_results<-paste(round(roc$auc,2), " [", round(roc$ci[1],2), ", ",round(roc$ci[3],2), "]",sep="")
    anu_map_auc[m,d]<-roc_results
}
}
rm(anu_data)
anu_map_auc
#train               test
#age_only          0.72 [0.71, 0.73]  0.71 [0.7, 0.73]
#UKBDRS_APOE_LASSO  0.79 [0.78, 0.8]  0.78 [0.76, 0.8]
#UKBDRS_LASSO      0.75 [0.74, 0.76] 0.74 [0.72, 0.75]
#ANU_ADRI          0.58 [0.57, 0.59]  0.58 [0.56, 0.6]



#### ANU ADRI (CVHS) ####

models <- c("age_only", "UKBDRS_APOE_LASSO", "UKBDRS_LASSO","ANU_ADRI")

#define dataframe to store results
anu_cvhs_auc=data.frame(matrix(ncol = 2, nrow = 4))
rownames(anu_cvhs_auc)<-models
colnames(anu_cvhs_auc)<-datasets

#subset data into drs age range
anu_data <- subset(df_test, Age_at_recruitment_0_0>=62)
summary(as.factor(anu_data$dataset))
#test  train 
#18979 75370

#compute auc for each model in both train and test data
#store auc and conf interval in anu_cvhs_auc table
for (d in datasets){
  for (m in models[1:3]){
    data <- subset(anu_data, dataset==d)
    print(paste0('computing AUC for ', m, ' in ', d))
    # AUC part
    roc <- pROC::roc(data[, c("dementia_BIN_TOTAL")], data[, paste(m, "predicted_prob", sep="_")], plot=TRUE, smooth = FALSE, ci=TRUE)
    print(roc$auc)
    print(roc$ci)
    roc_results<-paste(round(roc$auc,2), " [", round(roc$ci[1],2), ", ",round(roc$ci[3],2), "]",sep="")
    anu_cvhs_auc[m,d]<-roc_results
  }
}

# for anu_adri only
for (d in datasets){
  for (m in models[4]){
    data <- subset(anu_data, dataset==d)
    print(paste0('AUC for anu-adri in ', d))  
    roc <- pROC::roc(data[, c("dementia_BIN_TOTAL")], data[, c("ANU_ADRI")], plot=TRUE, smooth = FALSE, ci=TRUE)
    print(roc$auc)
    print(roc$ci)
    roc_results<-paste(round(roc$auc,2), " [", round(roc$ci[1],2), ", ",round(roc$ci[3],2), "]",sep="")
    anu_cvhs_auc[m,d]<-roc_results
  }
}
rm(anu_data)
anu_cvhs_auc
#train               test
#age_only          0.63 [0.62, 0.64] 0.63 [0.61, 0.65]
#UKBDRS_APOE_LASSO 0.73 [0.72, 0.74]  0.72 [0.7, 0.74]
#UKBDRS_LASSO      0.67 [0.66, 0.68] 0.66 [0.64, 0.69]
#ANU_ADRI           0.59 [0.57, 0.6] 0.58 [0.56, 0.61]
