#load required packages
library(ggplot2)
library(caret)
library(tidyverse)
library(broom)
library(rcompanion)
library(data.table)
library(sjPlot)
library(pROC)
library(CalibrationCurves)
#library(epicalc)

library(tidyverse)
library(rms)
library(Hmisc)
library(knitr)
library(broom)
library(pander)
library(ggbeeswarm)
library(gridExtra)
library(grid)
library(sjPlot)
library(sjmisc)
library(sjlabelled)
library(viridis)

load(file="../../raw_data/train_data_outliers_removed.rda")
load(file="../../raw_data/test_data_outliers_removed.rda")


#add versions of the UKB DRS models which include statins and aspirin, then test the auc diffs btwn those without
UKBDRS_LASSO_meds  <-  paste("dementia_BIN_TOTAL ~  Age_when_attended_assesment_centre_0_0 +  Sex + education_years +
                            Diabetes_BIN_FINAL_0_0  +  current_history_depression + stroke_TIA_BIN_FINAL + 
                            family_history_of_dementia + Townsend_deprivation_modelvar + Antihypertensive_meds_0_0 +
                            statins_0_0 + Aspirin_0_0")

UKBDRS_APOE_LASSO_meds <-   paste("dementia_BIN_TOTAL ~  Age_when_attended_assesment_centre_0_0 +  Sex + education_years +
                            Diabetes_BIN_FINAL_0_0  +  current_history_depression + stroke_TIA_BIN_FINAL + 
                            family_history_of_dementia + Townsend_deprivation_modelvar + Antihypertensive_meds_0_0 +
                            statins_0_0 + Aspirin_0_0 + APOE_genotype_bin")

models <- c("UKBDRS_LASSO_meds", "UKBDRS_APOE_LASSO_meds")

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

datasets <- c("test")

#dataframe for storing/organizing auc results for each model
auctable<-data.frame(matrix(ncol=3))
names(auctable)<-c("Model","train","test")

models <- c("UKBDRS_LASSO", "UKBDRS_APOE_LASSO","UKBDRS_LASSO_meds", "UKBDRS_APOE_LASSO_meds")

#cycle through train/test for each model, store auc
for (m in models){
  df_auc<-data.frame(matrix(ncol=3))
  names(df_auc)<-c("Model","train","test")
  df_auc$Model<-m
  for (d in c("train","test")){
    data <- subset(df_test, dataset==d)
    print(paste0('computing AUC for ', m, ' in ', d))
    # AUC part
    roc <- pROC::roc(data[, c("dementia_BIN_TOTAL")], data[, paste(m, "predicted_prob", sep="_")], plot=TRUE, smooth = FALSE, ci=TRUE)
    print(roc$auc)
    print(roc$ci)
    df_auc[d]<-paste(round(roc$auc[1],5), " [", round(roc$ci[1],5), ",", round(roc$ci[3],5),"]",sep="")
    print('========================================================')
    print('========================================================')
  }
  auctable<-rbind(auctable, df_auc)
}

#aucs in auctable look very similar. 
#do formal model comparison to be sure
UKBDRS_APOE_LASSO  <-pROC::roc(test.data$dementia_BIN_TOTAL, test.data$UKBDRS_APOE_LASSO_predicted_prob, plot=TRUE, smooth = FALSE, ci=TRUE)
UKBDRS_APOE_LASSO_meds  <-pROC::roc(test.data$dementia_BIN_TOTAL, test.data$UKBDRS_APOE_LASSO_meds_predicted_prob, plot=TRUE, smooth = FALSE, ci=TRUE)
UKBDRS_LASSO  <-pROC::roc(test.data$dementia_BIN_TOTAL, test.data$UKBDRS_LASSO_predicted_prob, plot=TRUE, smooth = FALSE, ci=TRUE)
UKBDRS_LASSO_meds  <-pROC::roc(test.data$dementia_BIN_TOTAL, test.data$UKBDRS_LASSO_meds_predicted_prob, plot=TRUE, smooth = FALSE, ci=TRUE)

#run all comparisons
all_tests <- combn(list(UKBDRS_APOE_LASSO, UKBDRS_LASSO, UKBDRS_APOE_LASSO_meds, UKBDRS_LASSO_meds),
                   FUN = function(x, ...) roc.test(x[[1]], x[[2]]),
                   m = 2,
                   simplify = FALSE, 
                   reuse.auc = TRUE, 
                   method = "delong", 
                   na.rm = TRUE
)

#create list of model names to be compared, using naming convention in paper
comparison_names <-combn(list("Model1", "Model2", "Model1meds","Model2meds"), 
                         m = 2, 
                         FUN = paste, 
                         simplify = TRUE, 
                         collapse = "_"
)

comparison_names<-data.frame(comparison_names)
#separate the combined comparison names into two new columns, for easy sorting
comparison_names_reorg <- comparison_names %>% separate(comparison_names, c("Score1","Score2"), sep="_", remove = FALSE)

# clean up results
all_tests <- setNames(all_tests, comparison_names)
tidy_results <- lapply(all_tests, broom::tidy)

# convert lists to df
comparison_list <- as.data.frame(comparison_names)
lstData <- Map(as.data.frame, tidy_results)
AUC_comparisons <- rbindlist(lstData) #fill=TRUE)

# add a column to serve as a key variable
AUC_comparisons$Number <- 1:nrow(AUC_comparisons) 
comparison_names_reorg$Number <- 1:nrow(comparison_names_reorg) 

# merge based on key variable
merged_results <- merge(AUC_comparisons, comparison_names_reorg, by.c="Number", all.x=TRUE)
merged_results <- dplyr::filter(merged_results, grepl('Model', comparison_names))

#Model 1 is no diff than model 1 with meds
#Model 2 is no diff than model  2 with meds