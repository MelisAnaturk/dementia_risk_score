# Created by Melis Anaturk (Feb 2020)

# This r script computes the discriminative ability and calibration of each model
# similar to Licher et al (2018)

# 1. Assess discriminative ability
# 2. Assess calibration
# 3. Cut-offs for UKB-DRS at specific sensitivities/specificities

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


#----- 1. ASSESS DISCRIMINATIVE ABILITY
# report the AUC and 95% confidence intervals for each risk model (using predicted probabilities)

load(file="../../raw_data/train_data_outliers_removed.rda")
load(file="../../raw_data/test_data_outliers_removed.rda")

#NB. Anu-adri has to be excluded from calibration calculations
models <- c("age_only", "UKBDRS_LASSO", "UKBDRS_LASSO_MAN", "UKBDRS_APOE_LASSO", "UKBDRS_APOE_LASSO_MAN", "CAIDE", "DRS")

test.data$dataset <- "test"
train.data$dataset <- "train"

df_test <- rbind(test.data, train.data) 

datasets <- c("test")

savepath = "../results/"

#dataframe for storing/organizing auc results for each model
df_sitable5<-data.frame(matrix(ncol=3))
names(df_sitable5)<-c("Model","train","test")

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
    df_auc[d]<-paste(round(roc$auc[1],2), " [", round(roc$ci[1],2), ",", round(roc$ci[2],2),"]",sep="")
    print('========================================================')
    print('========================================================')
  }
  df_sitable5<-rbind(df_sitable5, df_auc)
}

# for anu_adri only
df_auc<-data.frame(matrix(ncol=3))
names(df_auc)<-c("Model","train","test")
df_auc$Model<-"ANU_ADRI"
for (d in c("train","test")){
  data <- subset(df_test, dataset==d)
  print(paste0('AUC for anu-adri in ', d))  
  roc <- pROC::roc(data[, c("dementia_BIN_TOTAL")], data[, c("ANU_ADRI")], plot=TRUE, smooth = FALSE, ci=TRUE)
  print(roc$auc)
  print(roc$ci)
  df_auc[d]<-paste(round(roc$auc[1],2), " [", round(roc$ci[1],2), ",", round(roc$ci[2],2),"]",sep="")
}
df_sitable5<-rbind(df_sitable5, df_auc)


# pairwise comparisons of AUC

#----- 2. ASSESS CALIBRATION
# Calculate calibration metrics
# More info at https://academic.oup.com/jamia/article/27/4/621/5762806
# github page: https://github.com/easonfg/cali_tutorial

library(ResourceSelection)

### functions to compute spiegelhalter z and brier score
Spiegelhalter_z = function(y, prob){
  alpha = 0.05
  z_score = sum((y-prob)*(1-2*prob))/sqrt(sum(((1-2*prob)^2)*prob*(1-prob)))
  print(z_score)
  if (abs(z_score) > qnorm(1-alpha/2)){
    print('reject null. NOT calibrated')
  } else{
    print('fail to reject. calibrated')
  }
  cat('z score: ', z_score, '\n')
  cat('p value: ', 1-pnorm(abs(z_score)), '\n')
  return(z_score)
}

rescale_Brier = function(B, y){  # takes in brier score and vector of outcomes (0/1)
  Bmax  <- mean(y) * (1-mean(y))^2 + (1-mean(y)) * mean(y)^2
  Bscaled <- 1 - B/Bmax
  return(Bscaled)
}

# create empty df tos save results
df_calibration_table3 <- data.frame(matrix(nrow = 0,ncol=7))

# change model headings
names(df_calibration_table3) <- c("Model", "Intercept", "Slope", "Chi-squared", "Brier_Score", "Spiegelhalter_z_test", "p.value")

# for loop to populate df_calibration_table3
for (m in models){
  for (d in datasets){
    data <- subset(df_test, dataset==d)

    data$y <- ifelse(data[,"dementia_BIN_TOTAL"]==1,1,0)
    
    vec <-val.prob(data[, paste(m, "predicted_prob", sep="_")], data$y, g=10, pl=TRUE, smooth=TRUE, logistic.cal=FALSE, lim=c(0,0.4))
    print(vec)
    dev.off()
    
    print(round(vec[17],2))
    
    # Brier score
    print('Brier score')
    print(round(vec[11],2))
    
    print('Rescaled brier score')
    print(rescale_Brier(vec[17], data$y))
    #
    print('Spiegelhalter z test')  #intercept_slope
    Spiegelhalter_z(data$y, data[, paste(m, "predicted_prob", sep="_")])
    
    df_calib_stats <- data.frame(c(print(m), round(vec[12],3), round(vec[13],3), formatC(vec[11],format = "e", digits = 2),round(vec[17],3),formatC(vec[18],format = "e", digits = 2)))
    df_calib_stats  <- data.table::transpose(df_calib_stats)
    
    names(df_calib_stats) <- c("Model", "Intercept", "Slope", "Brier_Score", "Spiegelhalter_z_test", "p.value")
    
    df_calibration_table3 <- rbind(df_calibration_table3,df_calib_stats)
    
    print(paste0("Intercept: ",round(vec[12],2)))
    print(paste0("Slope: ",round(vec[13],2)))
    
    # HOSLEM TEST
    #print(paste0("HL test for ", m, " in ", d))
    #hl <- hoslem.test(data$y, data[, paste(m, "predicted_prob", sep="_")], g=10) 
    #print(hl)
    #print('========================================================')
    #print('========================================================')
    
    
    #pdf(file=paste0(savepath,"calibration_plot_for_",m, "_",d,"set_intercept_slope.pdf"))
    #pdf(file=paste0(savepath,"calibration_plot_for_",m, "_",d,"set_LATEST.pdf"))
    #plot <- val.prob(data[, paste(m, "predicted_prob", sep="_")], data$y, g=10, pl=TRUE, smooth=TRUE, logistic.cal=FALSE, lim=c(0,0.4))
    #plot <- val.prob(data[, paste(m, "predicted_prob_intercept_slope_update", sep="_")], data$y, g=10, pl=TRUE, smooth=TRUE, logistic.cal=FALSE, lim=c(0,1))
    #print(plot)
    #dev.off()
  }}

#---- Model comparison
# adding variables to the baseline model
#---- comparing the model to age only
#---- create ROC objects
# ROC curve
age_only    <-pROC::roc(test.data$dementia_BIN_TOTAL, test.data$age_only_predicted_prob, plot=TRUE, smooth = FALSE, ci=TRUE)
UKBDRS_APOE_LASSO  <-pROC::roc(test.data$dementia_BIN_TOTAL, test.data$UKBDRS_APOE_LASSO_predicted_prob, plot=TRUE, smooth = FALSE, ci=TRUE)
UKBDRS_APOE_LASSO_MAN  <-pROC::roc(test.data$dementia_BIN_TOTAL, test.data$UKBDRS_APOE_LASSO_MAN_predicted_prob, plot=TRUE, smooth = FALSE, ci=TRUE)
UKBDRS_LASSO  <-pROC::roc(test.data$dementia_BIN_TOTAL, test.data$UKBDRS_LASSO_predicted_prob, plot=TRUE, smooth = FALSE, ci=TRUE)
UKBDRS_LASSO_MAN  <-pROC::roc(test.data$dementia_BIN_TOTAL, test.data$UKBDRS_LASSO_MAN_predicted_prob, plot=TRUE, smooth = FALSE, ci=TRUE)

#FRS    <-pROC::roc(test.data$dementia_BIN_TOTAL, test.data$FRS_predicted_prob, plot=TRUE, smooth = FALSE, ci=TRUE)
CAIDE  <-pROC::roc(test.data$dementia_BIN_TOTAL, test.data$CAIDE_predicted_prob, plot=TRUE, smooth = FALSE, ci=TRUE)
DRS    <-pROC::roc(test.data$dementia_BIN_TOTAL, test.data$DRS_predicted_prob, plot=TRUE, smooth = FALSE, ci=TRUE)

#run all comparisons
all_tests <- combn(list(age_only, UKBDRS_LASSO, UKBDRS_LASSO_MAN, UKBDRS_APOE_LASSO, UKBDRS_APOE_LASSO_MAN,
                        CAIDE, DRS),
                   FUN = function(x, ...) roc.test(x[[1]], x[[2]]),
                   m = 2,
                   simplify = FALSE, 
                   reuse.auc = TRUE, 
                   method = "delong", 
                   na.rm = TRUE
)

# create list of names
#tests_names <-combn(list("age_only", "UKBDRS_LASSO", "UKBDRS_LASSO_MAN", "UKBDRS_APOE_LASSO", "UKBDRS_APOE_LASSO_MAN",
#                         "CAIDE", "DRS"), 
#                    m = 2, 
#                    FUN = paste, 
#                    simplify = TRUE, 
#                    collapse = "_"
#)


#create list of model names to be compared, using naming convention in paper
comparison_names <-combn(list("ageonly", "Model2", "Model4", "Model1", "Model3",
                         "CAIDE", "DRS"), 
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

# Now to ANU-ADRI (due to some missing data on this risk score)
# exclude ppl with missing ANU-ADRI
anu.test.data <- subset(test.data, complete.cases(ANU_ADRI))
ANU_ADRI  <-pROC::roc(anu.test.data$dementia_BIN_TOTAL, anu.test.data$ANU_ADRI, plot=TRUE, smooth = FALSE, ci=TRUE)
UKBDRS_APOE_LASSO_anu  <-pROC::roc(anu.test.data$dementia_BIN_TOTAL, anu.test.data$UKBDRS_APOE_LASSO_predicted_prob, plot=TRUE, smooth = FALSE, ci=TRUE)
UKBDRS_APOE_LASSO_MAN_anu  <-pROC::roc(anu.test.data$dementia_BIN_TOTAL, anu.test.data$UKBDRS_APOE_LASSO_MAN_predicted_prob, plot=TRUE, smooth = FALSE, ci=TRUE)
UKBDRS_LASSO_anu  <-pROC::roc(anu.test.data$dementia_BIN_TOTAL, anu.test.data$UKBDRS_LASSO_predicted_prob, plot=TRUE, smooth = FALSE, ci=TRUE)
UKBDRS_LASSO_MAN_anu  <-pROC::roc(anu.test.data$dementia_BIN_TOTAL, anu.test.data$UKBDRS_LASSO_MAN_predicted_prob, plot=TRUE, smooth = FALSE, ci=TRUE)

# plot ROC curves
library(extrafont)
g2 <- ggroc(list(Age_only=age_only, UKBDRS_Model1=UKBDRS_APOE_LASSO, UKBDRS_Model2=UKBDRS_LASSO, UKBDRS_Model3=UKBDRS_APOE_LASSO_MAN, UKBDRS_Model4=UKBDRS_LASSO_MAN, CAIDE = CAIDE, DRS = DRS, ANU_ADRI = ANU_ADRI))
plot <- g2 + theme_minimal()  +  theme(legend.title = element_blank(), panel.grid.major = element_blank(), 
                                                                      panel.grid.minor = element_blank(),
                                                                      panel.background = element_rect(colour = "black", size=1), text = element_text(size=14, family="LM Roman 10")) 
plot <- g2 + theme_minimal()  +  theme(legend.title = element_blank(), panel.grid.major = element_blank(), 
                                       panel.grid.minor = element_blank(),
                                       panel.background = element_rect(colour = "black", size=1), text = element_text(size=14)) 

ggsave(paste0(savepath,"roc_plotted_all.pdf"), width = 10, height = 10, units = "cm")

ggsave(paste0(savepath,"roc_plotted_all.pdf"), 
       plot = plot, 
       device = "pdf", 
       dpi = 320)

ggsave(paste0(savepath,"roc_plotted_all.png"), 
       plot = plot, 
       device = "png", 
       dpi = 320)


#run all comparisons
all_tests <- combn(list(UKBDRS_LASSO, UKBDRS_LASSO_MAN, UKBDRS_APOE_LASSO, UKBDRS_APOE_LASSO_MAN,
                        ANU_ADRI),
                   FUN = function(x, ...) roc.test(x[[1]], x[[2]]),
                   m = 2,
                   simplify = FALSE, 
                   reuse.auc = TRUE, 
                   method = "delong", 
                   na.rm = TRUE
)


# create list of names
#tests_names <-combn(list("UKBDRS_LASSO", "UKBDRS_LASSO_MAN", "UKBDRS_APOE_LASSO", "UKBDRS_APOE_LASSO_MAN",
#                         "ANU_ADRI"), 
#                    m = 2, 
#                    FUN = paste, 
#                    simplify = TRUE, 
#                    collapse = "_"
#)

#create list of model names to be compared, using naming convention in paper
comparison_names <-combn(list("Model2", "Model4", "Model1", "Model3",
                              "ANUADRI"), 
                         m = 2, 
                         FUN = paste, 
                         simplify = TRUE, 
                         collapse = "_"
)
comparison_names<-data.frame(comparison_names)
#separate the combined comparison names into two new columns, for easy sorting
comparison_names_reorg <- comparison_names %>% separate(comparison_names, c("Score1","Score2"), sep="_", remove = FALSE)

# clean up results
all_tests <- setNames(all_tests, comparison_names$comparison_names)
tidy_results <- lapply(all_tests, broom::tidy)

# convert lists to df
comparison_list <- as.data.frame(comparison_names)
lstData <- Map(as.data.frame, tidy_results)
AUC_comparisons <- rbindlist(lstData, fill=TRUE) #fill=TRUE)

# add a column to serve as a key variable
AUC_comparisons$Number <- 1:nrow(AUC_comparisons) 
comparison_names_reorg$Number <- 1:nrow(comparison_names_reorg) 

# merge based on key variable
merged_results2 <- merge(AUC_comparisons, comparison_names_reorg, by.c="Number", all.x=TRUE)
merged_results <- rbind(merged_results,merged_results2, fill=TRUE)


# sort according to p-value
attach(merged_results)
merged_results <- merged_results[order(p.value),]
detach(merged_results)

# correct for multiple comparisons
merged_results$FDR_BH = p.adjust(merged_results$p.value, method = "BH")
merged_results

# filter data to include only significant
AUC_comparisons_corrected <- dplyr::filter(merged_results, FDR_BH <= 0.05)

#si table 6 only needs some of this, select necessary columns for easy transfer
df_sitable6<-data.frame(cbind(AUC_comparisons_corrected$Score1, AUC_comparisons_corrected$Score2,
                   round(AUC_comparisons_corrected$estimate1,2), round(AUC_comparisons_corrected$estimate2,2),
                   round(AUC_comparisons_corrected$statistic,2), AUC_comparisons_corrected$p.value,
                   AUC_comparisons_corrected$FDR_BH))
names(df_sitable6)<-c("Risk Score 1","Risk Score 2","AUC 1","AUC 2","Z","p","pcorr")

#----- 3. RISK SCORE CUT-OFF
# getting threshold at specific sensitivities/specificities" https://stackoverflow.com/questions/33125558/get-optimal-threshold-with-at-least-75-sensitivity-with-proc-in-r


models <- c("UKBDRS_LASSO", "UKBDRS_LASSO_MAN", "UKBDRS_APOE_LASSO", "UKBDRS_APOE_LASSO_MAN") 
models <- c("UKBDRS_APOE_LASSO_MAN") #,"UKBDRS_APOE_LASSO_MAN")
for (m in models){for (d in c("train")){
  print(sprintf('-----------Reporting AUC for %s for %s model---------------', d, m))
  data <- subset(df_test, dataset==d)
  print(paste0('computing AUC for ', m, ' in ', d))
  roc <- pROC::roc(data[, c("dementia_BIN_TOTAL")], data[, paste(m, "predicted_prob", sep="_")], plot=TRUE, smooth = FALSE, ci=TRUE)
  print(roc$auc)
  print(roc$ci)
  print('potential cut-off based on 80% sensitivity')
  perform<-data.frame(pROC::coords(roc, x = "all", ret=c("threshold", "specificity", "sensitivity","npv", "ppv"), transpose = FALSE))
  perform <- round(perform[perform$sensitivity >= .8 & perform$sensitivity <= 0.811, ],3)
  print(perform)
  print('potential cut-off based on 85% sensitivity')
  perform<-data.frame(pROC::coords(roc, x = "all", ret=c("threshold", "specificity", "sensitivity","npv", "ppv"), transpose = FALSE))
  perform <- round(perform[perform$sensitivity >= .85 & perform$sensitivity <= 0.851, ],3)
  print(perform)
  print('potential cut-off based on 90% sensitivity')
  perform<-data.frame(pROC::coords(roc, x = "all", ret=c("threshold", "specificity", "sensitivity","npv", "ppv"), transpose = FALSE))
  perform <- round(perform[perform$sensitivity >= .9 & perform$sensitivity <= 0.901, ],3)
  print(perform)
  print('potential cut-off based on 95% sensitivity')
  perform<-data.frame(pROC::coords(roc, x = "all", ret=c("threshold", "specificity", "sensitivity","npv", "ppv"), transpose = FALSE))
  perform <- round(perform[perform$sensitivity >= .949 & perform$sensitivity <= 0.952, ],3)
  print(perform)
  print('----------------------------------------------------------------')
  print('potential cut-off based on 80% specificity')
  perform<-data.frame(pROC::coords(roc, x = "all", ret=c("threshold", "specificity", "sensitivity","npv", "ppv"), transpose = FALSE))
  perform <- round(perform[perform$specificity >= .79 & perform$specificity <= 0.801, ],3)
  print(perform)
  print('potential cut-off based on 85% specificity')
  perform<-data.frame(pROC::coords(roc, x = "all", ret=c("threshold", "specificity", "sensitivity","npv", "ppv"), transpose = FALSE))
  perform <- round(perform[perform$specificity >= .84 & perform$specificity <= 0.851, ],3)
  print(perform)
  print('potential cut-off based on 90% specificity')
  perform<-data.frame(pROC::coords(roc, x = "all", ret=c("threshold", "specificity", "sensitivity","npv", "ppv"), transpose = FALSE))
  perform <- round(perform[perform$specificity >= .9 & perform$specificity <= 0.901, ],3)
  print(perform)
  print('potential cut-off based on 95% specificity')
  perform<-data.frame(pROC::coords(roc, x = "all", ret=c("threshold", "specificity", "sensitivity","npv", "ppv"), transpose = FALSE))
  perform <- round(perform[perform$specificity >= .949 & perform$specificity <= 0.952, ],3)
  print(perform)
  print('----------------------------------------------------------------')
}}

