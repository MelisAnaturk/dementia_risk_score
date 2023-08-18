# UK Biobank Analyses
## Overview
This repository contains code used in producing the UKB-DRS risk score by Anatürk, Patel, and colleagues. This project uses data from the UK Biobank to develop a novel risk score for dementia prediction.  The score is mainly developed through a series of R scripts which derive diagnoses for dementia and other relevant conditions before deriving risk score formulae for the UKB-DRS and other previously developed risk scores.

## Contents
- ```scripts/``` contains R scripts used for data manipulation, deriving risk scores, and assessing performance
- ```results/``` contains output plots, lasso fit outputs, as well as the UKB-DRS calculator in excel format
- ```models/``` contains final lasso fitting in .rda format
All files in ```results/``` and ```models/``` are created via scripts in ```scripts/```

## Requirements
For optimal performance, we recommend a system with >=16GB of RAM and >=4 CPU cores. Development was carried out on a Ubuntu 18.04.6 LTS (Bionic Beaver) system. R version 3.6.3 was used. See scripts/session_info_R for specific package versions. For dementia ascertain python scripts, python version 3.6.5 was used.

## 1. Dementia ascertain
### 1.1 Mining GP prescription and clinical events data 
Dementia cases were identified through primary care and secondary care records, death reports and self-report (at baseline only, to exclude pre-existing cases). Participants were classified as a dementia case if they had a record of a primary or secondary diagnosis of dementia in their health/death records or had been prescribed a common 'dementia drugs' (e.g., memantine).

The script ```extract_diagnoses_from_primary_care.py``` searches through the gp clinical events records (i.e., ```gp_clinical.csv```) and identifies participants with read codes (e.g., bnf, dmd, read_v2) corresponding to several diseases of interest (e.g., dementia, stroke). The output is a csv file of all participants who have atleast one of the listed read codes of interest ("participants_with_dementia.csv"). These participants are then assigned a value of "1" and merged back to the original dataframe in the series of r scripts described below.
Similarly, ```extract_prescriptions_from_primary_care.py``` mines the database containing prescription records (i.e. ```gp_scripts.csv```). The list of read codes for each disease/medication class of interest was created using [all_lkps_maps_v2](https://biobank.ndph.ox.ac.uk/showcase/refer.cgi?id=592).

### 1.2 Self-report prescription data
To identify individuals taking hormone replacement therapies (HRTs), anti-hypertensive medications and other treatments of interest, we have grouped all of the drugs recorded in [data-field 20003](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=20003) according to their ATC codes [Wu et al. 2019, Supplementary Figure 1]. 

## 2. Main analysis
The main analysis is performed using a series of R scripts. A brief description of each script is provided below:
1. ```1_derive_diagnoses.R```: Generates variables corresponding to an individual's status for several diseases of interest (e.g,., "1" for dementia, "0" for no dementia). It also derives corresponding ```dates_of_diagnosis``` and ```time_to_diagnosis``` variables. Next, we exclude people based on prespecified criteria, such as anyone reporting a diagnosis of dementia at baseline or anyone who developed dementia within a year of enrolling into the UKB study. Finally, we also generate variables to reflect medication use (e.g., antidepressants). 

2. ```2_anu_adri.R```: Recodes all demographic and genetic variables of interest to compute the ANU-ADRI.

3. ```3_caide.R```: Recodes all relevant variables to compute the CAIDE, with the predicted probability of developing dementia (according to the CAIDE) also calculated. 

4. ```4_bdsi_retired.R```: Originally used to calculate the BDSI but this script has now been retired as one of the core components cannot (e.g. money problems) be calculated in the biobank due to a lack of information.
 
5. ```5_framingham.R```: Calculates several versions of the framingham risk score as well as predicted probability of developing CVD.

6. ```6_recode_variables.R```: Here, we're recoding demographic and other variables of interest prior to our main analyses.

7. ```7_DRS.R```: Computing the DRS risk score.
 
9. ```9_derive_hypertension_cholesterol_Dx.R```: Derive some extra diagnoses

10. ```10_cox_regression_LASSO.R```: Run LASSO for feature selection, identify most predictive factors

11. ```11_fit_cox.R```: Fit a cox model, derive coefficients

12. ```12_fit_crr.R```: Fit a competing risk model, derive coefficients

13. ```13_discrimination_and_calibration.R``` : This r script computes performance, via discriminative ability and calibration, of each model.

14. ```14_compute_agespecific_auc.R```: Test performance in age specific subsets

15. ```15_sensitivity_time_to_diagnosis.R```: Test performance in shorter time windows

16. ```16_include_vs_exclude.R```: Compare included vs excluded subjects

 
