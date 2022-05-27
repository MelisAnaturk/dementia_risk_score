# UK Biobank Analyses
## Overview
This repository contains code used in producing the UKB-DRS risk score by Anatürk, Patel, and colleagues. This project uses data from the UK Biobank to develop a novel risk score for dementia prediction.  The score is mainly developed through a series of R scripts which derive diagnoses for dementia and other relevant conditions before deriving risk score formulae for the UKB-DRS and other previously developed risk scores.
## 1. Initial data cleaning
Fsl's funpack (https://git.fmrib.ox.ac.uk/fsl/funpack) was used to initially clean the UKB data (e.g., demographics), e.g.:

```funpack --tsv_sep "," -cfg myconfig.cfg ukb43749_clean.csv ukb43749.csv```

Where myconfig.cfg is a configuration file including the following:

```
trust_types
noisy
noisy
noisy
variable_file     fmrib/variables_parentvalues.tsv
datacoding_file   fmrib/datacodings_navalues.tsv
datacoding_file   fmrib/datacodings_recoding.tsv
log_file          clean_log.txt
unknown_vars_file unknowns.tsv
description_file  descriptions.tsv
summary_file      summary.tsv
plugin_file       fmrib
loader            FMRIB_internal_info.txt FMRIBImaging
```
I've used ```funpack``` iteratively with various classes of variables. A log of the ```funpack``` commands I've used is included in ```fsl_funpack_log.txt```

## 2. Dementia ascertain
### 2.1 Mining GP prescription and clinical events data 
Dementia cases were identified through primary care and secondary care records, death reports and self-report (at baseline only, to exclude pre-existing cases). Participants were classified as a dementia case if they had a record of a primary or secondary diagnosis of dementia in their health/death records or had been prescribed a common 'dementia drugs' (e.g., memantine).

The script ```extract_diagnoses_from_primary_care.py``` searches through the gp clinical events records (i.e., ```gp_clinical.csv```) and identifies participants with read codes (e.g., bnf, dmd, read_v2) corresponding to several diseases of interest (e.g., dementia, stroke). The output is a csv file of all participants who have atleast one of the listed read codes of interest ("participants_with_dementia.csv"). These participants are then assigned a value of "1" and merged back to the original dataframe in the series of r scripts described below.
Similarly, ```extract_prescriptions_from_primary_care.py``` mines the database containing prescription records (i.e. ```gp_scripts.csv```). The list of read codes for each disease/medication class of interest was created using [all_lkps_maps_v2](https://biobank.ndph.ox.ac.uk/showcase/refer.cgi?id=592), with the final list contained in ```read_codes.csv```.

### 2.2 Self-report prescription data
To identify individuals taking hormone replacement therapies (HRTs), anti-hypertensive medications and other treatments of interest, we have grouped all of the drugs recorded in [data-field 20003](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=20003) according to their ATC codes [Wu et al. 2019, Supplementary Figure 1]. 

## 3. Main analysis
The main analysis is performed using a series of R scripts. A brief description of each script is provided below:
1. ```1_derive_diagnoses.R```: Generates variables corresponding to an individual's status for several diseases of interest (e.g,., "1" for dementia, "0" for no dementia). It also derives corresponding ```dates_of_diagnosis``` and ```time_to_diagnosis``` variables. Next, we exclude people based on prespecified criteria, such as anyone reporting a diagnosis of dementia at baseline or anyone who developed dementia within a year of enrolling into the UKB study. Finally, we also generate variables to reflect medication use (e.g., antidepressants). 

2. ```2_anu_adri.R```: Recodes all demographic and genetic variables of interest to compute the ANU-ADRI.

3. ```3_caide.R```: Recodes all relevant variables to compute the CAIDE, with the predicted probability of developing dementia (according to the CAIDE) also calculated. The version of CAIDE that includes APOE is used in our main analyses as it generally better at distinguishing between patients/controls than the version of this risk score without APOE information.

4. ```4_bdsi_retired.R```: Originally used to calculate the BDSI but this script has now been retired as one of the core components cannot (e.g. money problems) be calculated in the biobank due to a lack of information.
 
5. ```5_framingham.R```: Calculates several versions of the framingham risk score as well as predicted probability of developing CVD.

6. ```6_recode_variables.R```: Here, we're recoding demographic and other variables of interest prior to our main analyses.

7. ```7_DRS.R```: Computing the DRS.
 
9. ```9_logistic_regression_LASSO.R```: This is where the main analysis is performed. First LASSO regression is performed for feature selection, followed by logistic regression to calculate the beta-weights to be used in the UKB-DRS.
10. ```10_discrimination_and_calibration.R``` : This r script computes the discriminative ability and calibration of each model.

 
