# Dementia risk project
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

## 2. Dementia ascertain
### 2.1 Mining gp prescription and clinical events data 
Dementia cases were identified through HES records, death reports, primary care recodes and self-report (at baseline only, to exclude for pre-existing cases). Participants were classified as a dementia cases if they had a record of a primary or secondary diagnosis of dementia in their health/death records or had been prescribed a common 'dementia drugs' (e.g., memantine).

I recommend first familiarizing yourself with the [Primary care data](https://biobank.ndph.ox.ac.uk/showcase/ukb/docs/primary_care_data.pdf).
The script ```extract_diagnoses_from_primary_care.py``` searches through the gp clinical events records (i.e., ```gp_clinical.csv```) and identifies participants with read codes (e.g., bnf, dmd, read_v2) corresponding to several diseases of interest (e.g., dementia, stroke). Similarly, ```extract_prescriptions_from_primary_care.py``` mines the database containing prescription records (i.e. ```gp_scripts.csv```). The list of read codes for each disease/medication class of interest was created using [all_lkps_maps_v2](https://biobank.ndph.ox.ac.uk/showcase/refer.cgi?id=592). 

### 2.2 Self-report prescription data
To identify individuals taking hormone replacement therapies (HRTs), anti-hypertensive medications and other treatments of interest, I have grouped all of the drugs recorded in [data-field 20003](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=20003) according to their ATC codes [Wu et al. 2019, Supplementary Figure 1]. The file ```Wu_et_al._ATC_med_codes.xlsx``` contains our classifications for HRT and other medications.

## 3. Dementia risk prediction
These scripts preprocess and run logistic regression analyses in the UK biobank cohort. 
