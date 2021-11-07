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

## 2. Mining gp prescription and clinical events data 
Dementia cases were identified through HES records, death reports, primary care recodes and self-report (at baseline only, to exclude for pre-existing cases).

I recommend first familiarizing yourself with the [Primary care data](https://biobank.ndph.ox.ac.uk/showcase/ukb/docs/primary_care_data.pdf).
The script ```extract_diagnoses_from_primary_care.py``` searches through the gp clinical events records (i.e., ```gp_clinical.csv```) and identifies participants with read codes (e.g., bnf, dmd, read_v2) aligning to several diseases of interest (e.g., dementia, stroke). Similarly, ```extract_prescriptions_from_primary_care.py``` mines the data for prescription information for dementia. 

## 3. Dementia risk prediction
These scripts preprocess and run logistic regression analyses in the UK biobank cohort. 
