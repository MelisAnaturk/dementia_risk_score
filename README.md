# Dementia risk project

## 1. Initial data cleaning
I've used fsl's funpack here to clean the UKB data, e.g.:

```funpack --tsv_sep "," -cfg myconfig.cfg ukb43749_clean.csv ukb43749.csv```

## 2. Mining gp prescription and clinical events data 
```extract_diagnoses_from_primary_care.py``` and ```extract_prescriptions_from_primary_care.py```

## 3. R scripts
These scripts preprocess and run logistic regression analyses in the UK biobank cohort. 
