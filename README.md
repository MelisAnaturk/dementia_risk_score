# Dementia risk project

## 1. FSL's funpack
I've used fsl's funpack here to preprocess the UKB data.

```funpack --tsv_sep "," -cfg myconfig.cfg ukb43749_clean.csv ukb43749.csv```

## 2. R scripts
These scripts preprocess and run logistic regression analyses in the UK biobank cohort. 
