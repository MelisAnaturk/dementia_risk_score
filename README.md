# UK Biobank Analyses
N.B. I'll forward the UKB data via Oxfiles and the WHII data will be accessible via the DPUK portal.

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

I recommend first familiarizing yourself with [this pdf describing the primary care data](https://biobank.ndph.ox.ac.uk/showcase/ukb/docs/primary_care_data.pdf).
The script ```extract_diagnoses_from_primary_care.py``` searches through the gp clinical events records (i.e., ```gp_clinical.csv```) and identifies participants with read codes (e.g., bnf, dmd, read_v2) corresponding to several diseases of interest (e.g., dementia, stroke). The output is a csv file of all participants who have atleast one of the listed read codes of interest ("participants_with_dementia.csv"). These participants are then assigned a value of "1" and merged back to the original dataframe in the series of r scripts described below.
Similarly, ```extract_prescriptions_from_primary_care.py``` mines the database containing prescription records (i.e. ```gp_scripts.csv```). The list of read codes for each disease/medication class of interest was created using [all_lkps_maps_v2](https://biobank.ndph.ox.ac.uk/showcase/refer.cgi?id=592), with the final list contained in ```read_codes.csv```.

### 2.2 Self-report prescription data
To identify individuals taking hormone replacement therapies (HRTs), anti-hypertensive medications and other treatments of interest, I have grouped all of the drugs recorded in [data-field 20003](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=20003) according to their ATC codes [Wu et al. 2019, Supplementary Figure 1]. The file ```Wu_et_al._ATC_med_codes.xlsx``` contains our classifications for HRT and other medications.

## 3. Main analysis
The main analysis is performed using a series of R scripts. A brief description of each script is provided below:
1. ```1_derive_diagnoses.R```: Generates variables corresponding to an individual's status for several diseases of interest (e.g,., "1" for dementia, "0" for no dementia). It also derives corresponding ```dates_of_diagnosis``` and ```time_to_diagnosis``` variables. Next, we exclude people based on prespecified criteria, such as anyone reporting a diagnosis of dementia at baseline or anyone who developed dementia within a year of enrolling into the UKB study. Finally, we also generate variables to reflect medication use (e.g., antidepressants). 
***Please note:*** *I've used a few ```for loops``` in the above script which isn't best practice when writing code in R (It was just a quick solution!) and if you do plan to make these scripts publicly available, please do rewrite them as functions that can be applied e.g. using ```lapply```.*

2. ```2_anu_adri.R```: Recodes all demographic and genetic variables of interest to compute the ANU-ADRI.

3. ```3_caide.R```: Recodes all relevant variables to compute the CAIDE, with the predicted probability of developing dementia (according to the CAIDE) also calculated. The version of CAIDE that includes APOE is used in our main analyses as it generally better at distinguishing between patients/controls than the version of this risk score without APOE information.

4. ```4_bdsi_retired.R```: Originally used to calculate the BDSI but this script has now been retired as one of the core components cannot (e.g. money problems) be calculated in the biobank due to a lack of information.
 
5. ```5_framingham.R```: Calculates several versions of the framingham risk score as well as predicted probability of developing CVD.

6. ```6_recode_variables.R```: Here, we're recoding demographic and other variables of interest prior to our main analyses.

7. ```7_DRS.R```: Computing the DRS.
 
9. ```9_logistic_regression_LASSO.R```: This is where the main analysis is performed. First LASSO regression is performed for feature selection, followed by logistic regression to calculate the beta-weights to be used in the UKB-DRS.
10. ```10_discrimination_and_calibration.R``` : This r script computes the discriminative ability and calibration of each model.

# Whitehall Analyses
## 1. Accessing data
As we've applied for WHII data through the DPUK portal, you will first need to set up a few things up on your laptop/phone.

Mark Newbury will send you a username and a QR code once you've been added to the project. Use this username (please use lower case) to reset a password at https://portal.dementiasplatform.uk/Account/ResetPasswordRequest.

Should you experience problems logging in when following the instructions below, please try resetting your password at the following link instead of the link above: https://portal.dpuk.ukserp.ac.uk/RequestNewPassword.  

Instructions to connect are as follows:
1. Navigate to https://portal.dpuk.ukserp.ac.uk
2. For the Data Portal VDI connection, please download the latest version of **VMware Horizon Client**, which is towards the top right of the screen
3. Please also download the **FreeOTP** app on your mobile phone/tablet – available on iOS or Android.
4. Using your mobile device, open the FreeOTP and use the scan QR code function to scan the attached QR code on your PC screen. The app should recognise your user account and add this as a line below – this is now setup for use, and you do not need the QR code again.
5. Returning to https://portal.dpuk.ukserp.ac.uk , you can proceed to login with your username and password. 

To do this:
1. Select ```Launch Platform```
2. Now installed as above, allow VMware to launch, and accept the disclaimer that appears.
3. Your username should be present, if not, simply type this in the username field.
4. For the passcode box, use your mobile authentication app and tap the ```UKSeRP Mobile Authenticator line```, which will provide you with a code.
5. Type this code in the passcode box, and select login.
6. If successful, a second password prompt will appear – simply type your password and select login.
7. Double click on **DPUK Floating Desktop**.
8. The desktop will launch as a new window on your own PC.

You should then be able to access the **0346** study folder on the S Drive. The S Drive 0346 folder is to be used for shared work on the 0346 project. You will also find an individual user folder for your username on the P Drive. **I've added all scripts into the shared work folder**

 Any problems with logging in, contact Mark Newbury (m.s.newbury@swansea.ac.uk)
 
 ***Quick note about DPUK portal:*** 
*You can copy text into the desktop accessed through the portal but it's more difficult to export things out of the portal. For textual and numerical information such as model results, I've just had to hand write these into the manuscript. For plots, you need to go through a formal process whereby your request to export this information needs to be approved through a formal committee. Worth asking Mark Newbury about this if required.*

## 2. Main analysis
The scripts used on the Whitehall data are adapted from the UKB analyses described above. If you have any questions it's best just to drop me an email at melis.anaturk.14@ucl.ac.uk. The **important thing** to highlight here is that the "age" variable provided by DPUK is a **categorical variable**, i.e. age is coded in age bins. I've therefore used the median across all age bins as Mika has previous done in his papers.
