# UK Biobank Analyses
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
Just a note that I've used ```funpack``` iteratively with various categories of variables taken from the larger dataframe. A log of the various commands I've used is included in ```fsl_funpack_log.txt```

## 2. Dementia ascertain and identifying participants taking specific medications 
### 2.1 Mining gp prescription and clinical events data 
Dementia cases were identified through HES records, death reports, primary care recodes and self-report (at baseline only, to exclude for pre-existing cases). Participants were classified as a dementia cases if they had a record of a primary or secondary diagnosis of dementia in their health/death records or had been prescribed a common 'dementia drugs' (e.g., memantine).

I recommend first familiarizing yourself with the [Primary care data](https://biobank.ndph.ox.ac.uk/showcase/ukb/docs/primary_care_data.pdf).
The script ```extract_diagnoses_from_primary_care.py``` searches through the gp clinical events records (i.e., ```gp_clinical.csv```) and identifies participants with read codes (e.g., bnf, dmd, read_v2) corresponding to several diseases of interest (e.g., dementia, stroke). The output is a csv file of all participants who have atleast one of the listed read codes of interest ("participants_with_dementia.csv"). These participants are then assigned a value of "1" and merged back to the original dataframe in the series of r scripts described below.
Similarly, ```extract_prescriptions_from_primary_care.py``` mines the database containing prescription records (i.e. ```gp_scripts.csv```). The list of read codes for each disease/medication class of interest was created using [all_lkps_maps_v2](https://biobank.ndph.ox.ac.uk/showcase/refer.cgi?id=592), with the final list contained with ```read_codes.csv```.

### 2.2 Self-report prescription data
To identify individuals taking hormone replacement therapies (HRTs), anti-hypertensive medications and other treatments of interest, I have grouped all of the drugs recorded in [data-field 20003](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=20003) according to their ATC codes [Wu et al. 2019, Supplementary Figure 1]. The file ```Wu_et_al._ATC_med_codes.xlsx``` contains our classifications for HRT and other medications.

## 3. Main analysis
The main analysis is undertaken in a series of R scripts. A brief description of each script is provided below:
```test```

# Whitehall Analyses
## 1. Accessing data
As we've applied for WHII data through the DPUK portal, you will first need to set a few things up on your laptop.

Mark Newbury will send you a username and a QR code once you've been added to the project. Use this username (please use lower case) to reset a password at
https://portal.dementiasplatform.uk/Account/ResetPasswordRequest.

Should you experience problems logging in when following the instructions below, please try resetting your password at the following link instead of the link above: https://portal.dpuk.ukserp.ac.uk/RequestNewPassword.  

Instructions to connect are as follows:
1. Navigate to https://portal.dpuk.ukserp.ac.uk
2. For the Data Portal VDI connection, please download the latest version of VMware Horizon Client, which is towards the top right of the screen
3. Please also download the FreeOTP app on your mobile phone/tablet – available on iOS or Android.
4. Using your mobile device, open the FreeOTP and use the scan QR code function to scan the attached QR code on your PC screen. The app should recognise your user account and add this as a line below – this is now setup for use, and you do not need the QR code again.
5. Returning to https://portal.dpuk.ukserp.ac.uk , you can proceed to login with your username and password. 

To do this:
5.1 Select Launch Platform
5.2 Now installed as above, allow VMware to launch, and accept the disclaimer that appears.
5.3 Your username should be present, if not, simply type this in the username field.
5.4 For the passcode box, use your mobile authentication app and tap the UKSeRP Mobile Authenticator line, which will provide you with a code. 5.5 Type this code in the passcode box, and select login.
5.6 If successful, a second password prompt will appear – simply type your password and select login.
5.7 Double click on DPUK Floating Desktop.
5.8 The desktop will launch as a new window on your own PC.

You should then be able to access the 0346 study folder on the S Drive. The S Drive 0346 folder is to be used for shared work on the 0346 project. You will also find an individual user folder for your username on the P Drive.

 Any problems with logging in, contact Mark Newbury (m.s.newbury@swansea.ac.uk)
 
 ***Quick note about DPUK portal*** 
You can copy text into the desktop accessed through the portal but it's more difficult to export things out of the portal. For textual and numerical information such as model results, I've just had to hand write these into the manuscript. For plots, you need to go through a formal process whereby your request to export this information needs to be approved through a formal committee. Worth asking Mark Newbury about this if required.
