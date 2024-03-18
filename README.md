# Green space exposure and blood DNA methylation at birth and in childhood – a multi-cohort study

In this GitHub, you will find scripts for the analyses performed in this study. 

- The cohort-specific epigenome-wide association studies (EWAS):  You will find the scripts used in local and through Data Aggregation Through Anonymous Summary-statistics from Harmonised Individual-level Databases (DataSHIELD).

    - LC-LIT-ATH-EWAS-GreenSpaces-Cordblood_NoDataShield_code_v2.R
    - LC-LIT-ATH-EWAS-GreenSpaces-Childhood_NoDataShield_code_v2.R
    - LC-LIT-ATH_EWAS-GreenSpaces-DataShield.R

- Example in the main model of the Quality control (QC) of the cohort-specific results using the `EASIER R package`.
    -  QC_GS_CordBlood_M3.R
    - QC_GS_Child_Blood_M3.R
      
- Example in the main model of the fixed-effects inverse variance-weight meta-analyses using the using the `EASIER R package`
    - MA_GS_CordBlood_M3.R
    - MA_GS_ChildBlood_M3.R
    
- Example in the main model of the leave-one-out analysis, in which we re-ran the main analysis repeatedly with one of the cohorts removed each time
    - LOO_GS_CordBlood.R
    - LOO_GS_ChildBlood.R

-	Example in the in cord blood main model of the Differentially methylated regions (DMRs) using `DMRcate R package` and `Enmix-combp R package` on the meta-analysed results.
    - DMR_GS_CordBlood_M3.R

