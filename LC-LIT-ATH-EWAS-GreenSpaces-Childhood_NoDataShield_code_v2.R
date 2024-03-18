###########################################################################################
# EWAS Green Spaces and Child Blood- No DataSHIELD code (using a ExpressionSet)  
# Sofía Aguilar Lacasaña (sofia.aguilar@isglobal.org)
# 16/06/2022
#V2
###########################################################################################

############################################################################################################################################################
# The following R code will allow you to complete all the EWAS requested in the LifeCycle-LongITools-ATHLETE- EWAS of green spaces and Child blood 
#DNA methylation analysis plan.
# The code also produces files summarising the variables included in the EWAS.
# The code can be adapted to the needs of each cohort.
# There are two inputs required for this analysis:

# 1) A matrix with the methylation beta values (0-1)
# Each column is a sample and each row is a probe on the array (450K or EPIC)

# 2) A dataframe containing all the exposures + covariates + cell type proportions
# Each row is a sample(individual) and each column is a different variable. 
# Details on which variables should be included in this dataframe and how to code them are provided in the analysis plan. 
#Throughout the code the variables are named as indicated in the analysis plan. We recommend that you rename your variables as indicated. If not, you will need to manually modify the variable names in the code. 

### This code gets the methylation data and the metadata from an ExpressionSet.  If you have the matrix with the methylation data and the metadata 
#in different files, please create an ExpressionSet. In case you have any questions on how to do this, please contact us (sofia.aguilar@isglobal.org; mariona.bustamante@isglobal.org).


##########################################################################################
### STEP 0: Before starting 
##########################################################################################

##########################################
# Exclusion Criteria
### The following children will be excluded from the study:
# - Twins will be excluded. 
# - For non-twin siblings, cohorts will include only one child per mother, based on completeness of data and, if equal, randomly.

###########################################
### Change to your own working directory
setwd("/home/isglobal.lan/saguilar/data/WS_HELIX/HELIX_analyses/LCATH_EWAS_GreenSpace_SA/results/") 

###########################################
### Set initial parameters
cohort<-"HELIX" #define cohort name
Date<-Sys.Date()#define date

################################
################################
#PART 1: Prepare ExpressionSets
################################
################################

##########################################################################################
### STEP 1: Load required packages 
##########################################################################################
##(if these are not already installed, you will have to install them as the first step) --> install.packages("tidyverse") 
library(tidyverse)#data management
library(dplyr)#ata management
library(Biobase)# to be able to access and modify data in the ExpressionSet
library(EnvStats)#calculate the iqr
library(finalfit)#descriptive tables 
library(GGally)# Correlation plots
library(limma) #EWAS
library(qqman) # to get QQ plots


##########################################################################################
### STEP 2: Load the methylation data and the metadata
##########################################################################################

### Change to your own working directory where you can find the ExpressionSet
setwd("/home/isglobal.lan/saguilar/data/DS_ExpressionSet/Results")

### Load the RData with the ExpressionSet object
load("HELIX_Methyl_Child_Blood_8y_450K_20210318.Rdata")# Change the name for your ExpressionSet object

ls()#to see the R object loaded

### ExpressionSet summarized view
CreatedExpressionSet
#ExpressionSet (storageMode: lockedEnvironment)
#assayData: 480444 features, 1132 samples
#  element names: exprs
#protocolData: none
#phenoData
#  sampleNames: sample1 sample2 ... sample1132 (1132 total)
#  varLabels: id X ... gwas_pc10_eur (40 total)
#  varMetadata: labelDescription
#featureData: none
#experimentData: use 'experimentData(object)'
#Annotation: IlluminaHumanMethylation450kanno.ilmn12.hg19

###########################
### 2.1) Methylation data
###########################

### To extract the beta matrix from the ExpressionSet you should use exprs() function as you will see explained in section 8. 
# Each column is a sample and each row is a probe on the array (450k or EPIC). 

#beta_matrix<-exprs(CreatedExpressionSet) 

##################
### 2.2) Metadata 
##################

### Extract the metadata from the ExpressionSet 
# which includes some important covariates + cell type proportions

metadataEset<-pData(CreatedExpressionSet)
colnames(metadataEset)
# [1] "id"            "X"             "id_methyl"     "sex_methyl"
# [5] "age_methyl"    "anc_methyl"    "cohort_methyl" "zbmi_methyl"
# [9] "gwas_pc1"      "gwas_pc2"      "gwas_pc3"      "gwas_pc4"
#[13] "gwas_pc5"      "gwas_pc6"      "gwas_pc7"      "gwas_pc8"
#[17] "gwas_pc9"      "gwas_pc10"     "CD8T_H"        "CD4T_H"
#[21] "NK_H"          "Bcell_H"       "Mono_H"        "Gran_H"
#[25] "CD8T_S"        "CD4T_S"        "NK_S"          "Bcell_S"
#[29] "Mono_S"        "Neu_S"         "gwas_pc1_eur"  "gwas_pc2_eur"
#[33] "gwas_pc3_eur"  "gwas_pc4_eur"  "gwas_pc5_eur"  "gwas_pc6_eur"
#[37] "gwas_pc7_eur"  "gwas_pc8_eur"  "gwas_pc9_eur"  "gwas_pc10_eur"

##############################################
### STEP 3: Create cumulative Variables 
#############################################

# For the calculation of these variables, we will first calculate the average between these four time-periods:
#the age of DNA methylation assessment:
# - pregnancy
# -	infancy (>0 to <3y) 
# -	early childhood (=>3 to <=7y)
# - late childhood (>7 to <=10y)

#The cumulative variable will correspond to the average of these time-periods. Each cohort should include as many time-point measurements as
#possible while minimising the loss of sample size. Minimum of 2 time-points (pregnancy time-point and one postnatal time-point).

#NOTE: Cumulative variable is up to the age of DNA methylation. E.g. If the methylation age is at 5y, exposure up to this age will be included 
#in this variable. Therefore, in this case, the cumulative variable will correspond to the average  between pregnancy, infancy (1-2y) and early childhood (3-5y).

#Example Code in Helix: LC-ATH-Cumulative_variables_code.v1.R

##########################################################################################
### STEP 4: Load the dataframe with the exposures + additional covariates 
##########################################################################################

#### Read dataframe 
# Each row is a sample(individual) and each column is a different variable.
# Ensure all traits and covariates have been derived as specified in the analysis plan and that exclusions have been made.

df.exp.cov<-read.delim("/home/isglobal.lan/saguilar/data/WS_HELIX/HELIX_analyses/LCATH_EWAS_GreenSpace_SA/db/HELIX.pheno.total.EWASGS_Lifelong_N1132.28072022.txt",sep=",") ##change to your file name and directory
row.names(df.exp.cov) <- df.exp.cov$HelixID #to merge with the metadata from the ExpressionSet
colnames(df.exp.cov)
#[1] "h_cohort"      "e3_sex"        "agebirth_m_y"    "edu_m_0" 
#[5] "preg_smk"     "areases_tert_" "ndvi300_lf"     "ndvi100_lf"
#[9] "pm25_lf"



#################################################################################################
### STEP 5: Merge the metadata from the ExpressionSet and the dataframe with additional variables
#################################################################################################

### Merge both dataframes by the ID (it might be or not the same ID as the methylation dataset). 
dd<-merge(metadataEset,df.exp.cov,by.x="id",by.y="HelixID",all.x=TRUE)
# Now we have all variables needed for the analysis for children with methylation data

colnames(dd)
#[1] "id"            "id_methyl"     "sex_methyl"    "age_methyl"
# [5] "anc_methyl"    "cohort_methyl" "zbmi_methyl"   "gwas_pc1"
# [9] "gwas_pc2"      "gwas_pc3"      "gwas_pc4"      "gwas_pc5"
#[13] "gwas_pc6"      "gwas_pc7"      "gwas_pc8"      "gwas_pc9"
#[17] "gwas_pc10"     "CD8T_H"        "CD4T_H"        "NK_H"
#[21] "Bcell_H"       "Mono_H"        "Gran_H"        "CD8T_S"
#[25] "CD4T_S"        "NK_S"          "Bcell_S"       "Mono_S"
#[29] "Neu_S"         "gwas_pc1_eur"  "gwas_pc2_eur"  "gwas_pc3_eur"
#[33] "gwas_pc4_eur"  "gwas_pc5_eur"  "gwas_pc6_eur"  "gwas_pc7_eur"
#[37] "gwas_pc8_eur"  "gwas_pc9_eur"  "gwas_pc10_eur" "h_cohort"
#[41] "e3_sex"        "agebirth_m_y"  "edu_m_0"        "preg_smk"
#[45] "hs_age_years"  "areases_tert_" "ndvi300_lf"    "ndvi100_lf" 
#[49] "pm25_lf"


dd<-dd[,c(1,3:45,52,67:79,91:103,128:140)]
#########################################################################################################
### STEP 6: Include updated metadata (exposures+covariates+ cell type proportions) in the ExpressionSet
#########################################################################################################

###########
#IMPORTANT! 
###########

##################################################################
#6.1) Check if the samples in the merged dataframe are in the same order than in the metadata from the ExpressionSet

# If you do not order samples as they are in the ExpressionSet, you could incorrectly assign the values of the 
#variables to the samples and therefore also to the methylation. 
 
table(ifelse(metadataEset$id==dd$id,"Matched","--NOT MATCHED--"))

#--NOT MATCHED--         Matched
#            405             727

### Samples are not in the same order. We need to order in the same way

####################
#6.2) Order samples
dd.ord<-dd[order(match(dd$id,metadataEset$id)),]
table(ifelse(metadataEset$id==dd.ord$id,"Matched","--NOT MATCHED--"))

#Matched
#   1132

# Now they are in the same order

### We need the id_methyl in rows to link this new pheno+metadata with the methylation data from the ExpressionSet
rownames(dd.ord)<-dd.ord$id_methyl
                                 
#dd.ord<-rename(dd.ord,agebirth_m_y=h_age,preg_smk=e3_asmokyn_p,edu_m_0=h_edumc)                          

##############################################################
#6.3) Include the new ordered phenodata in the ExpressionSet 

pData(CreatedExpressionSet)<-dd.ord

##########################################
### STEP 7: Major ancestry group subset
##########################################

### Define group
#As ethnic groups are tested separatelly (if >100 samples), please indicate the major ancestry group you are working with

ancestry<-"EUR" 

##########################################################################
#7.1) Subset the ExpressionSet by the major ancestry group (>100 samples)

#NOTE: If you have closed the R session, re-extract the ordered phenodata from the expressionSet. If not, continue with dataframe created in STEP 6. 
#dd.ord<-pData(CreatedExpressionSet)

table(dd.ord$anc_methyl)
# 1   6
#1009  133

### I will subset children of European ancestry (group 1)
dd.EUR<-dd.ord[dd.ord$anc_methyl == 1,]
dim(dd.EUR)
#################################################
#7.2) Remove samples with NA in ancestry 
table(is.na(dd.EUR$anc_methyl))
# FALSE  
#  1009    
# There are NO NAS. 
#In case there are NAs, please delete them
#
#dd.EUR<-dd.EUR[(!is.na(dd.EUR$anc_methyl)),]
#table(is.na(dd.EUR$anc_methyl))


####################################
#7.3) Create the new ExpressionSet

### Now we create the ExpressionSet for European ancestry samples, without NAs in ancestry and without NAs in green spaces vars 

#European samples
samplesEUR<-rownames(dd.EUR)

#Subset the ExpressionSet with european samples
newmset.EUR<-CreatedExpressionSet[,(sampleNames(CreatedExpressionSet) %in% samplesEUR)]
newmset.EUR
#ExpressionSet (storageMode: lockedEnvironment)
#assayData: 480444 features, 1009 samples
#  element names: exprs
#protocolData: none
#phenoData
# sampleNames: sample1 sample2 ... sample1009 (1009 total)
#  varLabels: id X ... h_pm25_lf (133 total)
#  varMetadata: labelDescription
#featureData: none
#experimentData: use 'experimentData(object)'
#Annotation: IlluminaHumanMethylation450kanno.ilmn12.hg19

##########################################################################################
### STEP 8: WINSORIZE OUTLIERS
##########################################################################################

### Winsorize methylation beta values to remove potential outliers
### The functions are obtained from the ewaff package and edited to fit our needs (winsorize.pct = 0.005)
### https://github.com/perishky/ewaff/blob/master/R/handle-outliers.r

#######################################################
# 8.1. Extract methylation data from the ExpressionSet

###Extract the methylation data from the ExpressionSet created in STEP 6. 
beta_matrix<-exprs(newmset.EUR)

### Each column is a sample and each row is a probe on the array (450k or EPIC). 

###################################################
#7.2. Winsorize outliers 

### Functions to winsorize data

winsorize <- function(methylation,pct=winsorize.pct) {
  quantiles <- matrixStats::rowQuantiles(methylation, probs=c(pct,1-pct), na.rm=T)
  low <- quantiles[,1]
  upper <- quantiles[,2]
  
  outliers.lower <- rowSums(methylation < low, na.rm=T)
  outliers.upper <- rowSums(methylation > upper, na.rm=T)
  
  idx <- which(methylation < low, arr.ind=T)
  methylation[idx] <- low[idx[,1]]
  
  idx <- which(methylation > upper, arr.ind=T)
  methylation[idx] <- upper[idx[,1]]
  
  n <- rowSums(!is.na(methylation))
  log <- data.frame(outliers.lower, outliers.upper, n)
  
  return(list(methylation=methylation, log=log))
}

### Replace outliers using winsorizing
replace.outliers <- winsorize(beta_matrix, 0.005)
new.methylation <- replace.outliers$methylation
outlier.log <- replace.outliers$log

### Save the log
setwd("/home/isglobal.lan/saguilar/data/WS_HELIX/HELIX_analyses/LCATH_EWAS_GreenSpace_SA/results/HELIX/results.20220616/") #change to your own working directory
write.csv(outlier.log, paste0(cohort,".",ancestry,"_Methylation_Trim_Winsorize0.005_Log.csv"))

### It is important to include the modified beta matrix in the ExpressionSet
exprs(newmset.EUR)<-new.methylation

###################################################
#8.3. Save new ExpressionSet 

### Save this new ExpressionSet with the major ancestry group samples (in this case European ancestry), all the exposure variable, 
#covariables and cell type proportions.
setwd("/home/isglobal.lan/saguilar/data/WS_HELIX/HELIX_analyses/LCATH_EWAS_GreenSpace_SA/db/")#Change to your own working directory
save(newmset.EUR,file=paste0("ExpressionSet_EWAS.GreenSpaces_",cohort,".",ancestry,"_",Date,".Rdata"))

#NOTE: Please, repeat STEP 7 and STEP 8 for all major acestry groups with >100 samples in your cohort. 

#Now you should have the expressionset for each majority ethnic group (>100 samples) created and saved for further analysis (NEXT STEPS).


#################################
#################################
#PART II: Analyses
#################################
#################################

#NOTE: This example code will be run with Europeans as the major ancestry group. You will have to run the code and adapt it for each of 
#your major ancestry groups with >100 individuals.

### In case you have restarted the R Session, please:

#1.Load the ExpressionSet for your major ancestry group. We will start with Europeans. 

setwd("/home/isglobal.lan/saguilar/data/WS_HELIX/HELIX_analyses/LCATH_EWAS_GreenSpace_SA/db/")
#Change to your own working directory where you can find your ExpressionSet with European Samples

load("ExpressionSet_EWAS.GreenSpaces_HELIX.EUR_2022-06-16.Rdata") # Change the file to your own

#2.Set initial parameters again
cohort<-"HELIX" #define cohort name
ancestry<-"EUR" #define major ancestry group you are going to work with
Date<- Sys.Date()#define date 

### Now, choose the working directory where you want to create new directory to save all results for European Samples

setwd("/home/isglobal.lan/saguilar/data/WS_HELIX/HELIX_analyses/LCATH_EWAS_GreenSpace_SA/results/20220616/")
#Change to the directory where you want to save results

### Create a new directory
dir.create("EUR")

### Create an object with the path of this new directory 
res.dir<- "/home/isglobal.lan/saguilar/data/WS_HELIX/HELIX_analyses/LCATH_EWAS_GreenSpace_SA/results/20220616//EUR" 

### Go to the new directory 
setwd(res.dir)

#######################################################
# STEP 9. Check variables: summaries and plots 
#######################################################

### Go to the European working directory created just above.
setwd(res.dir)
 
### Create the directory to save all Descriptive information: 
dir.create("Descriptive.info")

### Go to the Descriptive.info directory 
setwd(paste0(res.dir,"/Descriptive.info"))

### IMPORTANT!! Check the analysis plan to make sure that you are coding the variables correctly 

### Now create a new directory for the Descriptive.plots 
dir.create("Descriptive.Plots")
setwd(paste0(res.dir,"/Descriptive.info/Descriptive.Plots"))

##################
# 9.1) Exposures
#################

#######################################################################
### Residential Sourrounding greenness: NDVI 300  Cumulative (numeric)
class(pData(newmset.EUR)$ndvi300_lf)
pData(newmset.EUR)$ndvi300_lf <-as.numeric(pData(newmset.EUR)$ndvi300_lf)
summary(pData(newmset.EUR)$ndvi300_lf)

### Histogram
jpeg(paste0("Histogram_",cohort,".",ancestry,"_ndvi300_lf_",Date,".jpg"))
hist(pData(newmset.EUR)$ndvi300_lf,main="Histogram ndvi300 cumulative") 
dev.off()

########################################################################
### Residential Sourrounding greenness: NDVI 100 Cumulative (numeric) 
class(pData(newmset.EUR)$ndvi100_lf)#ndvi100_lf
pData(newmset.EUR)$ndvi100_lf <-as.numeric(pData(newmset.EUR)$ndvi100_lf)
summary(pData(newmset.EUR)$ndvi100_lf)

### Histogram
jpeg(paste0("Histogram_",cohort,".",ancestry,"_ndvi100_lf_",Date,".jpg"))
hist(pData(newmset.EUR)$ndvi100_lf,main="Histogram ndvi100 cumulative") 
dev.off()

#########################
# 9.2) Covariates
#########################

###############################################################################
### Maternal education (Categorical). 3 Categories: 1 = high;2 = Medium;3= Low
class(pData(newmset.EUR)$edu_m_0
)
pData(newmset.EUR)$edu_m_0<-as.factor(pData(newmset.EUR)$edu_m_0)
summary(pData(newmset.EUR)$edu_m_0)


####################################################################
### Neighbourhood socio-economic status (Categorical).3 categories:
#1= 1st tertile (low deprivated);2= 2nd tertile (medium deprivated);3= 3rd tertile (high deprivated)
class(pData(newmset.EUR)$areases_tert)
pData(newmset.EUR)$areases_tert<-as.factor(pData(newmset.EUR)$areases_tert)
summary(pData(newmset.EUR)$areases_tert)

################################################
### Maternal age at delivery in years (numeric) 
class(pData(newmset.EUR)$agebirth_m_y)
pData(newmset.EUR)$agebirth_m_y<-as.numeric(pData(newmset.EUR)$agebirth_m_y)
summary(pData(newmset.EUR)$agebirth_m_y)

#####################################################################################################
### Maternal smoking pregnancy (Categorical).ANY smoking. Two categories:0= No; 1= Yes
class(pData(newmset.EUR)$preg_smk)
pData(newmset.EUR)$preg_smk<-as.factor(pData(newmset.EUR)$preg_smk)
summary(pData(newmset.EUR)$preg_smk)

###################################################################################
### Child sex (Categorical). Two categories: 1 = Male; 2 = Female
class(pData(newmset.EUR)$sex_methyl)
pData(newmset.EUR)$sex_methyl<-as.factor(pData(newmset.EUR)$sex_methyl)
summary(pData(newmset.EUR)$sex_methyl)

############################
# Pollution variables
############################

###################
### PM25 (Numeric)
class(pData(newmset.EUR)$pm25_lf)
pData(newmset.EUR)$pm25_lf<-as.numeric(pData(newmset.EUR)$pm25_lf)
summary(pData(newmset.EUR)$pm25_lf)

### Histogram
jpeg(paste0("Histogram_",cohort,".",ancestry,"_PM25_lf_",Date,".jpg"))
hist(pData(newmset.EUR)$pm25_lf,main="PM25 cumulative (air pollution)") 
dev.off()

###########################
# Reproductive variables
###########################

####################################
### Child Z-BMI (numeric)
class(pData(newmset.EUR)$zbmi_methyl)
pData(newmset.EUR)$zbmi_methyl<-as.numeric(pData(newmset.EUR)$zbmi_methyl)
summary(pData(newmset.EUR)$zbmi_methyl)

### Histogram
jpeg(paste0("Histogram_",cohort,".",ancestry,"_zbmi_",Date,".jpg"))
hist(pData(newmset.EUR)$zbmi_methyl,main="Histogram Child Z-BMI") 
dev.off()

############
# Celltypes
############

#HOUSEMAN 
#########
### Bcell_H
class(pData(newmset.EUR)$Bcell_H)
pData(newmset.EUR)$Bcell_H<-as.numeric(pData(newmset.EUR)$Bcell_H)
summary(pData(newmset.EUR)$Bcell_H)

#########
### CD4T_H
class(pData(newmset.EUR)$CD4T_H)
pData(newmset.EUR)$CD4T_H<-as.numeric(pData(newmset.EUR)$CD4T_H)
summary(pData(newmset.EUR)$CD4T_H)

#########
### CD8T_H
class(pData(newmset.EUR)$CD8T_H)
pData(newmset.EUR)$CD8T_H<-as.numeric(pData(newmset.EUR)$CD8T_H)
summary(pData(newmset.EUR)$CD8T_H)

#########
### Gran_H
class(pData(newmset.EUR)$Gran_H)
pData(newmset.EUR)$Gran_H<-as.numeric(pData(newmset.EUR)$Gran_H)
summary(pData(newmset.EUR)$Gran_H)

#########
### Mono_H
class(pData(newmset.EUR)$Mono_H)
pData(newmset.EUR)$Mono_H<-as.numeric(pData(newmset.EUR)$Mono_H)
summary(pData(newmset.EUR)$Mono_H)

#SALAS
#########
### NK_S
class(pData(newmset.EUR)$NK_S)
pData(newmset.EUR)$NK_S<-as.numeric(pData(newmset.EUR)$NK_S)
summary(pData(newmset.EUR)$NK_S)

## Bcell_S
class(pData(newmset.EUR)$Bcell_S)
pData(newmset.EUR)$Bcell_S<-as.numeric(pData(newmset.EUR)$Bcell_S)
summary(pData(newmset.EUR)$Bcell_S)

#########
### CD4T_S
class(pData(newmset.EUR)$CD4T_S)
pData(newmset.EUR)$CD4T_S<-as.numeric(pData(newmset.EUR)$CD4T_S)
summary(pData(newmset.EUR)$CD4T_S)

#########
### CD8T_S
class(pData(newmset.EUR)$CD8T_S)
pData(newmset.EUR)$CD8T_S<-as.numeric(pData(newmset.EUR)$CD8T_S)
summary(pData(newmset.EUR)$CD8T_S)

#########
### Mono_S
class(pData(newmset.EUR)$Mono_S)
pData(newmset.EUR)$Mono_S<-as.numeric(pData(newmset.EUR)$Mono_S)
summary(pData(newmset.EUR)$Mono_S)

### Now we have checked that variables in the ExpressionSet are coded correctly.  

###########################
#STEP 10: ANALYSES: MODEL 1
###########################

#M1:Child blood methylation = green spaces cumulative + maternal age + maternal smoking pregnancy + childs sex + child's age + 
#childs ancestry within major ancestry group (optional) + batch (optional) + cohort (optional) + selection variables (optional)

#NOTE: this model is not adjusted by maternal education, neighbourhood SES, celltypes, PM25 preg and child z-bmi
                                                
#IMPORTANT! Before running the analyses, please make sure that you have done STEP 9) Check variables: summaries and plots. In this way, 
#you make sure that the variables are coded as indicated in the analysis plan.

###########################################
#10.1) Select variables needed for Model 1
###########################################

dd.EUR.M1<-pData(newmset.EUR)[,c("ndvi300_lf","ndvi100_lf", "agebirth_m_y","preg_smk","sex_methyl","age_methyl","gwas_pc1_eur","gwas_pc2_eur","gwas_pc3_eur","gwas_pc4_eur","gwas_pc5_eur","gwas_pc6_eur","gwas_pc7_eur","gwas_pc8_eur","gwas_pc9_eur","gwas_pc10_eur")]


###########################
#10.2) Complete cases M1
###########################

############################################################################
#10.2.1) Create a new dataframe with complete cases from variables in Model 1

dd.EUR.M1.cc<-dd.EUR.M1[complete.cases(dd.EUR.M1),]
dim(dd.EUR.M1.cc)
#[1] 845  16

##################################
#10.2.2) Complete cases M1 samples

samplescmp.M1<-rownames(dd.EUR.M1.cc)

################################################################
#10.2.3) Subset the ExpressionSet with Complete cases M1 samples

newmset.EUR.M1<-newmset.EUR[,(sampleNames(newmset.EUR) %in% samplescmp.M1)]

##########################################################################################
#10.3)Create variables with surrounding residential greenness(NDVI Variable) transformation  
##########################################################################################

### The effect of Residential Surrounding Greenness variable (NDVI) will be reported by IQR, thus we will transform NDVI300 cumulative and 
#NDVI100 cumulative variables.

################################
#10.3.1) NDVI300 IQR cumulative

### Calculate the IQR of NDVI300
ndvi300.iqr<-iqr(dd.EUR.M1.cc$ndvi300_lf,na.rm=TRUE)

### The formula to create this variable corresponds to the original variable value divided by the interquartile value of the variable 
dd.EUR.M1.cc$ndvi300_lf.iqr<-dd.EUR.M1.cc$ndvi300_lf/ndvi300.iqr
summary(dd.EUR.M1.cc$ndvi300_lf.iqr)

###############################
#10.3.2) NDVI100 IQR cumulative

### Calculate the IQR of NDVI100
ndvi100.iqr<-iqr(dd.EUR.M1.cc$ndvi100_lf,na.rm=TRUE)

### The formula to create this variable corresponds to the original variable value divided by the interquartile value of the variable 
dd.EUR.M1.cc$ndvi100_lf.iqr<-dd.EUR.M1.cc$ndvi100_lf/ndvi100.iqr
summary(dd.EUR.M1.cc$ndvi100_lf.iqr)

#######################################################################################################################
#10.4)Include the new dataframe with complete cases M1 and the created NDVI variables transformed in the ExpressionSet
#######################################################################################################################

#Include the new dataframe with the complete case samples for M5 and the created NDVI variables transformed in the ExpressionSet subsetted in section 10.2.3.

###################
#10.4.1) Check order 

# Before including the new dataframe with the complete case samples for M1 and the created NDVI variables transformed in the ExpressionSet,
#check if the samples are in the same order. If you do not order samples as they are in the ExpressionSet, you could incorrectly assign the values 
#of the variables to the samples and therefore also to the methylation. 

table(ifelse(rownames(dd.EUR.M1.cc)==sampleNames(newmset.EUR.M1),"Matched","--NOT MATCHED--"))
#Matched
#    845

#In our case, samples are ordered in the same way. In case that your samples are not ordered, please, see STEP 6 for an example of how to do this. 

#############################################################
#10.4.2) Include the new dataframe in the ExpressionSet

pData(newmset.EUR.M1)<-dd.EUR.M1.cc

####################################
#13.5) Descriptive analyses for M1
####################################

#10.5.1) Descriptive tables
#10.5.2) Correlation between exposures
#10.5.3) Associations between green spaces vars and covariates

### Go to the Descriptive.info directory 
setwd(paste0(res.dir,"/Descriptive.info/"))

###  Create the Descriptive.Tables directory
dir.create("Descriptive.Tables")
### Create the Correlations directory
dir.create("Correlations")
### Create the GS vs Covariates directory
dir.create("GSvsCovariates")

###########################
#10.5.1) Descriptive tables

### Go to the Descriptive.Tables directory

setwd(paste0(res.dir,"/Descriptive.info/Descriptive.Tables"))

#a) Mean (SD)
##############

descriptive.vars.M1<-colnames(dd.EUR.M1.cc)#Select variables for the descriptives

dd.EUR.M1.cc %>%
  dplyr::mutate(
    ndvi300_lf = ff_label(ndvi300_lf, "NDVI300 cumulative"),ndvi300_lf.iqr = ff_label(ndvi300_lf.iqr, "NDVI300 cumulative IQR"),
    ndvi100_lf = ff_label(ndvi100_lf, "NDVI100 cumulative"),ndvi100_lf.iqr = ff_label(ndvi100_lf.iqr, "NDVI100 cumulative IQR"),
    agebirth_m_y = ff_label(agebirth_m_y, "Maternal Age at delivery"),preg_smk = ff_label(preg_smk, "Maternal smoking preg"), 
    sex_methyl = ff_label(sex_methyl, "Child sex"),age_methyl = ff_label(age_methyl, "Child age")
    )%>%
  summary_factorlist(dependent=NULL,descriptive.vars.M1,cont="mean",column=TRUE,na_include=TRUE,total_col=TRUE)->DescriptiveTable.Mean.M1

DescriptiveTable.Mean.M1<-DescriptiveTable.Mean.M1[,c("label","levels","Total")]
colnames(DescriptiveTable.Mean.M1)<-c("Variables","levels","Mean (SD) or n(%)")

#Median (IQR)
#############

dd.EUR.M1.cc %>%
  dplyr::mutate(
    ndvi300_lf = ff_label(ndvi300_lf, "NDVI300 cumulative"),ndvi300_lf.iqr = ff_label(ndvi300_lf.iqr, "NDVI300 cumulative IQR"),
    ndvi100_lf = ff_label(ndvi100_lf, "NDVI100 cumulative"),ndvi100_lf.iqr = ff_label(ndvi100_lf.iqr, "NDVI100 cumulative IQR"),
    agebirth_m_y = ff_label(agebirth_m_y, "Maternal Age at delivery"),preg_smk = ff_label(preg_smk, "Maternal smoking preg"),
    sex_methyl = ff_label(sex_methyl, "Child sex"),age_methyl = ff_label(age_methyl, "Child age")
    )%>%
  summary_factorlist(dependent=NULL,descriptive.vars.M1,cont="median",column=TRUE,na_include=TRUE,total_col=TRUE)->DescriptiveTable.Median.M1

DescriptiveTable.Median.M1<-DescriptiveTable.Median.M1[,c("label","levels","Total")]
colnames(DescriptiveTable.Median.M1)<-c("Variables","levels","Mean (SD) or n(%)")

### Create a single descriptive table with all variables

DescriptiveTable.final.M1<-cbind(DescriptiveTable.Mean.M1,DescriptiveTable.Median.M1)
DescriptiveTable.final.M1<-DescriptiveTable.final.M1[,c(1:3,5:6)]

write.table(
  DescriptiveTable.final.M1, file=paste0("M1.Descriptive.table.",cohort,".",ancestry,".cc.",Date,".txt"),col.names=T, row.names=F, quote=F, sep="\t")

######################################
#10.5.2) Correlation between exposures

### Go to the Correlations directory
setwd(paste0(res.dir,"/Descriptive.info/Correlations"))

### a)Correlation between NDVI300 and NDVI100 (Continuous vs continuous)
#########################################################################

cor.ndvi300.100<-dd.EUR.M1.cc[,c("ndvi300_lf","ndvi100_lf")]
#Create plot
plot.cor.ndvi300.100<-ggpairs(cor.ndvi300.100, title="Correlation between NDVI300 and NDVI100") 

#Define a path 
file.create(paste0("M1.CorrelationNDVI300vsNDVI100.",cohort,".",ancestry,"_",Date,".pdf"))

path1.M1<-file.path(paste0("M1.CorrelationNDVI300vsNDVI100.",cohort,".",ancestry,"_",Date,".pdf"))

pdf(path1.M1)
print(plot.cor.ndvi300.100)
dev.off()

###############################################################
#10.5.3) Associations between Green Spaces vars and Covariates

### Go to the GSvsCovariates directory
setwd(paste0(res.dir,"/Descriptive.info/GSvsCovariates"))

############################################################################
##a) Association beweeen Continuous Green spaces variables  and Covariables

GSexposures<-c("ndvi300_lf.iqr","ndvi100_lf.iqr")

Covariates<-c("agebirth_m_y","preg_smk","sex_methyl","age_methyl","gwas_pc1_eur","gwas_pc2_eur","gwas_pc3_eur","gwas_pc4_eur","gwas_pc5_eur",
              "gwas_pc6_eur","gwas_pc7_eur","gwas_pc8_eur","gwas_pc9_eur","gwas_pc10_eur")#Specify variables M1

coefs.M1=as.data.frame(matrix(NA,nrow=1,ncol=10))
colnames(coefs.M1)=c("DateTime","GreenSpaceVar","Covariates","Model","BetaNoStand", "SE", "P","R2","R2adj","N")

count=1
for(j in 1:length(GSexposures)) {
  print(GSexposures[j])
  results=NULL
  for(i in 1:length(Covariates)) {
    print(Covariates[i])
    ff <- paste0(GSexposures[j],"~",Covariates[i],sep="") 
    fit <- lm(ff,data=dd.EUR.M1.cc)
    s <-summary(fit)
    res <- s$coefficients
    R2<-s$r.squared
    R2.adj<- s$adj.r.squared
    
    #Save analysis details and coefficients
    coefs.M1[count,1]=gsub(" ","_",Sys.time())
    coefs.M1[count,2]=GSexposures[j]
    coefs.M1[count,3]=Covariates[i]
    coefs.M1[count,4]=ff
    coefs.M1[count,5:7]=res[2,c(1,2,4)]
    coefs.M1[count,8]= R2
    coefs.M1[count,9]=R2.adj
    coefs.M1[count,10]= s$df[2]
    
    count=count+1
    
  }
}

#Multiple testing correction by FDR 

### ndvi300 lf IQR vs Covariates
coefs.M1.subset.ndvi300<-subset(coefs.M1, coefs.M1$GreenSpaceVar == "ndvi300_lf.iqr")
coefs.M1.subset.ndvi300$Padj<-  p.adjust(coefs.M1.subset.ndvi300$P, method="fdr")

### ndvi100 lf IQR vs Covariates
coefs.M1.subset.ndvi100<-subset(coefs.M1, coefs.M1$GreenSpaceVar == "ndvi100_lf.iqr")
coefs.M1.subset.ndvi100$Padj<-  p.adjust(coefs.M1.subset.ndvi100$P, method="fdr")

### Create one table with all results
coefs.M1.all<-rbind(coefs.M1.subset.ndvi300,coefs.M1.subset.ndvi100)

write.csv(coefs.M1.all,paste0("M1.GScontVars.Covariates_",cohort,".",ancestry,"_",Date,".csv"), row.names = TRUE)


##############################################################################
#10.6) EWAS GREEN SPACES cumulative USING ROBUST LINEAR REGRESSIONS WITH LIMMA
##############################################################################

### Before starting with the EWAS analysis...

### Create a folder for the EWAS results 

setwd(res.dir)
dir.create("EWAS.Results")

### Go to the EWAS results directory
setwd(paste0(res.dir,"/EWAS.Results"))

#For the analysis use ExpressionSet subset created in section 10.1.3. 

#newmset.EUR.M1

# 10.6.1) Exposure NDVI300 cumulative
# 10.6.2) Exposure NDVI100 cumulative 

###################################
# 10.6.1) Exposure NDVI300 cumulative

################
# lf_N300_M1
################

# a) Run EWAS

### Design the model
design.lf_N300_M1<-model.matrix(~ndvi300_lf.iqr + agebirth_m_y + preg_smk + sex_methyl + age_methyl+ gwas_pc1_eur + gwas_pc2_eur + gwas_pc3_eur +
                                  gwas_pc4_eur + gwas_pc5_eur + gwas_pc6_eur + gwas_pc7_eur + gwas_pc8_eur + gwas_pc9_eur + 
                                  gwas_pc10_eur,data=pData(newmset.EUR.M1))

### Run Robust linear regression model with limma
fit.lf_N300_M1<-limma::lmFit(newmset.EUR.M1,design=design.lf_N300_M1, method="robust")
fit.lf_N300_M1<-limma::eBayes(fit.lf_N300_M1)
lf_N300_M1<-limma::topTable(fit.lf_N300_M1, coef =2, number = Inf,sort.by="none",confint=TRUE)
lf_N300_M1$SE <- (sqrt(fit.lf_N300_M1$s2.post) *fit.lf_N300_M1$stdev.unscaled)[, 2]
head(lf_N300_M1)

###  Export results

write.table(lf_N300_M1, paste0("EWAS_GS_",cohort,".",ancestry,"_lf_N300_M1_",Date,".txt"), na="NA") 
            
# b) Calculate lambda
lambda.lf_N300_M1<- qchisq(median(lf_N300_M1$P.Value,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
lambda.lf_N300_M1

# c) QQ plot 

#Create QQplot directory 
setwd(paste0(res.dir,"/EWAS.Results"))
dir.create("QQplots")

### Go to QQplots directory
setwd(paste0(res.dir,"/EWAS.Results/QQPlots"))

pvals<-lf_N300_M1$P.Value
jpeg(paste0("QQPlot_",cohort,".",ancestry,"_lf_N300_M1_",Date,".jpg"))
qq(pvals,main=paste0("QQPlot_",cohort,".",ancestry,"_lf_N300_M1")) 
dev.off()

###################################
# 10.6.2) Exposure NDVI100 cumulative
            
################
# lf_N100_M1
################


# a) Run EWAS

### Design the model
design.lf_N100_M1<-model.matrix(~ndvi100_lf.iqr + agebirth_m_y + preg_smk + sex_methyl + age_methyl + gwas_pc1_eur + gwas_pc2_eur +
                                  gwas_pc3_eur + gwas_pc4_eur + gwas_pc5_eur + gwas_pc6_eur + gwas_pc7_eur + gwas_pc8_eur + gwas_pc9_eur +
                                  gwas_pc10_eur, data=pData(newmset.EUR.M1))

### Run Robust linear regression model with limma
fit.lf_N100_M1<-limma::lmFit(newmset.EUR.M1,design=design.lf_N100_M1, method="robust")
fit.lf_N100_M1<-limma::eBayes(fit.lf_N100_M1)
lf_N100_M1<-limma::topTable(fit.lf_N100_M1, coef =2, number = Inf,sort.by="none",confint=TRUE)
lf_N100_M1$SE <- (sqrt(fit.lf_N100_M1$s2.post) *fit.lf_N100_M1$stdev.unscaled)[, 2]
head(lf_N100_M1)

###  Export results
setwd(paste0(res.dir,"/EWAS.Results/"))

write.table(lf_N100_M1, paste0("EWAS_GS_",cohort,".",ancestry,"_lf_N100_M1_",Date,".txt"), na="NA") 
            
# b) Calculate lambda
lambda.lf_N100_M1<- qchisq(median(lf_N100_M1$P.Value,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
lambda.lf_N100_M1

# c) QQ plot 

### Go to QQplots directory
setwd(paste0(res.dir,"/EWAS.Results/QQPlots"))

pvals<-lf_N100_M1$P.Value
jpeg(paste0("QQPlot_",cohort,".",ancestry,"_lf_N100_M1_",Date,".jpg"))
qq(pvals,main=paste0("QQPlot_",cohort,".",ancestry,"_lf_N100_M1")) 
dev.off()

              
####################
#10.7) LAMBDAS M1
####################

### Create lambdas table M1

### Create Lambdas directory
setwd(paste0(res.dir,"/EWAS.Results"))
dir.create("Lambdas")

setwd(paste0(res.dir,"/EWAS.Results/Lambdas"))
                                                          
lambdas.table<-rbind(lambda.lf_N300_M1,lambda.lf_N100_M1)
colnames(lambdas.table)<-"Lambdas M1"

write.table(lambdas.table, paste0("Lambdas_",cohort,".",ancestry,"_M1_",Date,".txt"), na="NA")

###########################
#STEP 11: ANALYSES: MODEL 2
###########################

#M2:Child blood methylation = green spaces cumulative + maternal age + maternal smoking pregnancy + childs sex + child's age + 
#childs ancestry within major ancestry group (optional) + batch (optional) + cohort (optional) + selection variables (optional) + maternal education +
#neighbourhood socio-economic status 

#NOTE: this model is not adjusted by celltypes, PM25 preg and child z-bmi
                                                
#IMPORTANT! Before running the analyses, please make sure that you have done STEP 9) Check variables: summaries and plots. In this way, you make sure that the variables are coded as indicated in the analysis plan.

###########################################
#11.1) Select variables needed for Model 2
###########################################

dd.EUR.M2<-pData(newmset.EUR)[,c("ndvi300_lf","ndvi100_lf","agebirth_m_y","preg_smk","sex_methyl","age_methyl","gwas_pc1_eur","gwas_pc2_eur",
                                 "gwas_pc3_eur","gwas_pc4_eur","gwas_pc5_eur","gwas_pc6_eur","gwas_pc7_eur","gwas_pc8_eur","gwas_pc9_eur","gwas_pc10_eur",
                                 "edu_m_0","areases_tert")]


###########################
#11.2) Complete cases M2
###########################

############################################################################
#11.2.1) Create a new dataframe with complete cases from variables in Model 2

dd.EUR.M2.cc<-dd.EUR.M2[complete.cases(dd.EUR.M2),]

##################################
#11.2.2) Complete cases M2 samples

samplescmp.M2<-rownames(dd.EUR.M2.cc)

################################################################
#11.2.3) Subset the ExpressionSet with Complete cases M2 samples

newmset.EUR.M2<-newmset.EUR[,(sampleNames(newmset.EUR) %in% samplescmp.M2)]

##########################################################################################
#11.3)Create variables with surrounding residential greenness(NDVI Variable) transformation  
##########################################################################################

### The effect of Residential Sorrounding Greenness variable (NDVI) will be reported by IQR, thus we will transform NDVI300 cumulative and
#NDVI100 cumulative variables.

################################
#11.3.1) NDVI300 IQR cumulative

### Calculate the IQR of NDVI300
ndvi300.iqr<-iqr(dd.EUR.M2.cc$ndvi300_lf,na.rm=TRUE)

### The formula to create this variable corresponds to the original variable value divided by the interquartile value of the variable 
dd.EUR.M2.cc$ndvi300_lf.iqr<-dd.EUR.M2.cc$ndvi300_lf/ndvi300.iqr
summary(dd.EUR.M2.cc$ndvi300_lf.iqr)

###############################
#11.3.2) NDVI100 IQR cumulative

### Calculate the IQR of NDVI100
ndvi100.iqr<-iqr(dd.EUR.M2.cc$ndvi100_lf,na.rm=TRUE)

### The formula to create this variable corresponds to the original variable value divided by the interquartile value of the variable 
dd.EUR.M2.cc$ndvi100_lf.iqr<-dd.EUR.M2.cc$ndvi100_lf/ndvi100.iqr
summary(dd.EUR.M2.cc$ndvi100_lf.iqr)

#######################################################################################################################
#11.4)Include the new dataframe with complete cases M2 and the created NDVI variables transformed in the ExpressionSet
#######################################################################################################################

#Include the new dataframe with the complete case samples for M4 and the created NDVI variables transformed in the ExpressionSet subsetted in section 11.2.3.

###################
#11.4.1) Check order 

# Before including the new dataframe with the complete case samples for M2 and the created NDVI variables transformed in the ExpressionSet, 
#check if the samples are in the same order. If you do not order samples as they are in the ExpressionSet, you could incorrectly assign the values 
#of the variables to the samples and therefore also to the methylation. 

table(ifelse(rownames(dd.EUR.M2.cc)==sampleNames(newmset.EUR.M2),"Matched","--NOT MATCHED--"))
#Matched
#    691

#In our case, samples are ordered in the same way. In case that your samples are not ordered, please, see STEP 6 for an example of how to do this. 

#############################################################
#11.4.2) Include the new dataframe in the ExpressionSet

pData(newmset.EUR.M2)<-dd.EUR.M2.cc

####################################
#11.5) Descriptive analyses for M2
###################################


#NOTE: We will not run descriptive analyses for M2 as we understand that the N will be the same as in M3 as model 3 is model 2 + 
#celltypes proportions which should not have NAs. 

##############################################################################
#11.6) EWAS GREEN SPACES cumulative USING ROBUST LINEAR REGRESSIONS WITH LIMMA
##############################################################################

### Before starting with the EWAS analysis...

### Go to the EWAS results directory
setwd(paste0(res.dir,"/EWAS.Results"))

#For the analysis use ExpressionSet subset created in section 11.1.3. 

#newmset.EUR.M2

# 11.6.1) Exposure NDVI300 cumulative
# 11.6.2) Exposure NDVI100 cumulative 

###################################
# 11.6.1) Exposure NDVI300 cumulative

################
# lf_N300_M2
################
dd.EUR.M2<-pData(newmset.EUR)[,c("ndvi300_lf","ndvi100_lf","agebirth_m_y","preg_smk","sex_methyl","age_methyl","gwas_pc1_eur","gwas_pc2_eur",
                                 "gwas_pc3_eur","gwas_pc4_eur","gwas_pc5_eur","gwas_pc6_eur","gwas_pc7_eur","gwas_pc8_eur","gwas_pc9_eur","gwas_pc10_eur",
                                 "edu_m_0","areases_tert")]
# a) Run EWAS

### Design the model
design.lf_N300_M2<-model.matrix(~ndvi300_lf.iqr + agebirth_m_y + preg_smk + sex_methyl + age_methyl + gwas_pc1_eur + gwas_pc2_eur + gwas_pc3_eur +
                                  gwas_pc4_eur + gwas_pc5_eur + gwas_pc6_eur + gwas_pc7_eur + gwas_pc8_eur + gwas_pc9_eur + gwas_pc10_eur +
                                  edu_m_0 + areases_tert, data=pData(newmset.EUR.M2))

### Run Robust linear regression model with limma
fit.lf_N300_M2<-limma::lmFit(newmset.EUR.M2,design=design.lf_N300_M2, method="robust")
fit.lf_N300_M2<-limma::eBayes(fit.lf_N300_M2)
lf_N300_M2<-limma::topTable(fit.lf_N300_M2, coef =2, number = Inf,sort.by="none",confint=TRUE)
lf_N300_M2$SE <- (sqrt(fit.lf_N300_M2$s2.post) *fit.lf_N300_M2$stdev.unscaled)[, 2]
head(lf_N300_M2)

###  Export results

write.table(lf_N300_M2, paste0("EWAS_GS_",cohort,".",ancestry,"_lf_N300_M2_",Date,".txt"), na="NA") 
            
# b) Calculate lambda
lambda.lf_N300_M2<- qchisq(median(lf_N300_M2$P.Value,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
lambda.lf_N300_M2

# c) QQ plot 

### Go to QQplots directory
setwd(paste0(res.dir,"/EWAS.Results/QQPlots"))

pvals<-lf_N300_M2$P.Value
jpeg(paste0("QQPlot_",cohort,".",ancestry,"_lf_N300_M2_",Date,".jpg")) #type error: was data. Changed to Date
qq(pvals,main=paste0("QQPlot_",cohort,".",ancestry,"_lf_N300_M2")) 
dev.off()

###################################
# 11.6.2) Exposure NDVI100 cumulative
            
################
# lf_N100_M2
################


# a) Run EWAS

### Design the model
design.lf_N100_M2<-model.matrix(~ndvi100_lf.iqr + agebirth_m_y + preg_smk + sex_methyl + age_methyl + gwas_pc1_eur + gwas_pc2_eur + gwas_pc3_eur +
                                  gwas_pc4_eur + gwas_pc5_eur + gwas_pc6_eur + gwas_pc7_eur + gwas_pc8_eur + gwas_pc9_eur + gwas_pc10_eur + 
                                  edu_m_0 + areases_tert, data=pData(newmset.EUR.M2))

### Run Robust linear regression model with limma
fit.lf_N100_M2<-limma::lmFit(newmset.EUR.M2,design=design.lf_N100_M2, method="robust")
fit.lf_N100_M2<-limma::eBayes(fit.lf_N100_M2)
lf_N100_M2<-limma::topTable(fit.lf_N100_M2, coef =2, number = Inf,sort.by="none",confint=TRUE)
lf_N100_M2$SE <- (sqrt(fit.lf_N100_M2$s2.post) *fit.lf_N100_M2$stdev.unscaled)[, 2]
head(lf_N100_M2)

###  Export results
setwd(paste0(res.dir,"/EWAS.Results"))
write.table(lf_N100_M2, paste0("EWAS_GS_",cohort,".",ancestry,"_lf_N100_M2",Date,".txt"), na="NA") 
            
# b) Calculate lambda
lambda.lf_N100_M2<- qchisq(median(lf_N100_M2$P.Value,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
lambda.lf_N100_M2

# c) QQ plot 

### Go to QQplots directory
setwd(paste0(res.dir,"/EWAS.Results/QQPlots"))

pvals<-lf_N100_M2$P.Value
jpeg(paste0("QQPlot_",cohort,".",ancestry,"_lf_N100_M2_",Date,".jpg"))
qq(pvals,main=paste0("QQPlot_",cohort,".",ancestry,"_lf_N100_M2")) 
dev.off()
 
####################
#11.7) LAMBDAS M2
####################

### Create lambdas table M2

setwd(paste0(res.dir,"/EWAS.Results/Lambdas"))
                                                          
lambdas.table<-rbind(lambda.lf_N300_M2,lambda.lf_N100_M2)
colnames(lambdas.table)<-"Lambdas M2"

write.table(lambdas.table, paste0("Lambdas_",cohort,".",ancestry,"_M2",Date,".txt"), na="NA")


###########################
#STEP 12: ANALYSES: MODEL 3
###########################

#M3:Child blood methylation = green spaces cumulative + maternal age + maternal smoking pregnancy + childs sex + child's age + 
#childs ancestry within major ancestry group (optional) + batch (optional) + cohort (optional) + selection variables (optional)  +
#maternal education + neighbourhood socio-economic status + blood cellular composition

#NOTE: this model is not adjusted by PM25 preg and child z-bmi
                                                
#IMPORTANT! Before running the analyses, please make sure that you have done STEP 9) Check variables: summaries and plots. In this way, you make sure that the variables are coded as indicated in the analysis plan.

###########################################
#12.1) Select variables needed for Model 3
###########################################

dd.EUR.M3<-pData(newmset.EUR)[,c("ndvi300_lf","ndvi100_lf", "agebirth_m_y","preg_smk","sex_methyl","age_methyl","gwas_pc1_eur","gwas_pc2_eur",
                                 "gwas_pc3_eur","gwas_pc4_eur","gwas_pc5_eur","gwas_pc6_eur","gwas_pc7_eur","gwas_pc8_eur","gwas_pc9_eur","gwas_pc10_eur",
                                 "edu_m_0","areases_tert","CD8T_H","CD4T_H","NK_H","Bcell_H","Mono_H","Gran_H")]

#añadir la cohorte seria importante!!
###########################
#12.2) Complete cases M3
###########################

############################################################################
#12.2.1) Create a new dataframe with complete cases from variables in Model 3

dd.EUR.M3.cc<-dd.EUR.M3[complete.cases(dd.EUR.M3),]

##################################
#12.2.2) Complete cases M3 samples

samplescmp.M3<-rownames(dd.EUR.M3.cc)

################################################################
#12.2.3) Subset the ExpressionSet with Complete cases M3 samples

newmset.EUR.M3<-newmset.EUR[,(sampleNames(newmset.EUR) %in% samplescmp.M3)]

##########################################################################################
#12.3)Create variables with sorrounding residential greenness(NDVI Variable) transformation  
##########################################################################################

### The effect of Residential Sorrounding Greenness variable (NDVI) will be reported by IQR, thus we will transform NDVI300 cumulative and
#NDVI100 cumulative variables.

################################
#12.3.1) NDVI300 IQR cumulative

### Calculate the IQR of NDVI300
ndvi300.iqr<-iqr(dd.EUR.M3.cc$ndvi300_lf,na.rm=TRUE)

### The formula to create this variable corresponds to the original variable value divided by the interquartile value of the variable 
dd.EUR.M3.cc$ndvi300_lf.iqr<-dd.EUR.M3.cc$ndvi300_lf/ndvi300.iqr
summary(dd.EUR.M3.cc$ndvi300_lf.iqr)

###############################
#12.3.2) NDVI100 IQR cumulative

### Calculate the IQR of NDVI100
ndvi100.iqr<-iqr(dd.EUR.M3.cc$ndvi100_lf,na.rm=TRUE)

### The formula to create this variable corresponds to the original variable value divided by the interquartile value of the variable 
dd.EUR.M3.cc$ndvi100_lf.iqr<-dd.EUR.M3.cc$ndvi100_lf/ndvi100.iqr
summary(dd.EUR.M3.cc$ndvi100_lf.iqr)

#######################################################################################################################
#12.4)Include the new dataframe with complete cases M3 and the created NDVI variables transformed in the ExpressionSet
#######################################################################################################################

#Include the new dataframe with the complete case samples for M3 and the created NDVI variables transformed in the ExpressionSet subsetted in section 12.2.3.

###################
#12.4.1) Check order 

# Before including the new dataframe with the complete case samples for M3 and the created NDVI variables transformed in the ExpressionSet,
#check if the samples are in the same order. If you do not order samples as they are in the ExpressionSet, you could incorrectly assign the values of 
#the variables to the samples and therefore also to the methylation. 

table(ifelse(rownames(dd.EUR.M3.cc)==sampleNames(newmset.EUR.M3),"Matched","--NOT MATCHED--"))
#Matched
#    727

#In our case, samples are ordered in the same way. In case that your samples are not ordered, please, see STEP 6 for an example of how to do this. 

#############################################################
#12.4.2) Include the new dataframe in the ExpressionSet

pData(newmset.EUR.M3)<-dd.EUR.M3.cc

####################################
#12.5) Descriptive analyses for M3
####################################

#12.5.1) Descriptive tables
#12.5.2) Correlation between exposures
#12.5.3) Associations between green spaces vars and covariates


###########################
#12.5.1) Descriptive tables

### Go to the Descriptive.Tables directory 
setwd(paste0(res.dir,"/Descriptive.info/Descriptive.Tables"))

#a) Mean (SD)
##############

descriptive.vars.M3<-colnames(dd.EUR.M3.cc)#Select variables for the descriptives

dd.EUR.M3.cc %>%
  dplyr::mutate(
    ndvi300_lf = ff_label(ndvi300_lf, "NDVI300 cumulative"),ndvi300_lf.iqr = ff_label(ndvi300_lf.iqr, "NDVI300 cumulative IQR"),
    ndvi100_lf = ff_label(ndvi100_lf, "NDVI100 cumulative"),ndvi100_lf.iqr = ff_label(ndvi100_lf.iqr, "NDVI100 cumulative IQR"),
    agebirth_m_y = ff_label(agebirth_m_y, "Maternal Age at delivery"),edu_m_0 = ff_label(edu_m_0, "Maternal education"),
    areases_tert = ff_label(areases_tert, "Neighbourhood socio-economic status"),preg_smk = ff_label(preg_smk, "Maternal smoking preg"),
    sex_methyl = ff_label(sex_methyl, "Child sex"),age_methyl = ff_label(age_methyl, "Child age")
    )%>%
  summary_factorlist(dependent=NULL,descriptive.vars.M3,cont="mean",column=TRUE,na_include=TRUE,total_col=TRUE)->DescriptiveTable.Mean.M3

DescriptiveTable.Mean.M3<-DescriptiveTable.Mean.M3[,c("label","levels","Total")]
colnames(DescriptiveTable.Mean.M3)<-c("Variables","levels","Mean (SD) or n(%)")

#Median (IQR)
#############

dd.EUR.M3.cc %>%
  dplyr::mutate(
    ndvi300_lf = ff_label(ndvi300_lf, "NDVI300 cumulative"),ndvi300_lf.iqr = ff_label(ndvi300_lf.iqr, "NDVI300 cumulative IQR"),
    ndvi100_lf = ff_label(ndvi100_lf, "NDVI100 cumulative"),ndvi100_lf.iqr = ff_label(ndvi100_lf.iqr, "NDVI100 cumulative IQR"),
    agebirth_m_y = ff_label(agebirth_m_y, "Maternal Age at delivery"),edu_m_0 = ff_label(edu_m_0, "Maternal education"),
    areases_tert = ff_label(areases_tert, "Neighbourhood socio-economic status"),preg_smk = ff_label(preg_smk, "Maternal smoking preg"),
    sex_methyl = ff_label(sex_methyl, "Child sex"),age_methyl = ff_label(age_methyl, "Child age")
    )%>%
  summary_factorlist(dependent=NULL,descriptive.vars.M3,cont="median",column=TRUE,na_include=TRUE,total_col=TRUE)->DescriptiveTable.Median.M3

DescriptiveTable.Median.M3<-DescriptiveTable.Median.M3[,c("label","levels","Total")]
colnames(DescriptiveTable.Median.M3)<-c("Variables","levels","Mean (SD) or n(%)")

### Create a single descriptive table with all variables

DescriptiveTable.final.M3<-cbind(DescriptiveTable.Mean.M3,DescriptiveTable.Median.M3)
DescriptiveTable.final.M3<-DescriptiveTable.final.M3[,c(1:3,5:6)]
             
write.table(
  DescriptiveTable.final.M3, file=paste0("M3.Descriptive.table.",cohort,".",ancestry,".cc.",Date,".txt"),col.names=T, row.names=F, quote=F, sep="\t")

######################################
#12.5.2) Correlation between exposures

### Go to the Correlations directory 
setwd(paste0(res.dir,"/Descriptive.info/Correlations"))

### a)Correlation between NDVI300 and NDVI100 (Continuous vs continuous)
#########################################################################

cor.ndvi300.100<-dd.EUR.M3.cc[,c("ndvi300_lf","ndvi100_lf")]
#Create plot
plot.cor.ndvi300.100<-ggpairs(cor.ndvi300.100, title="Correlation between NDVI300 and NDVI100") 

#Define a path 
file.create(paste0("M3.CorrelationNDVI300vsNDVI100.",cohort,".",ancestry,"_",Date,".pdf"))

path1.M3<-file.path(paste0("M3.CorrelationNDVI300vsNDVI100.",cohort,".",ancestry,"_",Date,".pdf"))

pdf(path1.M3)
print(plot.cor.ndvi300.100)
dev.off()

###############################################################
#12.5.3) Associations between Green Spaces vars and Covariates

### Go to the GSvsCovariates directory
setwd(paste0(res.dir,"/Descriptive.info/GSvsCovariates"))

############################################################################
##a) Association beweeen Continuous Green spaces variables  and Covariables

GSexposures<-c("ndvi300_lf.iqr","ndvi100_lf.iqr")

Covariates<-c("edu_m_0","areases_tert", "agebirth_m_y","preg_smk","sex_methyl","age_methyl","CD8T_H","CD4T_H","NK_H","Bcell_H","Mono_H","Gran_H",
              "gwas_pc1_eur","gwas_pc2_eur","gwas_pc3_eur","gwas_pc4_eur","gwas_pc5_eur","gwas_pc6_eur","gwas_pc7_eur","gwas_pc8_eur","gwas_pc9_eur",
              "gwas_pc10_eur")#Specify variables M3

coefs.M3=as.data.frame(matrix(NA,nrow=1,ncol=10))
colnames(coefs.M3)=c("DateTime","GreenSpaceVar","Covariates","Model","BetaNoStand", "SE", "P","R2","R2adj","N")

count=1
for(j in 1:length(GSexposures)) {
  print(GSexposures[j])
  results=NULL
  for(i in 1:length(Covariates)) {
    print(Covariates[i])
    ff <- paste0(GSexposures[j],"~",Covariates[i],sep="") 
    fit <- lm(ff,data=dd.EUR.M3.cc)
    s <-summary(fit)
    res <- s$coefficients
    R2<-s$r.squared
    R2.adj<- s$adj.r.squared
    
    #Save analysis details and coefficients
    coefs.M3[count,1]=gsub(" ","_",Sys.time())
    coefs.M3[count,2]=GSexposures[j]
    coefs.M3[count,3]=Covariates[i]
    coefs.M3[count,4]=ff
    coefs.M3[count,5:7]=res[2,c(1,2,4)]
    coefs.M3[count,8]= R2
    coefs.M3[count,9]=R2.adj
    coefs.M3[count,10]= s$df[2]
    
    count=count+1
    
  }
}

#Multiple testing correction by FDR 

### ndvi300 lf IQR vs Covariates
coefs.M3.subset.ndvi300<-subset(coefs.M3, coefs.M3$GreenSpaceVar == "ndvi300_lf.iqr")
coefs.M3.subset.ndvi300$Padj<-  p.adjust(coefs.M3.subset.ndvi300$P, method="fdr")

### ndvi100 lf IQR vs Covariates
coefs.M3.subset.ndvi100<-subset(coefs.M3, coefs.M3$GreenSpaceVar == "ndvi100_lf.iqr")
coefs.M3.subset.ndvi100$Padj<-  p.adjust(coefs.M3.subset.ndvi100$P, method="fdr")

### Create one table with all results
coefs.M3.all<-rbind(coefs.M3.subset.ndvi300,coefs.M3.subset.ndvi100)

write.csv(coefs.M3.all,paste0("M3.GScontVars.Covariates_",cohort,".",ancestry,"_",Date,".csv"), row.names = TRUE)


##############################################################################
#12.6) EWAS GREEN SPACES cumulative USING ROBUST LINEAR REGRESSIONS WITH LIMMA
##############################################################################

### Before starting with the EWAS analysis...

### Go to the EWAS results directory
setwd(paste0(res.dir,"/EWAS.Results"))

#For the analysis use ExpressionSet subset created in section 12.1.3. 

#newmset.EUR.M3

# 12.6.1) Exposure NDVI300 cumulative
# 12.6.2) Exposure NDVI100 cumulative 

###################################
# 12.6.1) Exposure NDVI300 cumulative

################
# lf_N300_M3
################

# a) Run EWAS

### Design the model
design.lf_N300_M3<-model.matrix(~ndvi300_lf.iqr + agebirth_m_y + preg_smk + sex_methyl + age_methyl + gwas_pc1_eur + gwas_pc2_eur + gwas_pc3_eur +
                                  gwas_pc4_eur + gwas_pc5_eur + gwas_pc6_eur + gwas_pc7_eur + gwas_pc8_eur + gwas_pc9_eur + gwas_pc10_eur + edu_m_0+
                                  areases_tert + CD8T_H + CD4T_H + NK_H + Bcell_H + Mono_H + Gran_H, data=pData(newmset.EUR.M3))

### Run Robust linear regression model with limma
fit.lf_N300_M3<-limma::lmFit(newmset.EUR.M3,design=design.lf_N300_M3, method="robust")
fit.lf_N300_M3<-limma::eBayes(fit.lf_N300_M3)
lf_N300_M3<-limma::topTable(fit.lf_N300_M3, coef =2, number = Inf,sort.by="none",confint=TRUE)
lf_N300_M3$SE <- (sqrt(fit.lf_N300_M3$s2.post) *fit.lf_N300_M3$stdev.unscaled)[, 2]
head(lf_N300_M3)

###  Export results

write.table(lf_N300_M3, paste0("EWAS_GS_",cohort,".",ancestry,"_lf_N300_M3_",Date,".txt"), na="NA") 
            
# b) Calculate lambda
lambda.lf_N300_M3<- qchisq(median(lf_N300_M3$P.Value,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
lambda.lf_N300_M3

# c) QQ plot 

### Go to QQplots directory
setwd(paste0(res.dir,"/EWAS.Results/QQPlots"))

pvals<-lf_N300_M3$P.Value
jpeg(paste0("QQPlot_",cohort,".",ancestry,"_lf_N300_M3_",Date,".jpg"))
qq(pvals,main=paste0("QQPlot_",cohort,".",ancestry,"_lf_N300_M3")) 
dev.off()

###################################
# 12.6.2) Exposure NDVI100 cumulative
            
################
# lf_N100_M3
################


# a) Run EWAS

### Design the model
design.lf_N100_M3<-model.matrix(~ndvi100_lf.iqr + agebirth_m_y + preg_smk + sex_methyl + age_methyl + gwas_pc1_eur + gwas_pc2_eur + gwas_pc3_eur +
                                  gwas_pc4_eur + gwas_pc5_eur + gwas_pc6_eur + gwas_pc7_eur + gwas_pc8_eur + gwas_pc9_eur + gwas_pc10_eur + edu_m_0 +
                                  areases_tert + CD8T_H + CD4T_H + NK_H + Bcell_H + Mono_H + Gran_H  , data=pData(newmset.EUR.M3))

### Run Robust linear regression model with limma
fit.lf_N100_M3<-limma::lmFit(newmset.EUR.M3,design=design.lf_N100_M3, method="robust")
fit.lf_N100_M3<-limma::eBayes(fit.lf_N100_M3)
lf_N100_M3<-limma::topTable(fit.lf_N100_M3, coef =2, number = Inf,sort.by="none",confint=TRUE)
lf_N100_M3$SE <- (sqrt(fit.lf_N100_M3$s2.post) *fit.lf_N100_M3$stdev.unscaled)[, 2]
head(lf_N100_M3)

###  Export results
setwd(paste0(res.dir,"/EWAS.Results"))
write.table(lf_N100_M3, paste0("EWAS_GS_",cohort,".",ancestry,"_lf_N100_M3",Date,".txt"), na="NA") 
            
# b) Calculate lambda
lambda.lf_N100_M3<- qchisq(median(lf_N100_M3$P.Value,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
lambda.lf_N100_M3

# c) QQ plot 

### Go to QQplots directory
setwd(paste0(res.dir,"/EWAS.Results/QQPlots"))

pvals<-lf_N100_M3$P.Value
jpeg(paste0("QQPlot_",cohort,".",ancestry,"_lf_N100_M3_",Date,".jpg"))
qq(pvals,main=paste0("QQPlot_",cohort,".",ancestry,"_lf_N100_M3")) 
dev.off()

              
####################
#12.7) LAMBDAS M3
####################

### Create lambdas table M3

setwd(paste0(res.dir,"/EWAS.Results/Lambdas"))
                                                          
lambdas.table<-rbind(lambda.lf_N300_M3,lambda.lf_N100_M3)
colnames(lambdas.table)<-"Lambdas M3"

write.table(lambdas.table, paste0("Lambdas_",cohort,".",ancestry,"_M3_",Date,".txt"), na="NA")


###########################
#STEP 13: ANALYSES: MODEL 4
###########################

#M4:Child blood methylation = green spaces cumulative + maternal age + maternal smoking pregnancy + childs sex + child's age +
#childs ancestry within major ancestry group (optional) + batch (optional) + cohort (optional) + selection variables (optional) + 
#maternal education + neighbourhood socio-economic status + blood cellular composition + PM2.5 cumulative 
                                                
#NOTE: this model is not adjusted by child z-bmi

#IMPORTANT! Before running the analyses, please make sure that you have done STEP 9) Check variables: summaries and plots. 
#In this way, you make sure that the variables are coded as indicated in the analysis plan.

###########################################
#13.1) Select variables needed for Model 4
###########################################

dd.EUR.M4<-pData(newmset.EUR)[,c("ndvi300_lf","ndvi100_lf", "agebirth_m_y","preg_smk","sex_methyl","age_methyl","gwas_pc1_eur","gwas_pc2_eur",
                                 "gwas_pc3_eur","gwas_pc4_eur","gwas_pc5_eur","gwas_pc6_eur","gwas_pc7_eur","gwas_pc8_eur","gwas_pc9_eur","gwas_pc10_eur",
                                 "edu_m_0" ,"areases_tert","CD8T_H","CD4T_H","NK_H","Bcell_H","Mono_H","Gran_H","pm25_lf")]


###########################
#13.2) Complete cases M4
###########################

############################################################################
#13.2.1) Create a new dataframe with complete cases from variables in Model 4

dd.EUR.M4.cc<-dd.EUR.M4[complete.cases(dd.EUR.M4),]

##################################
#13.2.2) Complete cases M4 samples

samplescmp.M4<-rownames(dd.EUR.M4.cc)

################################################################
#13.2.3) Subset the ExpressionSet with Complete cases M4 samples

newmset.EUR.M4<-newmset.EUR[,(sampleNames(newmset.EUR) %in% samplescmp.M4)]

##########################################################################################
#13.3)Create variables with sorrounding residential greenness(NDVI Variable) transformation  
##########################################################################################

### The effect of Residential Sorrounding Greenness variable (NDVI) will be reported by IQR, thus we will transform NDVI300 cumulative and
#NDVI100 cumulative variables.

################################
#13.3.1) NDVI300 IQR cumulative

### Calculate the IQR of NDVI300
ndvi300.iqr<-iqr(dd.EUR.M4.cc$ndvi300_lf,na.rm=TRUE)

### The formula to create this variable corresponds to the original variable value divided by the interquartile value of the variable 
dd.EUR.M4.cc$ndvi300_lf.iqr<-dd.EUR.M4.cc$ndvi300_lf/ndvi300.iqr
summary(dd.EUR.M4.cc$ndvi300_lf.iqr)

###############################
#13.3.2) NDVI100 IQR cumulative

### Calculate the IQR of NDVI100
ndvi100.iqr<-iqr(dd.EUR.M4.cc$ndvi100_lf,na.rm=TRUE)

### The formula to create this variable corresponds to the original variable value divided by the interquartile value of the variable 
dd.EUR.M4.cc$ndvi100_lf.iqr<-dd.EUR.M4.cc$ndvi100_lf/ndvi100.iqr
summary(dd.EUR.M4.cc$ndvi100_lf.iqr)

#######################################################################################################################
#13.4)Include the new dataframe with complete cases M4 and the created NDVI variables transformed in the ExpressionSet
#######################################################################################################################

#Include the new dataframe with the complete case samples for M4 and the created NDVI variables transformed in the ExpressionSet subsetted in section 13.2.3.

###################
#13.4.1) Check order 

# Before including the new dataframe with the complete case samples for M4 and the created NDVI variables transformed in the ExpressionSet,
#check if the samples are in the same order. If you do not order samples as they are in the ExpressionSet, you could incorrectly assign the values of
#the variables to the samples and therefore also to the methylation. 

table(ifelse(rownames(dd.EUR.M4.cc)==sampleNames(newmset.EUR.M4),"Matched","--NOT MATCHED--"))
#Matched
#    707

#In our case, samples are ordered in the same way. In case that your samples are not ordered, please, see STEP 6 for an example of how to do this. 

#############################################################
#13.4.2) Include the new dataframe in the ExpressionSet

pData(newmset.EUR.M4)<-dd.EUR.M4.cc

####################################
#13.5) Descriptive analyses for M4
####################################

#13.5.1) Descriptive tables
#13.5.2) Correlation between exposures
#13.5.3) Associations between green spaces vars and covariates


###########################
#13.5.1) Descriptive tables

### Go to the Descriptive.Tables directory 
setwd(paste0(res.dir,"/Descriptive.info/Descriptive.Tables"))

#a) Mean (SD)
##############

descriptive.vars.M4<-colnames(dd.EUR.M4.cc)#Select variables for the descriptives

dd.EUR.M4.cc %>%
  dplyr::mutate(
    ndvi300_lf = ff_label(ndvi300_lf, "NDVI300 cumulative"),ndvi300_lf.iqr = ff_label(ndvi300_lf.iqr, "NDVI300 cumulative IQR"),
    ndvi100_lf = ff_label(ndvi100_lf, "NDVI100 cumulative"),ndvi100_lf.iqr = ff_label(ndvi100_lf.iqr, "NDVI100 cumulative IQR"),
    agebirth_m_y = ff_label(agebirth_m_y, "Maternal Age at delivery"),edu_m_0 = ff_label(edu_m_0, "Maternal education"),
    areases_tert = ff_label(areases_tert, "Neighbourhood socio-economic status"),preg_smk= ff_label(preg_smk, "Maternal smoking preg"),
    sex_methyl = ff_label(sex_methyl, "Child sex"),age_methyl = ff_label(age_methyl, "Child age"),pm25_lf = ff_label(pm25_lf, "PM2.5 levels cumulative")
    )%>%
  summary_factorlist(dependent=NULL,descriptive.vars.M4,cont="mean",column=TRUE,na_include=TRUE,total_col=TRUE)->DescriptiveTable.Mean.M4

DescriptiveTable.Mean.M4<-DescriptiveTable.Mean.M4[,c("label","levels","Total")]
colnames(DescriptiveTable.Mean.M4)<-c("Variables","levels","Mean (SD) or n(%)")

#Median (IQR)
#############

dd.EUR.M4.cc %>%
  dplyr::mutate(
    ndvi300_lf = ff_label(ndvi300_lf, "NDVI300 cumulative"),ndvi300_lf.iqr = ff_label(ndvi300_lf.iqr, "NDVI300 cumulative IQR"),
    ndvi100_lf = ff_label(ndvi100_lf, "NDVI100 cumulative"),ndvi100_lf.iqr = ff_label(ndvi100_lf.iqr, "NDVI100 cumulative IQR"),
    agebirth_m_y = ff_label(agebirth_m_y, "Maternal Age at delivery"),edu_m_0 = ff_label(edu_m_0, "Maternal education"),
    areases_tert = ff_label(areases_tert, "Neighbourhood socio-economic status"),preg_smk= ff_label(preg_smk, "Maternal smoking preg"),
    sex_methyl = ff_label(sex_methyl, "Child sex"),age_methyl = ff_label(age_methyl, "Child age"),pm25_lf = ff_label(pm25_lf, "PM2.5 levels cumulative")
    )%>%
  summary_factorlist(dependent=NULL,descriptive.vars.M4,cont="median",column=TRUE,na_include=TRUE,total_col=TRUE)->DescriptiveTable.Median.M4

DescriptiveTable.Median.M4<-DescriptiveTable.Median.M4[,c("label","levels","Total")]
colnames(DescriptiveTable.Median.M4)<-c("Variables","levels","Mean (SD) or n(%)")

### Create a single descriptive table with all variables

DescriptiveTable.final.M4<-cbind(DescriptiveTable.Mean.M4,DescriptiveTable.Median.M4)
DescriptiveTable.final.M4<-DescriptiveTable.final.M4[,c(1:3,5:6)]
             
write.table(
  DescriptiveTable.final.M4, file=paste0("M4.Descriptive.table.",cohort,".",ancestry,".cc.",Date,".txt"),col.names=T, row.names=F, quote=F, sep="\t")

######################################
#13.5.2) Correlation between exposures

### Go to the Correlations directory 
setwd(paste0(res.dir,"/Descriptive.info/Correlations"))

### a)Correlation between NDVI300 and NDVI100 (Continuous vs continuous)
#########################################################################

cor.ndvi300.100<-dd.EUR.M4.cc[,c("ndvi300_lf","ndvi100_lf")]
#Create plot
plot.cor.ndvi300.100<-ggpairs(cor.ndvi300.100, title="Correlation between NDVI300 and NDVI100") 

#Define a path 
file.create(paste0("M4.CorrelationNDVI300vsNDVI100.",cohort,".",ancestry,"_",Date,".pdf"))

path1.M4<-file.path(paste0("M4.CorrelationNDVI300vsNDVI100.",cohort,".",ancestry,"_",Date,".pdf"))

pdf(path1.M4)
print(plot.cor.ndvi300.100)
dev.off()

###############################################################
#13.5.3) Associations between Green Spaces vars and Covariates

### Go to the GSvsCovariates directory
setwd(paste0(res.dir,"/Descriptive.info/GSvsCovariates"))

############################################################################
##a) Association beweeen Continuous Green spaces variables  and Covariables

GSexposures<-c("ndvi300_lf.iqr","ndvi100_lf.iqr")

Covariates<-c("edu_m_0","areases_tert", "agebirth_m_y","preg_smk","sex_methyl","age_methyl","pm25_lf","CD8T_H","CD4T_H","NK_H","Bcell_H","Mono_H",
              "Gran_H","gwas_pc1_eur","gwas_pc2_eur","gwas_pc3_eur","gwas_pc4_eur","gwas_pc5_eur","gwas_pc6_eur","gwas_pc7_eur","gwas_pc8_eur","gwas_pc9_eur",
              "gwas_pc10_eur")#Specify variables M4

coefs.M4=as.data.frame(matrix(NA,nrow=1,ncol=10))
colnames(coefs.M4)=c("DateTime","GreenSpaceVar","Covariates","Model","BetaNoStand", "SE", "P","R2","R2adj","N")

count=1
for(j in 1:length(GSexposures)) {
  print(GSexposures[j])
  results=NULL
  for(i in 1:length(Covariates)) {
    print(Covariates[i])
    ff <- paste0(GSexposures[j],"~",Covariates[i],sep="") 
    fit <- lm(ff,data=dd.EUR.M4.cc)
    s <-summary(fit)
    res <- s$coefficients
    R2<-s$r.squared
    R2.adj<- s$adj.r.squared
    
    #Save analysis details and coefficients
    coefs.M4[count,1]=gsub(" ","_",Sys.time())
    coefs.M4[count,2]=GSexposures[j]
    coefs.M4[count,3]=Covariates[i]
    coefs.M4[count,4]=ff
    coefs.M4[count,5:7]=res[2,c(1,2,4)]
    coefs.M4[count,8]= R2
    coefs.M4[count,9]=R2.adj
    coefs.M4[count,10]= s$df[2]
    
    count=count+1
    
  }
}

#Multiple testing correction by FDR 

### ndvi300 lf IQR vs Covariates
coefs.M4.subset.ndvi300<-subset(coefs.M4, coefs.M4$GreenSpaceVar == "ndvi300_lf.iqr")
coefs.M4.subset.ndvi300$Padj<-  p.adjust(coefs.M4.subset.ndvi300$P, method="fdr")

### ndvi100 lf IQR vs Covariates
coefs.M4.subset.ndvi100<-subset(coefs.M4, coefs.M4$GreenSpaceVar == "ndvi100_lf.iqr")
coefs.M4.subset.ndvi100$Padj<-  p.adjust(coefs.M4.subset.ndvi100$P, method="fdr")

### Create one table with all results
coefs.M4.all<-rbind(coefs.M4.subset.ndvi300,coefs.M4.subset.ndvi100)

write.csv(coefs.M4.all,paste0("M4.GScontVars.Covariates_",cohort,".",ancestry,"_",Date,".csv"), row.names = TRUE)


##############################################################################
#13.6) EWAS GREEN SPACES CUMULATIVE USING ROBUST LINEAR REGRESSIONS WITH LIMMA
##############################################################################

### Before starting with the EWAS analysis...

### Go to the EWAS results directory
setwd(paste0(res.dir,"/EWAS.Results"))

#For the analysis use ExpressionSet subset created in section 13.1.3. 

#newmset.EUR.M4

# 13.6.1) Exposure NDVI300 cumulative
# 13.6.2) Exposure NDVI100 cumulative 

###################################
# 13.6.1) Exposure NDVI300 cumulative

################
# lf_N300_M4
################

# a) Run EWAS

### Design the model
design.lf_N300_M4<-model.matrix(~ndvi300_lf.iqr + agebirth_m_y + preg_smk + sex_methyl + age_methyl + gwas_pc1_eur + gwas_pc2_eur + gwas_pc3_eur +
                                  gwas_pc4_eur + gwas_pc5_eur + gwas_pc6_eur + gwas_pc7_eur + gwas_pc8_eur + gwas_pc9_eur + gwas_pc10_eur + edu_m_0+ 
                                  areases_tert + CD8T_H + CD4T_H + NK_H + Bcell_H + Mono_H + Gran_H + pm25_lf, data=pData(newmset.EUR.M4))

### Run Robust linear regression model with limma
fit.lf_N300_M4<-limma::lmFit(newmset.EUR.M4,design=design.lf_N300_M4, method="robust")
fit.lf_N300_M4<-limma::eBayes(fit.lf_N300_M4)
lf_N300_M4<-limma::topTable(fit.lf_N300_M4, coef =2, number = Inf,sort.by="none",confint=TRUE)
lf_N300_M4$SE <- (sqrt(fit.lf_N300_M4$s2.post) *fit.lf_N300_M4$stdev.unscaled)[, 2]
head(lf_N300_M4)

###  Export results
setwd(paste0(res.dir,"/EWAS.Results"))
write.table(lf_N300_M4, paste0("EWAS_GS_",cohort,".",ancestry,"_lf_N300_M4_",Date,".txt"), na="NA") 
            
# b) Calculate lambda
lambda.lf_N300_M4<- qchisq(median(lf_N300_M4$P.Value,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
lambda.lf_N300_M4

# c) QQ plot 

### Go to QQplots directory
setwd(paste0(res.dir,"/EWAS.Results/QQPlots"))

pvals<-lf_N300_M4$P.Value
jpeg(paste0("QQPlot_",cohort,".",ancestry,"_lf_N300_M4_",Date,".jpg"))
qq(pvals,main=paste0("QQPlot_",cohort,".",ancestry,"_lf_N300_M4")) 
dev.off()

###################################
# 13.6.2) Exposure NDVI100 cumulative
            
################
# lf_N100_M4
################


# a) Run EWAS

### Design the model
design.lf_N100_M4<-model.matrix(~ndvi100_lf.iqr +  agebirth_m_y + preg_smk  + sex_methyl + age_methyl+ gwas_pc1_eur + gwas_pc2_eur + gwas_pc3_eur +
                                  gwas_pc4_eur + gwas_pc5_eur + gwas_pc6_eur + gwas_pc7_eur + gwas_pc8_eur + gwas_pc9_eur + gwas_pc10_eur + edu_m_0 +
                                  areases_tert + CD8T_H + CD4T_H + NK_H + Bcell_H + Mono_H + Gran_H +  pm25_lf, data=pData(newmset.EUR.M4))

### Run Robust linear regression model with limma
fit.lf_N100_M4<-limma::lmFit(newmset.EUR.M4,design=design.lf_N100_M4, method="robust")
fit.lf_N100_M4<-limma::eBayes(fit.lf_N100_M4)
lf_N100_M4<-limma::topTable(fit.lf_N100_M4, coef =2, number = Inf,sort.by="none",confint=TRUE)
lf_N100_M4$SE <- (sqrt(fit.lf_N100_M4$s2.post) *fit.lf_N100_M4$stdev.unscaled)[, 2]
head(lf_N100_M4)

###  Export results
setwd(paste0(res.dir,"/EWAS.Results"))
write.table(lf_N100_M4, paste0("EWAS_GS_",cohort,".",ancestry,"_lf_N100_M4_",Date,".txt"), na="NA") 
            
# b) Calculate lambda
lambda.lf_N100_M4<- qchisq(median(lf_N100_M4$P.Value,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
lambda.lf_N100_M4

# c) QQ plot 

### Go to QQplots directory
setwd(paste0(res.dir,"/EWAS.Results/QQPlots"))

pvals<-lf_N100_M4$P.Value
jpeg(paste0("QQPlot_",cohort,".",ancestry,"_lf_N100_M4_",Date,".jpg"))
qq(pvals,main=paste0("QQPlot_",cohort,".",ancestry,"_lf_N100_M4")) 
dev.off()

              
####################
#13.7) LAMBDAS M4
####################

### Create lambdas table M4

setwd(paste0(res.dir,"/EWAS.Results/Lambdas"))
                                                          
lambdas.table<-rbind(lambda.lf_N300_M4,lambda.lf_N100_M4)
colnames(lambdas.table)<-"Lambdas M4"

write.table(lambdas.table, paste0("Lambdas_",cohort,".",ancestry,"_M4_",Date,".txt"), na="NA")

#################
#IMPORTANT
#################

#Repeat from STEP 9 to STEP 13 for each of the major ancestry groups in your cohort with more than >100 samples.

############################
#THANK YOU!
############################
