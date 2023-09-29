###########################################################################################
# EWAS Green Spaces and Cord Blood- No DataSHIELD code (using a ExpressionSet)  
# Sofía Aguilar Lacasaña (sofia.aguilar@isglobal.org)
# 14/06/2022 
#V2 
###########################################################################################

############################################################################################################################################################
# The following R code will allow you to complete all the EWAS requested in the LifeCycle-LongITools-ATHLETE-EWAS of green spaces and Cord blood DNA 
#methylation analysis plan.
# The code also produces files summarising the variables included in the EWAS.
# The code can be adapted to the needs of each cohort.
# There are two inputs required for this analysis:

# 1) A matrix with the methylation beta values (0-1)
# Each column is a sample and each row is a probe on the array (450K or EPIC)

# 2) A dataframe containing all the exposures + covariates + cell type proportions
# Each row is a sample(individual) and each column is a different variable. 
# Details on which variables should be included in this dataframe and how to code them are provided in the analysis plan.  Throughout the code the variables
# are named as indicated in the analysis plan. We recommend that you rename your variables as indicated. If not, you will need to manually modify the variable names in the code. 

### This code gets the methylation data and the metadata from an ExpressionSet.  If you have the matrix with the methylation data and the metadata in 
#different files, please create an ExpressionSet. In case you have any questions on how to do this, please contact us 
# (sofia.aguilar@isglobal.org; mariona.bustamante@isglobal.org).

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
setwd("/home/isglobal.lan/saguilar/data/WS_INMA/Methylation_INMA/PACE/green_spaces_SA/db/") 

###########################################
### Set initial parameters
cohort<-"INMA" #define cohort name
Date<- Sys.Date()#define date 

#################################
#################################
#PART I: Prepare ExpressionSets
#################################
#################################

##########################################################################################
### STEP 1: Load required packages 
##########################################################################################
##(if these are not already installed, you will have to install them as the first step) --> install.packages("data.table") 

library(data.table) # to process 
library(Biobase)# to be able to access and modify data in the ExpressionSet
library(EnvStats)#calculate the iqr
library(finalfit)#descriptive tables
library(GGally)#correlation plots
library(stats) #glm
library(limma) #EWAS 
library(qqman) # to get QQ plots

##########################################################################################
### STEP 2: Load the methylation data and the metadata
##########################################################################################

### Change to your own working directory where you can find the ExpressionSet
setwd("/home/isglobal.lan/saguilar/data/WS_INMA/Methylation_INMA/450K_blood/DS_ExpressionSet/Results/Results_22122021/")

### Load the ExpressionSet object
load("INMA_Methyl_Child_Blood_0y_450K_20211222.Rdata")# Change the name for your ExpressionSet object

ls()#to see the R object loaded

### ExpressionSet summarized view
CreatedExpressionSet
# ExpressionSet (storageMode: lockedEnvironment)
# assayData: 476946 features, 361 samples
# element names: exprs
# protocolData: none
# phenoData
# sampleNames: 3998888002_R01C02 3998888002_R02C01 ...
#  7927554155_R03C01 (361 total)
# varLabels: id_methyl id ... nRBC (15 total)
# varMetadata: labelDescription
# featureData: none
# experimentData: use 'experimentData(object)'
# Annotation: IlluminaHumanMethylation450kanno.ilmn12.hg19

###########################
### 2.1) Methylation data
###########################

### To extract the beta matrix from the ExpressionSet you should use exprs() function as you will see explained in section 7. 
# Each column is a sample and each row is a probe on the array (450k or EPIC). 

#beta_matrix<-exprs(CreatedExpressionSet) 

##################
### 2.2) Metadata 
##################

### Extract the metadata from the ExpressionSet 
# which includes some important covariates + cell type proportions

metadataEset<-pData(CreatedExpressionSet)
colnames(metadataEset)
# [1] "id_methyl"     "id"            "sex_methyl"    "age_methyl"
# [5] "anc_methyl"    "cohort_methyl" "bw_methyl"     "ga_methyl"
# [9] "Bcell"         "CD4T"          "CD8T"          "Gran"
# [13] "Mono"          "NK"            "nRBC"


##########################################################################################
### STEP 3: Load the dataframe with the exposures + additional covariates 
##########################################################################################

#### Read dataframe 
# Each row is a sample(individual) and each column is a different variable.
# Ensure all traits and covariates have been derived as specified in the analysis plan and that exclusions have been made.

setwd("/home/isglobal.lan/saguilar/data/WS_INMA/Methylation_INMA/PACE/green_spaces_SA/db/") # Change to your own working directory
df.exp.cov<-read.table(file="INMA.exposures_covariables_20210823.csv",sep=",",header=TRUE) # Change the file name to your own
dim(df.exp.cov)
colnames(df.exp.cov)
#[1] "child_id"          "ndvi300_preg"      "ndvi100_preg"
#[4] "greenyn300_preg"   "agebirth_m_y"      "preg_smk"
#[7] "edu_m_0"            "areases_tert_preg" "pm25_preg"

###############################################################################################################
### STEP 4: Merge the metadata from the ExpressionSet and the dataframe with exposures + additional variables
###############################################################################################################

### Merge both dataframes by the ID (it might be or not the same ID as the methylation dataset). 
dd<-merge(metadataEset,df.exp.cov,by.x="id",by.y="child_id",all.x=TRUE)
dim(dd)

# Now we have all variables needed for the analysis for samples with methylation data
colnames(dd)
# [1] "id"                "id_methyl"         "sex_methyl"
# [4] "age_methyl"        "anc_methyl"        "cohort_methyl"
# [7] "bw_methyl"         "ga_methyl"         "Bcell"
#[10] "CD4T"              "CD8T"              "Gran"
#[13] "Mono"              "NK"                "nRBC"
#[16] "ndvi300_preg"      "ndvi100_preg"      "greenyn300_preg"
#[19] "agebirth_m_y"      "preg_smk"          "edu_m_0"
#[22] "areases_tert_preg" "pm25_preg"         


#########################################################################################################
### STEP 5: Include updated metadata (exposures+covariates+ cell type proportions) in the ExpressionSet
#########################################################################################################

###########
#IMPORTANT! 
###########

##################################################################
#5.1) Check if the samples in the merged dataframe are in the same order than in the metadata from the ExpressionSet

# If you do not order samples as they are in the ExpressionSet, you could incorrectly assign the values of the 
#variables to the samples and therefore also to the methylation. 

table(ifelse(metadataEset$id==dd$id,"Matched","--NOT MATCHED--"))
#--NOT MATCHED--         Matched
#            382               3
#

### Samples are not in the same order. We need to order in the same way

####################
#5.2) Order samples
dd.ord<-dd[order(match(dd$id,metadataEset$id)),]
table(ifelse(metadataEset$id==dd.ord$id,"Matched","--NOT MATCHED--"))
#Matched
#    385

# Now they are in the same order

### We need the id_methyl in rows to link this new pheno+metadata with the methylation data from the ExpressionSet
rownames(dd.ord)<-dd.ord$id_methyl

##############################################################
#5.3) Include the new ordered phenodata in the ExpressionSet 
pData(CreatedExpressionSet)<-dd.ord

##########################################
### STEP 6: Major ancestry group subset
##########################################

### Define group
#As ethnic groups are tested separatelly (if >100 samples), please indicate the major ancestry group you are working with

ancestry<-"EUR" 
##########################################################################
#6.1) Subset the ExpressionSet by the major ancestry group (>100 samples)

#NOTE: If you have closed the R session, re-extract the ordered phenodata from the expressionSet. If not, continue with dataframe created in STEP 5. 
#dd.ord<-pData(CreatedExpressionSet)

table(dd.ord$anc_methyl)
# 1   6
#361  22

### I will subset children of European ancestry (group 1)
dd.EUR<-dd.ord[dd.ord$anc_methyl == 1,]
dim(dd.EUR)
#################################################
#6.2) Remove samples with NA in ancestry 
table(is.na(dd.EUR$anc_methyl))
# FALSE  TRUE
#  361    2
# There are two NAS. 
dd.EUR<-dd.EUR[(!is.na(dd.EUR$anc_methyl)),]
dim(dd.EUR)
#[1] 361  81
table(is.na(dd.EUR$anc_methyl))
#FALSE
#  361

####################################
#6.4. Create the new ExpressionSet

### Now we create the ExpressionSet for European ancestry samples, without NAs in ancestry and without NAs in green spaces vars 

#European samples
samplesEUR<-rownames(dd.EUR)

#Subset the ExpressionSet with european samples
newmset.EUR<-CreatedExpressionSet[,(sampleNames(CreatedExpressionSet) %in% samplesEUR)]
newmset.EUR
# ExpressionSet (storageMode: lockedEnvironment)
# assayData: 476946 features, 312 samples
# element names: exprs
# protocolData: none
# phenoData
# sampleNames: 3998888002_R01C02 3998888002_R02C01 ...
#  7927554154_R06C02 (312 total)
# varLabels: id id_methyl ... no2_preg (81 total)
# varMetadata: labelDescription
# featureData: none
# experimentData: use 'experimentData(object)'
# Annotation: IlluminaHumanMethylation450kanno.ilmn12.hg19

##########################################################################################
### STEP 7: WINSORIZE OUTLIERS
##########################################################################################

### Winsorize methylation beta values to remove potential outliers
### The functions are obtained from the ewaff package and edited to fit our needs (winsorize.pct = 0.005)
### https://github.com/perishky/ewaff/blob/master/R/handle-outliers.r

#######################################################
# 7.1. Extract methylation data from the ExpressionSet

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
setwd("/home/isglobal.lan/saguilar/data/WS_INMA/Methylation_INMA/PACE/green_spaces_SA/results/Results.20220103/") #change to your own working directory
write.csv(outlier.log, paste0(cohort,".",ancestry,"_Methylation_Trim_Winsorize0.005_Log.csv"))

### It is important to include the modified beta matrix in the ExpressionSet
exprs(newmset.EUR)<-new.methylation

###################################################
#7.3. Save new ExpressionSet 

### Save this new ExpressionSet with the major ancestry group samples (in this case European ancestry), all the exposure variable,
# covariables and cell type proportions.
setwd("/home/isglobal.lan/saguilar/data/WS_INMA/Methylation_INMA/PACE/green_spaces_SA/db/")#Change to your own working directory
save(newmset.EUR,file=paste0("ExpressionSet_EWAS.GreenSpaces_",cohort,".",ancestry,"_",Date,".Rdata"))

#NOTE: Please, repeat STEP 6 and STEP 7 for all major acestry groups with >100 samples in your cohort. 

#Now you should have the expressionset for each majority ethnic group (>100 samples) created and saved for further analysis (NEXT STEPS)


#####################################################################################################################################################


#################################
#################################
#PART II: Analyses
#################################
#################################

#NOTE: This example code will be run with Europeans as the major ancestry group. You will have to run the code and adapt it for each of your major ancestry
# groups with >100 individuals.

### In case you restart the R Session, please:

#1.Load the ExpressionSet for your major ancestry group. We will start with Europeans. 

setwd("/home/isglobal.lan/saguilar/data/WS_INMA/Methylation_INMA/PACE/green_spaces_SA/db/")#Change to your own working directory where you can find your 
# ExpressionSet with European Samples
load("ExpressionSet_EWAS.GreenSpaces_INMA.EUR_2021-12-22.Rdata") # Change the file name to your own

#2.Set initial parameters again
cohort<-"INMA" #define cohort name
ancestry<-"EUR" #define major ancestry group you are going to work with
Date<- Sys.Date()#define date 

### Now, choose the working directory where you want to create the new directory to save all results for European Samples

setwd("/home/isglobal.lan/saguilar/data/WS_INMA/Methylation_INMA/PACE/green_spaces_SA/results/INMA/Results.20220614/")
#Change to the directory where you want to save results

### Create a new directory
dir.create("EUR")

### Create an object with the path of this new directory
res.dir<- "/home/isglobal.lan/saguilar/data/WS_INMA/Methylation_INMA/PACE/green_spaces_SA/results/INMA/Results.20220614/EUR/" 

### Go to the new directory
setwd(res.dir)

#######################################################
# STEP 8. Check variables: summaries and plots 
#######################################################

### Make sure you are in the European working directory created just above. 
setwd(res.dir)
### Create the directory to save all Descriptive information: 

dir.create("Descriptive.info")

### Go to the new directory
setwd(paste0(res.dir,"/Descriptive.info"))

### IMPORTANT!! Check the analysis plan to make sure that you are coding the variables correctly 

### Now create a new directory for the Descriptive.plots 
dir.create("Descriptive.Plots")
setwd(paste0(res.dir,"/Descriptive.info/Descriptive.Plots"))

##################
# 8.1) Exposures
#################

#######################################################################
### Residential Sourrounding greenness: NDVI 300  pregnancy (numeric)
class(pData(newmset.EUR)$ndvi300_preg)
pData(newmset.EUR)$ndvi300_preg<-as.numeric(pData(newmset.EUR)$ndvi300_preg)
summary(pData(newmset.EUR)$ndvi300_preg)

### Histogram
jpeg(paste0("Histogram_",cohort,".",ancestry,"_ndvi300_preg_",Date,".jpg"))
hist(pData(newmset.EUR)$ndvi300_preg,main="Histogram ndvi300 preg") 
dev.off()

########################################################################
### Residential Sourrounding greenness: NDVI 100 pregnancy (numeric) 
class(pData(newmset.EUR)$ndvi100_preg)
pData(newmset.EUR)$ndvi100_preg<-as.numeric(pData(newmset.EUR)$ndvi100_preg)
summary(pData(newmset.EUR)$ndvi100_preg)

### Histogram
jpeg(paste0("Histogram_",cohort,".",ancestry,"_ndvi100_preg_",Date,".jpg"))
hist(pData(newmset.EUR)$ndvi100_preg,main="Histogram ndvi100 preg") 
dev.off()

###############################################################
### Green access pregnancy(Greenyn 300 pregnancy) (Categorical)
# Two categories:
# 0=no
# 1=yes
class(pData(newmset.EUR)$greenyn300_preg)
pData(newmset.EUR)$greenyn300_preg<-as.factor(pData(newmset.EUR)$greenyn300_preg)
summary(pData(newmset.EUR)$greenyn300_preg)

### BarPlot
jpeg(paste0("BarPlot_",cohort,".",ancestry,"_Greenyn300_preg_",Date,".jpg"))
plot(pData(newmset.EUR)$greenyn300_preg,main="BarPlot Greenyn300 preg")
dev.off()

#########################
# 8.2) Covariates
#########################

###############################################################################
### Maternal education (Categorical). 3 Categories: 1 = high;2 = Medium;3= Low
class(pData(newmset.EUR)$edu_m_0)
pData(newmset.EUR)$edu_m_0<-as.factor(pData(newmset.EUR)$edu_m_0)
summary(pData(newmset.EUR)$edu_m_0)

####################################################################
### Neighbourhood socio-economic status (Categorical).3 categories:
#1= 1st tertile (low deprivated);2= 2nd tertile (medium deprivated);3= 3rd tertile (high deprivated)
class(pData(newmset.EUR)$areases_tert_preg)
pData(newmset.EUR)$areases_tert_preg<-as.factor(pData(newmset.EUR)$areases_tert_preg)
summary(pData(newmset.EUR)$areases_tert_preg)

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
class(pData(newmset.EUR)$pm25_preg)
pData(newmset.EUR)$pm25_preg<-as.numeric(pData(newmset.EUR)$pm25_preg)
summary(pData(newmset.EUR)$pm25_preg)

### Histogram
jpeg(paste0("Histogram_",cohort,".",ancestry,"_PM25_preg_",Date,".jpg"))
hist(pData(newmset.EUR)$pm25_preg,main="PM25 pregnancy (air pollution)") 
dev.off()

###########################
# Reproductive variables
###########################

####################################
### Gestational age in days (numeric)
class(pData(newmset.EUR)$ga_methyl)
pData(newmset.EUR)$ga_methyl<-as.numeric(pData(newmset.EUR)$ga_methyl)
summary(pData(newmset.EUR)$ga_methyl)

### Histogram
jpeg(paste0("Histogram_",cohort,".",ancestry,"_GestationalAge_",Date,".jpg"))
hist(pData(newmset.EUR)$ga_methyl,main="Histogram Gestational age") 
dev.off()

############################################
### Birth weight at birth in grams(numeric)
class(pData(newmset.EUR)$bw_methyl)
pData(newmset.EUR)$bw_methyl<-as.numeric(pData(newmset.EUR)$bw_methyl)
summary(pData(newmset.EUR)$bw_methyl)

### Histogram
jpeg(paste0("Histogram_",cohort,".",ancestry,"_BirthWeight_",Date,".jpg"))
hist(pData(newmset.EUR)$bw_methyl,main="Histogram BirthWeight") 
dev.off()

############
# Celltypes
############

#########
### Bcell
class(pData(newmset.EUR)$Bcell)
pData(newmset.EUR)$Bcell<-as.numeric(pData(newmset.EUR)$Bcell)
summary(pData(newmset.EUR)$Bcell)

#########
### CD4T
class(pData(newmset.EUR)$CD4T)
pData(newmset.EUR)$CD4T<-as.numeric(pData(newmset.EUR)$CD4T)
summary(pData(newmset.EUR)$CD4T)

#########
### CD8T
class(pData(newmset.EUR)$CD8T)
pData(newmset.EUR)$CD8T<-as.numeric(pData(newmset.EUR)$CD8T)
summary(pData(newmset.EUR)$CD8T)

#########
### Gran
class(pData(newmset.EUR)$Gran)
pData(newmset.EUR)$Gran<-as.numeric(pData(newmset.EUR)$Gran)
summary(pData(newmset.EUR)$Gran)

#########
### Mono
class(pData(newmset.EUR)$Mono)
pData(newmset.EUR)$Mono<-as.numeric(pData(newmset.EUR)$Mono)
summary(pData(newmset.EUR)$Mono)

#########
### NK
class(pData(newmset.EUR)$NK)
pData(newmset.EUR)$NK<-as.numeric(pData(newmset.EUR)$NK)
summary(pData(newmset.EUR)$NK)

#########
### nRBC
class(pData(newmset.EUR)$nRBC)
pData(newmset.EUR)$nRBC<-as.numeric(pData(newmset.EUR)$nRBC)
summary(pData(newmset.EUR)$nRBC)

### Now we have checked that variables in the ExpressionSet are coded correctly.  


###########################
#STEP 9: ANALYSES - MODEL 1
###########################

#M1: Cord blood methylation ~ Green Spaces pregnancy + maternal age + maternal smoking pregnancy + child's sex + 
#childs ancestry within major ancestry group (optional) + batch (optional) + cohort (optional) + selection variable (optional)

#NOTE: this model is not adjusted by maternal education, neighbourhood SES, celltypes, PM25 preg, gestational age and birth weight

#IMPORTANT! Before running the analyses, please make sure that you have done STEP 8) Check variables: summaries and plots. In this way,
#you make sure that the variables are coded as indicated in the analysis plan.

##################################################
#9.1) Select variables needed for this model
##################################################

dd.EUR.M1<-pData(newmset.EUR)[,c("ndvi300_preg","ndvi100_preg","greenyn300_preg","agebirth_m_y","preg_smk","sex_methyl")]

###########################
#9.2) Complete cases 
###########################

############################################################################
#9.2.1) Create a new dataframe with complete cases from variables in Model 1

dd.EUR.M1.cc<-dd.EUR.M1[complete.cases(dd.EUR.M1),]

##################################
#9.2.2) Complete cases M1 samples

samplescmp.M1<-rownames(dd.EUR.M1.cc)

################################################################
#9.2.3) Subset the ExpressionSet with Complete cases M1 samples

newmset.EUR.M1<-newmset.EUR[,(sampleNames(newmset.EUR) %in% samplescmp.M1)]

##########################################################################################
#9.3)Create variables with sorrounding residential greenness(NDVI Variable) transformation  
##########################################################################################

### The effect of Residential Sorrounding Greenness variable (NDVI) will be reported by IQR, thus we will transform NDVI300 pregnancy and
# NDVI100 pregnancy variables.

########################
#9.3.1) NDVI300 IQR preg

### Calculate the IQR of NDVI300
ndvi300.iqr<-iqr(dd.EUR.M1.cc$ndvi300_preg,na.rm=TRUE)

### The formula to create this variable corresponds to the original variable value divided by the interquartile value of the variable 
dd.EUR.M1.cc$ndvi300_preg.iqr<-dd.EUR.M1.cc$ndvi300_preg/ndvi300.iqr
summary(dd.EUR.M1.cc$ndvi300_preg.iqr)

#########################
#9.3.2) NDVI100 IQR preg

### Calculate the IQR of NDVI100
ndvi100.iqr<-iqr(dd.EUR.M1.cc$ndvi100_preg,na.rm=TRUE)

### The formula to create this variable corresponds to the original variable value divided by the interquartile value of the variable 
dd.EUR.M1.cc$ndvi100_preg.iqr<-dd.EUR.M1.cc$ndvi100_preg/ndvi100.iqr
summary(dd.EUR.M1.cc$ndvi100_preg.iqr)

#######################################################################################################################
#9.4)Include the new dataframe with complete cases and the created NDVI variables transformed in the ExpressionSet
#######################################################################################################################

#Include the new dataframe with the complete case samples for M1 and the created NDVI variables transformed in the ExpressionSet subsetted in section 9.2.3.

###################
#9.4.1) Check order 

# Before including the new dataframe with the complete case samples for M1 and the created NDVI variables transformed in the ExpressionSet,
#check if the samples are in the same order. If you do not order samples as they are in the ExpressionSet, you could incorrectly assign the values 
#of the variables to the samples and therefore also to the methylation. 

table(ifelse(rownames(dd.EUR.M1.cc)==sampleNames(newmset.EUR.M1),"Matched","--NOT MATCHED--"))
#Matched
#    357

#In our case, samples are ordered in the same way. In case that your samples are not ordered, please, see STEP 5 for an example of how to do this. 

#############################################################
#9.4.2) Include the new dataframe in the ExpressionSet
pData(newmset.EUR.M1)<-dd.EUR.M1.cc

####################################
#9.5) Descriptive analyses for M1 
####################################

#9.5.1) Descriptive tables
#9.5.2) Correlation between exposures
#9.5.3) Associations between green spaces vars and covariates

### Go to the Descriptive.info directory 
setwd(paste0(res.dir,"/Descriptive.info/"))

###  Create the Descriptive.Tables directory
dir.create("Descriptive.Tables")
### Create the Correlations directory
dir.create("Correlations")
### Create the GS vs Covariates directory
dir.create("GSvsCovariates")

###########################
#9.5.1) Descriptive tables

### Go to the Descriptive.Tables directory
setwd(paste0(res.dir,"/Descriptive.info/Descriptive.Tables"))


#a) Mean (SD)
##############

descriptive.vars.M1<-colnames(dd.EUR.M1.cc)#Select variables for the descriptives

dd.EUR.M1.cc %>%
  dplyr::mutate(
    ndvi300_preg = ff_label(ndvi300_preg, "NDVI300 preg"),ndvi300_preg.iqr = ff_label(ndvi300_preg.iqr, "NDVI300 preg IQR"),
    ndvi100_preg = ff_label(ndvi100_preg, "NDVI100 preg"),ndvi100_preg.iqr = ff_label(ndvi100_preg.iqr, "NDVI100 preg IQR"),
    greenyn300_preg = ff_label(greenyn300_preg, "Green Access"),agebirth_m_y = ff_label(agebirth_m_y, "Maternal Age at delivery"),
    preg_smk = ff_label(preg_smk, "Maternal smoking preg"), 
    sex_methyl = ff_label(sex_methyl, "Child sex")
    )%>%
  summary_factorlist(dependent=NULL,descriptive.vars.M1,cont="mean",column=TRUE,na_include=TRUE,total_col=TRUE)->DescriptiveTable.Mean.M1

DescriptiveTable.Mean.M1<-DescriptiveTable.Mean.M1[,c("label","levels","Total")]
colnames(DescriptiveTable.Mean.M1)<-c("Variables","levels","Mean (SD) or n(%)")

#b) Median (IQR)
#############

dd.EUR.M1.cc %>%
  dplyr::mutate(
    ndvi300_preg = ff_label(ndvi300_preg, "NDVI300 preg"),ndvi300_preg.iqr = ff_label(ndvi300_preg.iqr, "NDVI300 preg IQR"),
    ndvi100_preg = ff_label(ndvi100_preg, "NDVI100 preg"),ndvi100_preg.iqr = ff_label(ndvi100_preg.iqr, "NDVI100 preg IQR"),
    greenyn300_preg = ff_label(greenyn300_preg, "Green Access"),agebirth_m_y = ff_label(agebirth_m_y, "Maternal Age at delivery"),
    preg_smk = ff_label(preg_smk, "Maternal smoking preg"), 
    sex_methyl = ff_label(sex_methyl, "Child sex")
    )%>%
  summary_factorlist(dependent=NULL,descriptive.vars.M1,cont="median",column=TRUE,na_include=TRUE,total_col=TRUE)->DescriptiveTable.Median.M1

DescriptiveTable.Median.M1<-DescriptiveTable.Median.M1[,c("label","levels","Total")]
colnames(DescriptiveTable.Median.M1)<-c("Variables","levels","Median (IQR) or n(%)")


### Create a single descriptive table with all variables

DescriptiveTable.final.M1<-cbind(DescriptiveTable.Mean.M1,DescriptiveTable.Median.M1)
DescriptiveTable.final.M1<-DescriptiveTable.final.M1[,c(1:3,5:6)]
             
write.table(
  DescriptiveTable.final.M1, file=paste0("M1.Descriptive.table.",cohort,".",ancestry,".cc.",Date,".txt"),col.names=T, row.names=F, quote=F, sep="\t")

######################################
#9.5.2) Correlation between exposures

### Go to the Correlations directory
setwd(paste0(res.dir,"/Descriptive.info/Correlations"))


### a)Correlation between NDVI300 and NDVI100 (Continuous vs continuous)
#########################################################################

cor.ndvi300.100<-dd.EUR.M1.cc[,c("ndvi300_preg","ndvi100_preg")]
#Create plot
plot.cor.ndvi300.100<-ggpairs(cor.ndvi300.100, title="Correlation between NDVI300 and NDVI100") 


#Define a path 
file.create(paste0("M1.CorrelationNDVI300vsNDVI100.",cohort,".",ancestry,"_",Date,".pdf"))

path1.M1<-file.path(paste0("M1.CorrelationNDVI300vsNDVI100.",cohort,".",ancestry,"_",Date,".pdf"))

#change path to the file name created just above.

pdf(path1.M1)
print(plot.cor.ndvi300.100)
dev.off()

###b)Correlation between Green Access and NDVI300(Categorical vs Continuous)
################################################################################

#Create plot
plot.cor.ndvi300.GreenAccess<-ggplot(dd.EUR.M1.cc,aes(x=greenyn300_preg, y = ndvi300_preg.iqr)) + ggtitle ("Correlation between Green access vs NDVI300") + geom_boxplot()
  
#Define a path 
file.create(paste0("M1.BoxPlot.GreenAccessvsNDVI300.",cohort,".",ancestry,"_",Date,".pdf"))

path2.M1<-file.path(paste0("M1.BoxPlot.GreenAccessvsNDVI300.",cohort,".",ancestry,"_",Date,".pdf"))

pdf(path2.M1)
print(plot.cor.ndvi300.GreenAccess)
dev.off()
  
###c)Correlation between Green Access and NDVI100(Categorical vs Continuous)
################################################################################

#Create Plot
plot.cor.ndvi100.GreenAccess<-ggplot(dd.EUR.M1.cc,aes(x=greenyn300_preg, y = ndvi100_preg.iqr)) + ggtitle ("Correlation between Green access vs NDVI100") + geom_boxplot()
  
#Define a path 
file.create(paste0("M1.BoxPlot.GreenAccessvsNDVI100.",cohort,".",ancestry,"_",Date,".pdf"))

path3.M1<-file.path(paste0("M1.BoxPlot.GreenAccessvsNDVI100.",cohort,".",ancestry,"_",Date,".pdf")) 

pdf(path3.M1)
print(plot.cor.ndvi100.GreenAccess)
dev.off()
  

###############################################################
#9.5.3) Associations between Green Spaces vars and Covariates

### Go to the GSvsCovariates directory
setwd(paste0(res.dir,"/Descriptive.info/GSvsCovariates"))

############################################################################
##a) Association beweeen Continuous Green spaces variables  and Covariables

GSexposures<-c("ndvi300_preg.iqr","ndvi100_preg.iqr")

Covariates<-c("agebirth_m_y","preg_smk","sex_methyl")#Specify variables M1 

coefs.M1.cont=as.data.frame(matrix(NA,nrow=1,ncol=10))
colnames(coefs.M1.cont)=c("DateTime","GreenSpaceVar","Covariates","Model","BetaNoStand", "SE", "P","R2","R2adj","N")

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
    coefs.M1.cont[count,1]=gsub(" ","_",Sys.time())
    coefs.M1.cont[count,2]=GSexposures[j]
    coefs.M1.cont[count,3]=Covariates[i]
    coefs.M1.cont[count,4]=ff
    coefs.M1.cont[count,5:7]=res[2,c(1,2,4)]
    coefs.M1.cont[count,8]= R2
    coefs.M1.cont[count,9]=R2.adj
    coefs.M1.cont[count,10]= s$df[2]
    
    count=count+1
    
  }
}

#Multiple testing correction by FDR 

### ndvi300 pregnancy IQR vs Covariates
coefs.M1.subset.ndvi300<-subset(coefs.M1.cont, coefs.M1.cont$GreenSpaceVar == "ndvi300_preg.iqr")
coefs.M1.subset.ndvi300$Padj<-  p.adjust(coefs.M1.subset.ndvi300$P, method="fdr")

### ndvi100 pregnancy IQR vs Covariates
coefs.M1.subset.ndvi100<-subset(coefs.M1.cont, coefs.M1.cont$GreenSpaceVar == "ndvi100_preg.iqr")
coefs.M1.subset.ndvi100$Padj<-  p.adjust(coefs.M1.subset.ndvi100$P, method="fdr")

### Create one table with all results
coefs.M1.cont.all<-rbind(coefs.M1.subset.ndvi300,coefs.M1.subset.ndvi100)

write.csv(coefs.M1.cont.all,paste0("M1.GScontVars.Covariates_",cohort,".",ancestry,"_",Date,".csv"), row.names = TRUE)

#############################################################################
##b) Association beweeen Categorical Green spaces variable  and Covariables 
GSexposures<-"greenyn300_preg"

Covariates<-c("agebirth_m_y","preg_smk","sex_methyl")#Specify variables M1 

coefs.M1.cat=as.data.frame(matrix(NA,nrow=1,ncol=9))
colnames(coefs.M1.cat)=c("DateTime","GreenSpaceVar","Covariates","Model","BetaNoStand", "SE", "P","R2","N")

count=1
for(j in 1:length(GSexposures)) {
  print(GSexposures[j])
  results=NULL
  for(i in 1:length(Covariates)) {
    print(Covariates[i])
    ff <- paste0(GSexposures[j]," ~ ",Covariates[i],sep="") 
    fit <- glm(ff,data=dd.EUR.M1.cc,family="binomial",na.action=na.omit)
    s <-summary(fit)
    res <- s$coefficients
    R2<- with(s,1 - (deviance/null.deviance))
        
    #Save analysis details and coefficients
    coefs.M1.cat[count,1]=gsub(" ","_",Sys.time())
    coefs.M1.cat[count,2]=GSexposures[j]
    coefs.M1.cat[count,3]=Covariates[i]
    coefs.M1.cat[count,4]=ff
    coefs.M1.cat[count,5:7]=res[2,c(1,2,4)]
    coefs.M1.cat[count,8]= R2
    coefs.M1.cat[count,9]= s$df[2]
    
    count=count+1
    
  }
}

write.csv(coefs.M1.cat,paste0("M1.GScatVar.Covariates_",cohort,".",ancestry,"_",Date,".csv"), row.names = TRUE)

##############################################################################
#9.6) EWAS GREEN SPACES PREGNANCY USING ROBUST LINEAR REGRESSIONS WITH LIMMA
##############################################################################

### Before starting with the EWAS analysis...

### Create a folder for the EWAS results 

setwd(res.dir)
dir.create("EWAS.Results")

### Go to the EWAS results directory
setwd(paste0(res.dir,"/EWAS.Results"))

#For the analysis use ExpressionSet subset created in section 9.1.3. 

#newmset.EUR.M1

# 9.6.1) Exposure NDVI300 Pregnancy 
# 9.6.2) Exposure NDVI100 Pregnancy 
# 9.6.3) Exposure Green Access Pregnancy 

###################################
# 9.6.1) Exposure NDVI300 Pregnancy 

################
# Preg_N300_M1
################

# a) Run EWAS

### Design the model
design.Preg_N300_M1<-model.matrix(~ndvi300_preg.iqr + agebirth_m_y + preg_smk + sex_methyl, data=pData(newmset.EUR.M1))

### Run Robust linear regression model with limma
fit.Preg_N300_M1<-limma::lmFit(newmset.EUR.M1,design=design.Preg_N300_M1, method="robust")
fit.Preg_N300_M1<-limma::eBayes(fit.Preg_N300_M1)
Preg_N300_M1<-limma::topTable(fit.Preg_N300_M1, coef =2, number = Inf,sort.by="none",confint=TRUE)
Preg_N300_M1$SE <- (sqrt(fit.Preg_N300_M1$s2.post) *fit.Preg_N300_M1$stdev.unscaled)[, 2]
head(Preg_N300_M1)
#                 logFC    AveExpr        t      P.Value adj.P.Val         B
#cg13296371 0.017341479 0.27668253 4.995655 9.952075e-07 0.4746602 2.7206509
#cg25003922 0.015592547 0.81803102 4.724886 3.536803e-06 0.5358988 1.4758272
#cg06050372 0.002000149 0.03092754 4.686020 4.223528e-06 0.5358988 1.3048712
#cg21149357 0.077089491 0.40737168 4.672343 4.494419e-06 0.5358988 1.2493529
#cg06071701 0.003817574 0.05547216 4.618135 5.742062e-06 0.5477307 1.0004611
#cg20940607 0.003449854 0.05744883 4.530151 8.504761e-06 0.6760520 0.6322633
#SE
#cg13296371 0.0032392393
#cg25003922 0.0032082155
#cg06050372 0.0009906259
#cg21149357 0.0038793332
#cg06071701 0.0005316320
#cg20940607 0.0059676090

###  Export results
setwd(paste0(res.dir,"/EWAS.Results"))
write.table(Preg_N300_M1, paste0("EWAS_GS_",cohort,".",ancestry,"_Preg_N300_M1_",Date,".txt"), na="NA") 
            
# b) Calculate lambda
lambda.Preg_N300_M1<- qchisq(median(Preg_N300_M1$P.Value,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
lambda.Preg_N300_M1

# c) QQ plot 

#Create QQplot directory 
setwd(paste0(res.dir,"/EWAS.Results"))
dir.create("QQplots")

### Go to QQplots directory
setwd(paste0(res.dir,"/EWAS.Results/QQPlots"))

pvals<-Preg_N300_M1$P.Value
jpeg(paste0("QQPlot_",cohort,".",ancestry,"_Preg_N300_M1_",Date,".jpg"))
qq(pvals,main=paste0("QQPlot_",cohort,".",ancestry,"_Preg_N300_M1")) 
dev.off()

###################################
# 9.6.2) Exposure NDVI100 Pregnancy 
            
################
# Preg_N100_M1
################

# a) Run EWAS

### Design the model
design.Preg_N100_M1<-model.matrix(~ndvi100_preg.iqr + agebirth_m_y + preg_smk + sex_methyl, data=pData(newmset.EUR.M1))

### Run Robust linear regression model with limma
fit.Preg_N100_M1<-limma::lmFit(newmset.EUR.M1,design=design.Preg_N100_M1, method="robust")
fit.Preg_N100_M1<-limma::eBayes(fit.Preg_N100_M1)
Preg_N100_M1<-limma::topTable(fit.Preg_N100_M1, coef =2, number = Inf,sort.by="none",confint=TRUE)
Preg_N100_M1$SE <- (sqrt(fit.Preg_N100_M1$s2.post) *fit.Preg_N100_M1$stdev.unscaled)[, 2]
head(Preg_N100_M1)

###  Export results
setwd(paste0(res.dir,"/EWAS.Results"))

write.table(Preg_N100_M1, paste0("EWAS_GS_",cohort,".",ancestry,"_Preg_N100_M1_",Date,".txt"), na="NA") 
            
# b) Calculate lambda
lambda.Preg_N100_M1<- qchisq(median(Preg_N100_M1$P.Value,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
lambda.Preg_N100_M1

# c) QQ plot 

### Go to QQplots directory
setwd(paste0(res.dir,"/EWAS.Results/QQPlots"))

pvalsM1<-Preg_N100_M1$P.Value
jpeg(paste0("QQPlot_",cohort,".",ancestry,"_Preg_N100_M1_",Date,".jpg"))
qq(pvalsM1,main=paste0("QQPlot_",cohort,".",ancestry,"_Preg_N100_M1")) 
dev.off()   

#########################################
# 9.6.3) Exposure Green Access Pregnancy                                                                   
                                
################
# Preg_AcGS_M1
################

# a) Run EWAS

### Design the model
design.Preg_AcGS_M1<-model.matrix(~greenyn300_preg + agebirth_m_y + preg_smk + sex_methyl, data=pData(newmset.EUR.M1))

### Run Robust linear regression model with limma
fit.Preg_AcGS_M1<-limma::lmFit(newmset.EUR.M1,design=design.Preg_AcGS_M1, method="robust")
fit.Preg_AcGS_M1<-limma::eBayes(fit.Preg_AcGS_M1)
Preg_AcGS_M1<-limma::topTable(fit.Preg_AcGS_M1, coef =2, number = Inf,sort.by="none",confint=TRUE)
Preg_AcGS_M1$SE <- (sqrt(fit.Preg_AcGS_M1$s2.post) *fit.Preg_AcGS_M1$stdev.unscaled)[, 2]
head(Preg_AcGS_M1)

###  Export results
setwd(paste0(res.dir,"/EWAS.Results"))
write.table(Preg_AcGS_M1, paste0("EWAS_GS_",cohort,".",ancestry,"_Preg_AcGS_M1_",Date,".txt"), na="NA") 
            
# b) Calculate lambda
lambda.Preg_AcGS_M1<- qchisq(median(Preg_AcGS_M1$P.Value,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
lambda.Preg_AcGS_M1

# c) QQ plot 

setwd(paste0(res.dir,"/EWAS.Results/QQPlots"))

pvalsM1<-Preg_AcGS_M1$P.Value
jpeg(paste0("QQPlot_",cohort,".",ancestry,"_Preg_AcGS_M1_",Date,".jpg"))
qq(pvalsM1,main=paste0("QQPlot_",cohort,".",ancestry,"_Preg_AcGS_M1")) 
dev.off()
              
####################
#9.7) LAMBDAS M1
####################

### Create lambdas table M1

### Create Lambdas directory
setwd(paste0(res.dir,"/EWAS.Results"))

dir.create("Lambdas")

### Go to Lambdas directory
setwd(paste0(res.dir,"/EWAS.Results/Lambdas"))
                                                    
lambdas.table<-rbind(lambda.Preg_N300_M1,lambda.Preg_N100_M1,lambda.Preg_AcGS_M1)
colnames(lambdas.table)<-"Lambdas M1"

write.table(lambdas.table, paste0("Lambdas_",cohort,".",ancestry,"_M1_",Date,".txt"), na="NA")

###########################
#STEP 10: ANALYSES- MODEL 2
###########################

#M2:Cord blood methylation = green spaces pregnancy + maternal age + maternal smoking pregnancy + child’s sex + child’s ancestry 
#within major ancestry group (optional) + batch (optional) + cohort (optional) + selection variables (optional) + maternal education
#+ neighbourhood socio-economic status 

#NOTE: this model is not adjusted by celltypes, PM25 preg, gestational age and birth weight
                                                
#IMPORTANT! Before running the analyses, please make sure that you have done STEP 8) Check variables: summaries and plots. In this way, 
# you make sure that the variables are coded as indicated in the analysis plan.

###########################################
#10.1) Select variables needed for Model 2
###########################################

dd.EUR.M2<-pData(newmset.EUR)[,c("ndvi300_preg","ndvi100_preg","greenyn300_preg","agebirth_m_y","preg_smk","sex_methyl","edu_m_0",
                                 "areases_tert_preg")]

###########################
#10.2) Complete cases M2
###########################

############################################################################
#10.2.1) Create a new dataframe with complete cases from variables in Model 2

dd.EUR.M2.cc<-dd.EUR.M2[complete.cases(dd.EUR.M2),]

##################################
#10.2.2) Complete cases M2 samples

samplescmp.M2<-rownames(dd.EUR.M2.cc)

################################################################
#10.2.3) Subset the ExpressionSet with Complete cases M2 samples

newmset.EUR.M2<-newmset.EUR[,(sampleNames(newmset.EUR) %in% samplescmp.M2)]

##########################################################################################
#10.3)Create variables with sorrounding residential greenness(NDVI Variable) transformation  
##########################################################################################

### The effect of Residential Sorrounding Greenness variable (NDVI) will be reported by IQR, thus we will transform NDVI300 pregnancy 
# and NDVI100 pregnancy variables.

##########################
#10.3.1) NDVI300 IQR preg

### Calculate the IQR of NDVI300
ndvi300.iqr<-iqr(dd.EUR.M2.cc$ndvi300_preg,na.rm=TRUE)

### The formula to create this variable corresponds to the original variable value divided by the interquartile value of the variable 
dd.EUR.M2.cc$ndvi300_preg.iqr<-dd.EUR.M2.cc$ndvi300_preg/ndvi300.iqr
summary(dd.EUR.M2.cc$ndvi300_preg.iqr)

#########################
#10.3.2) NDVI100 IQR preg

### Calculate the IQR of NDVI100
ndvi100.iqr<-iqr(dd.EUR.M2.cc$ndvi100_preg,na.rm=TRUE)

### The formula to create this variable corresponds to the original variable value divided by the interquartile value of the variable 
dd.EUR.M2.cc$ndvi100_preg.iqr<-dd.EUR.M2.cc$ndvi100_preg/ndvi100.iqr
summary(dd.EUR.M2.cc$ndvi100_preg.iqr)

#######################################################################################################################
#10.4)Include the new dataframe with complete cases M2 and the created NDVI variables transformed in the ExpressionSet
#######################################################################################################################

#Include the new dataframe with the complete case samples for M2 and the created NDVI variables transformed in the ExpressionSet subsetted in section 10.2.3.

#####################
#10.4.1) Check order 

# Before including the new dataframe with the complete case samples for M2 and the created NDVI variables transformed in the ExpressionSet,
# check if the samples are in the same order. If you do not order samples as they are in the ExpressionSet, you could incorrectly assign the values
# of the variables to the samples and therefore also to the methylation. 

table(ifelse(rownames(dd.EUR.M2.cc)==sampleNames(newmset.EUR.M2),"Matched","--NOT MATCHED--"))
#Matched
#    357

#In our case, samples are ordered in the same way. In case that your samples are not ordered, please, see STEP 5 for an example of how to do this. 

#############################################################
#10.4.2) Include the new dataframe in the ExpressionSet

pData(newmset.EUR.M2)<-dd.EUR.M2.cc

####################################
#10.5) Descriptive analyses for M2
####################################

#NOTE: We will not run descriptive analyses for M2 as we understand that the N will be the same as in M3 as model 3 is model 2 + 
#celltypes proportions which should not have NAs. 

##############################################################################
#10.6) EWAS GREEN SPACES PREGNANCY USING ROBUST LINEAR REGRESSIONS WITH LIMMA
##############################################################################

### Before starting with the EWAS analysis...

### Go to the EWAS results directory
setwd(paste0(res.dir,"/EWAS.Results"))

#For the analysis use ExpressionSet subset created in section 10.1.3. 

#newmset.EUR.M2

# 10.6.1) Exposure NDVI300 Pregnancy 
# 10.6.2) Exposure NDVI100 Pregnancy 
# 10.6.3) Exposure Green Access Pregnancy 

###################################
# 10.6.1) Exposure NDVI300 Pregnancy 

################
# Preg_N300_M2
################

# a) Run EWAS

### Design the model
design.Preg_N300_M2<-model.matrix(~ndvi300_preg.iqr + agebirth_m_y + preg_smk + sex_methyl + edu_m_0 + areases_tert_preg, data=pData(newmset.EUR.M2))

### Run Robust linear regression model with limma
fit.Preg_N300_M2<-limma::lmFit(newmset.EUR.M2,design=design.Preg_N300_M2, method="robust")
fit.Preg_N300_M2<-limma::eBayes(fit.Preg_N300_M2)
Preg_N300_M2<-limma::topTable(fit.Preg_N300_M2, coef =2, number = Inf,sort.by="none",confint=TRUE)
Preg_N300_M2$SE <- (sqrt(fit.Preg_N300_M2$s2.post) *fit.Preg_N300_M2$stdev.unscaled)[, 2]
head(Preg_N300_M2)

###  Export results
setwd(paste0(res.dir,"/EWAS.Results"))
write.table(Preg_N300_M2, paste0("EWAS_GS_",cohort,".",ancestry,"_Preg_N300_M2_",Date,".txt"), na="NA") 
            
# b) Calculate lambda
lambda.Preg_N300_M2<- qchisq(median(Preg_N300_M2$P.Value,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
lambda.Preg_N300_M2

# c) QQ plot 

### Go to QQplots directory
setwd(paste0(res.dir,"/EWAS.Results/QQPlots"))

pvals<-Preg_N300_M2$P.Value
jpeg(paste0("QQPlot_",cohort,".",ancestry,"_Preg_N300_M2_",Date,".jpg"))
qq(pvals,main=paste0("QQPlot_",cohort,".",ancestry,"_Preg_N300_M2")) 
dev.off()

###################################
# 13.6.2) Exposure NDVI100 Pregnancy 
            
################
# Preg_N100_M2
################

# a) Run EWAS

### Design the model
design.Preg_N100_M2<-model.matrix(~ndvi100_preg.iqr  +  agebirth_m_y + preg_smk + sex_methyl + edu_m_0 + areases_tert_preg, data=pData(newmset.EUR.M2))

### Run Robust linear regression model with limma
fit.Preg_N100_M2<-limma::lmFit(newmset.EUR.M2,design=design.Preg_N100_M2, method="robust")
fit.Preg_N100_M2<-limma::eBayes(fit.Preg_N100_M2)
Preg_N100_M2<-limma::topTable(fit.Preg_N100_M2, coef =2, number = Inf,sort.by="none",confint=TRUE)
Preg_N100_M2$SE <- (sqrt(fit.Preg_N100_M2$s2.post) *fit.Preg_N100_M2$stdev.unscaled)[, 2]
head(Preg_N100_M2)

###  Export results
setwd(paste0(res.dir,"/EWAS.Results"))
write.table(Preg_N100_M2, paste0("EWAS_GS_",cohort,".",ancestry,"_Preg_N100_M2_",Date,".txt"), na="NA") 
            
# b) Calculate lambda
lambda.Preg_N100_M2<- qchisq(median(Preg_N100_M2$P.Value,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
lambda.Preg_N100_M2

# c) QQ plot 

setwd(paste0(res.dir,"/EWAS.Results/QQPlots"))

pvalsM2<-Preg_N100_M2$P.Value
jpeg(paste0("QQPlot_",cohort,".",ancestry,"_Preg_N100_M2_",Date,".jpg"))
qq(pvalsM2,main=paste0("QQPlot_",cohort,".",ancestry,"_Preg_N100_M2")) 
dev.off()   

#########################################
# 10.6.3) Exposure Green Access Pregnancy                                                                   
                                
################
# Preg_AcGS_M2
################

# a) Run EWAS

### Design the model

design.Preg_AcGS_M2<-model.matrix(~greenyn300_preg +  agebirth_m_y + preg_smk + sex_methyl + edu_m_0 + areases_tert_preg, data=pData(newmset.EUR.M2))

### Run Robust linear regression model with limma
fit.Preg_AcGS_M2<-limma::lmFit(newmset.EUR.M2,design=design.Preg_AcGS_M2, method="robust")
fit.Preg_AcGS_M2<-limma::eBayes(fit.Preg_AcGS_M2)
Preg_AcGS_M2<-limma::topTable(fit.Preg_AcGS_M2, coef =2, number = Inf,sort.by="none",confint=TRUE)
Preg_AcGS_M2$SE <- (sqrt(fit.Preg_AcGS_M2$s2.post) *fit.Preg_AcGS_M2$stdev.unscaled)[, 2]
head(Preg_AcGS_M2)

###  Export results
setwd(paste0(res.dir,"/EWAS.Results"))
write.table(Preg_AcGS_M2, paste0("EWAS_GS_",cohort,".",ancestry,"_Preg_AcGS_M2_",Date,".txt"), na="NA") 
            
# b) Calculate lambda
lambda.Preg_AcGS_M2<- qchisq(median(Preg_AcGS_M2$P.Value,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
lambda.Preg_AcGS_M2

# c) QQ plot 

setwd(paste0(res.dir,"/EWAS.Results/QQPlots"))

pvalsM2<-Preg_AcGS_M2$P.Value
jpeg(paste0("QQPlot_",cohort,".",ancestry,"_Preg_AcGS_M2_",Date,".jpg"))
qq(pvalsM2,main=paste0("QQPlot_",cohort,".",ancestry,"_Preg_AcGS_M2")) 
dev.off()
              
####################
#10.7) LAMBDAS M2
####################

### Create lambdas table M2

setwd(paste0(res.dir,"/EWAS.Results/Lambdas"))                                                            
lambdas.table<-rbind(lambda.Preg_N300_M2,lambda.Preg_N100_M2,lambda.Preg_AcGS_M2)
colnames(lambdas.table)<-"Lambdas M2"

write.table(lambdas.table, paste0("Lambdas_",cohort,".",ancestry,"_M2_",Date,".txt"), na="NA")

###########################
#STEP 11: ANALYSES- MODEL 3
###########################


#M3:Cord blood methylation = green spaces pregnancy + maternal age + maternal smoking pregnancy + child’s sex + child’s
#ancestry within major ancestry group (optional) + batch (optional) + cohort (optional) + selection variables (optional) + maternal education +
#neighbourhood socio-economic status + blood cellular composition

#NOTE: this model is not adjusted by PM25 preg, gestational age and birth weight
                                                 
#IMPORTANT! Before running the analyses, please make sure that you have done STEP 8) Check variables: summaries and plots. In this way,
#you make sure that the variables are coded as indicated in the analysis plan.

###########################################
#11.1) Select variables needed for Model 3
###########################################

dd.EUR.M3<-pData(newmset.EUR)[,c("ndvi300_preg","ndvi100_preg","greenyn300_preg","agebirth_m_y","preg_smk","sex_methyl","edu_m_0","areases_tert_preg","Bcell","CD4T","CD8T","Gran","Mono","NK","nRBC")]

###########################
#11.2) Complete cases M3
###########################

############################################################################
#11.2.1) Create a new dataframe with complete cases from variables in Model 3

dd.EUR.M3.cc<-dd.EUR.M3[complete.cases(dd.EUR.M3),]

##################################
#11.2.2) Complete cases M3 samples

samplescmp.M3<-rownames(dd.EUR.M3.cc)

################################################################
#11.2.3) Subset the ExpressionSet with Complete cases M3 samples

newmset.EUR.M3<-newmset.EUR[,(sampleNames(newmset.EUR) %in% samplescmp.M3)]

##########################################################################################
#11.3)Create variables with sorrounding residential greenness(NDVI Variable) transformation  
##########################################################################################

### The effect of Residential Sorrounding Greenness variable (NDVI) will be reported by IQR, thus we will transform NDVI300 pregnancy 
# and NDVI100 pregnancy variables.

########################
#11.3.1) NDVI300 IQR preg

### Calculate the IQR of NDVI300
ndvi300.iqr<-iqr(dd.EUR.M3.cc$ndvi300_preg,na.rm=TRUE)

### The formula to create this variable corresponds to the original variable value divided by the interquartile value of the variable 
dd.EUR.M3.cc$ndvi300_preg.iqr<-dd.EUR.M3.cc$ndvi300_preg/ndvi300.iqr
summary(dd.EUR.M3.cc$ndvi300_preg.iqr)

#########################
#11.3.2) NDVI100 IQR preg

### Calculate the IQR of NDVI100
ndvi100.iqr<-iqr(dd.EUR.M3.cc$ndvi100_preg,na.rm=TRUE)

### The formula to create this variable corresponds to the original variable value divided by the interquartile value of the variable 
dd.EUR.M3.cc$ndvi100_preg.iqr<-dd.EUR.M3.cc$ndvi100_preg/ndvi100.iqr
summary(dd.EUR.M3.cc$ndvi100_preg.iqr)

#######################################################################################################################
#11.4)Include the new dataframe with complete cases M3 and the created NDVI variables transformed in the ExpressionSet
#######################################################################################################################

#Include the new dataframe with the complete case samples for M3 and the created NDVI variables transformed in the ExpressionSet subsetted in section 11.2.3.

###################
#11.4.1) Check order 

# Before including the new dataframe with the complete case samples for M3 and the created NDVI variables transformed in the ExpressionSet, 
# check if the samples are in the same order. If you do not order samples as they are in the ExpressionSet, you could incorrectly assign the values
#of the variables to the samples and therefore also to the methylation. 

table(ifelse(rownames(dd.EUR.M3.cc)==sampleNames(newmset.EUR.M3),"Matched","--NOT MATCHED--"))
#Matched
#    357

#In our case, samples are ordered in the same way. In case that your samples are not ordered, please, see STEP 5 for an example of how to do this. 

#############################################################
#11.4.2) Include the new dataframe in the ExpressionSet

pData(newmset.EUR.M3)<-dd.EUR.M3.cc

####################################
#11.5) Descriptive analyses for M3
####################################

#11.5.1) Descriptive tables
#11.5.2) Correlation between exposures
#11.5.3) Associations between green spaces vars and covariates


###########################
#11.5.1) Descriptive tables

### Go to the Descriptive.Tables directory

setwd(paste0(res.dir,"/Descriptive.info/Descriptive.Tables"))

#a) Mean (SD)
##############

descriptive.vars.M3<-colnames(dd.EUR.M3.cc)#Select variables for the descriptives

dd.EUR.M3.cc %>%
  dplyr::mutate(
    ndvi300_preg = ff_label(ndvi300_preg, "NDVI300 preg"),ndvi300_preg.iqr = ff_label(ndvi300_preg.iqr, "NDVI300 preg IQR"),
    ndvi100_preg = ff_label(ndvi100_preg, "NDVI100 preg"),ndvi100_preg.iqr = ff_label(ndvi100_preg.iqr, "NDVI100 preg IQR"),
    greenyn300_preg = ff_label(greenyn300_preg, "Green Access"),agebirth_m_y = ff_label(agebirth_m_y, "Maternal Age at delivery"),
   edu_m_0 = ff_label(edu_m_0, "Maternal education"),
    areases_tert_preg = ff_label(areases_tert_preg, "Neighbourhood socio-economic status"),preg_smk = ff_label(preg_smk, "Maternal smoking preg"),
    sex_methyl = ff_label(sex_methyl, "Child sex")
    )%>%
  summary_factorlist(dependent=NULL,descriptive.vars.M3,cont="mean",column=TRUE,na_include=TRUE,total_col=TRUE)->DescriptiveTable.Mean.M3

DescriptiveTable.Mean.M3<-DescriptiveTable.Mean.M3[,c("label","levels","Total")]
colnames(DescriptiveTable.Mean.M3)<-c("Variables","levels","Mean (SD) or n(%)")

#Median (IQR)
#############

dd.EUR.M3.cc %>%
  dplyr::mutate(
    ndvi300_preg = ff_label(ndvi300_preg, "NDVI300 preg"),ndvi300_preg.iqr = ff_label(ndvi300_preg.iqr, "NDVI300 preg IQR"),
    ndvi100_preg = ff_label(ndvi100_preg, "NDVI100 preg"),ndvi100_preg.iqr = ff_label(ndvi100_preg.iqr, "NDVI100 preg IQR"),
    greenyn300_preg = ff_label(greenyn300_preg, "Green Access"),agebirth_m_y = ff_label(agebirth_m_y, "Maternal Age at delivery"),
    edu_m_0 = ff_label(edu_m_0, "Maternal education"),areases_tert_preg = ff_label(areases_tert_preg, "Neighbourhood socio-economic status"),
    preg_smk = ff_label(preg_smk, "Maternal smoking preg"),sex_methyl = ff_label(sex_methyl, "Child sex")
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
#11.5.2) Correlation between exposures

### Go to the Correlations directory
setwd(paste0(res.dir,"/Descriptive.info/Correlations"))

### a)Correlation between NDVI300 and NDVI100 (Continuous vs continuous)
#########################################################################

cor.ndvi300.100<-dd.EUR.M3.cc[,c("ndvi300_preg","ndvi100_preg")]
#Create plot
plot.cor.ndvi300.100<-ggpairs(cor.ndvi300.100, title="Correlation between NDVI300 and NDVI100") 

#Define a path 
file.create(paste0("M3.CorrelationNDVI300vsNDVI100.",cohort,".",ancestry,"_",Date,".pdf"))

path1.M3<-file.path(paste0("M3.CorrelationNDVI300vsNDVI100.",cohort,".",ancestry,"_",Date,".pdf"))

pdf(path1.M3)
print(plot.cor.ndvi300.100)
dev.off()

###b)Correlation between Green Access and NDVI300(Categorical vs Continuous)
################################################################################

#Ceate Plot
plot.cor.ndvi300.GreenAccess<-ggplot(dd.EUR.M3.cc,aes(x=greenyn300_preg, y = ndvi300_preg.iqr)) +
  ggtitle ("Correlation between Green access vs NDVI300") + geom_boxplot()
  
#Define a path 
file.create(paste0("M3.BoxPlot.GreenAccessvsNDVI300.",cohort,".",ancestry,"_",Date,".pdf"))

path2.M3<-file.path(paste0("M3.BoxPlot.GreenAccessvsNDVI300.",cohort,".",ancestry,"_",Date,".pdf"))

pdf(path2.M3)
print(plot.cor.ndvi300.GreenAccess)
dev.off()
  
###c)Correlation between Green Access and NDVI100(Categorical vs Continuous)
################################################################################

#Create Plot
plot.cor.ndvi100.GreenAccess<-ggplot(dd.EUR.M3.cc,aes(x=greenyn300_preg, y = ndvi100_preg.iqr))+ 
ggtitle ("Correlation between Green access vs NDVI100") + geom_boxplot()
  
#Define a path 
file.create(paste0("M3.BoxPlot.GreenAccessvsNDVI100.",cohort,".",ancestry,"_",Date,".pdf"))

path3.M3<-(paste0("M3.BoxPlot.GreenAccessvsNDVI100.",cohort,".",ancestry,"_",Date,".pdf")) 

pdf(path3.M3)
print(plot.cor.ndvi100.GreenAccess)
dev.off()
  

###############################################################
#11.5.3) Associations between Green Spaces vars and Covariates

### Go to the GSvsCovariates directory
setwd(paste0(res.dir,"/Descriptive.info/GSvsCovariates"))

############################################################################
##a) Association beweeen Continuous Green spaces variables  and Covariables

GSexposures<-c("ndvi300_preg.iqr","ndvi100_preg.iqr")

Covariates<-c("agebirth_m_y","preg_smk","sex_methyl","edu_m_0","areases_tert_preg",
              "Bcell","CD4T","CD8T","Gran","Mono","NK","nRBC")#Specify variables M3 

coefs.M3.cont=as.data.frame(matrix(NA,nrow=1,ncol=10))
colnames(coefs.M3.cont)=c("DateTime","GreenSpaceVar","Covariates","Model","BetaNoStand", "SE", "P","R2","R2adj","N")

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
    coefs.M3.cont[count,1]=gsub(" ","_",Sys.time())
    coefs.M3.cont[count,2]=GSexposures[j]
    coefs.M3.cont[count,3]=Covariates[i]
    coefs.M3.cont[count,4]=ff
    coefs.M3.cont[count,5:7]=res[2,c(1,2,4)]
    coefs.M3.cont[count,8]= R2
    coefs.M3.cont[count,9]=R2.adj
    coefs.M3.cont[count,10]= s$df[2]
    
    count=count+1
    
  }
}

#Multiple testing correction by FDR 

### ndvi300 pregnancy IQR vs Covariates
coefs.M3.subset.ndvi300<-subset(coefs.M3.cont, coefs.M3.cont$GreenSpaceVar == "ndvi300_preg.iqr")
coefs.M3.subset.ndvi300$Padj<-  p.adjust(coefs.M3.subset.ndvi300$P, method="fdr")

### ndvi100 pregnancy IQR vs Covariates
coefs.M3.subset.ndvi100<-subset(coefs.M3.cont, coefs.M3.cont$GreenSpaceVar == "ndvi100_preg.iqr")
coefs.M3.subset.ndvi100$Padj<-  p.adjust(coefs.M3.subset.ndvi100$P, method="fdr")

### Create one table with all results
coefs.M3.cont.all<-rbind(coefs.M3.subset.ndvi300,coefs.M3.subset.ndvi100)

write.csv(coefs.M3.cont.all,paste0("M3.GScontVars.Covariates_",cohort,".",ancestry,"_",Date,".csv"), row.names = TRUE)

#############################################################################
##b) Association beweeen Categorical Green spaces variable  and Covariables 
GSexposures<-"greenyn300_preg"

Covariates<-c("agebirth_m_y","preg_smk","sex_methyl","edu_m_0","areases_tert_preg",
              "Bcell","CD4T","CD8T","Gran","Mono","NK","nRBC")#Specify variables M3 

coefs.M3.cat=as.data.frame(matrix(NA,nrow=1,ncol=9))
colnames(coefs.M3.cat)=c("DateTime","GreenSpaceVar","Covariates","Model","BetaNoStand", "SE", "P","R2","N")

count=1
for(j in 1:length(GSexposures)) {
  print(GSexposures[j])
  results=NULL
  for(i in 1:length(Covariates)) {
    print(Covariates[i])
    ff <- paste0(GSexposures[j]," ~ ",Covariates[i],sep="") 
    fit <- glm(ff,data=dd.EUR.M3.cc,family="binomial",na.action=na.omit)
    s <-summary(fit)
    res <- s$coefficients
    R2<- with(s,1 - (deviance/null.deviance))
        
    #Save analysis details and coefficients
    coefs.M3.cat[count,1]=gsub(" ","_",Sys.time())
    coefs.M3.cat[count,2]=GSexposures[j]
    coefs.M3.cat[count,3]=Covariates[i]
    coefs.M3.cat[count,4]=ff
    coefs.M3.cat[count,5:7]=res[2,c(1,2,4)]
    coefs.M3.cat[count,8]= R2
    coefs.M3.cat[count,9]= s$df[2]
    
    count=count+1
    
  }
}

write.csv(coefs.M3.cat,paste0("M3.GScatVar.Covariates_",cohort,".",ancestry,"_",Date,".csv"), row.names = TRUE)

##############################################################################
#11.6) EWAS GREEN SPACES PREGNANCY USING ROBUST LINEAR REGRESSIONS WITH LIMMA
##############################################################################

### Before starting with the EWAS analysis...

### Go to  the EWAS results directory
setwd(paste0(res.dir,"/EWAS.Results"))

#For the analysis use ExpressionSet subset created in section 13.1.3. 

#newmset.EUR.M3

# 11.6.1) Exposure NDVI300 Pregnancy 
# 11.6.2) Exposure NDVI100 Pregnancy 
# 11.6.3) Exposure Green Access Pregnancy 

###################################
# 11.6.1) Exposure NDVI300 Pregnancy 

################
# Preg_N300_M3
################

# a) Run EWAS

### Design the model

t <- proc.time() # Inicia el cronómetro

design.Preg_N300_M3<-model.matrix(~ndvi300_preg.iqr + agebirth_m_y + preg_smk + sex_methyl + edu_m_0 + areases_tert_preg 
+ Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, data=pData(newmset.EUR.M3))

# NUESTRO CODIGO
### Run Robust linear regression model with limma
fit.Preg_N300_M3<-limma::lmFit(newmset.EUR.M3,design=design.Preg_N300_M3, method="robust")
fit.Preg_N300_M3<-limma::eBayes(fit.Preg_N300_M3)
Preg_N300_M3<-limma::topTable(fit.Preg_N300_M3, coef =2, number = Inf,sort.by="none",confint=TRUE)
Preg_N300_M3$SE <- (sqrt(fit.Preg_N300_M3$s2.post) *fit.Preg_N300_M3$stdev.unscaled)[, 2]
head(Preg_N300_M3)
proc.time()-t 

###  Export results
setwd(paste0(res.dir,"/EWAS.Results"))
write.table(Preg_N300_M3, paste0("EWAS_GS_",cohort,".",ancestry,"_Preg_N300_M3_",Date,".txt"), na="NA") 
            
# b) Calculate lambda
lambda.Preg_N300_M3<- qchisq(median(Preg_N300_M3$P.Value,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
lambda.Preg_N300_M3

# c) QQ plot 

### Go to QQplots directory
setwd(paste0(res.dir,"/EWAS.Results/QQPlots"))

pvals<-Preg_N300_M3$P.Value
jpeg(paste0("QQPlot_",cohort,".",ancestry,"_Preg_N300_M3_",Date,".jpg"))
qq(pvals,main=paste0("QQPlot_",cohort,".",ancestry,"_Preg_N300_M3")) 
dev.off()

###################################
# 11.6.2) Exposure NDVI100 Pregnancy 
            
################
# Preg_N100_M3
################

# a) Run EWAS

### Design the model
design.Preg_N100_M3<-model.matrix(~ndvi100_preg.iqr  +  agebirth_m_y + preg_smk + sex_methyl + edu_m_0 + areases_tert_preg + Bcell  + CD4T + CD8T + Gran + Mono + NK + nRBC, data=pData(newmset.EUR.M3))

### Run Robust linear regression model with limma
fit.Preg_N100_M3<-limma::lmFit(newmset.EUR.M3,design=design.Preg_N100_M3, method="robust")
fit.Preg_N100_M3<-limma::eBayes(fit.Preg_N100_M3)
Preg_N100_M3<-limma::topTable(fit.Preg_N100_M3, coef =2, number = Inf,sort.by="none",confint=TRUE)
Preg_N100_M3$SE <- (sqrt(fit.Preg_N100_M3$s2.post) *fit.Preg_N100_M3$stdev.unscaled)[, 2]
head(Preg_N100_M3)

###  Export results
setwd(paste0(res.dir,"/EWAS.Results"))

write.table(Preg_N100_M3, paste0("EWAS_GS_",cohort,".",ancestry,"_Preg_N100_M3_",Date,".txt"), na="NA") 
            
# b) Calculate lambda
lambda.Preg_N100_M3<- qchisq(median(Preg_N100_M3$P.Value,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
lambda.Preg_N100_M3

# c) QQ plot 

setwd(paste0(res.dir,"/EWAS.Results/QQPlots"))

pvalsM3<-Preg_N100_M3$P.Value
jpeg(paste0("QQPlot_",cohort,".",ancestry,"_Preg_N100_M3_",Date,".jpg"))
qq(pvalsM3,main=paste0("QQPlot_",cohort,".",ancestry,"_Preg_N100_M3")) 
dev.off()   

#########################################
# 13.6.3) Exposure Green Access Pregnancy                                                                   
                                
################
# Preg_AcGS_M3
################

# a) Run EWAS

### Design the model

design.Preg_AcGS_M3<-model.matrix(~greenyn300_preg  +  agebirth_m_y + preg_smk + sex_methyl+ edu_m_0 + areases_tert_preg 
                                  + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, data=pData(newmset.EUR.M3))

### Run Robust linear regression model with limma
fit.Preg_AcGS_M3<-limma::lmFit(newmset.EUR.M3,design=design.Preg_AcGS_M3, method="robust")
fit.Preg_AcGS_M3<-limma::eBayes(fit.Preg_AcGS_M3)
Preg_AcGS_M3<-limma::topTable(fit.Preg_AcGS_M3, coef =2, number = Inf,sort.by="none",confint=TRUE)
Preg_AcGS_M3$SE <- (sqrt(fit.Preg_AcGS_M3$s2.post) *fit.Preg_AcGS_M3$stdev.unscaled)[, 2]
head(Preg_AcGS_M3)

###  Export results
setwd(paste0(res.dir,"/EWAS.Results"))
write.table(Preg_AcGS_M3, paste0("EWAS_GS_",cohort,".",ancestry,"_Preg_AcGS_M3_",Date,".txt"), na="NA") 
            
# b) Calculate lambda
lambda.Preg_AcGS_M3<- qchisq(median(Preg_AcGS_M3$P.Value,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
lambda.Preg_AcGS_M3

# c) QQ plot 
setwd(paste0(res.dir,"/EWAS.Results/QQPlots"))

pvalsM3<-Preg_AcGS_M3$P.Value
jpeg(paste0("QQPlot_",cohort,".",ancestry,"_Preg_AcGS_M3_",Date,".jpg"))
qq(pvalsM3,main=paste0("QQPlot_",cohort,".",ancestry,"_Preg_AcGS_M3")) 
dev.off()
              
####################
#11.7) LAMBDAS M3
####################

### Create lambdas table M3
setwd(paste0(res.dir,"/EWAS.Results/Lambdas"))
                                                           
lambdas.table<-rbind(lambda.Preg_N300_M3,lambda.Preg_N100_M3,lambda.Preg_AcGS_M3)
colnames(lambdas.table)<-"Lambdas M3"

write.table(lambdas.table, paste0("Lambdas_",cohort,".",ancestry,"_M3_",Date,".txt"), na="NA")

###########################
#STEP 12: ANALYSES- MODEL 4
###########################

#M4:Cord blood methylation = green spaces pregnancy + maternal education + neighbourhood socio-economic status + maternal age +
#maternal smoking pregnancy + child’s sex + child’s ancestry within major ancestry group (optional) + batch (optional) + cohort (optional) + 
#selection variables (optional) + blood cellular composition + PM2.5 

#NOTE: this model is not adjusted by gestational age and birth weight
                                                
#IMPORTANT! Before running the analyses, please make sure that you have done STEP 8) Check variables: summaries and plots. In this way,
#you make sure that the variables are coded as indicated in the analysis plan.

###########################################
#12.1) Select variables needed for Model 4
###########################################

dd.EUR.M4<-pData(newmset.EUR)[,c("ndvi300_preg","ndvi100_preg","greenyn300_preg","agebirth_m_y","preg_smk","sex_methyl","edu_m_0",
                                 "areases_tert_preg","Bcell","CD4T","CD8T","Gran","Mono","NK","nRBC","pm25_preg")]

###########################
#12.2) Complete cases M4
###########################

############################################################################
#12.2.1) Create a new dataframe with complete cases from variables in Model 4

dd.EUR.M4.cc<-dd.EUR.M4[complete.cases(dd.EUR.M4),]

##################################
#12.2.2) Complete cases M4 samples

samplescmp.M4<-rownames(dd.EUR.M4.cc)

################################################################
#12.2.3) Subset the ExpressionSet with Complete cases M4 samples

newmset.EUR.M4<-newmset.EUR[,(sampleNames(newmset.EUR) %in% samplescmp.M4)]

##########################################################################################
#12.3)Create variables with sorrounding residential greenness(NDVI Variable) transformation  
##########################################################################################

### The effect of Residential Sorrounding Greenness variable (NDVI) will be reported by IQR, thus we will transform NDVI300 pregnancy
# and NDVI100 pregnancy variables.

##########################
#12.3.1) NDVI300 IQR preg

### Calculate the IQR of NDVI300
ndvi300.iqr<-iqr(dd.EUR.M4.cc$ndvi300_preg,na.rm=TRUE)

### The formula to create this variable corresponds to the original variable value divided by the interquartile value of the variable 
dd.EUR.M4.cc$ndvi300_preg.iqr<-dd.EUR.M4.cc$ndvi300_preg/ndvi300.iqr
summary(dd.EUR.M4.cc$ndvi300_preg.iqr)

#########################
#12.3.2) NDVI100 IQR preg

### Calculate the IQR of NDVI100
ndvi100.iqr<-iqr(dd.EUR.M4.cc$ndvi100_preg,na.rm=TRUE)

### The formula to create this variable corresponds to the original variable value divided by the interquartile value of the variable 
dd.EUR.M4.cc$ndvi100_preg.iqr<-dd.EUR.M4.cc$ndvi100_preg/ndvi100.iqr
summary(dd.EUR.M4.cc$ndvi100_preg.iqr)

#######################################################################################################################
#12.4)Include the new dataframe with complete cases M4 and the created NDVI variables transformed in the ExpressionSet
#######################################################################################################################

#Include the new dataframe with the complete case samples for M4 and the created NDVI variables transformed in the ExpressionSet subsetted in section 12.2.3.

###################
#12.1) Check order 

# Before including the new dataframe with the complete case samples for M4 and the created NDVI variables transformed in the ExpressionSet,
# check if the samples are in the same order. If you do not order samples as they are in the ExpressionSet, you could incorrectly assign 
# the values of the variables to the samples and therefore also to the methylation. 

table(ifelse(rownames(dd.EUR.M4.cc)==sampleNames(newmset.EUR.M4),"Matched","--NOT MATCHED--"))
#Matched
#    357

#In our case, samples are ordered in the same way. In case that your samples are not ordered, please, see STEP 5 for an example of how to do this. 

#############################################################
#12.4.2) Include the new dataframe in the ExpressionSet

pData(newmset.EUR.M4)<-dd.EUR.M4.cc

####################################
#12.5) Descriptive analyses for M4
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

descriptive.vars.M4<-colnames(dd.EUR.M4.cc)#Select variables for the descriptives

dd.EUR.M4.cc %>%
  dplyr::mutate(
    ndvi300_preg = ff_label(ndvi300_preg, "NDVI300 preg"),ndvi300_preg.iqr = ff_label(ndvi300_preg.iqr, "NDVI300 preg IQR"),
    ndvi100_preg = ff_label(ndvi100_preg, "NDVI100 preg"),ndvi100_preg.iqr = ff_label(ndvi100_preg.iqr, "NDVI100 preg IQR"),
    greenyn300_preg = ff_label(greenyn300_preg, "Green Access"),agebirth_m_y = ff_label(agebirth_m_y, "Maternal Age at delivery"),
    edu_m_0 = ff_label(edu_m_0, "Maternal education"),areases_tert_preg = ff_label(areases_tert_preg, "Neighbourhood socio-economic status"),
    preg_smk = ff_label(preg_smk, "Maternal smoking preg"),sex_methyl = ff_label(sex_methyl, "Child sex"),pm25_preg = ff_label(pm25_preg, "PM2.5 levels preg")
    )%>%
  summary_factorlist(dependent=NULL,descriptive.vars.M4,cont="mean",column=TRUE,na_include=TRUE,total_col=TRUE)->DescriptiveTable.Mean.M4

DescriptiveTable.Mean.M4<-DescriptiveTable.Mean.M4[,c("label","levels","Total")]
colnames(DescriptiveTable.Mean.M4)<-c("Variables","levels","Mean (SD) or n(%)")

#Median (IQR)
#############

dd.EUR.M4.cc %>%
  dplyr::mutate(
    ndvi300_preg = ff_label(ndvi300_preg, "NDVI300 preg"),ndvi300_preg.iqr = ff_label(ndvi300_preg.iqr, "NDVI300 preg IQR"),
    ndvi100_preg = ff_label(ndvi100_preg, "NDVI100 preg"),ndvi100_preg.iqr = ff_label(ndvi100_preg.iqr, "NDVI100 preg IQR"),
    greenyn300_preg = ff_label(greenyn300_preg, "Green Access"),agebirth_m_y = ff_label(agebirth_m_y, "Maternal Age at delivery"),
    edu_m_0 = ff_label(edu_m_0, "Maternal education"),areases_tert_preg = ff_label(areases_tert_preg, "Neighbourhood socio-economic status"),
    preg_smk = ff_label(preg_smk, "Maternal smoking preg"),sex_methyl = ff_label(sex_methyl, "Child sex"),pm25_preg = ff_label(pm25_preg, "PM2.5 levels preg")
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
#12.5.2) Correlation between exposures

### Go to the Correlations directory
setwd(paste0(res.dir,"/Descriptive.info/Correlations"))

### a)Correlation between NDVI300 and NDVI100 (Continuous vs continuous)
#########################################################################

cor.ndvi300.100<-dd.EUR.M4.cc[,c("ndvi300_preg","ndvi100_preg")]
#Create plot
plot.cor.ndvi300.100<-ggpairs(cor.ndvi300.100, title="Correlation between NDVI300 and NDVI100") 

#Define a path 
file.create(paste0("M4.CorrelationNDVI300vsNDVI100.",cohort,".",ancestry,"_",Date,".pdf"))

path1.M4<-file.path(paste0("M4.CorrelationNDVI300vsNDVI100.",cohort,".",ancestry,"_",Date,".pdf"))

pdf(path1.M4)
print(plot.cor.ndvi300.100)
dev.off()

###b)Correlation between Green Access and NDVI300(Categorical vs Continuous)
################################################################################

#Create Plot
plot.cor.ndvi300.GreenAccess<-ggplot(dd.EUR.M4.cc,aes(x=greenyn300_preg, y = ndvi300_preg.iqr)) +
  ggtitle ("Correlation between Green access vs NDVI300") + geom_boxplot()
  
#Define a path 
file.create(paste0("M4.BoxPlot.GreenAccessvsNDVI300.",cohort,".",ancestry,"_",Date,".pdf"))

path2.M4<-file.path(paste0("M4.BoxPlot.GreenAccessvsNDVI300.",cohort,".",ancestry,"_",Date,".pdf"))
pdf(path2.M4)
print(plot.cor.ndvi300.GreenAccess)
dev.off()
  
###c)Correlation between Green Access and NDVI100(Categorical vs Continuous)
################################################################################

#Create Plot
plot.cor.ndvi100.GreenAccess<-ggplot(dd.EUR.M4.cc,aes(x=greenyn300_preg, y = ndvi100_preg.iqr)) +
  ggtitle ("Correlation between Green access vs NDVI100") + geom_boxplot()
  
#Define a path 
file.create(paste0("M4.BoxPlot.GreenAccessvsNDVI100.",cohort,".",ancestry,"_",Date,".pdf"))

path3.M4<-file.path(paste0("M4.BoxPlot.GreenAccessvsNDVI100.",cohort,".",ancestry,"_",Date,".pdf")) 

pdf(path3.M4)
print(plot.cor.ndvi100.GreenAccess)
dev.off()
  

###############################################################
#13.5.3) Associations between Green Spaces vars and Covariates

### Go to the GSvsCovariates directory
setwd(paste0(res.dir,"/Descriptive.info/GSvsCovariates"))

############################################################################
##a) Association beweeen Continuous Green spaces variables  and Covariables

GSexposures<-c("ndvi300_preg.iqr","ndvi100_preg.iqr")

Covariates<-c("agebirth_m_y","preg_smk","sex_methyl","edu_m_0","areases_tert_preg",
              "Bcell","CD4T","CD8T","Gran","Mono","NK","nRBC","pm25_preg")#specify variables M4

coefs.M4.cont=as.data.frame(matrix(NA,nrow=1,ncol=10))
colnames(coefs.M4.cont)=c("DateTime","GreenSpaceVar","Covariates","Model","BetaNoStand", "SE", "P","R2","R2adj","N")

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
    coefs.M4.cont[count,1]=gsub(" ","_",Sys.time())
    coefs.M4.cont[count,2]=GSexposures[j]
    coefs.M4.cont[count,3]=Covariates[i]
    coefs.M4.cont[count,4]=ff
    coefs.M4.cont[count,5:7]=res[2,c(1,2,4)]
    coefs.M4.cont[count,8]= R2
    coefs.M4.cont[count,9]=R2.adj
    coefs.M4.cont[count,10]= s$df[2]
    
    count=count+1
    
  }
}

#Multiple testing correction by FDR 

### ndvi300 pregnancy IQR vs Covariates
coefs.M4.subset.ndvi300<-subset(coefs.M4.cont, coefs.M4.cont$GreenSpaceVar == "ndvi300_preg.iqr")
coefs.M4.subset.ndvi300$Padj<-  p.adjust(coefs.M4.subset.ndvi300$P, method="fdr")

### ndvi100 pregnancy IQR vs Covariates
coefs.M4.subset.ndvi100<-subset(coefs.M4.cont, coefs.M4.cont$GreenSpaceVar == "ndvi100_preg.iqr")
coefs.M4.subset.ndvi100$Padj<-  p.adjust(coefs.M4.subset.ndvi100$P, method="fdr")

### Create one table with all results
coefs.M4.cont.all<-rbind(coefs.M4.subset.ndvi300,coefs.M4.subset.ndvi100)

write.csv(coefs.M4.cont.all,paste0("M4.GScontVars.Covariates_",cohort,".",ancestry,"_",Date,".csv"), row.names = TRUE)

#############################################################################
##b) Association beweeen Categorical Green spaces variable  and Covariables 
GSexposures<-"greenyn300_preg"

Covariates<-c("agebirth_m_y","preg_smk","sex_methyl","edu_m_0","areases_tert_preg",
              "Bcell","CD4T","CD8T","Gran","Mono","NK","nRBC","pm25_preg")#Specify variables M4 +

coefs.M4.cat=as.data.frame(matrix(NA,nrow=1,ncol=9))
colnames(coefs.M4.cat)=c("DateTime","GreenSpaceVar","Covariates","Model","BetaNoStand", "SE", "P","R2","N")

count=1
for(j in 1:length(GSexposures)) {
  print(GSexposures[j])
  results=NULL
  for(i in 1:length(Covariates)) {
    print(Covariates[i])
    ff <- paste0(GSexposures[j]," ~ ",Covariates[i],sep="") 
    fit <- glm(ff,data=dd.EUR.M4.cc,family="binomial",na.action=na.omit)
    s <-summary(fit)
    res <- s$coefficients
    R2<- with(s,1 - (deviance/null.deviance))
        
    #Save analysis details and coefficients
    coefs.M4.cat[count,1]=gsub(" ","_",Sys.time())
    coefs.M4.cat[count,2]=GSexposures[j]
    coefs.M4.cat[count,3]=Covariates[i]
    coefs.M4.cat[count,4]=ff
    coefs.M4.cat[count,5:7]=res[2,c(1,2,4)]
    coefs.M4.cat[count,8]= R2
    coefs.M4.cat[count,9]= s$df[2]
    
    count=count+1
    
  }
}

write.csv(coefs.M4.cat,paste0("M4.GScatVar.Covariates_",cohort,".",ancestry,"_",Date,".csv"), row.names = TRUE)

##############################################################################
#12.6) EWAS GREEN SPACES PREGNANCY USING ROBUST LINEAR REGRESSIONS WITH LIMMA
##############################################################################

### Before starting with the EWAS analysis...

### Go to  the EWAS results directory
setwd(paste0(res.dir,"/EWAS.Results"))

#For the analysis use ExpressionSet subset created in section 12.1.3. 

#newmset.EUR.M4

# 12.6.1) Exposure NDVI300 Pregnancy 
# 12.6.2) Exposure NDVI100 Pregnancy 
# 12.6.3) Exposure Green Access Pregnancy 

###################################
# 12.6.1) Exposure NDVI300 Pregnancy 

################
# Preg_N300_M4
################

# a) Run EWAS

### Design the model
design.Preg_N300_M4<-model.matrix(~ndvi300_preg.iqr + edu_m_0 + areases_tert_preg +  agebirth_m_y + preg_smk + sex_methyl+
                                    Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + pm25_preg, data=pData(newmset.EUR.M4))

### Run Robust linear regression model with limma
fit.Preg_N300_M4<-limma::lmFit(newmset.EUR.M4,design=design.Preg_N300_M4, method="robust")
fit.Preg_N300_M4<-limma::eBayes(fit.Preg_N300_M4)
Preg_N300_M4<-limma::topTable(fit.Preg_N300_M4, coef =2, number = Inf,sort.by="none",confint=TRUE)
Preg_N300_M4$SE <- (sqrt(fit.Preg_N300_M4$s2.post) *fit.Preg_N300_M4$stdev.unscaled)[, 2]
head(Preg_N300_M4)

###  Export results

setwd(paste0(res.dir,"/EWAS.Results"))
write.table(Preg_N300_M4, paste0("EWAS_GS_",cohort,".",ancestry,"_Preg_N300_M4_",Date,".txt"), na="NA") 
            
# b) Calculate lambda
lambda.Preg_N300_M4<- qchisq(median(Preg_N300_M4$P.Value,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
lambda.Preg_N300_M4

# c) QQ plot 

### Go to QQplots directory

### Go to  the EWAS results directory
setwd(paste0(res.dir,"/EWAS.Results/QQPlots"))

pvals<-Preg_N300_M4$P.Value
jpeg(paste0("QQPlot_",cohort,".",ancestry,"_Preg_N300_M4_",Date,".jpg"))
qq(pvals,main=paste0("QQPlot_",cohort,".",ancestry,"_Preg_N300_M4")) 
dev.off()

###################################
# 12.6.2) Exposure NDVI100 Pregnancy 
            
################
# Preg_N100_M4
################

# a) Run EWAS

### Design the model
design.Preg_N100_M4<-model.matrix(~ndvi100_preg.iqr + edu_m_0 + areases_tert_preg +  agebirth_m_y + preg_smk + sex_methyl+
                                    Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + pm25_preg, data=pData(newmset.EUR.M4))

### Run Robust linear regression model with limma
fit.Preg_N100_M4<-limma::lmFit(newmset.EUR.M4,design=design.Preg_N100_M4, method="robust")
fit.Preg_N100_M4<-limma::eBayes(fit.Preg_N100_M4)
Preg_N100_M4<-limma::topTable(fit.Preg_N100_M4, coef =2, number = Inf,sort.by="none",confint=TRUE)
Preg_N100_M4$SE <- (sqrt(fit.Preg_N100_M4$s2.post) *fit.Preg_N100_M4$stdev.unscaled)[, 2]
head(Preg_N100_M4)

###  Export results

setwd(paste0(res.dir,"/EWAS.Results"))

write.table(Preg_N100_M4, paste0("EWAS_GS_",cohort,".",ancestry,"_Preg_N100_M4_",Date,".txt"), na="NA") 
            
# b) Calculate lambda
lambda.Preg_N100_M4<- qchisq(median(Preg_N100_M4$P.Value,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
lambda.Preg_N100_M4

# c) QQ plot 

setwd(paste0(res.dir,"/EWAS.Results/QQPlots"))

pvalsM4<-Preg_N100_M4$P.Value
jpeg(paste0("QQPlot_",cohort,".",ancestry,"_Preg_N100_M4",Date,".jpg"))
qq(pvalsM4,main=paste0("QQPlot_",cohort,".",ancestry,"_Preg_N100_M4")) 
dev.off()   

#########################################
# 12.6.3) Exposure Green Access Pregnancy                                                                   
                                
################
# Preg_AcGS_M4
################

# a) Run EWAS

### Design the model

design.Preg_AcGS_M4<-model.matrix(~greenyn300_preg + edu_m_0 + areases_tert_preg +  agebirth_m_y + preg_smk + sex_methyl+ 
                                    Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + pm25_preg, data=pData(newmset.EUR.M4))

### Run Robust linear regression model with limma
fit.Preg_AcGS_M4<-limma::lmFit(newmset.EUR.M4,design=design.Preg_AcGS_M4, method="robust")
fit.Preg_AcGS_M4<-limma::eBayes(fit.Preg_AcGS_M4)
Preg_AcGS_M4<-limma::topTable(fit.Preg_AcGS_M4, coef =2, number = Inf,sort.by="none",confint=TRUE)
Preg_AcGS_M4$SE <- (sqrt(fit.Preg_AcGS_M4$s2.post) *fit.Preg_AcGS_M4$stdev.unscaled)[, 2]
head(Preg_AcGS_M4)

###  Export results
setwd(paste0(res.dir,"/EWAS.Results"))
write.table(Preg_AcGS_M4, paste0("EWAS_GS_",cohort,".",ancestry,"_Preg_AcGS_M4_",Date,".txt"), na="NA") 
            
# b) Calculate lambda
lambda.Preg_AcGS_M4<- qchisq(median(Preg_AcGS_M4$P.Value,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
lambda.Preg_AcGS_M4

# c) QQ plot 


### Go to  the EWAS results directory
setwd(paste0(res.dir,"/EWAS.Results/QQPlots"))

pvalsM4<-Preg_AcGS_M4$P.Value
jpeg(paste0("QQPlot_",cohort,".",ancestry,"_Preg_AcGS_M4_",Date,".jpg"))
qq(pvalsM4,main=paste0("QQPlot_",cohort,".",ancestry,"_Preg_AcGS_M4")) 
dev.off()
              
####################
#13.7) LAMBDAS M4
####################

### Create lambdas table M4

### Go to Lambdas directory

setwd(paste0(res.dir,"/EWAS.Results/Lambdas"))
                                                          
lambdas.table<-rbind(lambda.Preg_N300_M4,lambda.Preg_N100_M4,lambda.Preg_AcGS_M4)
colnames(lambdas.table)<-"Lambdas M4"

write.table(lambdas.table, paste0("Lambdas_",cohort,".",ancestry,"_M4_",Date,".txt"), na="NA")


###########################
#STEP 13: ANALYSES: MODEL 5
###########################

#M5:Cord blood methylation = green spaces pregnancy + maternal age + maternal smoking pregnancy + child’s sex + child’s ancestry within
# major ancestry group (optional) + batch (optional) + cohort (optional) + selection variables (optional) + maternal education + 
# neighbourhood socio-economic status+ blood cellular composition + PM2.5 + birth weight + gestational age
                                                
#IMPORTANT! Before running the analyses, please make sure that you have done STEP 8) Check variables: summaries and plots. In this way, 
# you make sure that the variables are coded as indicated in the analysis plan.

###########################################
#13.1) Select variables needed for Model 5
###########################################

dd.EUR.M5<-pData(newmset.EUR)[,c("ndvi300_preg","ndvi100_preg","greenyn300_preg","agebirth_m_y","preg_smk","sex_methyl","edu_m_0",
                                 "areases_tert_preg","Bcell","CD4T","CD8T","Gran","Mono","NK","nRBC","pm25_preg","bw_methyl","ga_methyl")]

###########################
#13.2) Complete cases M5
###########################

############################################################################
#13.2.1) Create a new dataframe with complete cases from variables in Model 5

dd.EUR.M5.cc<-dd.EUR.M5[complete.cases(dd.EUR.M5),]

##################################
#13.2.2) Complete cases M5 samples

samplescmp.M5<-rownames(dd.EUR.M5.cc)

################################################################
#13.2.3) Subset the ExpressionSet with Complete cases M5 samples

newmset.EUR.M5<-newmset.EUR[,(sampleNames(newmset.EUR) %in% samplescmp.M5)]

##########################################################################################
#13.3)Create variables with sorrounding residential greenness(NDVI Variable) transformation  
##########################################################################################

### The effect of Residential Sorrounding Greenness variable (NDVI) will be reported by IQR, thus we will transform NDVI300 pregnancy
# and NDVI100 pregnancy variables.

########################
#13.3.1) NDVI300 IQR preg

### Calculate the IQR of NDVI300
ndvi300.iqr<-iqr(dd.EUR.M5.cc$ndvi300_preg,na.rm=TRUE)

### The formula to create this variable corresponds to the original variable value divided by the interquartile value of the variable 
dd.EUR.M5.cc$ndvi300_preg.iqr<-dd.EUR.M5.cc$ndvi300_preg/ndvi300.iqr
summary(dd.EUR.M5.cc$ndvi300_preg.iqr)

#########################
#13.3.2) NDVI100 IQR preg

### Calculate the IQR of NDVI100
ndvi100.iqr<-iqr(dd.EUR.M5.cc$ndvi100_preg,na.rm=TRUE)

### The formula to create this variable corresponds to the original variable value divided by the interquartile value of the variable 
dd.EUR.M5.cc$ndvi100_preg.iqr<-dd.EUR.M5.cc$ndvi100_preg/ndvi100.iqr
summary(dd.EUR.M5.cc$ndvi100_preg.iqr)

#######################################################################################################################
#13.4)Include the new dataframe with complete cases M5 and the created NDVI variables transformed in the ExpressionSet
#######################################################################################################################

#Include the new dataframe with the complete case samples for M5 and the created NDVI variables transformed in the ExpressionSet subsetted in section 13.2.3.

###################
#13.4.1) Check order 

# Before including the new dataframe with the complete case samples for M5 and the created NDVI variables transformed in the ExpressionSet,
# check if the samples are in the same order. If you do not order samples as they are in the ExpressionSet, you could incorrectly assign the values 
# of the variables to the samples and therefore also to the methylation. 

table(ifelse(rownames(dd.EUR.M5.cc)==sampleNames(newmset.EUR.M5),"Matched","--NOT MATCHED--"))
#Matched
#    357

#In our case, samples are ordered in the same way. In case that your samples are not ordered, please, see STEP 5 for an example of how to do this. 

#############################################################
#13.4.2) Include the new dataframe in the ExpressionSet

pData(newmset.EUR.M5)<-dd.EUR.M5.cc

####################################
#13.5) Descriptive analyses for M5
###################################

#13.5.1) Descriptive tables
#13.5.2) Correlation between exposures
#13.5.3) Associations between green spaces vars and covariates


###########################
#13.5.1) Descriptive tables

### Go to the Descriptive.Tables directory
setwd(paste0(res.dir,"Descriptive.info/Descriptive.Tables"))

#a) Mean (SD)
##############

descriptive.vars.M5<-colnames(dd.EUR.M5.cc)#Select variables for the descriptives

dd.EUR.M5.cc %>%
  dplyr::mutate(
    ndvi300_preg = ff_label(ndvi300_preg, "NDVI300 preg"),ndvi300_preg.iqr = ff_label(ndvi300_preg.iqr, "NDVI300 preg IQR"),
    ndvi100_preg = ff_label(ndvi100_preg, "NDVI100 preg"),ndvi100_preg.iqr = ff_label(ndvi100_preg.iqr, "NDVI100 preg IQR"),
    greenyn300_preg = ff_label(greenyn300_preg, "Green Access"),agebirth_m_y = ff_label(agebirth_m_y, "Maternal Age at delivery"),
    edu_m_0 = ff_label(edu_m_0, "Maternal education"),areases_tert_preg = ff_label(areases_tert_preg, "Neighbourhood socio-economic status"),
    preg_smk = ff_label(preg_smk, "Maternal smoking preg"),sex_methyl = ff_label(sex_methyl, "Child sex"),pm25_preg = ff_label(pm25_preg, "PM2.5 levels preg"),
    bw_methyl = ff_label(bw_methyl, "Birth weight"),ga_methyl = ff_label(ga_methyl, "Gestational age")
    )%>%
  summary_factorlist(dependent=NULL,descriptive.vars.M5,cont="mean",column=TRUE,na_include=TRUE,total_col=TRUE)->DescriptiveTable.Mean.M5

DescriptiveTable.Mean.M5<-DescriptiveTable.Mean.M5[,c("label","levels","Total")]
colnames(DescriptiveTable.Mean.M5)<-c("Variables","levels","Mean (SD) or n(%)")

#Median (IQR)
#############

dd.EUR.M5.cc %>%
  dplyr::mutate(
    ndvi300_preg = ff_label(ndvi300_preg, "NDVI300 preg"),ndvi300_preg.iqr = ff_label(ndvi300_preg.iqr, "NDVI300 preg IQR"),
    ndvi100_preg = ff_label(ndvi100_preg, "NDVI100 preg"),ndvi100_preg.iqr = ff_label(ndvi100_preg.iqr, "NDVI100 preg IQR"),
    greenyn300_preg = ff_label(greenyn300_preg, "Green Access"),agebirth_m_y = ff_label(agebirth_m_y, "Maternal Age at delivery"),
    edu_m_0 = ff_label(edu_m_0, "Maternal education"),areases_tert_preg = ff_label(areases_tert_preg, "Neighbourhood socio-economic status"),
    preg_smk = ff_label(preg_smk, "Maternal smoking preg"),sex_methyl = ff_label(sex_methyl, "Child sex"),pm25_preg = ff_label(pm25_preg, "PM2.5 levels preg"),
    bw_methyl = ff_label(bw_methyl, "Birth weight"),ga_methyl = ff_label(ga_methyl, "Gestational age")
    )%>%
  summary_factorlist(dependent=NULL,descriptive.vars.M5,cont="median",column=TRUE,na_include=TRUE,total_col=TRUE)->DescriptiveTable.Median.M5

DescriptiveTable.Median.M5<-DescriptiveTable.Median.M5[,c("label","levels","Total")]
colnames(DescriptiveTable.Median.M5)<-c("Variables","levels","Mean (SD) or n(%)")

### Create a single descriptive table with all variables

DescriptiveTable.final.M5<-cbind(DescriptiveTable.Mean.M5,DescriptiveTable.Median.M5)
DescriptiveTable.final.M5<-DescriptiveTable.final.M5[,c(1:3,5:6)]
             
write.table(
  DescriptiveTable.final.M5, file=paste0("M5.Descriptive.table.",cohort,".",ancestry,".cc.",Date,".txt"),col.names=T, row.names=F, quote=F, sep="\t")

######################################
#13.5.2) Correlation between exposures

### Go to the Correlations directory
setwd(paste0(res.dir,"Descriptive.info/Correlations"))

### a)Correlation between NDVI300 and NDVI100 (Continuous vs continuous)
#########################################################################

cor.ndvi300.100<-dd.EUR.M5.cc[,c("ndvi300_preg","ndvi100_preg")]
#Create plot
plot.cor.ndvi300.100<-ggpairs(cor.ndvi300.100, title="Correlation between NDVI300 and NDVI100") 

#Define a path 
file.create(paste0("M5.CorrelationNDVI300vsNDVI100.",cohort,".",ancestry,"_",Date,".pdf"))

path1.M5<-file.path(paste0("M5.CorrelationNDVI300vsNDVI100.",cohort,".",ancestry,"_",Date,".pdf"))

pdf(path1.M5)
print(plot.cor.ndvi300.100)
dev.off()

###b)Correlation between Green Access and NDVI300(Categorical vs Continuous)
################################################################################

#Create Plot
plot.cor.ndvi300.GreenAccess<-ggplot(dd.EUR.M5.cc,aes(x=greenyn300_preg, y = ndvi300_preg.iqr)) + ggtitle ("Correlation between Green access vs NDVI300") +
  geom_boxplot()
  
#Define a path 
file.create(paste0("M5.BoxPlot.GreenAccessvsNDVI300.",cohort,".",ancestry,"_",Date,".pdf"))

path2.M5<-file.path(paste0("M5.BoxPlot.GreenAccessvsNDVI300.",cohort,".",ancestry,"_",Date,".pdf"))

pdf(path2.M5)
print(plot.cor.ndvi300.GreenAccess)
dev.off()
  
###c)Correlation between Green Access and NDVI100(Categorical vs Continuous)
################################################################################

#Create Plot
plot.cor.ndvi100.GreenAccess<-ggplot(dd.EUR.M5.cc,aes(x=greenyn300_preg, y = ndvi100_preg.iqr)) + ggtitle ("Correlation between Green access vs NDVI100") +
  geom_boxplot()
  
#Define a path 
file.create(paste0("M5.BoxPlot.GreenAccessvsNDVI100.",cohort,".",ancestry,"_",Date,".pdf"))

path3.M5<-file.path(paste0("M5.BoxPlot.GreenAccessvsNDVI100.",cohort,".",ancestry,"_",Date,".pdf"))

pdf(path3.M5)
print(plot.cor.ndvi100.GreenAccess)
dev.off()
  

###############################################################
#13.5.3) Associations between Green Spaces vars and Covariates

### Go to the GSvsCovariates directory
setwd(paste0(res.dir,"Descriptive.info/GSvsCovariates"))


############################################################################
##a) Association beweeen Continuous Green spaces variables  and Covariables

GSexposures<-c("ndvi300_preg.iqr","ndvi100_preg.iqr")

Covariates<-c("agebirth_m_y","preg_smk","sex_methyl","edu_m_0","areases_tert_preg",
              "Bcell","CD4T","CD8T","Gran","Mono","NK","nRBC","pm25_preg","bw_methyl","ga_methyl")#Specify variables M5 

coefs.M5.cont=as.data.frame(matrix(NA,nrow=1,ncol=10))
colnames(coefs.M5.cont)=c("DateTime","GreenSpaceVar","Covariates","Model","BetaNoStand", "SE", "P","R2","R2adj","N")

count=1
for(j in 1:length(GSexposures)) {
  print(GSexposures[j])
  results=NULL
  for(i in 1:length(Covariates)) {
    print(Covariates[i])
    ff <- paste0(GSexposures[j],"~",Covariates[i],sep="") 
    fit <- lm(ff,data=dd.EUR.M5.cc)
    s <-summary(fit)
    res <- s$coefficients
    R2<-s$r.squared
    R2.adj<- s$adj.r.squared
    
    #Save analysis details and coefficients
    coefs.M5.cont[count,1]=gsub(" ","_",Sys.time())
    coefs.M5.cont[count,2]=GSexposures[j]
    coefs.M5.cont[count,3]=Covariates[i]
    coefs.M5.cont[count,4]=ff
    coefs.M5.cont[count,5:7]=res[2,c(1,2,4)]
    coefs.M5.cont[count,8]= R2
    coefs.M5.cont[count,9]=R2.adj
    coefs.M5.cont[count,10]= s$df[2]
    
    count=count+1
    
  }
}

#Multiple testing correction by FDR 

### ndvi300 pregnancy IQR vs Covariates
coefs.M5.subset.ndvi300<-subset(coefs.M5.cont, coefs.M5.cont$GreenSpaceVar == "ndvi300_preg.iqr")
coefs.M5.subset.ndvi300$Padj<-  p.adjust(coefs.M5.subset.ndvi300$P, method="fdr")

### ndvi100 pregnancy IQR vs Covariates
coefs.M5.subset.ndvi100<-subset(coefs.M5.cont, coefs.M5.cont$GreenSpaceVar == "ndvi100_preg.iqr")
coefs.M5.subset.ndvi100$Padj<-  p.adjust(coefs.M5.subset.ndvi100$P, method="fdr")

### Create one table with all results
coefs.M5.cont.all<-rbind(coefs.M5.subset.ndvi300,coefs.M5.subset.ndvi100)

write.csv(coefs.M5.cont.all,paste0("M5.GScontVars.Covariates_",cohort,".",ancestry,"_",Date,".csv"), row.names = TRUE)

#############################################################################
##b) Association beweeen Categorical Green spaces variable  and Covariables 
GSexposures<-"greenyn300_preg"

Covariates<-c("agebirth_m_y","preg_smk","sex_methyl","edu_m_0","areases_tert_preg",
              "Bcell","CD4T","CD8T","Gran","Mono","NK","nRBC","pm25_preg","bw_methyl","ga_methyl")#Specify variables M5  

coefs.M5.cat=as.data.frame(matrix(NA,nrow=1,ncol=9))
colnames(coefs.M5.cat)=c("DateTime","GreenSpaceVar","Covariates","Model","BetaNoStand", "SE", "P","R2","N")

count=1
for(j in 1:length(GSexposures)) {
  print(GSexposures[j])
  results=NULL
  for(i in 1:length(Covariates)) {
    print(Covariates[i])
    ff <- paste0(GSexposures[j]," ~ ",Covariates[i],sep="") 
    fit <- glm(ff,data=dd.EUR.M5.cc,family="binomial",na.action=na.omit)
    s <-summary(fit)
    res <- s$coefficients
    R2<- with(s,1 - (deviance/null.deviance))
        
    #Save analysis details and coefficients
    coefs.M5.cat[count,1]=gsub(" ","_",Sys.time())
    coefs.M5.cat[count,2]=GSexposures[j]
    coefs.M5.cat[count,3]=Covariates[i]
    coefs.M5.cat[count,4]=ff
    coefs.M5.cat[count,5:7]=res[2,c(1,2,4)]
    coefs.M5.cat[count,8]= R2
    coefs.M5.cat[count,9]= s$df[2]
    
    count=count+1
    
  }
}

write.csv(coefs.M5.cat,paste0("M5.GScatVar.Covariates_",cohort,".",ancestry,"_",Date,".csv"), row.names = TRUE)

##############################################################################
#13.6) EWAS GREEN SPACES PREGNANCY USING ROBUST LINEAR REGRESSIONS WITH LIMMA
##############################################################################

### Before starting with the EWAS analysis...

### Go to the EWAS results directory
setwd(paste0(res.dir,"/EWAS.Results"))

#For the analysis use ExpressionSet subset created in section 13.1.3. 

#newmset.EUR.M5

# 13.6.1) Exposure NDVI300 Pregnancy 
# 13.6.2) Exposure NDVI100 Pregnancy 
# 13.6.3) Exposure Green Access Pregnancy 

###################################
# 13.6.1) Exposure NDVI300 Pregnancy 

################
# Preg_N300_M5
################

# a) Run EWAS

### Design the model
design.Preg_N300_M5<-model.matrix(~ndvi300_preg.iqr +  agebirth_m_y + preg_smk + sex_methyl + edu_m_0 + areases_tert_preg +
                                    Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + pm25_preg+bw_methyl + ga_methyl, data=pData(newmset.EUR.M5))

### Run Robust linear regression model with limma
fit.Preg_N300_M5<-limma::lmFit(newmset.EUR.M5,design=design.Preg_N300_M5, method="robust")
fit.Preg_N300_M5<-limma::eBayes(fit.Preg_N300_M5)
Preg_N300_M5<-limma::topTable(fit.Preg_N300_M5, coef =2, number = Inf,sort.by="none",confint=TRUE)
Preg_N300_M5$SE <- (sqrt(fit.Preg_N300_M5$s2.post) *fit.Preg_N300_M5$stdev.unscaled)[, 2]
head(Preg_N300_M5)

###  Export results
setwd(paste0(res.dir,"/EWAS.Results"))
write.table(Preg_N300_M5, paste0("EWAS_GS_",cohort,".",ancestry,"_Preg_N300_M5_",Date,".txt"), na="NA") 
            
# b) Calculate lambda
lambda.Preg_N300_M5<- qchisq(median(Preg_N300_M5$P.Value,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
lambda.Preg_N300_M5

# c) QQ plot 

### Go to QQplots directory
setwd(paste0(res.dir,"/EWAS.Results/QQPlots"))

pvals<-Preg_N300_M5$P.Value
jpeg(paste0("QQPlot_",cohort,".",ancestry,"_Preg_N300_M5_",Date,".jpg"))
qq(pvals,main=paste0("QQPlot_",cohort,".",ancestry,"_Preg_N300_M5")) 
dev.off()

###################################
# 13.6.2) Exposure NDVI100 Pregnancy 
            
################
# Preg_N100_M5
################

# a) Run EWAS

### Design the model
design.Preg_N100_M5<-model.matrix(~ndvi100_preg.iqr  +  agebirth_m_y + preg_smk + sex_methyl+ edu_m_0 + areases_tert_preg + 
                                    Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + pm25_preg+bw_methyl + ga_methyl, data=pData(newmset.EUR.M5))

### Run Robust linear regression model with limma
fit.Preg_N100_M5<-limma::lmFit(newmset.EUR.M5,design=design.Preg_N100_M5, method="robust")
fit.Preg_N100_M5<-limma::eBayes(fit.Preg_N100_M5)
Preg_N100_M5<-limma::topTable(fit.Preg_N100_M5, coef =2, number = Inf,sort.by="none",confint=TRUE)
Preg_N100_M5$SE <- (sqrt(fit.Preg_N100_M5$s2.post) *fit.Preg_N100_M5$stdev.unscaled)[, 2]
head(Preg_N100_M5)

###  Export results

setwd(paste0(res.dir,"/EWAS.Results"))
write.table(Preg_N100_M5, paste0("EWAS_GS_",cohort,".",ancestry,"_Preg_N100_M5_",Date,".txt"), na="NA") 
            
# b) Calculate lambda
lambda.Preg_N100_M5<- qchisq(median(Preg_N100_M5$P.Value,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
lambda.Preg_N100_M5

# c) QQ plot 

setwd(paste0(res.dir,"/EWAS.Results/QQPlots"))

pvalsM5<-Preg_N100_M5$P.Value
jpeg(paste0("QQPlot_",cohort,".",ancestry,"_Preg_N100_M5_",Date,".jpg"))
qq(pvalsM5,main=paste0("QQPlot_",cohort,".",ancestry,"_Preg_N100_M5")) 
dev.off()   

#########################################
# 13.6.3) Exposure Green Access Pregnancy                                                                   
                                
################
# Preg_AcGS_M5
################

# a) Run EWAS

### Design the model

design.Preg_AcGS_M5<-model.matrix(~greenyn300_preg +  agebirth_m_y + preg_smk + sex_methyl + edu_m_0 + areases_tert_preg +
                                    Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + pm25_preg+bw_methyl + ga_methyl, data=pData(newmset.EUR.M5))

### Run Robust linear regression model with limma
fit.Preg_AcGS_M5<-limma::lmFit(newmset.EUR.M5,design=design.Preg_AcGS_M5, method="robust")
fit.Preg_AcGS_M5<-limma::eBayes(fit.Preg_AcGS_M5)
Preg_AcGS_M5<-limma::topTable(fit.Preg_AcGS_M5, coef =2, number = Inf,sort.by="none",confint=TRUE)
Preg_AcGS_M5$SE <- (sqrt(fit.Preg_AcGS_M5$s2.post) *fit.Preg_AcGS_M5$stdev.unscaled)[, 2]
head(Preg_AcGS_M5)

###  Export results

setwd(paste0(res.dir,"/EWAS.Results"))

write.table(Preg_AcGS_M5, paste0("EWAS_GS_",cohort,".",ancestry,"_Preg_AcGS_M5_",Date,".txt"), na="NA") 
            
# b) Calculate lambda
lambda.Preg_AcGS_M5<- qchisq(median(Preg_AcGS_M5$P.Value,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
lambda.Preg_AcGS_M5

# c) QQ plot 

setwd(paste0(res.dir,"/EWAS.Results/QQPlots"))

pvalsM5<-Preg_AcGS_M5$P.Value
jpeg(paste0("QQPlot_",cohort,".",ancestry,"_Preg_AcGS_M5_",Date,".jpg"))
qq(pvalsM5,main=paste0("QQPlot_",cohort,".",ancestry,"_Preg_AcGS_M5")) 
dev.off()
              
####################
#13.7) LAMBDAS M5
####################

### Create lambdas table M5

setwd(paste0(res.dir,"/EWAS.Results/Lambdas"))
                                                            
lambdas.table<-rbind(lambda.Preg_N300_M5,lambda.Preg_N100_M5,lambda.Preg_AcGS_M5)
colnames(lambdas.table)<-"Lambdas M5"

write.table(lambdas.table, paste0("Lambdas_",cohort,".",ancestry,"_M5_",Date,".txt"), na="NA")



#################
#IMPORTANT
#################

#Repeat from STEP 8 to STEP 13 for each of the major ancestry groups in your cohort with more than >100 samples.

############################
#THANK YOU!
############################

