######################################################################
#LC_ATH_EWAS_GS_DatashieldCOde
#Update: 12.07.2022
######################################################################

##########################################################################################
### STEP 0: Before starting 
##########################################################################################
##########################################
# Exclusion Criteria
### The following children will be excluded from the study:
# - Twins will be excluded. 
# - For non-twin siblings, cohorts will include only one child per mother, based on completeness of data and, if equal, randomly.

###########################################

###########################################
### Change to your own working directory
setwd("~/DOCTORADO/TERCER_AÑO/EWAS_GREENSPACES_blood/DataSHIELD/Results_18112022")

###########################################
### Set initial parameters
cohort<-"INMA" #define cohort name Depends in the analysis taht I am doing  
Date<- Sys.Date()#define date 

#################################
#################################
#PART I: Prepare ExpressionSets
#################################
#################################

##########################################################################################
### STEP 1: Install, Load required packages and connect TO DataShield Server
##########################################################################################
#Install devtools
 install.packages("devtools")
 install.packages("DSOpal")
 install.packages("opalr")
 install.packages('dsBaseClient', repos=c(getOption('repos'), 'http://cran.obiba.org'), dependencies=TRUE)
 
devtools::install_github("isglobal-brge/dsOmicsClient",force=TRUE)
 
 install.packages("remotes")
 library(remotes)
 install_github("lifecycle-project/ds-helper")
 install.packages("DSMolgenisArmadillo")


#Load packages
library(opalr)
library(DSOpal)
library(dsBaseClient)
library(dsOmicsClient)
library(DSMolgenisArmadillo)
library(dsHelper)
library(dplyr)

#################
#INMA SERVER
#################

#Log in Opal

o <- opalr::opal.login(username = "username", password = '',
                       url = "https://opal.isglobal.org/repo")
opal.projects(o)

builder <- DSI::newDSLoginBuilder()

builder$append(server = 'INMA_cohort', url = "https://opal.isglobal.org/repo", 
               user = 'lifecycle_saguilar', password = '', profile = "rock-inma")


logindata <- builder$build()
conns <- DSI::datashield.login(logins = logindata, assign = TRUE, 
                               symbol = 'res')
datashield.tables(conns)
DSI::datashield.logout(conns)

##########################################################################################
### STEP 2: Prepare and load the methylation data and the metadata
##########################################################################################

########################################
### 2.1) Methylation data and Metadata
#########################################


#EXPRESSIONSET 0Y 
datashield.assign.resource(conns, symbol = "eSet_0y",
                           resource = list(INMA_cohort = "Athlete_Resources.INMA_Methyl_Child_Blood_0y_450K"))
ds.class('eSet_0y',datasources = conns)
datashield.assign.expr(conns, symbol = "methy_0y",
                       expr = quote(as.resource.object(eSet_0y)))

#Explore Eset
ds.class('methy_0y',datasources=conns) #class
ds.nFeatures("methy_0y",datasources = conns) #number of CpGs
ds.featureNames("methy_0y",datasources = conns) # name of CpGs
ds.nSamples("methy_0y",datasources = conns) # number of samples
# $INMA_cohort
# Samples 
#385 
# 
# attr(,"class")
# [1] "dsnSamples" "list"
ds.varLabels("methy_0y",datasources = conns) # column names (variables)
# $INMA_cohort
# [1] "id_methyl"     "id"            "sex_methyl"    "age_methyl"    "anc_methyl"    "cohort_methyl" "bw_methyl"    
# [8] "ga_methyl"     "Bcell"         "CD4T"          "CD8T"          "Gran"          "Mono"          "NK"           
# [15] "nRBC"     
##########################################
### STEP 3: Major ancestry group subset
##########################################

### Define group
#As ethnic groups are tested separatelly (if >100 samples), please indicate the major ancestry group you are working with

ancestry<-"EUR" 

#Subset the ExpressionSet by the major ancestry group (>100 samples)

ds.subsetExpressionSet(
  eSet='methy_0y',
  objective_variable ="anc_methyl",
  objective_value = "1",
  complete_cases = FALSE,
  newobj.name = 'methy_0y_EUR',
  datasources = conns)

ds.class('methy_0y_EUR',datasources=conns)
ds.varLabels('methy_0y_EUR',datasources = conns)
ds.nFeatures('methy_0y_EUR',datasources = conns)
# $INMA_cohort
# Features 
# 476946 

ds.nSamples('methy_0y_EUR',datasources = conns)
# $INMA_cohort
# Samples 
# 361

##########################################################################################
### STEP 4: WINSORIZE OUTLIERS
##########################################################################################
### Winsorize methylation beta values to remove potential outliers
### The functions are obtained from the ewaff package and edited to fit our needs (winsorize.pct = 0.005)

ds.removeOutliers('methy_0y_EUR', pct = 0.005, newobj.name = "methy_0y_EUR.winz", datasources =conns)

ds.class('methy_0y_EUR.winz',datasources=conns)
ds.varLabels('methy_0y_EUR.winz',datasources = conns)
ds.nFeatures('methy_0y_EUR.winz',datasources = conns)
ds.nSamples('methy_0y_EUR.winz',datasources = conns)
# $INMA_cohort
# Samples 
# 361 
# 
# attr(,"class")
# [1] "dsnSamples" "list" 
##########################################################################################
### STEP 5: Load the dataframe with the exposures + additional covariates 
##########################################################################################

#EXPOSURE AND COVARIABLES. Connect with the tables from opal or armadillo


#EXPOSURE AND COVARIABLES. Connect with the tables from opal
datashield.assign.table(conns, symbol = "Exposures_norep",
                        table = list(INMA_cohort = "lc_isglobal_core_2_1.2_1_core_1_1_non_rep_saguilar"))
datashield.assign.table(conns, symbol = "Exposures_yearly",
                        table = list(INMA_cohort = "lc_isglobal_core_2_1.2_1_core_1_1_yearly_rep_saguilar"))

#Explore the tables

#Dimension
ds.dim('Exposures_norep',datasources=conns)
ds.dim('Exposures_yearly',datasources=conns)

#Name of variables
ds.colnames('Exposures_norep',datasources=conns)

ds.colnames('Exposures_yearly',datasources=conns)


#Select from the nonrep dataframe the variables of interest 

# "child_id"            "cohort_id"           "coh_country" , "areases_tert_preg", "greenyn300_preg","ndvi100_preg"        "ndvi300_preg"
# "preg_smk", "agebirth_m_y","pm25_preg"

ds.dataFrameSubset(df.name = "Exposures_norep",
                   keep.cols = c(1:3,6,8,11,12,14,15,29,31,42), 
                   rm.cols = NULL,
                   keep.NAs = TRUE,
                   newobj = "Exposures_norep_new",
                   datasources = conns, 
                   notify.of.progress = FALSE)


ds.colnames("Exposures_norep_new",datasources = conns)
# $INMA_cohort
# [1] "child_id"          "cohort_id"         "coh_country"       "areases_tert_preg" "greenyn300_preg"   "ndvi100_preg"     
# [7] "ndvi300_preg"      "recruit_age"       "sex"               "preg_smk"          "agebirth_m_y"      "pm25_preg"
#Select from the yearly dataframe the variables of interest 

#For Green spaces preg - Cordblood

ds.colnames("Exposures_yearly",datasources = conns)

ds.dataFrameSubset(df.name = "Exposures_yearly",
                   keep.cols =  c(1,3,22), 
                   keep.NAs = TRUE,
                   newobj = "Exposures_yearly_new_Preg",
                   datasources = conns, #all servers are used
                   notify.of.progress = FALSE)

ds.colnames("Exposures_yearly_new_Preg",datasources = conns)
# $INMA_cohort
# [1] "age_years" "child_id"  "edu_m_" 

#Create table with variables per year
ds.reShape(data.name = "Exposures_yearly_new_Preg",
           timevar.name = "age_years",
           idvar.name = "child_id",
           direction = "wide",
           newobj = "Exposures_yearly_wide_Preg",
           datasources = conns)

ds.colnames('Exposures_yearly_wide_Preg',datasources = conns)
# $INMA_cohort
# [1] "child_id"  "edu_m_.5"  "edu_m_.3"  "edu_m_.7"  "edu_m_.8"  "edu_m_.0"  "edu_m_.2"  "edu_m_.6"  "edu_m_.11" "edu_m_.1" 
# [11] "edu_m_.9"  "edu_m_.12" "edu_m_.10" "edu_m_.4" 

ds.dataFrameSubset(df.name = "Exposures_yearly_wide_Preg",
                   keep.cols = c(1,6), 
                   keep.NAs = TRUE,
                   newobj = "Exposures_yearly_def",
                   datasources = conns, #all servers are used
                   notify.of.progress = FALSE)

ds.dim('Exposures_yearly_def',datasources = conns)
ds.colnames('Exposures_yearly_def',datasources = conns)
#[1] "child_id" "edu_m_.0"
#Merge to have a dataframe with all vars

ds.dim("Exposures_norep_new",datasources = conns)
ds.dim("Exposures_yearly_def",datasources = conns)

ds.merge(
  x.name ="Exposures_norep_new",
  y.name ="Exposures_yearly_def",
  by.x.names = "child_id",
  by.y.names = "child_id",
  all.x = TRUE,
  all.y = TRUE,
  newobj = "df.exp.cov",
  datasources = conns
)
ds.dim("df.exp.cov",datasources = conns)
# $`dimensions of df.exp.cov in INMA_cohort`
# [1] 2270   13

ds.colnames("df.exp.cov",datasources = conns)
# $INMA_cohort
# [1] "child_id"          "cohort_id"         "coh_country"       "areases_tert_preg" "greenyn300_preg"   "ndvi100_preg"     
# [7] "ndvi300_preg"      "recruit_age"       "sex"               "preg_smk"          "agebirth_m_y"      "pm25_preg"        
# [13] "edu_m_.0"    

###############################################################################################################
### STEP 6: Merge the metadata from the ExpressionSet and the dataframe with exposures + additional variables
###############################################################################################################
#Merge Exposure Variables and covariates inside the DataFrame of the eSet

#Add dataframe with exposures + additional variables
ds.addPhenoData2eSet('methy_0y_EUR.winz', 'df.exp.cov', identifier = "child_id", alternate_eset_id = "id",  newobj.name = 'methy_0y_allVars',datasources = conns)
ds.class('methy_0y_allVars',datasources=conns)
ds.varLabels('methy_0y_allVars',datasources = conns)
#Aggregated (varLabelsDS(methy_0y_allVars)) [===========================================] 100% / 0s
# $INMA_cohort
# [1] "id_methyl"         "id"                "sex_methyl"        "age_methyl"        "anc_methyl"        "cohort_methyl"     "bw_methyl"        
# [8] "ga_methyl"         "Bcell"             "CD4T"              "CD8T"              "Gran"              "Mono"              "NK"               
# [15] "nRBC"              "cohort_id"         "coh_country"       "areases_tert_preg" "greenyn300_preg"   "ndvi100_preg"      "ndvi300_preg"     
# [22] "recruit_age"       "sex"               "preg_smk"          "agebirth_m_y"      "pm25_preg"         "edu_m_.0"       

ds.nFeatures('methy_0y_allVars',datasources = conns)
# $INMA_cohort
# Features 
# 476946 
ds.nSamples('methy_0y_allVars',datasources = conns)
# $INMA_cohort
# Samples 
# 361 


#################################
#################################
#PART II: Analyses
#################################
#################################

cohort<-"INMA"
ancestry<-"EUR"
Date<- Sys.Date()

### Now, choose the working directory where you want to create the new directory to save all results for European Samples
setwd("~/DOCTORADO/TERCER_AÑO/EWAS_GREENSPACES_blood/DataSHIELD/Results_18112022")
#Change to the directory where you want to save results

### Create a new directory
dir.create("EUR")

### Create an object with the path of this new directory
res.dir<- "C:/Users/saguilar/Documents/DOCTORADO/TERCER_AÑO/EWAS_GREENSPACES_blood/DataSHIELD/Results_18112022/EUR" 

### Go to the new directory
setwd(res.dir)

#######################################################
# STEP 7. Check variables: summaries and plots 
#######################################################

##################
# 7.1) Exposures
#################

ds.pData('methy_0y_allVars','pheno',datasources = conns)

ls#######################################################################
### Residential Sourrounding greenness: NDVI 300  pregnancy (numeric)
ds.summary("pheno$ndvi300_preg",datasources = conns)

### Histogram
png("Histogram_INMA.EUR_ndvi300_preg.png", width=800, height=800)
ds.histogram(x="pheno$ndvi300_preg", datasources = conns)# editar
dev.off() 


########################################################################
### Residential Sourrounding greenness: NDVI 100 pregnancy (numeric) 
ds.summary("pheno$ndvi100_preg",datasources = conns)

### Histogram
png("Histogram_INMA.EUR_ndvi100_preg.png", width=800, height=800)
ds.histogram(x="pheno$ndvi100_preg", datasources = conns)# editar
dev.off() 

###############################################################
### Green access pregnancy(Greenyn 300 pregnancy) (Categorical)
# Two categories:
# 0=no
# 1=yes
ds.summary("pheno$greenyn300_preg",datasources = conns)


# #########################
# 8.2) Covariates
#########################

###############################################################################
### Maternal education (Categorical). 3 Categories: 1 = high;2 = Medium;3= Low
ds.summary("pheno$edu_m_.0",datasources = conns)
# 

####################################################################
### Neighbourhood socio-economic status (Categorical).3 categories:
#1= 1st tertile (low deprivated);2= 2nd tertile (medium deprivated);3= 3rd tertile (high deprivated)
ds.summary("pheno$areases_tert_preg",datasources = conns)

################################################
### Maternal age at delivery in years (numeric) 
ds.summary("pheno$agebirth_m_y",datasources = conns)

#####################################################################################################
### Maternal smoking pregnancy (Categorical).ANY smoking. Two categories:0= No; 1= Yes
ds.summary("pheno$preg_smk",datasources = conns)


###################################################################################
### Child sex (Categorical). Two categories: 1 = Male; 2 = Female
ds.summary("pheno$sex_methyl",datasources = conns)
# 


############################
# Pollution variables
############################

###################
### PM25 (Numeric)
ds.summary("pheno$pm25_preg",datasources = conns)

### Histogram
jpeg(paste0("Histogram_",cohort,".",ancestry,"_PM25_preg_",Date,".jpg"))
ds.histogram(x="pheno$pm25_preg", datasources = conns)# editar
dev.off()

###########################
# Reproductive variables
###########################

####################################
### Gestational age in days (numeric)
ds.summary("pheno$ga_methyl",datasources = conns)


### Histogram
jpeg(paste0("Histogram_",cohort,".",ancestry,"_GestationalAge_",Date,".jpg"))
ds.histogram(x="pheno$ga_methyl", datasources = conns)# editar
dev.off()

############################################
### Birth weight at birth in grams(numeric)
ds.summary("pheno$bw_methyl",datasources = conns)

### Histogram
jpeg(paste0("Histogram_",cohort,".",ancestry,"_BirthWeight_",Date,".jpg"))
ds.histogram(x="pheno$bw_methyl", datasources = conns)# editar
dev.off()

############
# Celltypes
############

#########
### Bcell
ds.summary("pheno$Bcell",datasources = conns)

#########
### CD4T
ds.summary("pheno$CD4T",datasources = conns)
# 

#########
### CD8T
ds.summary("pheno$CD8T",datasources = conns)

#########
### Gran
ds.summary("pheno$Gran",datasources = conns)

#########
### Mono
ds.summary("pheno$Mono",datasources = conns)

#########
### NK
ds.summary("pheno$NK",datasources = conns)


#########
### nRBC
ds.summary("pheno$nRBC",datasources = conns)



### Now we have checked that variables in the ExpressionSet are coded correctly.  

#EXAMPLE MAIN MODEL 
################################################################################################################
#ESET M3
#Select variables needed for model 3
#M3: Cord blood methylation ~ Green Spaces pregnancy + maternal age + maternal smoking pregnancy + child's sex + maternal education + SES neigh + celltype proportions

# First of all, I extract the phenodata with all the necessary variables from the expressionSet created above. This way, I have the subset of European individuals with
# methylation + variables of interest.

ds.pData("methy_0y_allVars", 'pheno.M3')

ds.colnames("pheno.M3",datasources = conns)


# 
#I create a subset of the variables I need for M3. I do not select the cell types because they are already in the Eset.
ds.dataFrameSubset(df.name = "pheno.M3",
                   keep.cols = c(1:6,18:21,24:25,27), 
                   keep.NAs = TRUE,
                   newobj = "dd.M3",
                   datasources = conns, #all servers are used
                   notify.of.progress = FALSE)

ds.dim("dd.M3",datasources = conns)
ds.colnames("dd.M3",datasources = conns)
# Aggregated (exists("dd.M3")) [=========================================================] 100% / 1s
# Aggregated (classDS("dd.M3")) [========================================================] 100% / 0s
# Aggregated (colnamesDS("dd.M3")) [=====================================================] 100% / 0s
# $INMA_cohort
# [1] "id_methyl"         "id"                "sex_methyl"        "age_methyl"        "anc_methyl"        "cohort_methyl"     "areases_tert_preg" "greenyn300_preg"   "ndvi100_preg"      "ndvi300_preg"     
# [11] "preg_smk"          "agebirth_m_y"      "edu_m_.0" 
#Now I do the complete cases of variables of model M3
ds.completeCases(x1 = "dd.M3", newobj ="dd.M3" , datasources = conns)
ds.dim("dd.M3",datasources = conns)
# $`dimensions of dd.M3 in INMA_cohort`
# [1] 357  13
# 
# $`dimensions of dd.M3 in combined studies`
# [1] 357  13

ds.colnames("dd.M3",datasources = conns)
# $INMA_cohort
# [1] "id_methyl"         "id"                "sex_methyl"        "age_methyl"        "anc_methyl"        "cohort_methyl"     "areases_tert_preg" "greenyn300_preg"   "ndvi100_preg"      "ndvi300_preg"     
# [11] "preg_smk"          "agebirth_m_y"      "edu_m_.0"         

#I am going to remove those variables that are already in the expressionSet from the beginning; otherwise, I get an error when merging because there are variables with the same name.

ds.dataFrameSubset(df.name = "dd.M3",
                   keep.cols = c(1,7:13), 
                   keep.NAs = TRUE,
                   newobj = "dd.M3",
                   datasources = conns, #all servers are used
                   notify.of.progress = FALSE)

ds.dim("dd.M3",datasources = conns)
ds.colnames("dd.M3",datasources = conns)
# $INMA_cohort
# [1] "id_methyl"         "areases_tert_preg" "greenyn300_preg"   "ndvi100_preg"      "ndvi300_preg"      "preg_smk"          "agebirth_m_y"      "edu_m_.0"  

#Now, I'm going to calculate NDVI IQR.

#Calculate iqr ndvi 300

# OBTAIN THE QUANTILES
quants.ndvi300 <- ds.quantileMean("dd.M3$ndvi300_preg")

#  CALCULATE THE IQR (25-75 IN THIS CASE)
iqr.ndvi300 <- quants.ndvi300[5] - quants.ndvi300[3] #0.1035744

ds.make(toAssign = "dd.M3$ndvi300_preg /0.1035744", newobj = 'ndvi300_preg.iqr')
ds.dataFrame(x = c("dd.M3", "ndvi300_preg.iqr"), newobj = "dd.M3")
ds.colnames("dd.M3",datasources = conns)
ds.dim("dd.M3",datasources = conns)

#Calculate iqr ndvi 100

# OBTAIN THE QUANTILES
quants.ndvi100 <- ds.quantileMean("dd.M3$ndvi100_preg")

#  CALCULATE THE IQR (25-75 IN THIS CASE)
iqr.ndvi100 <- quants.ndvi100[5] - quants.ndvi100[3] #0.0716195

ds.make(toAssign = "dd.M3$ndvi100_preg /0.0716195 ", newobj = 'ndvi100_preg.iqr')
ds.dataFrame(x = c("dd.M3", "ndvi100_preg.iqr"), newobj = "dd.M3")
ds.colnames("dd.M3",datasources = conns)
ds.dim("dd.M3",datasources = conns)


# Now, we include in the expressionSet this new dataframe with the calculated IQR variables. 
#Additionally, by setting complete.cases to TRUE, the expressionSet will only retain those with
#complete methylation and information for the rest of the variables (complete.cases).
ds.addPhenoData2eSet('methy_0y_EUR.winz', "dd.M3", complete_cases = TRUE, newobj.name = "methy_0y_M3", identifier = "id_methyl")

ds.varLabels("methy_0y_M3",datasources = conns)
# $INMA_cohort
# [1] "id_methyl"         "id"                "sex_methyl"        "age_methyl"        "anc_methyl"        "cohort_methyl"     "bw_methyl"         "ga_methyl"         "Bcell"             "CD4T"             
# [11] "CD8T"              "Gran"              "Mono"              "NK"                "nRBC"              "areases_tert_preg" "greenyn300_preg"   "ndvi100_preg"      "ndvi300_preg"      "preg_smk"         
# [21] "agebirth_m_y"      "edu_m_.0"          "ndvi300_preg.iqr"  "ndvi100_preg.iqr" 

ds.nSamples("methy_0y_M3",datasources = conns)
# Aggregated (nSamplesDS(methy_0y_M3)) [=================================================] 100% / 2s
# $INMA_cohort
# Samples 
# 357 

####################################
#8.2) Descriptive analyses for M3 
####################################
#8.2.1) Descriptive tables
#8.2.2) Correlation between exposures
#8.2.3) Associations between green spaces vars and covariates

### Go to the Descriptive.info directory 
setwd(paste0(res.dir,"/Descriptive.info/"))

###  Create the Descriptive.Tables directory

###########################
#8.2.1) Descriptive tables

### Go to the Descriptive.Tables directory
setwd(paste0(res.dir,"/Descriptive.info/Descriptive.Tables"))

#excract from the Eset M3, the df 
ds.pData("methy_0y_M3","dd.M3",datasources = conns)

DescriptiveTable.final.M3<-dh.getStats(df="dd.M3",vars=c("ndvi300_preg","ndvi100_preg","ndvi300_preg.iqr","ndvi100_preg.iqr","greenyn300_preg","agebirth_m_y","preg_smk","sex","anc_methyl","Bcell","CD4T",            
                                                         "CD8T","Gran","Mono","NK","nRBC","areases_tert_preg","edu_m_.0"),digits=2,conns=conns)
class(DescriptiveTable.final.M3)
DescriptiveTable.final.M3<-as.data.frame(DescriptiveTable.final.M3)
summary(DescriptiveTable.final.M3)
# Length Class  Mode
# categorical 10     tbl_df list
# continuous  15     tbl_df list

#There are two dataframes, one caterogical variables and one for continuous

write.table(
  DescriptiveTable.final.M3$categorical, file=paste0("M3.Descriptive.table.CategoricalVars.",cohort,".",ancestry,".cc.",Date,".txt"),col.names=T, row.names=F, quote=F, sep="\t")

write.table(
  DescriptiveTable.final.M3$continuous, file=paste0("M3.Descriptive.table.ContinuousVars.",cohort,".",ancestry,".cc.",Date,".txt"),col.names=T, row.names=F, quote=F, sep="\t")

######################################
#8.2.2) Correlation between exposures

### Go to the Correlations directory
setwd(paste0(res.dir,"/Descriptive.info/Correlations"))


### a)Correlation between NDVI300 and NDVI100 (Continuous vs continuous)
#########################################################################

#To get the correlation matrix of the exposures you can collect all of the exposures in one new dataframe using the ds.dataFrame function and then use the ds.cor to get all pairwise correlations of the variables in that new dataframe


#correlations 

correlations_ndvi300_ndvi100<-ds.cor(x="dd.M3$ndvi300_preg", y="dd.M3$ndvi300_preg",type= "split") 



##############################################################################
#9.3) EWAS GREEN SPACES PREGNANCY USING ROBUST LINEAR REGRESSIONS WITH LIMMA
##############################################################################

### Before starting with the EWAS analysis...

### Create a folder for the EWAS results 

setwd(res.dir)
dir.create("EWAS.Results")

### Go to the EWAS results directory
setwd(paste0(res.dir,"/EWAS.Results"))

# 9.3.1) Exposure NDVI300 Pregnancy 
# 9.3.2) Exposure NDVI100 Pregnancy 
# 9.3.3) Exposure Green Access Pregnancy 

###################################
# 9.3.1) Exposure NDVI300 Pregnancy 

################
# Preg_N300_M3
################
#"methy_0y_M3"


# a) Run EWAS

Preg_N300_M3 <- ds.limma(model = ~ ndvi300_preg.iqr +  agebirth_m_y + preg_smk + sex_methyl+ areases_tert_preg+ edu_m_.0 + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, Set = "methy_0y_M3",method="robust",sort.by = "none", datasources = conns)
lapply(Preg_N300_M3, head)
lapply(Preg_N300_M3, dim)
png("QQPlot_INMA.EUR_Preg_N300_M3.png", width=800, height=800)
qqplot(Preg_N300_M3$INMA_cohort$P.Value, ci = 0.95)
dev.off() 
#Lambda is inside the picture

write.table(Preg_N300_M3, paste0("EWAS_GS_",cohort,".",ancestry,"_Preg_N300_M3_",Date,".txt"), na="NA") 

###################################
# 8.3.2) Exposure NDVI100 Pregnancy 

################
# Preg_N100_M3
################

# a) Run EWAS
Preg_N100_M3 <- ds.limma(model = ~ ndvi100_preg.iqr +  agebirth_m_y + preg_smk + sex_methyl+ areases_tert_preg+ edu_m_.0 + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, Set = "methy_0y_M3",method="robust",sort.by = "none", datasources = conns)
lapply(Preg_N100_M3, head)
lapply(Preg_N100_M3, dim)
png("QQPlot_INMA.EUR_Preg_N100_M3.png", width=800, height=800)
qqplot(Preg_N100_M3$INMA_cohort$P.Value, ci = 0.95)
dev.off() 

write.table(Preg_N100_M3, paste0("EWAS_GS_",cohort,".",ancestry,"_Preg_N100_M3_",Date,".txt"), na="NA") 

#########################################
# 8.3.3) Exposure Green Access Pregnancy                                                                   

################
# Preg_AcGS_M3
################

# a) Run EWAS

Preg_AcGS_M3 <- ds.limma(model = ~ greenyn300_preg +  agebirth_m_y + preg_smk + sex_methyl+ areases_tert_preg+ edu_m_.0 + agebirth_m_y + preg_smk + sex_methyl+ areases_tert_preg+ edu_m_.0 + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, Set = "methy_0y_M3",method="robust",sort.by = "none", datasources = conns)
lapply(Preg_AcGS_M3, head)
lapply(Preg_AcGS_M3, dim)
png("QQPlot_INMA.EUR_Preg_AcGS_M3.png", width=800, height=800)
qqplot(Preg_AcGS_M3$INMA_cohort$P.Value, ci = 0.95)
dev.off() 

write.table(Preg_AcGS_M3, paste0("EWAS_GS_",cohort,".",ancestry,"_Preg_AcGS_M3_",Date,".txt"), na="NA") 


DSI::datashield.logout(conns)
opal.logout(o)




