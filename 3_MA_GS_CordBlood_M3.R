###########################################
# META-ANALYSIS M3 Green Spaces Cord blood 
##Sofía Aguilar Lacasaña
###########################################

## ######################################### ##
##  Meta-Analysis to use with EASIER package ##
## ######################################### ##


## -------------------------------------
##  Install EASIER Package Code
## -------------------------------------
##
##  Uncomment this code to install EASIER package
#
# # Install devtools
# install.packages("devtools")
#
# # Install required packages
# devtools::source_url("https://raw.githubusercontent.com/isglobal-brge/EASIER/HEAD/installer.R")

# # Install EASIER package
# devtools::install_github("isglobal-brge/EASIER@HEAD")

##  END -  Install EASIER Package Code
## -------------------------------------

# load package
library(EASIER)


########## ----------  VARIABLES DEFINED BY USER  ----------  ##########

# Set working directory to metaanalysis folder
setwd("/PROJECTES/INMA_OMICS/analyses/PACE/green_spaces_SA/results/MA/new/03052023")

# Files used in QC, needed in meta-analysis to plot ForestPlot
files<-c("/PROJECTES/INMA_OMICS/analyses/PACE/green_spaces_SA/results/QC/preQC/EWAS_GS_ALSPAC.EUR_preg_N300_M3_QCdf.txt",
"/PROJECTES/INMA_OMICS/analyses/PACE/green_spaces_SA/results/QC/preQC/EWAS_GS_BIB.EUR_preg_N300_M3_QCdf.txt",
"/PROJECTES/INMA_OMICS/analyses/PACE/green_spaces_SA/results/QC/preQC/EWAS_GS_EDEN.EUR_preg_N300_M3_QCdf.txt",
"/PROJECTES/INMA_OMICS/analyses/PACE/green_spaces_SA/results/QC/preQC/EWAS_GS_ENVIRONAGE.450K.EUR_preg_N300_M3_QCdf.txt",
"/PROJECTES/INMA_OMICS/analyses/PACE/green_spaces_SA/results/QC/preQC/EWAS_GS_ENVIRONAGE.EPIC.EUR_preg_N300_M3_QCdf.txt",
"/PROJECTES/INMA_OMICS/analyses/PACE/green_spaces_SA/results/QC/preQC/EWAS_GS_GenR.EUR_preg_N300_M3_QCdf.txt",
"/PROJECTES/INMA_OMICS/analyses/PACE/green_spaces_SA/results/QC/preQC/EWAS_GS_INMA.EUR_preg_N300_M3_QCdf.txt",
"/PROJECTES/INMA_OMICS/analyses/PACE/green_spaces_SA/results/QC/preQC/EWAS_GS_PFC.EUR_preg_N300_M3_QCdf.txt",
"/PROJECTES/INMA_OMICS/analyses/PACE/green_spaces_SA/results/QC/preQC/EWAS_GS_ALSPAC.EUR_preg_N100_M3_QCdf.txt",
"/PROJECTES/INMA_OMICS/analyses/PACE/green_spaces_SA/results/QC/preQC/EWAS_GS_BIB.EUR_preg_N100_M3_QCdf.txt",
"/PROJECTES/INMA_OMICS/analyses/PACE/green_spaces_SA/results/QC/preQC/EWAS_GS_EDEN.EUR_preg_N100_M3_QCdf.txt",
"/PROJECTES/INMA_OMICS/analyses/PACE/green_spaces_SA/results/QC/preQC/EWAS_GS_ENVIRONAGE.450K.EUR_preg_N100_M3_QCdf.txt",
"/PROJECTES/INMA_OMICS/analyses/PACE/green_spaces_SA/results/QC/preQC/EWAS_GS_ENVIRONAGE.EPIC.EUR_preg_N100_M3_QCdf.txt",
"/PROJECTES/INMA_OMICS/analyses/PACE/green_spaces_SA/results/QC/preQC/EWAS_GS_GenR.EUR_preg_N100_M3_QCdf.txt",
"/PROJECTES/INMA_OMICS/analyses/PACE/green_spaces_SA/results/QC/preQC/EWAS_GS_INMA.EUR_preg_N100_M3_QCdf.txt",
"/PROJECTES/INMA_OMICS/analyses/PACE/green_spaces_SA/results/QC/preQC/EWAS_GS_PFC.EUR_preg_N100_M3_QCdf.txt",
"/PROJECTES/INMA_OMICS/analyses/PACE/green_spaces_SA/results/QC/preQC/EWAS_GS_ALSPAC.EUR_preg_AcGS_M3_QCdf.txt",
"/PROJECTES/INMA_OMICS/analyses/PACE/green_spaces_SA/results/QC/preQC/EWAS_GS_BIB.EUR_preg_AcGS_M3_QCdf.txt",
"/PROJECTES/INMA_OMICS/analyses/PACE/green_spaces_SA/results/QC/preQC/EWAS_GS_EDEN.EUR_preg_AcGS_M3_QCdf.txt",
"/PROJECTES/INMA_OMICS/analyses/PACE/green_spaces_SA/results/QC/preQC/EWAS_GS_GenR.EUR_preg_AcGS_M3_QCdf.txt",
"/PROJECTES/INMA_OMICS/analyses/PACE/green_spaces_SA/results/QC/preQC/EWAS_GS_INMA.EUR_preg_AcGS_M3_QCdf.txt",
"/PROJECTES/INMA_OMICS/analyses/PACE/green_spaces_SA/results/QC/preQC/EWAS_GS_PFC.EUR_preg_AcGS_M3_QCdf.txt")

prefixes <- c('ALSPAC_M3_NDVI300','BIB_M3_NDVI300','EDEN_M3_NDVI300','ENV.450K_M3_NDVI300','ENV.EPIC_M3_NDVI300','GenR_M3_NDVI300','INMA_M3_NDVI300','PFC_M3_NDVI300','ALSPAC_M3_NDVI100',
              'BIB_M3_NDVI100','EDEN_M3_NDVI100','ENV.450K_M3_NDVI100','ENV.EPIC_M3_NDVI100','GenR_M3_NDVI100','INMA_M3_NDVI100','PFC_M3_NDVI100','ALSPAC_M3_AcGS','BIB_M3_AcGS','EDEN_M3_AcGS',
              'GenR_M3_AcGS','INMA_M3_AcGS','PFC_M3_AcGS')
# Samples in original files used in QC

N<-c(618,325,137,188,345,1171,357,172,618,325,137,188,345,1171,357,172,618,325,95,1171,357,172)
# Define data for each meta-analysis
metafiles <- list(
   'MetaM3_NDVI300' = c('ALSPAC_M3_NDVI300','EDEN_M3_NDVI300','ENV.450K_M3_NDVI300','ENV.EPIC_M3_NDVI300','GenR_M3_NDVI300','INMA_M3_NDVI300','PFC_M3_NDVI300'),
   'MetaM3_NDVI100' = c('ALSPAC_M3_NDVI100','EDEN_M3_NDVI100','ENV.450K_M3_NDVI100','ENV.EPIC_M3_NDVI100','GenR_M3_NDVI100','INMA_M3_NDVI100','PFC_M3_NDVI100'),
   'MetaM3_AcGS' = c('ALSPAC_M3_AcGS','EDEN_M3_AcGS','GenR_M3_AcGS','INMA_M3_AcGS','PFC_M3_AcGS'),
   'MetaM3_AcGS_excl.EDEN' = c('ALSPAC_M3_AcGS','GenR_M3_AcGS','INMA_M3_AcGS','PFC_M3_AcGS')
   )


# Array type, used in each meta-analysis : EPIC or 450K
artype <- c('450K','450K','450K','450K')

# Define maximum percent missing for each CpG
#     if pcenMissin = 0 only runs meta-analysis with all data
pcentMissing <- 0.5 # CpGs with precense lower than pcentMissing after GWAS meta-analysis will be deleted from the study.


# Paths with QCResults and path to store GWAMA results
results_folder <- '/PROJECTES/INMA_OMICS/analyses/PACE/green_spaces_SA/results/QC/QC_03052023/'
results_gwama <- '.'


# Venn diagrams ==> IMPORTANT : maximum 5 meta-analysis by venn diagram
venndiag_threshold <- 0.05
venn_diagrams <- list(c("MetaM3_NDVI300","MetaM3_NDVI100","MetaM3_AcGS_excl.EDEN"))



########## ----------  END VARIABLES DEFINED BY USER  ----------  ##########


# GWAMA binary path  (GWAMA IsGlobal Server installation)
gwama.dir <-"/PROJECTES/HELIX_OMICS/software/GWAMA/"

## Create directory for GWAMA configuration files and GWAMA_Results
if(!dir.exists(file.path(paste(results_gwama, "GWAMA", sep="/") )))
   suppressWarnings(dir.create(file.path( paste(results_gwama, "GWAMA", sep="/"))))

## Create directory for GWAMA_Results
outputfolder <- paste0(results_gwama, "/GWAMA_Results")
if(!dir.exists(file.path( outputfolder )))
   suppressWarnings(dir.create(file.path( outputfolder)))


# Create hapmap files for the different artypes that we cab use (450K and EPIC)
# map file is used in Manhattan plots
hapmapfile_450K <- paste(results_gwama,"GWAMA", "hapmap_450K.map" ,sep = "/")
generate_hapmap_file("450K", hapmapfile_450K)
hapmapfile_EPIC <- paste(results_gwama,"GWAMA", "hapmap_EPIC.map" ,sep = "/")
generate_hapmap_file("EPIC", hapmapfile_EPIC)



for( metf in 1:length(metafiles))
{

   list.lowCpGs <- NULL

   # Create folder for a meta-analysis in GWAMA folder, here we store the GWAMA input files for each meta-analysis,
   # We create one for complete meta-analysis
   if(!dir.exists(file.path( paste(results_gwama,"GWAMA", names(metafiles)[metf] ,sep="/") )))
      suppressWarnings(dir.create(file.path( paste(results_gwama,"GWAMA", names(metafiles)[metf], sep="/"))))
   # We create another for meta-analysis without filtered CpGs with low percentage (sufix _Filtr)
   if(!dir.exists(file.path( paste0(results_gwama,"/GWAMA/", names(metafiles)[metf],"_Filtr") )))
      suppressWarnings(dir.create(file.path( paste0(results_gwama,"/GWAMA/", names(metafiles)[metf],"_Filtr"))))

   # GWAMA File name base
   inputfolder <- paste0(results_gwama,"/GWAMA/",  names(metafiles)[metf])

   modelfiles <- unlist(metafiles[metf])

   runs <- c('Normal', 'lowcpgs') # Execution with all CpGs and without filtered CpGs
   lowCpGs = FALSE;
   outputfiles <- list()

   outputgwama <- paste(outputfolder,names(metafiles)[metf],sep = '/')

   for(j in 1:length(runs))
   {
      if(runs[j]=='lowcpgs') {
         lowCpGs = TRUE
         # Get low presence CpGs in order to exclude this from the new meta-analysis
         list.lowCpGs <- get_low_presence_CpGs(outputfiles[[j-1]], pcentMissing)
         inputfolder <- paste0(results_gwama,"/GWAMA/",  names(metafiles)[metf], "_Filtr")
         outputgwama <- paste0(outputgwama,"_Filtr")
      }

      # Create GWAMA files for each file in meta-analysis and execute GWAMA
      for ( i in 1:length(modelfiles) )
         create_GWAMA_files(file.path(results_folder,modelfiles[i]),  modelfiles[i], inputfolder, N[which(prefixes==modelfiles[i])], list.lowCpGs )


      # Get hapmapfile attending to current metaanalysis artype
      hapmapfile <- hapmapfile_450K
      if(artype[metf]=='EPIC'){
         hapmapfile <- hapmapfile_EPIC
      }

      #.Original.#outputfiles[[runs[j]]] <- execute_GWAMA_MetaAnalysis(prefixgwama, names(metafiles)[metf])
      outputfiles[[runs[j]]] <- run_GWAMA_MetaAnalysis(inputfolder, outputgwama, names(metafiles)[metf], gwama.dir, hapmapfile)

      # Post-metha-analysis QC --- >>> adds BN and FDR adjustment
      dataPost <- get_descriptives_postGWAMA(outputgwama, outputfiles[[runs[j]]], modelfiles, names(metafiles)[metf], artype[metf], N[which(prefixes %in% modelfiles)] )

      # Forest-Plot
      plot_ForestPlot( dataPost, metafiles[[metf]], runs[j], inputfolder, names(metafiles)[metf], files, outputgwama, nsignificatives = 30)

   }

}


# Venn_Diagrams for for meta-analysis with fixed effects

#for (i in 1:length(venn_diagrams)){
#   if(length(venn_diagrams[[i]])>1){
#      plot_venndiagram(venn_diagrams[[i]], qcpath = outputfolder, plotpath =  paste0(results_gwama, "/GWAMA_Results"), pattern = '_Fixed_Modif.out',bn='Bonferroni', fdr='FDR', venndiag_threshold)
#   }
#}

if(dir.exists(file.path( paste(results_gwama, "GWAMA", sep="/") )))
   unlink(file.path(results_gwama, "GWAMA"), recursive=TRUE)
