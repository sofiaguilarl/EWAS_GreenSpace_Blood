###########################################
# LEAVE-ONE-OUT ANALYSIS
#21.04.2023
#Sofía Aguilar Lacasaña
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

setwd("/PROJECTES/INMA_OMICS/analyses/PACE/green_spaces_SA/results/MA/LifeLong/")

# Files used in QC, needed in meta-analysis to plot ForestPlot
files<-c("/PROJECTES/INMA_OMICS/analyses/PACE/green_spaces_SA/results/QC/preQC/EWAS_GS_ALSPAC.EUR_life-long_N300_M3_QCdf.txt",
"/PROJECTES/INMA_OMICS/analyses/PACE/green_spaces_SA/results/QC/preQC/EWAS_GS_GenR.EUR_life-long_N300_M3_QCdf.txt",
"/PROJECTES/INMA_OMICS/analyses/PACE/green_spaces_SA/results/QC/preQC/EWAS_GS_HELIX.EUR_life-long_N300_M3_QCdf.txt",
"/PROJECTES/INMA_OMICS/analyses/PACE/green_spaces_SA/results/QC/preQC/EWAS_GS_ALSPAC.EUR_life-long_N100_M3_QCdf.txt",
"/PROJECTES/INMA_OMICS/analyses/PACE/green_spaces_SA/results/QC/preQC/EWAS_GS_GenR.EUR_life-long_N100_M3_QCdf.txt",
"/PROJECTES/INMA_OMICS/analyses/PACE/green_spaces_SA/results/QC/preQC/EWAS_GS_HELIX.EUR_life-long_N100_M3_QCdf.txt")

prefixes <- c('ALSPAC_M3_NDVI300_lf','GenR_M3_NDVI300_lf','HELIX_M3_NDVI300_lf','ALSPAC_M3_NDVI100_lf','GenR_M3_NDVI100_lf','HELIX_M3_NDVI100_lf')
# Samples in original files used in QC

N<-c(682,440,727,682,440,727)
# Define data for each meta-analysis
metafiles <- list(
   'Leave-one-M3_NDVI300_lf_excALSPAC' = c('GenR_M3_NDVI300_lf','HELIX_M3_NDVI300_lf'),
   'Leave-one-M3_NDVI300_lf_excGenR' = c('ALSPAC_M3_NDVI300_lf','HELIX_M3_NDVI300_lf'),
   'Leave-one-M3_NDVI300_lf_excHELIX' = c('ALSPAC_M3_NDVI300_lf','GenR_M3_NDVI300_lf'),
   'Leave-one-M3_NDVI100_lf_excALSPAC' = c('GenR_M3_NDVI100_lf','HELIX_M3_NDVI100_lf'),
   'Leave-one-M3_NDVI100_lf_excGenR' = c('ALSPAC_M3_NDVI100_lf','HELIX_M3_NDVI100_lf'),
   'Leave-one-M3_NDVI100_lf_excHELIX' = c('ALSPAC_M3_NDVI100_lf','GenR_M3_NDVI100_lf'))

# Array type, used in each meta-analysis : EPIC or 450K
artype <- c('450K','450K','450K','450K','450K','450K')

# Define maximum percent missing for each CpG
#     if pcenMissin = 0 only runs meta-analysis with all data
pcentMissing <- 0.5 # CpGs with precense lower than pcentMissing after GWAS meta-analysis will be deleted from the study.


# Paths with QCResults and path to store GWAMA results
results_folder <- '/PROJECTES/INMA_OMICS/analyses/PACE/green_spaces_SA/results/QC/QC_lifelong/'
results_gwama <- '.'


# Venn diagrams ==> IMPORTANT : maximum 5 meta-analysis by venn diagram
venndiag_threshold <- 0.05
venn_diagrams <- list(c('Leave-one-M3_NDVI300_lf_excALSPAC','Leave-one-M3_NDVI100_lf_excALSPAC'))
                      




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
