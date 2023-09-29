########################
#DMR comb-p
#Sofía Aguilar Lacasaña
#15.05.2023
########################


#install packages

#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

#BiocManager::install("ENmix")

#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

#BiocManager::install("DMRcate")

#load

library(ENmix)
library(data.table)
library(tidyverse)


#MA NDVI300
MA.NDVI300<-read.table("/PROJECTES/INMA_OMICS/analyses/PACE/green_spaces_SA/results/MA/Pregnancy/GWAMA_Results/MetaM3_NDVI300_Filtr/MetaM3_NDVI300_Fixed_Modif.out",header=TRUE, sep="\t")

MA.NDVI300.df<-MA.NDVI300[,c(1,28,29,10)]
chr<-MA.NDVI300.df$chr
chr_a<- str_split_fixed(chr,"r",n=2)
Chr<-chr_a[,2]
length(Chr)


MA.NDVI300.df$chr<-Chr
MA.NDVI300.df$start<-MA.NDVI300.df$pos
MA.NDVI300.df$end<-(MA.NDVI300.df$pos)+1

MA.NDVI300.df<-MA.NDVI300.df[,c(1:2,4:6)]
colnames(MA.NDVI300.df)<-c("probe","chr","p","start","end")

#comb-p () function 

setwd("/PROJECTES/INMA_OMICS/analyses/PACE/green_spaces_SA/results/DMR/combp/pregnancy/NDVI300/")

DMR1<-combp(MA.NDVI300.df,dist.cutoff=1000,bin.size=310,seed=0.05,region_plot=TRUE,mht_plot=TRUE,nCores=10,verbose=TRUE)

#MA NDVI100
MA.NDVI100<-read.table("/PROJECTES/INMA_OMICS/analyses/PACE/green_spaces_SA/results/MA/Pregnancy/GWAMA_Results/MetaM3_NDVI100_Filtr/MetaM3_NDVI100_Fixed_Modif.out",header=TRUE, sep="\t")

MA.NDVI100.df<-MA.NDVI100[,c(1,28,29,10)]
chr<-MA.NDVI100.df$chr
chr_a<- str_split_fixed(chr,"r",n=2)
Chr<-chr_a[,2]
length(Chr)


MA.NDVI100.df$chr<-Chr
MA.NDVI100.df$start<-MA.NDVI100.df$pos
MA.NDVI100.df$end<-(MA.NDVI100.df$pos)+1

MA.NDVI100.df<-MA.NDVI100.df[,c(1:2,4:6)]
colnames(MA.NDVI100.df)<-c("probe","chr","p","start","end")
#comb-p () function 

setwd("/PROJECTES/INMA_OMICS/analyses/PACE/green_spaces_SA/results/DMR/combp/pregnancy/NDVI100/")

DMR1<-combp(MA.NDVI100.df,dist.cutoff=1000,bin.size=310,seed=0.05,region_plot=TRUE,mht_plot=TRUE,nCores=10,verbose=TRUE)


#MA AcGS
MA.AcGS<-read.table("/PROJECTES/INMA_OMICS/analyses/PACE/green_spaces_SA/results/MA/Pregnancy/GWAMA_Results/MetaM3_AcGS_excl.EDEN_Filtr/MetaM3_AcGS_excl.EDEN_Fixed_Modif.out",header=TRUE, sep="\t")

MA.AcGS.df<-MA.AcGS[,c(1,25,26,10)]
chr<-MA.AcGS.df$chr
chr_a<- str_split_fixed(chr,"r",n=2)
Chr<-chr_a[,2]
length(Chr)

MA.AcGS.df$chr<-Chr
MA.AcGS.df$start<-MA.AcGS.df$pos
MA.AcGS.df$end<-(MA.AcGS.df$pos)+1

MA.AcGS.df<-MA.AcGS.df[,c(1:2,4:6)]
colnames(MA.AcGS.df)<-c("probe","chr","p","start","end")

#comb-p () function 

setwd("/PROJECTES/INMA_OMICS/analyses/PACE/green_spaces_SA/results/DMR/combp/pregnancy/AcGS/")

DMR1<-combp(MA.AcGS.df,dist.cutoff=1000,bin.size=310,seed=0.05,region_plot=TRUE,mht_plot=TRUE,nCores=10,verbose=TRUE)

