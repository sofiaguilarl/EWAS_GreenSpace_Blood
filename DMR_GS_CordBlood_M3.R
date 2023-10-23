###################################
#DMR (ENmix:combp and DMRcate)
#Sofía Aguilar Lacasaña
###################################

# Example with NDVI100 MA results 
#load

library(ENmix)
library(data.table)
library(tidyverse)
library(Biobase)
library(dplyr)
library(DMRcate)
library(minfi)
library(bumphunter)
library("TxDb.Hsapiens.UCSC.hg19.knownGene")

################
#NDVI100 PREG
################

##############
#ENmix:comb.p
##############

#MA NDVI100
MA.NDVI100<-read.table("/PROJECTES/INMA_OMICS/analyses/PACE/green_spaces_SA/results/MA/Pregnancy/GWAMA_Results/MetaM3_NDVI100_Filtr/MetaM3_NDVI100_Fixed_Modif.out",header=TRUE, sep="\t")
# 

MA.NDVI100.df<-MA.NDVI100[,c("CpGId","chr","pos","p.value")]
chr<-MA.NDVI100.df$chr
chr_a<- str_split_fixed(chr,"r",n=2)
Chr<-chr_a[,2]
length(Chr)

MA.NDVI100.df$chr<-Chr
MA.NDVI100.df$start<-MA.NDVI100.df$pos
MA.NDVI100.df$end<-(MA.NDVI100.df$pos)+1

MA.NDVI100.df<-MA.NDVI100.df[,c(1:2,4:6)] #select probe, chr, position, start position, end position
colnames(MA.NDVI100.df)<-c("probe","chr","p","start","end")

#comb-p () function 

setwd("/PROJECTES/INMA_OMICS/analyses/PACE/green_spaces_SA/results/DMR/03102023/")

DMR1<-combp(MA.NDVI100.df,dist.cutoff=1000,bin.size=310,seed=0.05,region_plot=TRUE,mht_plot=TRUE,nCores=10,verbose=TRUE)

DMR1<-read.table("resu_combp.csv",sep=",",header=TRUE)

setwd("/PROJECTES/INMA_OMICS/analyses/PACE/green_spaces_SA/results/DMR/03102023/")

write.table(DMR1,"combp_M3_NDVI100_pregnancy.txt",sep=",")

#rm(list=ls())

#############
#DMRcate
#############

#load cpg.annotate object with MA data

load("/PROJECTES/INMA_OMICS/analyses/PACE/green_spaces_SA/results/DMR/cpg.annotate/myannotation.NDVI100.preg.RData")

dmrcoutput <- dmrcate(myannotation, lambda=1000, C=2)
DMR_RESULTS<-extractRanges(dmrcoutput, genome = "hg19")
DMR_RESULTS_1<-as.data.frame(DMR_RESULTS)

setwd("/PROJECTES/INMA_OMICS/analyses/PACE/green_spaces_SA/results/DMR/03102023/")
write.table(DMR_RESULTS_1,"DMRcate_M3_NDVI100_pregnancy.txt",sep=",")

#rm(list=ls())

############
#Merge DMRs
############

genome <- annotateTranscripts(TxDb.Hsapiens.UCSC.hg19.knownGene)

#Match genes

# #Comb-p

combp.NDVI100<-read.table("/PROJECTES/INMA_OMICS/analyses/PACE/green_spaces_SA/results/DMR/03102023/combp_M3_NDVI100_pregnancy.txt",sep=",",header=TRUE)

#subset to significant DMRs
combp.NDVI100<- combp.NDVI100[combp.NDVI100$sidak <0.05,]

combp.NDVI100$chr <- paste0("chr", combp.NDVI100$chr)
combp.NDVI100_subset<-makeGRangesFromDataFrame(combp.NDVI100, keep.extra.columns=TRUE, start.field="start", end.field="end", seqnames.field="chr")
combp.NDVI100_genes<-matchGenes(combp.NDVI100_subset, genome, type = "any",  skipExons = TRUE, verbose = TRUE)

dim(combp.NDVI100)
dim(combp.NDVI100_genes)

combp.final.results<-cbind(combp.NDVI100,combp.NDVI100_genes)
combp.final.results_s<-combp.final.results[,c(1:9,12)]
colnames(combp.final.results_s)<-c("combp.chr","combp.start","combp.end","combp.p","combp.fdr","combp.sidak" , "combp.nprobe","combp.probe","combp.name","combp.region")

#DMRcate

DMRcate.NDVI100<-read.table("/PROJECTES/INMA_OMICS/analyses/PACE/green_spaces_SA/results/DMR/03102023/DMRcate_M3_NDVI100_pregnancy.txt",sep=",")

#subset to significant DMRs
DMRcate.NDVI100<-DMRcate.NDVI100[DMRcate.NDVI100$HMFDR < 0.05,]

DMRcate.NDVI100_subset<-makeGRangesFromDataFrame(DMRcate.NDVI100, keep.extra.columns=TRUE, start.field="start", end.field="end", seqnames.field="seqnames")
DMRcate.NDVI100_genes<-matchGenes(DMRcate.NDVI100_subset, genome, type = "any",  skipExons = TRUE, verbose = TRUE)
# 
DMRcate.final.results<-cbind(DMRcate.NDVI100,DMRcate.NDVI100_genes)
DMRcate.final.results_s<-DMRcate.final.results[,c(1:14,17)]
#View(DMRcate.final.results)
colnames(DMRcate.final.results_s)<-c("DMRcate.seqnames","DMRcate.start","DMRcate.end","DMRcate.width","DMRcate.strand","DMRcate.no.cpgs","DMRcate.min_smoothed_fdr", "DMRcate.Stouffer", "DMRcate.HMFDR",
"DMRcate.Fisher", "DMRcate.maxdiff","DMRcate.meandiff","DMRcate.overlapping.genes", "DMRcate.name","DMRcate.region")

#Merge combp and DMRcate results by Gene

OVERLAP_GENES<-merge(combp.final.results_s,DMRcate.final.results_s,by.x="combp.name",by.y="DMRcate.name")

setwd("/PROJECTES/INMA_OMICS/analyses/PACE/green_spaces_SA/results/DMR/03102023/")
write.table(OVERLAP_GENES,"DMR_NDVI100_preg_corrected_19102023.txt",sep=",")


