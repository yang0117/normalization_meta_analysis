#
# Step1_LoadData.R
#
# Script that:
#  1. loads the meth450k .idat files and
#  2. load the sample specific data
#  3. perform first pass QC
#  4. filter out problematic probes 
#

#---------------------------------------------------------------------------#
# Install commands
#---------------------------------------------------------------------------#

if(F){ # trip flag to True to re-install/update packages

  # note: may need to use developer version:
  # library(BiocInstaller); useDevel() 
  
  source("http://bioconductor.org/biocLite.R"); 
  biocLite("minfi"); biocLite("minfiData"); biocLite("IlluminaHumanMethylation450kmanifest")
  biocLite("methylumi"); 
 # biocLite("ClassDiscovery")
  
}


rm(list=ls())

#---------------------------------------------------------------------------#
# Section Flags
#---------------------------------------------------------------------------#

LoadDataFlag=T # load the data from the raw files

QCFlag=T # run various QC reports


#---------------------------------------------------------------------------#
# Libraries / Options
#---------------------------------------------------------------------------#

options(width=70)
require(minfi)
require(minfiData)
require(IlluminaHumanMethylation450kmanifest)
require(methylumi)
#require(ClassDiscovery)
#require(shinyMethyl)

# set the working directory here, if running on laptop and not Rstudio server
# setwd("C:/Users/Jingjing/Dropbox/XRenMeth14")



#---------------------------------------------------------------------------#
# Load the Data
#---------------------------------------------------------------------------#

if(LoadDataFlag){

#  baseDir=c("~/Dropbox/Work/Consulting/XRen/ArsenicMethFinal/Data")
  baseDir=c("Data/iDats")
  targets <- read.450k.sheet(baseDir)
  targets
  
  # remove cancer patients and JQ samples for the time
  # being, as they may negatively impact the SWAN normalization
  RGset <- read.450k.exp(targets = targets[c(1:20,25:64),])
  
  # look at pheno data:
  pd=pData(RGset)
  
  #  load the supplemental phenotype data stored in SamplePheno.csv
  
  rawPheno=read.table("Data/SamplePheno.csv",sep=",",header=T,skip=1)
  
  mtch=match(pd[,1],rawPheno[,1])
  
  pd$Family.ID=rawPheno[mtch,2]
  pd$Member.ID=rawPheno[mtch,3]
  pd$Arsenic=rawPheno[mtch,4]
  
  save.image("RData/Step1.RData")
  
}



if((!LoadDataFlag)&QCFlag) load("RData/Step1.RData")


#---------------------------------------------------------------------------#
# Run QC 
#   note: this also poulates MSet.raw, Mset.norm, and MSet.swan
#---------------------------------------------------------------------------#


if(QCFlag){
  
  # generate QC report
  qcReport(RGset, sampNames = pd$Sample_Name, 
           sampGroups = pd$Sample_Group, pdf = "Figures/Step1_QC/qcReport.pdf")
  
  
  # try some normalization methods.
  
  MSet.raw <- preprocessRaw(RGset)
  
 # MSet.norm <- preprocessIllumina(RGset, bg.correct = TRUE,
  #                                normalize = "controls", reference = 2)
  # note: reference = 2 is an arbitrary choice!
  
  MSet.swan <- preprocessSWAN(RGset, MSet.raw)
  
  pdf("Figures/Step1_QC/BetasByType.pdf")
  par(mfrow=c(1,2))
  plotBetasByType(MSet.raw[,1], main = "Raw")
  plotBetasByType(MSet.swan[,1], main = "SWAN")
  dev.off()
  
  # "A simple way to quickly check if a sample failed is to look at the log median
  # intensity in both the methylated and unmethylated channels. When plotting
  # the U channel against the M channel, it has been observed that good samples
  # cluster together, while failed samples tend to separate and to have lower median
  # intensities.
  # But wait, so far we only have green and red intensities. We need the methylated
  # and unmethylated signals; the function preprocessRaw is done for that. It takes
  # as input a RGChannelSet, convert the red and green intensities to methylated
  # and unmethylated signals according to the probe design stored in the manifest
  # object, and returns the converted signals in a new object of class MethylSet."
  # 

  qc <- getQC(MSet.raw)
  pdf("Figures/Step1_QC/getQCMsetraw_Yang.pdf")
  plotQC(qc)
  dev.off()
  
  qc.swan <- getQC(MSet.swan)
  pdf("Figures/Step1_QC/getQCMsetswan_Yang.pdf")
  plotQC(qc.swan)
  dev.off()
  
  
  pdf("Figures/Step1_QC/densPlotMSetraw_Yang.pdf")
  densityPlot(MSet.raw)
  dev.off()
  
  pdf("Figures/Step1_QC/densPlotMSetswam_Yang.pdf") 
  densityPlot(MSet.swan)
  dev.off()
    
  save.image("RData/Step1.RData")
  
}

if(!QCFlag) load("RData/Step1.RData")


#---------------------------------------------------------------------------#
# Filter the Data
#---------------------------------------------------------------------------#


# load the .csv files with IMA-snp and cross-reactive probeset lists
badsnp=read.csv("FilterFiles/IMA-snp.csv",header = T)
crossreact=read.csv("FilterFiles/cross-reactive.csv",header = T)

# need to grab the autosomal meth sites..
manifest <- getManifest(RGset)
manifest

probeI=getProbeInfo(manifest,type = "I")[,1]
probeII=getProbeInfo(manifest,type = "II")[,1]
probecontrol=getProbeInfo(manifest,type = "Control")[,1]
probesnpI=getProbeInfo(manifest,type = "SnpI")[,1]
probesnpII=getProbeInfo(manifest,type = "SnpII")[,1]
allprobe=c(probeI,probeII,probecontrol,probesnpI,probesnpII)

probenames=cbind(allprobe,seq(1:length(allprobe)))
colnames(probenames)=c("TargetID","Index")

#
# Filter by bad SNP and cross-reactive
#
badsnp$snpind=rep("snp",length.out=dim(badsnp)[1])
crossreact$crtind=rep("crt",length.out=dim(crossreact)[1])

badsnp$TargetID=as.character(badsnp[,1])
crossreact$TargetID=as.character(crossreact[,1])

tempsnp=merge(badsnp,probenames, by = "TargetID")
tempcrt=merge(crossreact,probenames, by = "TargetID")

badSNPDX=as.numeric(tempsnp$Index)
crossreactDX=as.numeric(tempcrt$Index)

#
# Filter by detection p-value
#

# "minfi provides several functions and diagnostic plots to assess quality of the
# methylation samples. As a starting point, we suggest to look at the function
# detectionP() which identifies failed positions defined as both the methylated
# and unmethylated channel reporting background signal levels:"

detP <- detectionP(RGset)
failed <- detP > 0.01  
head(failed, n=3)

#To see the fraction of failed positions per sample:
colMeans(failed)

# let's look at histogram of failed detections across probesets
hist(rowMeans(failed))
prbFail=rowMeans(failed)
hist(prbFail[prbFail>=(1/60)])

# of the assays, identify which ones had a detection p-value of greater than 0.01 
# for at least 1 of the 72 samples - maybe too stringent?
# perhaps we will change this..

DetectFailDX=which(prbFail>=(1/60)) 


#
# Filter by autosomal location
#

# use methylumi package to get ratioSet then get genome information

ratioSet <- ratioConvert(MSet.swan, what = "beta", keepCN = F)
ratioSet

gRatioSet <- mapToGenome(ratioSet, mergeManifest = TRUE)
gRatioSet

ratioSet.raw <- ratioConvert(MSet.raw, what = "beta", keepCN = F)
ratioSet.raw

gRatioSet.raw <- mapToGenome(ratioSet.raw, mergeManifest = TRUE)
gRatioSet.raw

showClass("GenomicRatioSet")

# Note that the GenomicRatioSet extends the class SummarizedExperiment. Here
# are the main accessors functions to access the data:
#   > getBeta(gRatioSet)
# > getM(gRatioSet)
# > getCN(gRatioSet)
# > sampleNames <- sampleNames(gRatioSet)
# > probeNames <- featureNames(gRatioSet)
# > pheno <- pData(gRatioSet)

# To extract the genomic locations:
gRanges <- rowData(gRatioSet)
head(gRanges, n= 3)

# identify X and Y linked assays:
XchromDX=which(gRanges@seqnames=="chrX")
YchromDX=which(gRanges@seqnames=="chrY")

#
# Get data matrices that we will operate on later with clustering and so forth..
#

XB=getBeta(gRatioSet)
XB.raw=getBeta(gRatioSet.raw)

colnames(XB)=pd[,1]
colnames(XB.raw)=pd[,1]


#
# Here is the filter call that I used previously:
#

keepDX=1:(dim(XB)[1])

summary(DetectFailDX)
summary(badSNPDX)
summary(crossreactDX)

# NAdx - identify probes which are missing data
NAdx=which(rowSums(is.na(XB.raw))!=0)

keepDX=setdiff(keepDX,c(XchromDX,YchromDX,DetectFailDX,badSNPDX,NAdx,crossreactDX))################ this is the Index of the final dataset after filtering.

length(keepDX)

# 
# Now, set up the data matrices
#

XB=XB[keepDX,]
XB.raw=XB.raw[keepDX,]


#-----------------------------------#
# save workspace
#-----------------------------------#


save(file="RData/Step1_Xmats.RData",XB,XB.raw,pd,keepDX)



