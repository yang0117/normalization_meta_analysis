#man_fid
rm(list=ls())
library(oligo)
library(GEOquery)
getGEOSuppFiles("GSE38792")
list.files("GSE38792")
untar("GSE38792/GSE38792_RAW.tar", exdir = "GSE38792/CEL")
list.files("GSE38792/CEL")
celfiles <- list.files("GSE38792/CEL", full = TRUE)
rawData <- read.celfiles(celfiles)
rawData
class(rawData)
getClass("GeneFeatureSet")
exprs(rawData)[1:4,1:3]
max(exprs(rawData))

filename <- sampleNames(rawData)
pData(rawData)$filename <- filename
sampleNames <- sub(".*_", "", filename)
sampleNames <- sub(".CEL.gz$", "", sampleNames)
sampleNames(rawData) <- sampleNames
pData(rawData)$group <- ifelse(grepl("^OSA", sampleNames(rawData)),
                               "OSA", "Control")
normData <- oligo::rma(rawData)
pData(rawData)

boxplot(rawData)

probeInfo = oligo:::getFidProbeset(rawData)
str(probeInfo)
normData <- rma(rawData)
backgroundCorrect(rawData)
