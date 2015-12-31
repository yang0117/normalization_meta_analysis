#
# Step1A.LoadTGEN.R
# 
# load the example TGEN data
#

# clear environment
rm(list=ls())
source("http://www.bioconductor.org/biocLite.R")

#install and load package
if(!require(pd.genomewidesnp.6)){
  biocLite("pd.genomewidesnp.6")
  biocLite("genomewidesnp6cdf")
  require(pd.genomewidesnp.6)
}

if(!require(oligo)){
  biocLite("oligo")
  require(oligo)
}
library(oligo)


# old code
#install.packages("genomewidesnp6cdf")
#biocLite("makecdfenv")

# require(affy)
# require(limma)
# require(annotate)
# require(hgu95a.db)
# require(makecdfenv)

#make.cdf.package("GenomeWideSNP_6.cdf", species = "Homo_sapiens")

#get the path of each .cel file
tgen_path = "./SampleRawData/TGEN/CEL_FILES_AFFY6"
file_vec <- list.celfiles(tgen_path,full.name=T)
#load data
dat <- read.celfiles(file_vec,pkgname = "pd.genomewidesnp.6")
dat <- read.celfiles(file_vec)
dat
test1 <- phenoData(dat)

library(convert)
as(test1, "ExpressionSet")
dat1 <- ReadAffy(filenames = file_vec) 
#get expression set
dat2 <- justRMA(filenames = file_vec)
probeInfo = oligo:::getFidProbeset(dat)
str(probeInfo)
idx = probeInfo[["fid"]]

showClass("SnpCnvFeatureSet")

setMethod("rma", "SnpCnvFeatureSet",
          function(object, background=TRUE, normalize=TRUE, subset=NULL){
            sql <- paste("SELECT man_fsetid, fid FROM pmfeatureCNV")
            probeInfo <- dbGetQuery(db(object), sql)
            probeInfo <- probeInfo[order(probeInfo[["man_fsetid"]]),]
            pms <- exprs(object)[probeInfo[["fid"]],,drop=FALSE]
            exprs <- basicRMA(pms, probeInfo[["man_fsetid"]], normalize, background)
            out <- new("ExpressionSet",
                       phenoData = phenoData(object),
                       annotation = annotation(object),
                       experimentData = experimentData(object),
                       exprs = exprs)
            return(out)
          })

filename <- sampleNames(dat)
pData(dat)$filename <- filename
sampleNames <- sub(".*_", "", filename)
sampleNames <- sub(".CEL.gz$", "", sampleNames)
sampleNames(dat) <- sampleNames
pData(dat)

exp_dat <- oligo::rma(dat,normalize=F)
test1 <- phenoData(dat)
test2 <- pData(test1)
browseVignettes("convert")

pData(dat)
#check the express
dat_exprs <- exprs(dat)
dim(dat_exprs)
dat_exprs[1:5,1:5]
#look at the probe intensities across the samples
boxplot(dat)
