#
# Step1B.LoadTARC.R
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
tgen_path = "./Rcode/SampleRawData/TARC/CEL files/"
file_vec <- list.celfiles(tgen_path,full.name=T)
#load data
dat <- read.celfiles(file_vec)
dat
pData(dat)
#check the express
dat_exprs <- exprs(dat)
dim(dat_exprs)
dat_exprs[1:5,1:5]
#look at the probe intensities across the samples
boxplot(dat)

