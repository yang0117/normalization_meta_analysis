#
# Step1A.LoadTGEN.R
# 
# load the example TGEN data
#

# clear environment
rm(list=ls())

#install and load package
source("http://www.bioconductor.org/biocLite.R")

if(!require(oligo)){
  biocLite("oligo")
  require(oligo)
}
browseVignettes("oligo")

# if(!require(pd.genomewidesnp.6)){
#   biocLite("pd.genomewidesnp.6")
#   biocLite("genomewidesnp6cdf")
#   require(pd.genomewidesnp.6)
# }

# old code
#install.packages("genomewidesnp6cdf")
#biocLite("makecdfenv")

# require(limma)
# require(annotate)
# require(hgu95a.db)
# require(makecdfenv)

#make.cdf.package("GenomeWideSNP_6.cdf", species = "Homo_sapiens")

#get the path of each .cel file
tgen_path = "./SampleRawData/TGEN/CEL_FILES_AFFY6"
#tgen_path = "~/genedata/TGEN/CEL_FILES_AFFY6"
file_vec <- list.celfiles(tgen_path,full.name=T)
#load data
#dat <- read.celfiles(file_vec,pkgname = "pd.genomewidesnp.6")
dat <- read.celfiles(file_vec)
dat
showClass("SnpCnvFeatureSet")

test1 <- phenoData(dat)

#get expression set
probeInfo = oligo:::getFidProbeset(dat)
str(probeInfo)

#detach("package:affy", unload=TRUE)
exp_dat <- rma(dat,normalize=F)
class(exp_dat)

pData(dat)
#check the express
dat_exprs <- exprs(dat)
dim(dat_exprs)
dat_exprs[1:5,1:5]
#look at the probe intensities across the samples
boxplot(dat)
