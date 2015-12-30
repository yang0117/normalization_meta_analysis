#load cel file
rm(list=ls())
source("https://bioconductor.org/biocLite.R")
biocLite("affy")
biocLite("pd.genomewidesnp.6")
biocLite("genomewidesnp6cdf")
biocLite("hgu95a.db")
biocLite("genomewidesnp6cdf")
install.packages("genomewidesnp6cdf")

require(pd.genomewidesnp.6)
require(affy)
require(limma)
require(hgu95a.db)
require(annotate)

rm(list=ls())

#read TARC
tarc_path = "./SampleRawData/TARC/CEL files/"
dat <- ReadAffy(celfile.path = tarc_path)
dat@cdfName <- "hgu95av2"
object.size(dat)
sampleNames(dat)
ptest1 <- phenoData(dat)

#read TGEN
tgen_path = "./SampleRawData/TGEN/CEL_FILES_AFFY6/"
dat <- ReadAffy(celfile.path = tgen_path)



