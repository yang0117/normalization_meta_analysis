#load cel data with oligo

rm(list=ls())
#load package
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("oligo")

require(oligo)

tgen_path = "./Rcode/SampleRawData/TGEN/CEL_FILES_AFFY6/"

file_vec <- list.celfiles(tgen_path,full.name=T)
dat <- read.celfiles(file_vec)
dat
exprs(dat)[1:4,]
