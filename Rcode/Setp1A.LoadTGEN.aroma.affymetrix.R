

#Setp1A.LoadTGEN.aroma.affymetrix

#clean memory
rm(list=ls())

#install and load package
source("https://www.bioconductor.org/biocLite.R")
if(!require(aroma.affymetrix)){
  biocLite("aroma.affymetrix",dependencies = T)
  require(aroma.affymetrix)
}

