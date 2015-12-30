#
# Step1C.LoadLOAD.R
# 
# load the example TGEN data
#

#load library
# clear environment
rm(list=ls())
source("http://www.bioconductor.org/biocLite.R")

#install and load package
if(!require(minfi)){
  biocLite("minfi")
}

if(!require(stringr)){
  install.packages("stringr")
  library(stringr)
}


#read path

path <- "./Rcode/SampleRawData/LOAD/IDAT/"
file_name <- list.files(path)
file_name
class(file_name)
file_name <- str_replace(file_name,".idat","")
file_name
file_name <- str_replace(file_name,"_Red","")
file_name
file_name <- str_replace(file_name,"_Grn","")
file_name
base_name <- unique(file_name)
base_name
base_name_path <- paste(path,base_name,sep = "")
base_name_path
#laod data

all_data <- read.450k(base_name_path, verbose = TRUE)
all_data
pData(all_data)

dim(getRed(all_data))
getGreen(all_data)
