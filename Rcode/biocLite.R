source("https://bioconductor.org/biocLite.R")
biocLite("BiocUpgrade")
biocLite("inSilicoMerging")
biocLite("inSilicoDb")
browseVignettes("inSilicoMerging")

library(inSilicoDb)
library(inSilicoMerging)
InSilicoLogin("rpackage_tester@insilicodb.com", "5c4d0b231e5cba4a0bc54783b385cc9a");
eset1 = getDataset("GSE19804", "GPL570", norm="FRMA", features = "gene", curation = 17470);
eset2 = getDataset("GSE10072", "GPL96",  norm="FRMA", features = "gene", curation = 17469);
esets = list(eset1,eset2);

table(pData(eset1)[,"Disease"]);
table(pData(eset2)[,"Disease"]);

library(inSilicoMerging);
eset_FRMA = merge(esets);
