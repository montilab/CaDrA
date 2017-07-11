# CaDrA
Candidate Drivers Analysis: Multi-Omic Search for Candidate Drivers of Functional Signatures

(1) Please install the following packages first before installing CaDrA:

cran.dependencies <- c("doParallel","plyr","ggplot2","gplots","gtable","R.cache","reshape2")

bioconduct.dependencies <- "Biobase"

install.packages(cran.dependencies)

source("https://bioconductor.org/biocLite.R")

biocLite(bioconduct.dependencies)

(2) Installing CaDrA

Once you have access to the source CaDrA package tar ball, and have installed the above dependencies without error, you may install CaDrA:

install.packages("<source_file>",repos=NULL,type="source")

require(CaDrA)
