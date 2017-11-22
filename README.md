# CaDrA
Candidate Drivers Analysis: Multi-Omic Search for Candidate Drivers of Functional Signatures

CaDrA is an R package that allows one to query a binary genomic dataset (where the rows are 1/0 vectors indicating the presence/absence of genomic features such as somatic mutations or copy number alteration events) with an associated sample ranking (where the samples are ranked by a certain phenotypic readout of interest such as protein expression, pathway activity etc.) in order to determine a set of features that, together (based on their union or 'logical OR'), provide the best score associated with the observed ranking.

For more information, please see the associated manuscript (Kartha et al. CaDrA: A computational framework for performing candidate driver analyses using binary genomic features)

In order to get CaDrA up and running in R, you'll need to make sure some dependencies are pre-installed:

(1) Please install the following packages first before installing CaDrA:

cran.dependencies <- c("doParallel","plyr","ggplot2","gplots","gtable","R.cache","reshape2")

bioconduct.dependencies <- "Biobase"

install.packages(cran.dependencies)

source("https://bioconductor.org/biocLite.R")

biocLite(bioconduct.dependencies)

(2) Installing CaDrA

Once you have access to the source CaDrA package tar ball (see available tag releases), and have installed the above dependencies without error, you may install CaDrA:

install.packages("<source_file>",repos=NULL,type="source")

require(CaDrA)
