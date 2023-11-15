
# Make sure BiocManager is installed
if (!require('BiocManager', quietly = TRUE))
  utils::install.packages('BiocManager')

# Load BiocManager package
library(BiocManager)

# Required bioconductor packages
bioconductor_pkgs <- c("SummarizedExperiment")

# Select only the packages that aren't currently installed.
install_bioconductor_lib <- bioconductor_pkgs[!bioconductor_pkgs %in% utils::installed.packages()]

# And finally we install the missing bioconductor packages
for(lib in install_bioconductor_lib) BiocManager::install(lib)

# Required R packages
cran_pkgs <- c('doParallel','ggplot2','gplots','graphics','grid','gtable','MASS','methods','misc3d','plyr','ppcor','R.cache','reshape2','stats')

# Select only the packages that aren't currently installed.
install_cran_lib <- cran_pkgs[!cran_pkgs %in% utils::installed.packages()]

# And finally we install the missing R packages, including their dependency.
for(lib in install_cran_lib) utils::install.packages(lib, dependencies = TRUE, repos='http://cran.rstudio.com/')





