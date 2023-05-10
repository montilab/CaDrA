
# Create a build argument
ARG R_VERSION
ARG R_VERSION=${R_VERSION:-4.3.0}

# Get shiny+tidyverse+devtools packages from rocker image
FROM rocker/shiny-verse:${R_VERSION}

# Set up the maintainer information
MAINTAINER Reina Chau (lilychau999@gmail.com)
    
# Define a system argument
ARG DEBIAN_FRONTEND=noninteractive

# Install system libraries of general use
RUN apt-get update && apt-get -y --no-install-recommends install \
    librsvg2-dev \
    libudunits2-dev \
    libv8-dev \
    libsodium-dev \
    libbz2-dev \
    liblzma-dev \
    tcl8.6-dev \
    tk8.6-dev 
    
# Install required bioconductor packages to run CaDrA
RUN R -e "BiocManager::install('GSVA')"
RUN R -e "BiocManager::install('SummarizedExperiment')"

# Install required R packages to run CaDrA
RUN R -e "install.packages('doParallel', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('gplots', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('misc3d', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('plyr', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('ppcor', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('R.cache', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('reshape2', dependencies=TRUE, repos='http://cran.rstudio.com/')"

# Install CaDrA package with all dependencies already installed previously
RUN R -e "devtools::install_github('montilab/CaDrA', dependencies=TRUE)"
