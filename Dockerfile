# Build according to a specified version of R
ARG R_VERSION
ARG R_VERSION=${R_VERSION:-4.3.0}

############# Build Stage: base ##################

# Get shiny+tidyverse packages from rocker image
FROM rocker/tidyverse:${R_VERSION} as base

# Install system libraries of general use
RUN apt-get update --allow-releaseinfo-change --fix-missing \
  && apt-get -y --no-install-recommends install \
  tk \
  libxtst6 \
  libxt6 \
  && apt clean autoclean \
  && apt autoremove --yes \
  && rm -rf /var/lib/{apt,dpkg,cache,log}/
 
# Build according to a specified release of CaDrA
ARG CADRA_BRANCH
ENV CADRA_BRANCH=${CADRA_BRANCH:-devel}

# Set working directory to install CaDrA
WORKDIR / 

# Clone CaDrA repo
RUN git clone https://github.com/montilab/CaDrA.git

# Set working directory to CaDrA
WORKDIR /CaDrA 

# Checkout the desired branch for the build
RUN git checkout ${CADRA_BRANCH} && git pull

# Install CaDrA denpendencies
RUN Rscript "/CaDrA/install_r_packages.R"

# Install all package dependencies
RUN Rscript "${PACKAGE_DIR}/install_r_packages.R"

# Install package
RUN Rscript -e 'pkg_dir <- Sys.getenv("PACKAGE_DIR"); devtools::install(pkg_dir, dependencies=TRUE)'

