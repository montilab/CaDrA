# Build according to a specified version of R
ARG R_VERSION
ARG R_VERSION=${R_VERSION:-4.3.0}

############# Build Stage: base ##################

# Get shiny+tidyverse packages from rocker image
FROM rocker/tidyverse:${R_VERSION} as base

# Install system libraries of general use
RUN apt-get update --allow-releaseinfo-change --fix-missing \
  && apt-get -y --no-install-recommends install \
  apt-utils \
  libxtst6 \
  libxt6 \
  tk \
  && apt clean autoclean \
  && apt autoremove --yes \
  && rm -rf /var/lib/{apt,dpkg,cache,log}/
 
# Create package directory 
ENV PACKAGE_DIR=/CaDrA 

# Make package as working directory
WORKDIR ${PACKAGE_DIR}

# Copy package code to Docker image
COPY . ${PACKAGE_DIR}

# Install all package dependencies
RUN Rscript "${PACKAGE_DIR}/install_r_packages.R"

# Install package
RUN Rscript -e 'pkg_dir <- Sys.getenv("PACKAGE_DIR"); devtools::install(pkg_dir, dependencies=TRUE)'

