# Get shiny+tidyverse+devtools packages from rocker image
FROM rocker/shiny-verse:4.0.3

# Set up the maintainer information
MAINTAINER Reina Chau (lilychau999@gmail.com)
    
# Set up a volume directory
VOLUME /srv/shiny-server/   

# Set up working directory to the app
WORKDIR /srv/shiny-server/

# Define a system argument
ARG DEBIAN_FRONTEND=noninteractive

# Install system libraries of general use
RUN apt-get update && apt-get install -y \
    libudunits2-dev \
    libv8-dev \
    libsodium-dev \
    python-dev \
    libbz2-dev \
    liblzma-dev
  
# Install the required bioconductor packages to run the app
RUN R -e "install.packages('BiocManager', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install('Biobase')"

# Install the required R packages to run the app
RUN R -e "install.packages('doParallel', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('plyr', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('ggplot2', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('grid', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('gplots', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('gtable', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('R.cache', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('reshape2', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('magrittr', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('purrr', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('stats', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('MASS', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('misc3d', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('ppcor', dependencies=TRUE, repos='http://cran.rstudio.com/')"

# Install additional packages for shiny applications
RUN R -e "install.packages('shinyjs', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('shinycssloaders', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('DT', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('plotly', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('heatmaply', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('RcolorBrewer', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('promises', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('future', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('cachem', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('ipc', dependencies=TRUE, repos='http://cran.rstudio.com/')"

# Install K2taxonomer and its dependencies
RUN R -e "devtools::install_github('montilab/CaDrA', dependencies=TRUE)"

# Install packages for xposome-api
RUN R -e "install.packages('unix', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('plumber', dependencies=TRUE, repos='http://cran.rstudio.com/')"

# Make the ShinyApp available at port 3838
EXPOSE 3838

# Copy configuration files to Docker image
COPY shiny-server.sh /usr/bin/shiny-server.sh

# Allow permission
RUN ["chmod", "+rwx", "/srv/shiny-server/"]
RUN ["chmod", "+x", "/usr/bin/shiny-server.sh"]

# Execute the app
CMD ["/usr/bin/shiny-server.sh"]
