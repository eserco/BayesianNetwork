FROM rocker/verse:latest

# environmental variable necesarry to run the container
ENV PASSWORD gyara2

RUN apt-get update && apt-get install -y \
  build-essential \
  libcurl4-gnutls-dev \
  libssl-dev 

# Expose the necessary ports rstudio runs on 8787 on the container
EXPOSE 8787
RUN mkdir -p /setup/R
COPY install_packages.R /setup/R 
RUN Rscript /setup/R/install_packages.R

#VOLUME /src /home/rstudio/scripts
