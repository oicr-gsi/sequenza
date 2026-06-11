FROM modulator:latest

MAINTAINER Fenglin Chen <f73chen@uwaterloo.ca>

# packages should already be set up in modulator:latest
USER root

# move in the yaml to build modulefiles from
COPY recipes/sequenza_recipe.yaml /modulator/code/gsi/recipe.yaml

# build the modules and set folder / file permissions
RUN ./build-local-code /modulator/code/gsi/recipe.yaml --initsh /usr/share/modules/init/sh --output /modules && \
	find /modules -type d -exec chmod 777 {} \; && \
	find /modules -type f -exec chmod 777 {} \;

# install required packages
RUN apt-get -m update && apt-get install -y gzip zip unzip

# add the user
RUN groupadd -r -g 1000 ubuntu && useradd -r -g ubuntu -u 1000 ubuntu
USER ubuntu

# copy the setup file to load the modules at startup
COPY .bashrc /home/ubuntu/.bashrc

# set environment paths for modules
#ENV BIOCONDUCTOR_ROOT="/modules/gsi/modulator/sw/Ubuntu18.04/bioconductor-3.8-rstats3.6"
#ENV RSTATS_ROOT="/modules/gsi/modulator/sw/Ubuntu18.04/rstats-3.6"
#ENV SEQUENZA_ROOT="/modules/gsi/modulator/sw/Ubuntu18.04/sequenza-2.1.2"
#ENV SEQUENZA_RES_ROOT="/modules/gsi/modulator/sw/Ubuntu18.04/sequenza-res-2.1.2"
#ENV SEQUENZA_SCRIPTS_ROOT="/modules/gsi/modulator/sw/Ubuntu18.04/sequenza-scripts-2.1.2"

#ENV PATH="/modules/gsi/modulator/sw/Ubuntu18.04/rstats-3.6/bin:/modules/gsi/modulator/sw/Ubuntu18.04/sequenza-scripts-2.1.2/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
#ENV MANPATH="/modules/gsi/modulator/sw/Ubuntu18.04/rstats-3.6/share/man"
#ENV LD_LIBRARY_PATH="/modules/gsi/modulator/sw/Ubuntu18.04/rstats-3.6/lib:/modules/gsi/modulator/sw/Ubuntu18.04/sequenza-2.1.2/lib:/modules/gsi/modulator/sw/Ubuntu18.04/sequenza-scripts-2.1.2/lib:/modules/gsi/modulator/sw/Ubuntu18.04/bioconductor-3.8-rstats3.6/lib"
#ENV R_LIBS_SITE="/modules/gsi/modulator/sw/Ubuntu18.04/rstats-3.6/lib/R/library:/modules/gsi/modulator/sw/Ubuntu18.04/sequenza-2.1.2/lib/R/library:/modules/gsi/modulator/sw/Ubuntu18.04/sequenza-scripts-2.1.2/lib/R/library:/modules/gsi/modulator/sw/Ubuntu18.04/bioconductor-3.8-rstats3.6/lib/R/library" 

CMD /bin/bash
