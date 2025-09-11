FROM rocker/shiny-verse:latest
LABEL maintainer="David Chisanga"
## -------------------------
## R Package Installation
## -------------------------
# Installing the R packages needed
# 'shiny' and its dependencies are already in the base image.
# We use a single RUN command for efficiency.
RUN R -e "install.packages(c('reticulate', 'pheatmap', 'DT'), \
    repos = 'https://cloud.r-project.org/')"

## -------------------------
## Miniconda Installation
## -------------------------
# Install Miniconda to manage Python environments.
# The rocker base is Debian, so we use apt-get and download the Linux installer.
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    bzip2 && \
    wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/miniconda && \
    rm ~/miniconda.sh && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Add conda to the system's PATH. This makes 'conda' command available directly.
ENV PATH="/opt/miniconda/bin:${PATH}"

## -------------------------
## Python Environment Setup
## -------------------------
# Define the conda environment name as a variable for easy reuse.
ARG CONDA_ENV=tcrdist3_env

# Create the conda environment and install the Python packages in a single step.
# This follows the logic of your R script.
RUN conda create -n ${CONDA_ENV} python=3.8 -c conda-forge -c bioconda -c defaults -y && \
    # Use 'conda run' to execute pip within the newly created environment
    conda run -n ${CONDA_ENV} pip install tcrdist3 && \
    # Clean up conda caches to minimize the final image size
    conda clean --all -f -y

## -------------------------
## Configure Reticulate
## -------------------------
# Set the RETICULATE_PYTHON environment variable. This tells the R 'reticulate'
# package exactly which Python executable to use by default, so you don't need
# to call use_condaenv() in your R scripts.
ENV RETICULATE_PYTHON="/opt/miniconda/envs/${CONDA_ENV}/bin/python"

RUN mkdir /root/tcrdist3_shiny 
RUN mkdir /root/results

COPY R /root/tcrdist3_shiny

EXPOSE 3838
CMD ["R", "-e", "shiny::runApp('/root/tcrdist3_shiny',host = '0.0.0.0', port=3838)"]