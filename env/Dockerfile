FROM mambaorg/micromamba:1.5.10

# NEED TO MAKE THE USER DECLARE THE PATH TO MAGMA
# wget https://ctg.cncr.nl/software/MAGMA/prog/magma_v1.10.zip
# unzip magma_v1.10.zip
# docker build -t neurobridge:1 . - loop in nf and check if already there if not install
# > Adapt SMR installation

USER root

RUN apt-get update && apt-get install -y --no-install-recommends \
    r-base r-base-dev \
    build-essential gfortran pkg-config \
    libcurl4-openssl-dev libssl-dev libxml2-dev \
    libgit2-dev \
    && rm -rf /var/lib/apt/lists/*

RUN R -q -e "install.packages(c('remotes','data.table','progress','tidyverse','readr','dplyr','coloc'), repos='https://cloud.r-project.org')" && \
    R -q -e "if(!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager', repos='https://cloud.r-project.org')" && \
    R -q -e "BiocManager::install(c('GenomicRanges','IRanges','S4Vectors','BiocParallel'), ask=FALSE, update=FALSE)" && \
    R -q -e "remotes::install_github('alexploner/cfdr.pleio', build_vignettes=FALSE)" && \
    R -q -e "remotes::install_github('zhenin/HDL', subdir='HDL', upgrade='never', build_vignettes=FALSE)" && \
    R -q -e "remotes::install_github('GenomicSEM/GenomicSEM', build_vignettes=FALSE)" && \
    R -q -e "remotes::install_github('josefin-werme/LAVA', build_vignettes=FALSE)" && \
    R -q -e "remotes::install_github('stephenslab/susieR', build_vignettes=FALSE)"

USER $MAMBA_USER
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba create -y -n neurobridge_v1 -f /tmp/environment.yml && micromamba clean -a -y
ENV PATH=/opt/conda/envs/neurobridge_v1/bin:$PATH
ENV CONDA_DEFAULT_ENV=neurobridge_v1

CMD ["bash"]
