FROM nciccbr/ccbr_ubuntu_22.04:v4

# build time variables
ARG BUILD_DATE="000000"
ENV BUILD_DATE=${BUILD_DATE}
ARG BUILD_TAG="000000"
ENV BUILD_TAG=${BUILD_TAG}
ARG REPONAME="000000"
ENV REPONAME=${REPONAME}

ARG R_VERSION=4.3.2
ENV R_VERSION=${R_VERSION}

SHELL ["/bin/bash", "-lc"]

# Install conda and give write permissions to conda folder
RUN echo 'export PATH=/opt2/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh" -O ~/miniforge3.sh && \
    /bin/bash ~/miniforge3.sh -b -p /opt2/conda && \
    rm ~/miniforge3.sh && chmod 777 -R /opt2/conda/
ENV PATH="/opt2/conda/bin:$PATH"

# Pin channels and update
RUN conda config --add channels conda-forge \
 && conda config --add channels bioconda \
 && conda config --set channel_priority strict

# install conda packages
RUN mamba install -y -c conda-forge \
    r-base=${R_VERSION} \
    r-devtools \
    r-ggplot2 \
    r-ggrepel r-viridis r-upsetr r-patchwork r-plotly \
    r-matrix r-mgcv r-survival \
    bioconductor-genomicranges \
    bioconductor-summarizedexperiment \
    bioconductor-delayedarray \
    bioconductor-s4arrays \
    bioconductor-annotationdbi \
    bioconductor-annotate \
    bioconductor-keggrest \
  && conda clean -afy

# install R package
COPY . /opt2/MOSuite
RUN R -e "devtools::install_local('/opt2/MOSuite', dependencies = TRUE, repos='http://cran.rstudio.com')"

# add mosuite exec to the path
RUN chmod -R +x /opt2/conda/lib/R/library/MOSuite/exec
ENV PATH="$PATH:/opt2/conda/lib/R/library/MOSuite/exec"
RUN mosuite --help

# copy example script & json to data
COPY ./inst/extdata/example_script.sh /data2/
COPY ./inst/extdata/json_args/ /data2/json_args/

# Save Dockerfile in the docker
COPY Dockerfile /opt2/Dockerfile_${REPONAME}.${BUILD_TAG}
RUN chmod a+r /opt2/Dockerfile_${REPONAME}.${BUILD_TAG}

# cleanup
WORKDIR /data2
RUN apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
