FROM nciccbr/ccbr_ubuntu_24.04:v1

# build time variables
ARG BUILD_DATE="000000"
ENV BUILD_DATE=${BUILD_DATE}
ARG BUILD_TAG="000000"
ENV BUILD_TAG=${BUILD_TAG}
ARG REPONAME="000000"
ENV REPONAME=${REPONAME}

ARG R_VERSION=4.4
ENV R_VERSION=${R_VERSION}

# install conda packages
RUN mamba install -y -c conda-forge \
    r-base=${R_VERSION} \
    r-argparse \
    r-amap \
    r-broom \
    r-cffr \
    r-colorspace \
    r-dendextend \
    r-devtools \
    "r-ggplot2 <4.0.0" \
    r-ggrepel \
    r-gridExtra \
    r-here \
    r-matrix \
    r-mgcv \
    r-pak \
    r-patchwork \
    r-plotly \
    r-rcolorbrewer \
    r-readr \
    r-rmarkdown \
    r-survival \
    r-upsetr \
    r-viridis \
    bioconductor-annotationdbi \
    bioconductor-annotate \
    bioconductor-complexheatmap \
    bioconductor-delayedarray \
    bioconductor-edger \
    bioconductor-genomicranges \
    bioconductor-keggrest \
    bioconductor-limma \
    bioconductor-s4arrays \
    bioconductor-summarizedexperiment \
    bioconductor-sva \
  && conda clean -afy

# install R package
COPY . /opt2/MOSuite
RUN R -e "devtools::install_local('/opt2/MOSuite', dependencies = TRUE, repos='http://cran.rstudio.com', upgrade='never')" && \
  R -e "library(MOSuite); devtools::test('/opt2/MOSuite')" && \
  R -s -e "readr::write_csv(tibble::as_tibble(installed.packages()), '/data2/r-packages_mosuite.csv')"

# add mosuite exec to the path
RUN chmod -R +x /opt2/conda/lib/R/library/MOSuite/exec
ENV PATH="$PATH:/opt2/conda/lib/R/library/MOSuite/exec"
RUN mosuite --help

# copy .Rprofile to R_HOME
ADD .github/.Rprofile /opt2/conda/lib/R/

# copy example script & json to data
COPY ./inst/extdata/example_script.sh /data2/
COPY ./inst/extdata/json_args/ /data2/json_args/

# cleanup
WORKDIR /data2
RUN apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
