ARG R_VERSION=4.5
FROM rocker/r-ver:${R_VERSION}

ARG R_VERSION
ENV R_VERSION=${R_VERSION}
ARG QUARTO_VERSION=1.10.19

# build time variables
ARG BUILD_DATE="000000"
ENV BUILD_DATE=${BUILD_DATE}
ARG BUILD_TAG="000000"
ENV BUILD_TAG=${BUILD_TAG}
ARG REPONAME="000000"
ENV REPONAME=${REPONAME}

RUN apt-get update && \
  apt-get install -y --no-install-recommends ca-certificates git pandoc && \
  apt-get clean && rm -rf /var/lib/apt/lists/*

RUN /rocker_scripts/install_quarto.sh "${QUARTO_VERSION}"

ENV PKG_INCLUDE_LINKINGTO=TRUE

# install R package
# --mount avoids baking a copy of the source into the image; rw lets the test
# suite write to the mount (discarded afterwards)
RUN --mount=type=bind,source=.,target=/opt/MOSuite,rw \
  Rscript -e 'install.packages("pak", repos = sprintf("https://r-lib.github.io/p/pak/stable/%s/%s/%s", .Platform$pkgType, R.Version()$os, R.Version()$arch))' && \
  Rscript -e 'pak::pkg_install(c("devtools", "readr", "tibble"))' && \
  Rscript -e "devtools::install_local('/opt/MOSuite', dependencies = TRUE, repos='http://cran.rstudio.com', upgrade='never')" && \
  Rscript -e "library(MOSuite); devtools::test('/opt/MOSuite')" && \
  Rscript -e "readr::write_csv(tibble::as_tibble(installed.packages()), '/data/r-packages_mosuite.csv')"

# add mosuite exec to the path
RUN Rscript -e 'executable <- system.file("exec", "mosuite", package = "MOSuite", mustWork = TRUE); Sys.chmod(executable, mode = "0755"); stopifnot(file.symlink(executable, "/usr/local/bin/mosuite"))' && \
  mosuite --help

# copy .Rprofile to R_HOME
RUN --mount=type=bind,source=.github/.Rprofile,target=/opt/mosuite.Rprofile \
  cat /opt/mosuite.Rprofile >> "${R_HOME}/etc/Rprofile.site"

# cleanup
WORKDIR /data2
RUN apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
CMD ["/bin/bash"]
