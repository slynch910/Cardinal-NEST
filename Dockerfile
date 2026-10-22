# Cardinal NEST - Docker image
#
# Build:  docker build -t cardinal-nest .
# Run:    docker run --rm -p 8080:8080 -v /path/to/data:/data cardinal-nest
# Open:   http://localhost:8080

ARG R_VERSION=4.6.1
FROM docker.io/rocker/r-ver:${R_VERSION}

ENV DEBIAN_FRONTEND=noninteractive
# Rootless builds (podman/buildah userns) cannot restore archived file modes.
ENV TAR_OPTIONS="--no-same-owner --no-same-permissions"

RUN apt-get update \
 && apt-get install -y --no-install-recommends \
      make g++ gfortran curl ca-certificates cmake pandoc git \
      libuv1-dev zlib1g-dev \
      libcurl4-openssl-dev libssl-dev libxml2-dev \
      libfontconfig1-dev libfreetype6-dev \
      libpng-dev libjpeg-dev libtiff-dev \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/cardinal-nest

# Restore pinned packages first so the layer is cached across app code changes.
COPY renv.lock ./
RUN . /etc/os-release \
 && export CRAN_REPO="https://packagemanager.posit.co/cran/__linux__/${VERSION_CODENAME}/latest" \
 && export RENV_PATHS_CACHE=/tmp/renv-cache \
 && Rscript -e "install.packages(c('renv', 'BiocManager'), repos = Sys.getenv('CRAN_REPO'))" \
 && Rscript -e "options(repos = c(CRAN = Sys.getenv('CRAN_REPO'))); bioc <- renv::lockfile_read('renv.lock')\$Bioconductor\$Version; repos <- BiocManager::repositories(version = bioc); renv::restore(lockfile = 'renv.lock', library = .libPaths()[1], repos = repos, prompt = FALSE)" \
 && Rscript -e "library(Cardinal); library(shiny); library(bslib); library(shinyjs); library(shinyWidgets); library(shinycssloaders)" \
 && rm -rf /tmp/renv-cache /tmp/Rtmp* /root/.cache/R

COPY DESCRIPTION NAMESPACE ./src/
COPY R ./src/R
COPY man ./src/man
COPY inst ./src/inst
RUN R CMD INSTALL --no-docs src \
 && rm -rf src \
 && Rscript -e "library(CardinalNEST)"

RUN useradd --create-home --uid 1000 --shell /bin/bash nest \
 && mkdir -p /opt/cardinal-nest/logs /data \
 && chown -R nest:nest /opt/cardinal-nest /data

USER nest

ENV PORT=8080
EXPOSE 8080
VOLUME ["/data"]

CMD ["Rscript", "-e", "CardinalNEST::runCardinalNEST(host = '0.0.0.0', port = as.integer(Sys.getenv('PORT', '8080')), launch.browser = FALSE)"]
