FROM bioconductor/bioconductor_docker:devel

WORKDIR /opt/mutseqr

# Dependency installation is cached until DESCRIPTION or the devel base image
# changes. The package source itself is mounted at runtime, so ordinary source
# edits never require rebuilding this image.
COPY DESCRIPTION ./DESCRIPTION

RUN Rscript -e 'options(Ncpus = max(1L, parallel::detectCores(logical = FALSE)), repos = BiocManager::repositories(version = "devel")); BiocManager::install(c("BiocCheck", "devtools", "rcmdcheck", "remotes"), ask = FALSE, update = FALSE); remotes::install_deps(".", dependencies = TRUE, upgrade = "never", repos = getOption("repos"))'

# R CMD check builds the package's PDF manual.
RUN apt-get update \
    && apt-get install --yes --no-install-recommends \
        texinfo \
        texlive-fonts-recommended \
        texlive-latex-base \
        texlive-latex-extra \
        texlive-latex-recommended \
    && rm -rf /var/lib/apt/lists/*
