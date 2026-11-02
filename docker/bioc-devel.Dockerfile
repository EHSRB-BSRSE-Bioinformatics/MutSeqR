FROM bioconductor/bioconductor_docker:devel

WORKDIR /opt/mutseqr

# Dependency installation is cached until DESCRIPTION or the devel base image
# changes. The package source itself is mounted at runtime, so ordinary source
# edits never require rebuilding this image.
COPY DESCRIPTION ./DESCRIPTION

RUN Rscript -e 'options(Ncpus = max(1L, parallel::detectCores(logical = FALSE)), repos = BiocManager::repositories(version = "devel")); BiocManager::install(c("BiocCheck", "devtools", "rcmdcheck", "remotes"), ask = FALSE, update = FALSE); remotes::install_deps(".", dependencies = TRUE, upgrade = "never", repos = getOption("repos"))'

# Generate documentation with exactly the version specified by the package.
RUN Rscript -e 'tools <- c("pkgdown", "rmarkdown"); missing <- tools[!vapply(tools, requireNamespace, logical(1), quietly = TRUE)]; if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org"); expected <- read.dcf("DESCRIPTION")[1, "Config/roxygen2/version"]; remotes::install_version("roxygen2", version = expected, repos = "https://cloud.r-project.org", upgrade = "never")'

# R CMD check builds the package's PDF manual.
RUN apt-get update \
    && apt-get install --yes --no-install-recommends \
        libgsl-dev \
        pandoc \
        texinfo \
        texlive-fonts-recommended \
        texlive-latex-base \
        texlive-latex-extra \
        texlive-latex-recommended \
    && rm -rf /var/lib/apt/lists/*

# Avoid the optional Inconsolata font while still building the indexed manual.
ENV R_RD4PDF=times,hyper
