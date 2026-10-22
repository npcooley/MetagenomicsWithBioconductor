FROM bioconductor/bioconductor_docker:RELEASE_3_23

WORKDIR /home/rstudio

COPY --chown=rstudio:rstudio . /home/rstudio/

RUN Rscript -e "options(repos = c(CRAN = 'https://cran.r-project.org')); BiocManager::install(ask=FALSE, version = '3.23')"

RUN Rscript -e "options(repos = BiocManager::repositories(version = '3.23')); devtools::install('.', dependencies=TRUE, build_vignettes=TRUE)"
