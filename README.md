Metagenomics with Bioconductor
================
Nicholas Cooley
2026-10-19

- [Introduction](#introduction)
- [Workshop information](#workshop-information)
- [Setup](#setup)
- [General](#general)

## Introduction

Bioconductor has tools for working with metagenomic data. This workshop
is designed to introduce those tools to data science students who have a
strong background in informatics generally and R, but may not have a
background in computational biology.

## Workshop information

**Learning goals.** By the end of this workshop you should be able to:

- describe how MGnify results are organised in a `MultiAssayExperiment`;
- explore taxonomic composition, diversity, and between-sample
  differences;
- explore functional profiles (GO-slim, GO terms, InterPro) and explain
  why their resolution matters;
- compare taxonomic and functional views of the same samples.

**Prerequisites.** Basic R; some familiarity with
`SummarizedExperiment`-style objects is helpful but not required.

**Packages used.** `MGnifyR`, `MultiAssayExperiment`,
`TreeSummarizedExperiment`, `mia`, `miaViz`, `scater`, `UpSetR`,
`ComplexHeatmap`, `vegan`, `ggplot2`.

## Setup

You can participate using Orchestra in your browser or R on your own
computer. You only need to choose one option.

### Option 1: Use Orchestra in your browser

An Orchestra environment is being prepared for this workshop. The
organisers will confirm when it is available and provide the template
name.

Once available:

1.  Go to [Orchestra workshop
    templates](https://app.orchestraplatform.org/templates) and sign in
    with your Google account.
2.  Find the metagenomics workshop template.
3.  Click **Launch** to start your workshop session.

The environment will include the required packages and data. You will
not need to install R or RStudio on your computer or run the local
installation instructions below.

### Option 2: Use R on your own computer

Please complete these steps before the workshop, as installing the
packages can take some time.

Use R 4.6.x and a recent version of RStudio. These instructions use
Bioconductor release 3.23, the same version as the single-cell workshop.

Run each section below in order in the R console.

#### 1. Install the workshop package and its dependencies

``` r
## Check R version
if (getRversion() < "4.6.0")
  stop(
    "This workshop needs R 4.6.0 or later (you have ",
    getRversion(), "). Please update R and re-run."
  )

## Install BiocManager
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

## Use Bioconductor release 3.23
BiocManager::install(version = "3.23", ask = FALSE)
options(repos = BiocManager::repositories(version = "3.23"))

## Install remotes and usethis
install.packages(c("remotes", "usethis"))

## Install the workshop package, dependencies and bundled data
remotes::install_github(
  "npcooley/MetagenomicsWithBioconductor@devel",
  dependencies = TRUE,
  upgrade = "always"
)
```

Here, `@devel` refers to the workshop repository’s branch, not the
development version of Bioconductor. The instructions above select
Bioconductor release 3.23.

#### 2. Download the workshop notebooks

``` r
usethis::use_course(
  "https://github.com/npcooley/MetagenomicsWithBioconductor/archive/refs/heads/devel.zip",
  destdir = "~"
)
```

Answer any prompts before continuing. The workshop notebook is
`vignettes/MetagenomicsWithBioconductor.Rmd` in the downloaded project.

#### 3. Check that the packages and data load

``` r
## Load the packages used in the workshop
pkgs <- c(
  "MGnifyR",
  "MultiAssayExperiment",
  "TreeSummarizedExperiment",
  "mia",
  "miaViz",
  "scater",
  "UpSetR",
  "ComplexHeatmap",
  "vegan",
  "ggplot2"
)
invisible(lapply(pkgs, library, character.only = TRUE))

## Load the bundled workshop data
data_path <- system.file(
  "extdata", "workshop_data.RData",
  package = "MetagenomicsWithBioconductor",
  mustWork = TRUE
)
load(data_path, verbose = TRUE)

## Inspect the available experiments
names(mse)

## Check the Bioconductor installation
BiocManager::valid()
```

The data-loading step should load objects called `tse` and `mse`. The
`names(mse)` command should list four experiments: `Taxonomic analysis`,
`go-slim`, `go-terms` and `interpro-identifiers`.

The final check should return `TRUE`. If you get an error or a report of
incompatible packages, share the output with the organisers.

The data are included in the workshop package. You do not need to run
the MGnify download examples. For the hands-on analysis, start at
**Exploring the data** in the [workshop
materials](https://npcooley.github.io/MetagenomicsWithBioconductor/articles/MetagenomicsWithBioconductor.html).

## General

A full end-to-end metagenomics experiment would include:

- sample collection
- sample preparation
- sequencing
- quality control
- assembly
- taxonomic and functional assignments
- informatics analysis

This workshop will cover only the last step. The brief discussion of
`MGnify` analysis pipeline versions refers to the second-to-last step.

The following code shows how the workshop data were retrieved. The data
are already included in the workshop package, so participants do not
need to run this code.

``` r
library(MGnifyR)
# increase R's default timeout limit of server calls, just because
options(timeout = 600)
# use a cache in R's temp directory
mg <- MgnifyClient(useCache = TRUE,
                   cacheDir = file.path(tempdir(),
                                        "MGnifyR_cache"))

# danish waste water collection study from copenhagen
study_id <- "MGYS00001135"

analysis_ids <- searchAnalysis(mg, "studies", study_id)

# if you haven't changed R's default timeout option, this may complain
analysis_meta_data <- getMetadata(x = mg,
                                  accession = analysis_ids)

# this is an idiomatic call to tapply; without specifying `FUN`
# tapply will return integer ids for the shared factor combinations,
# in this case, latitude and longitude
# it can contain NAs when NAs are present in a supplied INDEX
sample_locations <- tapply(X = analysis_meta_data,
                           INDEX = list(analysis_meta_data$`sample_geographic location (longitude)`,
                                        analysis_meta_data$`sample_geographic location (latitude)`))

# metadata can be messy, we know ahead of time that this data can contain NAs,
# but this is a form you would use to check if you didn't
# know that ahead of time
# any(is.na(analysis_meta_data$`sample_geographic location (latitude)`))
# any(is.na(analysis_meta_data$`sample_geographic location (longitude)`))

# we want to select the sampling location with the fewest samples, for simplicity
# table the sample location index, find the lowest value, get the name
# that in this case is just an integer
w1 <- as.integer(names(which.min(table(sample_locations))))

# which() manages NAs for us
analysis_subset <- analysis_meta_data[which(sample_locations == w1), ]

# getting a single tree summarized experiment
tse <- getResult(x = mg,
                 accession = analysis_subset[["analysis_accession"]][1],
                 get.func = FALSE)

# getting a multiassayexperiment containing the tses for several accessions
mse <- getResult(x = mg,
                 accession = analysis_subset[["analysis_accession"]])
```
