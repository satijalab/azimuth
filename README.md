
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Azimuth v0.1.0

<!-- badges: start -->

[![Lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://github.com/satijalab/azimuth)
<!-- badges: end -->

Azimuth is a Shiny app demonstrating a query-reference mapping algorithm
for single-cell data. The reference data accompanying the app and the
algorithms used are described in the publication “Integrated analysis of
multimodal single-cell data” (Y. Hao, S. Hao, et al., bioRxiv 2020).

We have made an instance of the app available for public use, [described
here](https://satijalab.org/azimuth).

All the analysis and visualization functionality available in the app -
and much more - is available in version 4 of the the [Seurat R
package](https://satijalab.org/seurat).

## Installation

You can install Azimuth from GitHub with:

``` r
if (!requireNamespace('remotes', quietly = TRUE) {
  install.packages('remotes')
}
remotes::install_github('satijalab/azimuth', ref = 'master')
```

## Running the app

The app is launched as:

``` r
Azimuth::AzimuthApp()
```

By default, the appropriate reference files are loaded into memory by
accessing a web URL. If you instead have a directory containing
reference files at `/path/to/reference`, specify it as:

``` r
Azimuth::AzimuthApp(reference = '/path/to/reference')
```

### Downloading the PBMC reference

You can download the PBMC reference files that would be automatically
loaded by default using the following command:

    wget -m -R 'index.html*' -P pbmc -nd https://seurat.nygenome.org/references/pbmc/

### Specifying options

You can set options by passing a parameter to the `AzimuthApp` function:

``` r
Azimuth::AzimuthApp(maxcells = 100000)
```

or setting the option in R (e.g. if it is not a parameter to the
`AzimuthApp` function):

``` r
options('Azimuth.map.pbcorthresh' = 0.5)
Azimuth::AzimuthApp()
```

## Docker

First, build the Docker image. Clone the repository and run the
following while in the root of the repository to build the image and tag
it with the name “azimuth”:

    docker build -t azimuth .

Next, launch a container based on the image `azimuth` with a bind mount
mounting the directory on the host containing the reference files
(e.g. `/path/to/reference`) as `/reference-data` in the container.

    docker run -it -p 3838:3838 -v /path/to/reference:/reference-data:ro azimuth

If port 3838 is already in use on the host or you wish to use a
different port, use `-p NNNN:3838` in the run command instead, to bind
port NNNN on the host to port 3838 on the container. The container runs
the command `R -e "Azimuth::AzimuthApp(reference = '/reference-data')"`
by default.

### Specifying options

You can set options by passing a parameter to the `AzimuthApp` function:

    docker run -it -p 3838:3838 -v /path/to/reference:/reference-data:ro azimuth R -e "Azimuth::AzimuthApp(reference = '/reference-data', max.cells = 100000)"

or setting the option in R (e.g. if it is not a parameter to the
`AzimuthApp` function):

    docker run -it -p 3838:3838 -v /path/to/reference:/reference-data:ro azimuth R -e "options('Azimuth.map.pbcorthresh' = 0.5)" -e "Azimuth::AzimuthApp(reference = '/reference-data')"

or just starting a shell in the container, from which you can launch an
interactive R session and set options as desired:

    docker run -it -p 3838:3838 -v /path/to/reference:/reference-data:ro azimuth /bin/bash

## Support

We do not actively support users running the app themselves and suggest
you use version 4 of the [Seurat package](https://satijalab.org/seurat)
to run the reference mapping workflow and related visualizations on your
local system. Please see the [Seurat mapping
vignette](https://satijalab.org/seurat/reference_mapping.html) for an
example of how to use Seurat for reference mapping. If you use the
instance of the app we are hosting on the web, you can download a Seurat
v4 R script once your analysis is complete that will guide you in
reproducing the analysis. You do not need Azimuth to reproduce the
analysis.

If you would like to help us improve the app, and you believe a dataset
meets the requirements and it is publicly available for us to use for
debugging but the app doesn’t work, please file a Github issue linking
to the dataset and describing the problem on [the issues
page](https://github.com/satijalab/azimuth/issues).
