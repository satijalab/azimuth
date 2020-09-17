
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SeuratMapper v0.0.0.9005

<!-- badges: start -->

[![Lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://github.com/mojaveazure/seurat-mapper)
<!-- badges: end -->

SeuratMapper is a Shiny app demonstrating a query-reference mapping
algorithm for single-cell data. The reference data accompanying the app
and the algorithms used are described [in our
publication](Seurat%204%20biorxiv%20link).

We have made an instance of the app available for public use, [available
here](app%20landing%20page%20on%20lab%20website).

All the analysis and visualization functionality available in the app -
and much more - is available in [version 4 of the the Seurat R
package](Seurat4%20landing%20page%20on%20lab%20website).

## Installation

You can install SeuratMapper from GitHub with:

``` r
if (!requireNamespace('remotes', quietly = TRUE) {
  install.packages('remotes')
}
remotes::install_github('mojaveazure/seurat-mapper', ref = 'master')
```

## Running the app

The app is launched as:

``` r
SeuratMapper::AzimuthApp()
```

By default, the appropriate reference files are loaded into memory by
accessing a web URL. If you instead have a directory containing
reference files at `/path/to/reference`, specify it as:

``` r
SeuratMapper::AzimuthApp(reference = '/path/to/reference')
```

### Specifying options

You can set options by passing a parameter to the `AzimuthApp` function:

``` r
SeuratMapper::AzimuthApp(max.cells = 100000)
```

or setting the option in R (e.g. if it is not a parameter to the
`AzimuthApp` function):

``` r
options('Azimuth.map.pbcorthresh' = 0.5)
SeuratMapper::AzimuthApp()
```

## Docker

First, build the Docker image. Clone the repository and run the
following while in the root of the repository to build the image and tag
it with the name “seurat-mapper”:

    docker build -t seurat-mapper .

Next, launch a container based on the image `seurat-mapper` with a bind
mount mounting the directory on the host containing the reference files
(e.g. `/path/to/reference`) as `/reference-data` in the container.

    docker run -it -p 3838:3838 -v /path/to/reference:/reference-data:ro seurat-mapper

If port 3838 is already in use on the host or you wish to use a
different port, use `-p NNNN:3838` in the run command instead, to bind
port NNNN on the host to port 3838 on the container. The container runs
the command `R -e "SeuratMapper::AzimuthApp(reference =
'/reference-data')"` by default.

### Specifying options

You can set options by passing a parameter to the `AzimuthApp` function:

    docker run -it -p 3838:3838 -v /path/to/reference:/reference-data:ro seurat-mapper R -e "SeuratMapper::AzimuthApp(reference = '/reference-data', max.cells = 100000)"

or setting the option in R (e.g. if it is not a parameter to the
`AzimuthApp` function):

    docker run -it -p 3838:3838 -v /path/to/reference:/reference-data:ro seurat-mapper R -e "options('Azimuth.map.pbcorthresh' = 0.5)" -e "SeuratMapper::AzimuthApp(reference = '/reference-data')"

or just starting a shell in the container, from which you can launch an
interactive R session and set options as desired:

    docker run -it -p 3838:3838 -v /path/to/reference:/reference-data:ro seurat-mapper /bin/bash

## Run asynchronous mapping score

If the app is launched in an environment with `future::plan('multicore',
workers = 2)`, the mapping score will continue to compute in the
background once the visualizations become available. Please note: This
feature is unavailable when running in Rstudio, and we have observed
that it may not work as intended on Mac OS.

## Support

We do not actively support users running the app themselves and suggest
you use [version 4 of the Seurat
package](Seurat4%20landing%20page%20on%20lab%20website) to run the
reference mapping workflow and related visualizations on your local
system. Please see the [Seurat mapping vignette](link%20to%20vignette)
for an example of how to use Seurat for reference mapping. If you use
the instance of the app we are hosting on the web, you can download a
Seurat v4 R script once your analysis is complete that will guide you in
reproducing the analysis. You do not need seurat-mapper to reproduce the
analysis.

If you would like to help us improve the app, and you believe a dataset
meets the requirements and it is publicly available for us to use for
debugging but the app doesn’t work, please file a Github issue linking
to the dataset and describing the problem on [the issues
page](app%20repository%20issue%20page).
