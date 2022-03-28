
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Azimuth v0.4.3

<!-- badges: start -->

[![Lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://github.com/satijalab/azimuth)
<!-- badges: end -->

Azimuth is a Shiny app demonstrating a query-reference mapping algorithm
for single-cell data. The reference data accompanying the app and the
algorithms used are described in the publication [“Integrated analysis
of multimodal single-cell data” (Y. Hao, S. Hao, et al., Cell
2021)](https://doi.org/10.1016/j.cell.2021.04.048).

We have made instances of the app available for public use, [described
here](https://azimuth.hubmapconsortium.org).

All the analysis and visualization functionality available in the app -
and much more - is available in version 4 of the the [Seurat R
package](https://satijalab.org/seurat).

## Installation

**Note**: you may need to update some packages prior to installing
Azimuth; from a fresh R session run:

``` r
update.packages(oldPkgs = c("withr", "rlang"))
```

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
reference files at `/path/to/reference` (directory must contain files
named `ref.Rds` and `idx.annoy`), specify it as:

``` r
Azimuth::AzimuthApp(reference = '/path/to/reference')
```

### Downloading the app reference files

You can download the reference files that would be automatically loaded
by default from Zenodo. Links are available on the Azimuth website
[here](https://azimuth.hubmapconsortium.org/references/).

### Specifying options

You can set options by passing a parameter to the `AzimuthApp` function.
Options in the `Azimuth.app` namespace (e.g. `max_cells` as shown in the
example below) can omit the “Azimuth.app.” prefix. Options in other
namespaces (e.g. `Azimuth.de.digits` as shown in the example below)
including non-Azimuth namespaces, must be specified using their full
name.

``` r
Azimuth::AzimuthApp(max_cells = 100000)


Azimuth::AzimuthApp('Azimuth.de.digits' = 5)
```

We also support reading options from a JSON-formatted config file.
Provide the path to the config file as the parameter `config` to
`AzimuthApp`. [Example config file](inst/resources/config.json). As
described above regarding setting options through parameters, the
“Azimuth.app.” prefix may be omitted.

``` r
Azimuth::AzimuthApp(config = 'config.json')
```

You can also set Azimuth or other options in R. (The full name must
always be specified, even for options in the Azimuth.app namespace.)

``` r
options('Azimuth.de.digits' = 5)
Azimuth::AzimuthApp()
```

Options can be set in any of these three ways simultaneously. Please
note that options set in R will be overwritten by the same option
specified in a config file, which are both overwritten by the same
option provided as a parameter to the `AzimuthApp` function.

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

### Rebuilding the Docker image more quickly in certain cases

The docker image takes about 20 minutes to build from scratch. To save
time, adding the argument `--build-arg SEURAT_VER=$(date +%s)` to the
`docker build` command will use cached layers of the image (if
available) and only reinstall Seurat and Azimuth (and not any of the
dependencies), which takes less than a minute. Alternatively, to only
reinstall Azimuth (and not Seurat or other dependencies) use the
argument `--build-arg AZIMUTH_VER=$(date +%s)`.

### Specifying options

You can set options by passing a parameter to the `AzimuthApp` function:

    docker run -it -p 3838:3838 -v /path/to/reference:/reference-data:ro azimuth R -e "Azimuth::AzimuthApp(reference = '/reference-data', max_cells = 100000)"

or providing the path to a config file (in this example, for
convenience, the config file is assumed to be in the reference directory
that is bind mounted to the container):

    docker run -it -p 3838:3838 -v /path/to/reference:/reference-data:ro azimuth R -e "Azimuth::AzimuthApp(config = '/reference-data/config.json', max_cells = 100000)"

or setting the option in R:

    docker run -it -p 3838:3838 -v /path/to/reference:/reference-data:ro azimuth R -e "options('Azimuth.map.pbcorthresh' = 0.5)" -e "Azimuth::AzimuthApp(reference = '/reference-data')"

or just starting a shell in the container, from which you can launch an
interactive R session and set options as desired:

    docker run -it -p 3838:3838 -v /path/to/reference:/reference-data:ro azimuth /bin/bash

## Support

We do not actively support users running the app themselves and suggest
you use version 4 of the [Seurat package](https://satijalab.org/seurat)
to run the reference mapping workflow and related visualizations on your
local system. Please see the [Seurat mapping
vignette](https://satijalab.org/seurat/articles/multimodal_reference_mapping.html)
for an example of how to use Seurat for reference mapping. If you use
the instance of the app we are hosting on the web, you can download a
Seurat v4 R script once your analysis is complete that will guide you in
reproducing the analysis. You do not need Azimuth to reproduce the
analysis.

If you would like to help us improve the app, and you believe a dataset
meets the requirements and it is publicly available for us to use for
debugging but the app doesn’t work, please file a Github issue linking
to the dataset and describing the problem on [the issues
page](https://github.com/satijalab/azimuth/issues).
