% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/azimuth.R
\name{AzimuthReference}
\alias{AzimuthReference}
\title{Create a Seurat object compatible with Azimuth.}
\usage{
AzimuthReference(
  object,
  refUMAP = "umap",
  refDR = "spca",
  refAssay = "SCT",
  dims = 1:50,
  k.param = 31,
  plotref = "umap",
  plot.metadata = NULL,
  ori.index = NULL,
  colormap = NULL,
  assays = NULL,
  metadata = NULL,
  reference.version = "0.0.0",
  verbose = FALSE
)
}
\arguments{
\item{object}{Seurat object}

\item{refUMAP}{Name of UMAP in reference to use for mapping}

\item{refDR}{Name of DimReduc in reference to use for mapping}

\item{refAssay}{Name of SCTAssay to use in reference}

\item{dims}{Dimensions to use in reference neighbor finding}

\item{k.param}{Defines k for the k-nearest neighbor algorithm}

\item{plotref}{Either the name of the DimReduc in the provided Seurat object
to use for the plotting reference or the DimReduc object itself.}

\item{plot.metadata}{A data.frame of discrete metadata fields for the cells
in the plotref.}

\item{ori.index}{Index of the cells used in mapping in the original object on
which UMAP was run. Only need to provide if UMAP was run on different set of
cells.}

\item{colormap}{A list of named and ordered vectors specifying the colors and levels
for the metadata. See \code{\link{CreateColorMap}} for help
generating your own.}

\item{assays}{Assays to retain for transfer}

\item{metadata}{Metadata to retain for transfer}

\item{reference.version}{Version of the Azimuth reference}

\item{verbose}{Display progress/messages}
}
\value{
Returns a Seurat object with AzimuthData stored in the tools slot for
use with Azimuth.
}
\description{
Create a Seurat object compatible with Azimuth.
}
