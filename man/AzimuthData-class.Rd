% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/azimuth.R
\docType{class}
\name{AzimuthData-class}
\alias{AzimuthData-class}
\alias{AzimuthData}
\title{AzimuthData}
\description{
The AzimuthData class is used to store reference info needed for Azimuth
}
\section{Slots}{

\describe{
\item{\code{plotref}}{DimReduc object containing UMAP for plotting and projection.
This should also contain the cell IDs in the misc slot}

\item{\code{colormap}}{Vector of id-color mapping for specifying the plots.}

\item{\code{seurat.version}}{Version of Seurat used in reference construction}

\item{\code{azimuth.version}}{Version of Azimuth used in reference construction}

\item{\code{reference.version}}{Version of the Azimuth reference}
}}

