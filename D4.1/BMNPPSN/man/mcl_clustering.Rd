% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clustering.R
\name{mcl_clustering}
\alias{mcl_clustering}
\title{Performs MCL clustering with inference of the inflation parameter}
\usage{
mcl_clustering(resDirectory, outDirectory)
}
\arguments{
\item{resDirectory}{Path to the input file}

\item{outDirectory}{Path to the output directory}
}
\value{
Stores the MCL results as an attribute of the igraph network and the updated igraph network
}
\description{
This function loads .sif igraph networks present in the input resDirectory,
after inflation parameter tuning for the MCL clustering (such that the number of clusters
with more than 10 genes is the max), it stores the clustering for the optimized inflation value
as a vertex attributes for the igraph object, the updated igraph object in a .sif format, and the
session as an Rdata object.
}
