% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_gencode_coords.R
\name{get_gencode_coords}
\alias{get_gencode_coords}
\title{Parse Gencode gtf files into a data.frame of ids and coordinates.}
\usage{
get_gencode_coords(gtf, feature = c("gene", "transcript", "exon"))
}
\arguments{
\item{gtf}{character string. Path to gtf file.}

\item{feature}{character string.}
}
\value{
data.frame
}
\description{
Parse Gencode gtf files into a data.frame of ids and coordinates.
}
