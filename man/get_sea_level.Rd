% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_sea_level.R
\name{get_sea_level}
\alias{get_sea_level}
\title{Get sea level estimate}
\usage{
get_sea_level(time_bp)
}
\arguments{
\item{time_bp}{the time of interest}
}
\value{
a vector of sea levels in meters from present level
}
\description{
This function returns the estimated sea level from Spratt et al. 2016, using
the long PC1. Sea levels are from contemporary sea level (note that the original
data are with reference to the sea level during the Holocene ~5k year ago).
}
\keyword{internal}
