% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/time_bp.R
\name{time_bp}
\alias{time_bp}
\title{Extract time in years before present from SpatRaster}
\usage{
time_bp(x)
}
\arguments{
\item{x}{a \code{terra::SpatRaster}}
}
\value{
a date in years BP (where negative numbers indicate a date in the past)
}
\description{
A wrapper around \code{terra::time}, which converts time into years before
present
}
