% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/location_slice_from_region_series.R
\name{location_slice_from_region_series}
\alias{location_slice_from_region_series}
\title{Extract local climate for one or more locations for a given time slice.}
\usage{
location_slice_from_region_series(
  x,
  time_bp = NULL,
  time_ce = NULL,
  coords = NULL,
  region_series,
  nn_interpol = TRUE,
  buffer = FALSE,
  directions = 8
)
}
\arguments{
\item{x}{a data.frame with columns x and y coordinates(see \code{coords} for
standard coordinate names, or
how to use custom ones), plus optional
columns \code{time_bp} or \code{time_ce} (depending on the units used) and
\code{name}. Alternatively, a vector of cell numbers.}

\item{time_bp}{used if no \code{time_bp} column is present in \code{x}: the dates in
years before present (negative
values represent time before present, i.e. 1950, positive values time in the future)
for each location.}

\item{time_ce}{time in years CE as an alternative to \code{time_bp}. Only one of
\code{time_bp} or \code{time_ce} should be used.}

\item{coords}{a vector of length two giving the names of the "x" and "y"
coordinates, as found in \code{data}. If left to NULL, the function will
try to guess the columns based on standard names \code{c("x", "y")}, \code{c("X","Y")},
\code{c("longitude", "latitude")}, or \code{c("lon", "lat")}}

\item{region_series}{a \code{\link[terra:SpatRaster-class]{terra::SpatRasterDataset}} obtained with \code{\link[=region_series]{region_series()}}}

\item{nn_interpol}{boolean determining whether nearest neighbour
interpolation is used to estimate climate for cells that lack such
information (i.e. they are under water or ice). By default, interpolation is only
performed from the first ring of nearest neighbours; if climate is not
available, NA will be returned for that location. The number of neighbours
can be changed with the argument \code{directions}. \code{nn_interpol} defaults to TRUE.}

\item{buffer}{boolean determining whether the variable will be returned
as the mean of a buffer around the focal cell. If set to TRUE, it overrides
\code{nn_interpol} (which provides the same estimates as \code{buffer} but only for
locations that are in cells with an NA). The buffer size is determined
by the argument \code{directions}. \code{buffer} defaults to FALSE.}

\item{directions}{character or matrix to indicate the directions in which
cells are considered connected when using \code{nn_interpol} or \code{buffer}.
The following character values are allowed: "rook" or "4" for the
horizontal and vertical neighbours; "bishop" to get the diagonal neighbours;
"queen" or "8" to get the vertical, horizontal and diagonal neighbours;
or "16" for knight and one-cell queen move neighbours. If directions
is a matrix it should have odd dimensions and have logical (or 0, 1) values.}
}
\value{
a data.frame with the climatic variables of interest.
}
\description{
This function extract local climate for a set of locations
at the appropriate times (selecting the closest time slice available for the
specific date associated with each location).
}
