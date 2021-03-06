% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/buildGRUmaps.R
\name{buildGRUmaps}
\alias{buildGRUmaps}
\title{Build Grouped Response Units in maps}
\usage{
buildGRUmaps(gruLoc, parsValues)
}
\arguments{
\item{gruLoc}{raster file that is comprised by numbers from 1 to the total number of \code{GRUs} that were defined.}

\item{parsValues}{data frame that has the values of the four parameters of each \code{GRU}. It must have equal number of
rows as number of GRU that were defined, and must have four columns which define the \code{alpha1}, \code{alpha2}, \code{d}
and \code{Smax} parameters.}
}
\value{
a list which consists of four vectors and four raster, each one of them has the values of a parameter spatialized according with
the GRU raster layer.
}
\description{
This function builds raster maps for each parameter based on a raster file where the location of the
Grouped Response Units (\code{GRUs}) are defined. This raster must have the same resolution as the forcing files
(i.e., for each cell that is planned to be simulated, there must be forcing time series and a cell assigned to a \code{GRU}).
}
\examples{
data(GRU)    
data(param)
gru_maps <- buildGRUmaps(GRU, param)


}
\author{
Nicolas Duque Gardeazabal <nduqueg@unal.edu.co> \cr
Pedro Felipe Arboleda Obando <pfarboledao@unal.edu.co> \cr
Carolina Vega Viviescas <cvegav@unal.edu.co> \cr
David Zamora <dazamoraa@unal.edu.co> \cr

Water Resources Engineering Research Group - GIREH
Universidad Nacional de Colombia - sede Bogota
}
