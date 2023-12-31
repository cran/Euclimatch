% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/climatch_plot.R
\name{climatch_plot}
\alias{climatch_plot}
\title{Plot or Create SpatRaster of Climatch Data}
\usage{
climatch_plot(climdat,
                     recipient,
                     source = NULL,
                     climatch = NULL,
                     provide_SpatRaster = FALSE,
                     xlim = terra::ext(recipient)[1:2],
                     ylim = terra::ext(recipient)[3:4],
                     plg = list(title = "Climatch", size=1),
                     xlab = expression(paste("Longitude (",degree,")")),
                     ylab = expression(paste("Latitude (",degree, ")")),
                     ...
                     )
}
\arguments{
\item{climdat}{A SpatRaster, RasterStack or RaterLayer of the climate data to extract.}

\item{recipient}{An object specifying location of where the recipient (i.e., target) region. Can be a SpatialPolygosDataFrame, SpatialPolygons, SpatVector.}

\item{source}{An object, like 'recipient', specifying the location of the source region.}

\item{climatch}{Vector of climatch values to use in creating SpatRaster of recipient.}

\item{provide_SpatRaster}{Logical. If TRUE then function returns SpatRaster object, if FALSE (default) return plot.}

\item{xlim}{Numeric, specify the limits of the x axis. Default is extent of x-axis from recipient SpatRaster.}

\item{ylim}{Numeric, specify the limits of the y axis. Default is extent of y-axis from recipient SpatRaster.}

\item{plg}{A list of parameters for specifying the legend. Default is "Climatch" see 'plot' in 'terra' for more documentation.}

\item{xlab}{Character for x axis label.}

\item{ylab}{Character for y axis label.}

\item{...}{Pass arguments to plot function.}
}
\value{
A plot of the climatch within the recipient region. A SpatRater if provide_SpatRaster is TRUE.
}
\description{
Create a plot or SpatRaster of climatch values within recipient region.
}
\examples{
r1 <- data.frame()
for(i in 1:100){r1 <- rbind(r1, runif(60))}
rclim1 <- terra::rast(as(r1, "matrix")) #Create the RasterLayer

# Dummy lon lat mimicking species occurrence records
spec_occ <- data.frame("lon" = 1:10, "lat" = 11:20)

# Create dummy polygons
x_coor <- c(1, 5,  10, 8, 3)
y_coor <- c(15, 20,  27, 30, 29)
dum_coor <- cbind(x_coor, y_coor)
dum_poly <- terra::vect(dum_coor, type = "polygon")

# Run and plot the climatch
climatch_plot(recipient = dum_poly, source = spec_occ, climdat = rclim1, xlab = "Lon", ylab = "Lat")
}
