
#' Plot or Create SpatRaster of Climatch Data
#'
#' Create a plot or SpatRaster of climatch values within recipient region.
#'
#' @param climdat A SpatRaster, RasterStack or RaterLayer of the climate data to extract.
#' @param recipient An object specifying location of where the recipient (i.e., target) region. Can be a SpatialPolygosDataFrame, SpatialPolygons, SpatVector.
#' @param source An object, like 'recipient', specifying the location of the source region.
#' @param climatch Vector of climatch values to use in creating SpatRaster of recipient.
#' @param provide_SpatRaster Logical. If TRUE then function returns SpatRaster object, if FALSE (default) return plot.
#' @param xlim Numeric, specify the limits of the x axis. Default is extent of x-axis from recipient SpatRaster.
#' @param ylim Numeric, specify the limits of the y axis. Default is extent of y-axis from recipient SpatRaster.
#' @param plg A list of parameters for specifying the legend. Default is "Climatch" see 'plot' in 'terra' for more documentation.
#' @param xlab Character for x axis label.
#' @param ylab Character for y axis label.
#' @param ... Pass arguments to plot function.
#'
#' @return A plot of the climatch within the recipient region. A SpatRater if provide_SpatRaster is TRUE.
#'
#' @import RColorBrewer
#' @import grDevices
#' @import terra
#'
#' @usage climatch_plot(climdat,
#'                      recipient,
#'                      source = NULL,
#'                      climatch = NULL,
#'                      provide_SpatRaster = FALSE,
#'                      xlim = terra::ext(recipient)[1:2],
#'                      ylim = terra::ext(recipient)[3:4],
#'                      plg = list(title = "Climatch", size=1),
#'                      xlab = expression(paste("Longitude (",degree,")")),
#'                      ylab = expression(paste("Latitude (",degree, ")")),
#'                      ...
#'                      )
#'
#' @examples
#'r1 <- data.frame()
#'for(i in 1:100){r1 <- rbind(r1, runif(60))}
#'rclim1 <- terra::rast(as(r1, "matrix")) #Create the RasterLayer
#'
#'# Dummy lon lat mimicking species occurrence records
#'species.occurr <- data.frame("lon" = 1:10, "lat" = 11:20)
#'
#'# Create dummy polygons
#'x.coor <- c(1, 5,  10, 8, 3)
#'y.coor <- c(15, 20,  27, 30, 29)
#'dummy_coordinates <- cbind(x.coor, y.coor)
#'dummy_polygon <- terra::vect(dummy_coordinates, type = "polygon")
#'
#'# Run and plot the climatch
#'climatch_plot(recipient = dummy_polygon, source = species.occurr, climdat = rclim1)
#' @export
climatch_plot <- function(climdat, recipient, source = NULL, climatch = NULL, provide_SpatRaster = FALSE,
                          xlim = terra::ext(recipient)[1:2], ylim = terra::ext(recipient)[3:4],
                          plg = list(title = "Climatch", size=1),
                          xlab = expression(paste("Longitude (",degree,")")),
                          ylab = expression(paste("Latitude (",degree, ")")),
                          ...){

  rec <- Euclimatch::extract_clim_data(climdat = climdat, locations = recipient, id = T) # Keep "cell" number for plotting

  if(is.null(source) == FALSE){
    sour <- Euclimatch::extract_clim_data(climdat = climdat, locations = source) # No need for "cell" in source
  }

  if(is.null(climatch)){
    gv <- apply(stats::na.omit(terra::values(climdat, dataframe = T)), 2, stats::var) # Calculate variance of variables

    # # Run the climate match
    cv <- Euclimatch::climatch_vec(recipient = rec, source = sour, globvar = gv) # climatch_vec runs for length globavar so no need to remove "cell" number column on the end
  }

  # If climatch is provided then assign to cv
  if(is.null(climatch) == FALSE){
    cv <- climatch
  }

  # Create an empty raster to fill with climatch data
  empty.ras <-rast(nrows=nrow(climdat), ncols=ncol(climdat),
                   crs = as.character(terra::crs(climdat)),
                   extent = terra::ext(climdat),
                   resolution = terra::res(climdat))

  # Create the raster with climatch
  climatch_raster <- terra::rasterize(x = terra::xyFromCell(empty.ras, cell=rec$cell),
                                      y = empty.ras, value = cv,
                                      source = source)

  # Logical condition if provide_raster == TRUE then return rast and nothing else
  if(provide_SpatRaster == TRUE){return(climatch_raster)}

  # Set the colours
  cols <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, 'RdYlBu')))(100)[c(floor(min(cv)*10)+1) : c(floor(max(cv)*10)+1)]

  # Create the plot
  plot(climatch_raster,
       xlim = xlim, ylim = ylim,
       col = cols,
       buffer = T,
       xlab = xlab,
       ylab = ylab,
       plg = plg,
       ...
  )

}
