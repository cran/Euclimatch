#' Extract Climate Data
#'
#' Extracts climate data from several types of inputs.
#'
#' @param climdat A SpatRaster, RasterStack or RaterLayer of the climate data to extract.
#' @param locations An object specifying location of where to extract the climate data from. Can be a SpatialPolygosDataFrame, SpatialPolygons, SpatVector, a single data.frame or a list of data.frames. If data.frame of list of data.frames, must provide only two columns with column names "lon" for longitude and "lat" for latitude.
#' @param id Choose to include cell numbers with climate data. 'TRUE' includes cell numbers, the default 'FALSE' does not include cell numbers. Used in the climatch_plot() function.
#'
#' @return A data.frame or list of data.frames of the extracted climate data.
#'
#' @import terra
#'
#' @usage extract_clim_data(climdat, locations, id = FALSE)
#'
#' @examples
#'# Create fake climate data in as a SpatRaster object
#'r1 <- data.frame()
#'for(i in 1:100){r1 <- rbind(r1, runif(60))}
#'rclim1 <- terra::rast(as(r1, "matrix")) #Create the RasterLayer
#'rclim2 <- c(rclim1, rclim1) # Create the stack
#'
#'# Dummy lon lat data i.e., species occurrences. Cols must be labelled "lon" and "lat"
#'# Cols must be labelled "lon" and "lat"
#'species.occurr <- data.frame("lon" = 1:10, "lat" = 11:20)
#'
#'# Create dummy polygons
#'x.coor <- c(1, 5,  10, 8, 3)
#'y.coor <- c(15, 20,  27, 30, 29)
#'dummy_coordinates <- cbind(x.coor, y.coor)
#'dummy_polygon <- terra::vect(dummy_coordinates, type = "polygon")
#'dummy_polygon2 <- rbind(dummy_polygon, dummy_polygon)
#'
#'# Extract the dummy data
#'# Extract dummy lon lat data
#'extract_clim_data(climdat = rclim2, locations = species.occurr)
#'# Extract dummy SpatVector with single polygon
#'extract_clim_data(climdat = rclim2, locations = dummy_polygon2)
#' @export

extract_clim_data <- function(climdat, locations, id = FALSE){

  # If the climdat is a 'RasterLayer' or 'RasterStack' convert to 'SpatRaster'
  if(class(climdat)[[1]] %in% c("RasterLayer", "RasterStack")){climdat <- terra::rast(climdat)}

  ##################################################
  # Extract climate data if 'locations' is a 'SpatialPolygonsDataFrame' or 'SpatialPolygons' object
  if(class(locations)[[1]] %in% c("SpatialPolygonsDataFrame", "SpatialPolygons")){

    n <- length(locations) # Number of SPDFs in locations

    # Extract climate data if 'locations' if they a 'SpatialPolygonsDataFrame' or "SpatialPolygons"
    if(n > 1){
      clim.dataframe.list <- list() # Create empty list

      # Loop for extracting climate data for each polygon and set them in list
      for(i in 1:n){
        #Extracts the climate data and puts them in data.frame and removes NAs
        if(id == TRUE){clim.dataframe.list[[i]] <- stats::na.omit(terra::extract(climdat, terra::vect(locations[i,]), cells = TRUE, ID = FALSE))} # [,-1] remove ID column
        if(id == FALSE){clim.dataframe.list[[i]] <- stats::na.omit(terra::extract(climdat, terra::vect(locations[i,]), ID = FALSE))}
      }
      return(clim.dataframe.list)
    }

    else{ # If there is only one polygon then extract and return a single data.frame
      clim.dataframe <- data.frame()
      if(id == TRUE){clim.dataframe <- stats::na.omit(terra::extract(climdat, terra::vect(locations), cells = TRUE, ID = FALSE))}
      if(id == FALSE){clim.dataframe <- stats::na.omit(terra::extract(climdat, terra::vect(locations), ID = FALSE))}
      return(clim.dataframe)
    }
  }

  ##################################################
  # Extract climate data if 'locations' is a 'SpatVector' object containing polygon(s)
  if(class(locations)[[1]] %in% c("SpatVector")){

    n <- length(locations) # Number of Polygons in locations

    # Extract climate data
    if(n > 1){
      clim.dataframe.list <- list() # Create empty list

      # Loop for extracting climate data for each polygon and set them in list
      for(i in 1:n){
        if(id == TRUE){clim.dataframe.list[[i]] <- stats::na.omit(terra::extract(climdat, locations[i,], cells = TRUE, ID = FALSE))}
        if(id == FALSE){clim.dataframe.list[[i]] <- stats::na.omit(terra::extract(climdat, locations[i,], ID = FALSE))}
      }
      return(clim.dataframe.list)
    }

    else{ # If there is only one polygon then extract and return a single data.frame
      clim.dataframe <- data.frame()
      if(id == TRUE){clim.dataframe <- stats::na.omit(terra::extract(climdat, locations, cells = TRUE, ID = FALSE))}
      if(id == FALSE){clim.dataframe <- stats::na.omit(terra::extract(climdat, locations, ID = FALSE))}
      return(clim.dataframe)
    }
  }

  ##################################################
  # Extract climate data if 'locations' is a 'list' of data.frames or a single "data.frame" with columns of longitude and latitude
  if(class(locations)[[1]] %in% c("list", "data.frame")){

    if(class(locations)[[1]] == "list"){
      clim.dataframe.list <- list() # Create empty list
      n <- length(locations)
      for(i in 1:n){
        locations.temp <- terra::vect(locations[[i]], geom=c("lon", "lat"), crs=as.character(terra::crs(climdat))) # convert lon lat to spatial points
        if(id == TRUE){clim.dataframe.list[[i]] <- stats::na.omit(terra::extract(climdat, locations.temp, cells = TRUE, ID = FALSE))} # Extract the climate data and apply to list
        if(id == FALSE){clim.dataframe.list[[i]] <- stats::na.omit(terra::extract(climdat, locations.temp, ID = FALSE))}
      }
      return(clim.dataframe.list)
    }

    else{
      clim.dataframe <- data.frame()
      locations.temp <- terra::vect(locations, geom=c("lon", "lat"), crs=as.character(terra::crs(climdat))) # convert long lat to spatial points
      if(id == TRUE){clim.dataframe <- stats::na.omit(terra::extract(climdat, locations.temp, cells = TRUE, ID = FALSE))} # Extract the climate data
      if(id == FALSE){clim.dataframe <- stats::na.omit(terra::extract(climdat, locations.temp, ID = FALSE))}
      return(clim.dataframe)
    }
  }
}
