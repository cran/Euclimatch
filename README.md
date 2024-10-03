#  Euclimatch

<!-- badges: start -->
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version-last-release/Euclimatch)](https://cran.r-project.org/package=Euclimatch)
[![cran checks](https://badges.cranchecks.info/summary/Euclimatch.svg)](https://cran.r-project.org/web/checks/check_results_Euclimatch.html)
[![Downloads last.mnth](https://cranlogs.r-pkg.org/badges/Euclimatch)](https://cran.r-project.org/package=Euclimatch)
[![metacran downloads](https://cranlogs.r-pkg.org/badges/grand-total/Euclimatch)](https://cran.r-project.org/package=Euclimatch)
<!--[![CRAN checks](https://cranchecks.info/badges/worst/Euclimatch)](https://cranchecks.info/pkgs/Euclimatch) -->
<!-- badges: end -->

## Introduction

<img align="right"  src="https://github.com/JustinHubbard/R.package.media.JH/blob/main/Euclimatch_sticker.jpg?raw=true" width="20%"/>

The `Euclimatch` package provides the Euclidean “Climatch” algorithm [1] in an R environment to deliver versatility in the use of user-defined historical or projected climate data (e.g., WorldClim [2], MERRAclim [3], Global Climate Models [4]) and geographic location data (e.g., species occurrence records typically based on longitude and latitude, ecoregions, watersheds, global administrative areas), in climate matching.

Climate matching is a method is used in biological risk assessment frameworks, such as horizon scanning and invasive species risk assessment tools (e.g., Freshwater Invasiveness Scoring Kit; FISK [5]), or independent research [6,7,8,9], to estimate non-native species survival in recipient (target) non-native regions, or the survival of species transported in trade pathways among regions.

The `Euclimatch` package also provides functions to assist in climate data extraction and visualizations of climate match data, and offers the use of parallelization to maximize processing speed of larger datasets. To further quicken processing, the `climatch_vec()` function, the engine of the package, which runs the “Climatch” algorithm, was coded in `C++` and integrated with the Rcpp package [10]. This package imports and relies on the `terra` package [11] for working with the spatial data, such as extraction, and `foreach` [12] and `doParallel` [13] for parallel computing.

## Euclimatch Functions
Function name   |  Description
|---   |  -----------|
`extract_clim_data()` |  Extracts the climate data of single or multiple locations
`climatch_vec()`      |  Runs “Climatch” algorithm, provides vector of climatch score (0-10) for each grid cell in the recipient region
`climatch_sum()`  |  Provides a summary climatch score of the percentage of grid cells within recipient region(s)
`climatch_plot()` |  Provides a plot of `climatch_vec()`, using `terra::plot()`, or a `SpatRaster` object that can then be used in a different visualization package e.g., `ggplot2`, `tmap`, `rasterVis`
`climatch_par()`  |  Runs `climatch_vec()` or `climatch_sum()` in parallel for faster computing
---

## Installation
For Windows operating systems a recent verion of Rtools is required to compile the C++ code. See https://cran.r-project.org/bin/windows/Rtools/
```
# Install from CRAN
install.packages(“Euclimatch”)

# Load the library
library(Euclimatch)
```
## Examples

Here, we provide several examples of the use of Euclimatch in climate matching.

We begin the workflow with loading other helpful packages, climate data as a SpatRaster or raster object (e.g., .tif), location data as SpatVector, SpatialPolygonsDataFrame, SpatialPolygons (e.g., .shp), or longitude and latitude points as a dataframe (e.g., .csv) or list of dataframes. Here, we use Freshwater Ecoregions of the World [14], downloaded at [15], for the recipient regions. For our source region, we use the species occurrence records of the Oscar (*Astronotus ocellatus*) drawn from gbif [16]. Select climatic variables were based on [17], though it is important to note other variable sets are also used [6,7,8,9].


```
# If 'terra', 'geodata' or 'dplyr' are  not installed, install then load climate and spatial data
#install.packages(c('terra', 'geodata','dplyr'))
library(terra) # Spatial data
library(geodata) # Download climate data
library(dplyr) # For %>% operator

# set the projection
proj <- c("+proj=lonlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") 

# Load location data:
# Freshwater ecoregions of the world GIS shapefile data: https://www.feow.org/download
feow <- terra::vect("~/feow.shp", crs = proj)
feow <- feow[1:426,] # Trim ecoregions to described ecoregions

# Load and clean species occurrence data e.g., for Oscar (Astronotus ocellatus) 
# from GBIF.org (24 June 2023) GBIF Occurrence Download https://doi.org/10.15468/dl.yzjcv7
ast_oce_occurrences <- read.csv("~/astronotus.ocellatus.csv", header=T)

# Remove occurrences with invalid precision, datum, count, and coordinate mismatch.
ast_oce_occur <- ast_oce_occurrences %>% filter(!grepl("COORDINATE_PRECISION_INVALID|GEODETIC_DATUM_INVALID|CONTINENT_COORDINATE_MISMATCH|INDIVIDUAL_COUNT_INVALID", "issue"))
# IMPORTANT # Retain only lon and lat columns
ast_oce_occur <- data.frame("lon" = ast_oce_occur[,"decimalLongitude"], "lat" = ast_oce_occur[,"decimalLatitude"]) 

# Climate data from worldclim - see also https://www.worldclim.org/
hist_clim <- worldclim_global(var = "bio", res = 10, version="2.1", path=tempdir())
terra::crs(hist_clim) <- proj
# Climate projections
CanESM_370_2070 <- cmip6_world(model = "CanESM5", ssp = "370", time = "2061-2080", var = "bio", res = 10, version="2.1", path=tempdir())
terra::crs(CanESM_370_2070) <- proj

# Lets select variables predictive of freshwater animal survival [17]
freshwater_var <- c(1,6,9,11)
hist_clim <- hist_clim[[freshwater_var]]
CanESM_370_2070 <- CanESM_370_2070[[freshwater_var]]
```

With data loaded, we can calculate the global variances of our climate data, which are used in the “climatch” algorithm [1]. First, we rename our variable columns so that they are the same across climate data objects if we are using multiple climate data sets. Then, compute the variances with `var()` and `apply()`. For climate matching with climate change projections and historical data sets we can use the mean values between variances of each climate data set.

```
#Rename variables
names(hist_clim) <- c("bioclim1", "bioclim6","bioclim9","bioclim11")
names(CanESM_370_2070) <- c("bioclim1", "bioclim6","bioclim9","bioclim11")

# IMPORTANT # Compute global variance of climate data to use in climate match
gv_hist <- apply(na.omit(terra::values(hist_clim, dataframe = T)), 2, var) 
```

Here, we can use `extract_clim_dat()` to extract climate data from points or polygons.

```
# Extract climate data
ast_oce_hist <- extract_clim_data(climdat = hist_clim, locations = ast_oce_occur) # Oscar occurrences, historical climate
gl_hist <- extract_clim_data(climdat = hist_clim, locations = feow[16,]) # Recipient Laurentian Great Lakes, historical climate
gl_fut <- extract_clim_data(climdat = CanESM_370_2070, locations = feow[16,]) # Recipient Laurentian Great Lakes, future climate
```

Now we can run our climate matching. The `climatch_vec()` function returns a vector of climatch scores for each grid cell in the recipient region.

```
# Both recipient and source historical
climatch_hist <- climatch_vec(recipient = gl_hist, source = ast_oce_hist, globvar = gv_hist)
# Recipient under climate change and source historical
climatch_fut <- climatch_vec(recipient = gl_fut, source = ast_oce_hist, globvar = gv_hist)
```

Plot the climate matches with the `climatch_plot()` function. Provide the original recipient region spatial object and the climate data SpatRaster or raster so the function can first create a SpatRaster of the climatch data.

```
par(mfrow =c(1,2)) # Plots in 1 row two columns
climatch_plot(recipient = feow[16,], climdat = hist_clim, climatch = climatch_hist)
plot(feow[16,], add=T) # Add outiline of our recipient region. 
climatch_plot(recipient = feow[16,], climdat = hist_clim, climatch = climatch_fut)
plot(feow[16,], add=T)
```
<p align="center">
	<img src="https://github.com/JustinHubbard/R.package.media.JH/blob/main/climatch_plot().png?raw=true" width="85%"/>
	<figcaption> Figure 1. Climate match of Oscar (*Astronotus ocellatus*) under historical (left) and global climate model CanESM SSP3-7.0 2070 (right) to the Laurentian Great Lakes.</figcaption>


If we want to utilize different packages, such as tmap, levelplot, ggplot or rasterVis for visualizations, we can use the `climatch_plot()` function but with the argument and command `provide_raster = TRUE` to return a SpatRaster without the plot.

```
library(rasterVis) # Load the plotting library
library(RColorBrewer)
# Create the SpatRaster
climatch_raster <- climatch_plot(recipient = feow[16,], climdat = hist_clim, climatch = climatch_hist, provide_SpatRaster = TRUE)

ex.fe <- terra::ext(feow[16,]) # Grab extent for defining xlim and ylim

plotTheme <- rasterTheme(region=rev(brewer.pal(10,"RdBu")))
cmatch_lplot <- levelplot(climatch_raster, margin=F, par.settings=plotTheme,
                          xlim = ex.fe[1:2], ylim = ex.fe[3:4], maxpixels = 2e20)
cmatch_lplot

```
<p align="center">
	<img src="https://github.com/JustinHubbard/R.package.media.JH/blob/main/levelplot.jpeg?raw=true" width="85%"/>
	<figcaption>Climate match of Oscar Oscar (_Astronotus ocellatus_) under historical climatic conditions to the Laurentian Great Lakes.</figcaption>
</p>


If global variance can be grabbed from the climate data and a source is provided then `climatch_plot()` can  act as a wrapper for `climate_vec()`, and create the plot/SpatRaster.

```
par(mfrow =c(1,1)) # Change back to 1 row 1 column of plotting
climatch_plot(recipient = feow[16,], source = ast_oce_hist, climdat = hist_clim)
plot(feow[16,], add=T)
```

To perform an assessment where climate match scores need to be summarized across regions, we can use the `climatch_sum()` function. There are currently two score types that can be specified with the `type` argument, the default “perc” command summarizes the climate match scores as the percentage of climatch scores within a recipient >= a threshold of 6 (default), that can be changed (e.g., 4,5,7) see [6,7]. “mean” provides the mean climatch score across the region e.g., in [9].

```
# Summarize the climatch score within the Laurentian Great Lakes as percent match >= 6
climatch_sum(recipient = gl_hist, source = ast_oce_hist, globvar = gv_hist)
climatch_sum(recipient = gl_fut, source = ast_oce_hist, globvar = gv_hist)
```

To determine the climate match across a combination of single or multiple recipient and source regions, we can use the same function but supply a list of data.frames for the regions. Once data are extracted, run the climate match with `climatch_sum()` and assign the values to the SpatVector using ‘$’. ‘tmap’ may give a warning with the SpatVector so we can first convert to a SpatialPolygonsDataFrames with `terra::as()` then plot. Plotting here is slower than the climate matching and will likely take several minutes

```
# Extract data for all ecoregions
feow_rec <- extract_clim_data(climdat = hist_clim, locations = feow) 
feow$feow_rec_cm <- climatch_sum(recipient = feow_rec, source = ast_oce_hist, globvar = gv_hist)
feow_spdf <- as(feow, 'Spatial') # tmap seems to prefer SPDF over SpatVector

# Create the plots
feow_spdf_plot <- tm_shape(feow_spdf) +
  tm_graticules(alpha =.5, lwd = .03, labels.size = 0.8, n.y = 5, n.x = 5) + # Plots the grid i.e., graticules
  tm_fill("feow_rec_cm", title = "Climatch (%)", # Fills the polygons
          palette = "seq",
          style = "cont", 
          breaks = seq(0, 100, by = 10),
          legend.show = T) +
  tm_layout(aes.palette = list(seq = "-RdYlBu"),legend.outside = TRUE) + # Defines colour palette and location of legend
  tm_borders("black", lwd = .005, alpha=1) # Places an outer box.

feow_spdf_plot
```
<p align="center">
	<img src="https://github.com/JustinHubbard/R.package.media.JH/blob/main/feow_spdf_plot.png?raw=true" width="85%"/>
	<figcaption>Climate match of Oscar (*Astronotus ocellatus*) under historical climatic conditions to the global freshwater ecoregions.</figcaption>
</p>


If there is a need to run hundreds or thousands of climate matches, e.g. for a horizon scanning or screening assessment of many species and depending on the resolution of the climate data, it may be desirable to run the climate matches in parallel across multiple CPUs to quicken the process. Parallel processing is available with the `climatch_par()` function. This function operates like the others and can use either the `climatch_sum()` or `climatch_vec()` functions. The `ncores` argument specifies the number of CPU cores to utilize (default is 1); be careful not to use too many cores otherwise your computer may crash. Typically, a maximum of 1 less than the number of cores in your computer is recommended.

As an example, we match multiple source ecoregions in the Neotropics to recipient ecoregions of Canada (any ecoregion that intersects). The same thing could be done with species’ occurrence records. In this example, we will parallelize each computation with three cores.

Then, we can visualize the results as a sort of “hot spot” analysis where we plot the number of matches above a certain percent climatch score for each recipient ecoregion. In a study of non-native freshwater fishes introduced to the Great Lakes, Howeth et al. [7] found that a climate match 71.7% was predictive of successful species establishment following introduction; therefore, we use 71.7% as our threshold to categorize high matches. We can also investigate how climate change may affect this by repeating the assessment with historical and projected climates. See Hubbard et al. [8] for details on these combinations when using regional sources (e.g., ecoregions, watersheds, GADM) instead of species occurrences.

```
# The vector of freshwater ecoregions ids that intersect with Canada. The ids from result of raster::intersect(canada, feow) (quite slow) where Canada is a polygon on of Canada
can_feow_id <- c("142", "118", "119", "116", "117", "114", "110", "113", "120", "109", "107", "108", "103", "104", "111", "105", "102", "106", "112", "115", "101")
canada_feow <- feow[feow$FEOW_ID %in% can_feow_id,]

# Create vector of ids for ecoregions in the neotopics https://www.feow.org/ecoregions/list?page=7
neotropic_feow_id <- as.character(c(201:216, 301:312))
neotropic_feow <- feow[feow$FEOW_ID %in% neotropic_feow_id,]

# Extract the data for historical 
feow_canada_hist <- extract_clim_data(climdat = hist_clim, locations = canada_feow)
feow_neotropic_hist <- extract_clim_data(climdat = hist_clim, locations = neotropic_feow) 

# Extract for climate change
feow_canada_fut <- extract_clim_data(climdat = CanESM_370_2070, locations = canada_feow)
feow_neotropic_fut <- extract_clim_data(climdat = CanESM_370_2070, locations = neotropic_feow) 

# Run the climate matches 3 ways with historical and climate projections see Hubbard et al.[8] for details on these combinations
feow_can_match_hist <- climatch_par(recipient = feow_canada_hist, source = feow_neotropic_hist, globvar = gv_hist, ncores = 3, type = "perc", threshold = 6)
feow_can_match_hist_fut <- climatch_par(recipient = feow_canada_fut, source = feow_neotropic_hist, globvar = gv_hist, ncores = 3, type = "perc", threshold = 6)
feow_can_match_fut_fut <- climatch_par(recipient = feow_canada_fut, source = feow_neotropic_fut, globvar = gv_hist, ncores = 3, type = "perc", threshold = 6)

# We can use 'apply' with a simple function to count the number of matches > 71.7% or whatever threshold we decide
# And assign to our canada_feow SpatVector
canada_feow$count_can_hist <- apply(feow_can_match_hist, 2, function(x){sum(x>71.7)})
canada_feow$count_can_hist_fut <- apply(feow_can_match_hist_fut, 2, function(x){sum(x>71.7)})
canada_feow$count_can_fut_fut <- apply(feow_can_match_fut_fut, 2, function(x){sum(x>71.7)})

# Convert to spdf
canada_feow <- as(canada_feow, 'Spatial')

# Now plot these values
canada_feow_plot <- tm_shape(canada_feow, bbox = c(-177, 39, -45, 83)) +
  tm_graticules(alpha =.5, lwd = .03, labels.size = 0.8, n.y = 5, n.x = 5) +
  tm_fill(c("count_can_hist", "count_can_hist_fut", "count_can_fut_fut"),  
          palette = "seq", legend.show = T) +
  tm_borders("black", lwd = .005, alpha=1) +
  tm_layout(panel.show = T,
            panel.labels = c("Both Historical", "Canada Future to Neotropic Historical", "Both Future"),
            panel.label.size = 0.8, 
            aes.palette = list(seq = "-RdYlBu"), 
            legend.outside = TRUE)

canada_feow_plot
```
<p align="center">
	<img src="https://github.com/JustinHubbard/R.package.media.JH/blob/main/canada_feow_plot.png?raw=true" width="85%"/>
	<figcaption> The number of Neotropical freshwater ecoregions with a climate match >71.7% to freeshwater ecoregions of Canada under historical climatic conditions and climate change projection CanESM SSP3-7.0 2070.</figcaption>
</p>

## How to cite this package

To cite package ‘Euclimatch’ in publications use:

  Hubbard JAG, Drake DAR, Mandrak NE (2023). “Euclimatch: Euclidean Climatch Algorithm in
  R.” R package version 1.0.0, <https://CRAN.R-project.org/package=Euclimatch>.

A BibTeX entry for LaTeX users is

  @Misc{,
    title = {Euclimatch: Euclidean Climatch Algorithm in R},
    author = {Justin A. G. Hubbard and D. Andrew R. Drake and Nicholas E. Mandrak},
    year = {2023},
    note = {R package version 1.0.1},
    url = {https://CRAN.R-project.org/package=Euclimatch},
  }


### References
[1] Crombie, J., Brown, L., Lizzio, J., and Hood, G. 2008. Climatch user manual. Aust. Gov. Bur. Rural Sci.: 16.

[2] Fick, S. E., & Hijmans, R. J. (2017). WorldClim 2: new 1-km spatial resolution climate surfaces for global land areas. International Journal of Climatology, 37(12), 4302–4315. https://doi.org/10.1002/joc.5086

[3] Vega, G.C., Pertierra, L.R., and Olalla-Tárraga, M.Á. 2017. MERRAclim, a high-resolution global dataset of remotely sensed bioclimatic variables for ecological modelling. Sci. Data 4(June). https://doi.org/10.1038/sdata.2017.78.

[4] https://worldclim.org/data/cmip6/cmip6_clim10m.html

[5] Copp, G.H. 2013. The Fish Invasiveness Screening Kit (FISK) for non-native freshwater fishes-A summary of current applications. Risk Anal. 33(8): 1394–1396. https://doi.10.1111/risa.12095.

[6] Bomford, M., Barry, S. C., & Lawrence, E. (2010). Predicting establishment success for introduced freshwater fishes: A role for climate matching. Biological Invasions, 12(8), 2559–2571. https://doi.org/10.1007/s10530-009-9665-3

[7] Howeth, J.G., et al. 2016. Predicting invasiveness of species in trade: climate match, trophic guild and fecundity influence establishment and impact of non-native freshwater fishes. Divers. Distrib. 22(2): 148–160. https://doi.org/10.1111/ddi.12391.

[8] Hubbard, J.A.G., Drake, D.A.R., and Mandrak, N.E. 2023. Estimating potential global sources and secondary spread of freshwater invasions under historical and future climates. Divers. Distrib. (00): 1–11. https://doi.10.1111/ddi.13695.

[9] Britton, J.R., et al. 2010. Non-native fishes and climate change: Predicting species responses to warming temperatures in a temperate region. Freshw. Biol. 55(5): 1130–1141. https://doi.10.1111/j.1365-2427.2010.02396.x.

[10] Eddelbuettel D, François R (2011). Rcpp: Seamless R and C++ Integration. Journal of Statistical Software, 40(8), 1–18.https://doi.org/10.18637/jss.

[11] Hijmans R (2023). terra: Spatial Data Analysis. R package version 1.7-39, https://CRAN.R-project.org/package=terra.

[12] Microsoft, Weston S (2022). foreach: Provides Foreach Looping Construct. R package version 1.5.2, https://CRAN.R-project.org/package=foreach.

[13] Corporation M, Weston S (2022). doParallel: Foreach Parallel Adaptor for the 'parallel' Package. R package version 1.0.17,
  https://CRAN.R-project.org/package=doParallel.

[14] Abell, R., et al. 2008. Freshwater ecoregions of the world: A new map of biogeo-graphic units for freshwater biodiversity conservation. Bioscience 58(5): 403–414. https://doi.org/10.1641/B580507.

[15] https://www.feow.org/download

[16] GBIF.org (24 June 2023) GBIF Occurrence Download https://doi.org/10.15468/dl.yzjcv7

[17] Bradie, J., Pietrobon, A., and Leung, B. 2015. Beyond species-specific assessments: an analysis and validation of environmental distance metrics for non-indigenous species risk assessment. Biol. Invasions 17(12): 3455–3465. https://doi.org/10.1007/s10530-015-0970-8.

