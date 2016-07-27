
library(dplyr)
library(testthat)
library(statmod)
library(RColorBrewer)
library(ggplot2)
library(gridExtra)
library(mgcv)
library(plyr)
library(raster)
library(reshape2)
library(rgdal)
library(rgeos)
library(maptools)
library(sp)
library(spdep)
library(vegan)
library(cluster)

#  We use two different projection systems here.  This is the test to create the
#  base resolution.
if (model.proj == '+init=epsg:4326') {
  #  lat/long
  base.rast <- raster(xmn = -98.6, xmx = -66.1, ncols = 391,
                      ymn = 36.5,  ymx = 49.75, nrows = 160,
                      crs = '+init=epsg:4326')
  numbered.rast <- setValues(base.rast, 1:ncell(base.rast))
}

if (model.proj == '+init=epsg:3175') {
  base.rast <- raster(xmn = -71000, xmx = 2297000, ncols = 296,
                      ymn = 58000,  ymx = 1498000, nrows = 180,
                      crs = '+init=epsg:3175')
  numbered.rast <- setValues(base.rast, 1:ncell(base.rast))
}

if (model.proj == '+init=epsg:3175') {
  # These are narrower boundaries than the base raster for reasons associated with
  # the broader PalEON project.
  xylimits <- c(-100000, 1050000, 600000, 1600000)
}
if (model.proj == '+init=epsg:4326') {
  xylimits <- c(-98, -83, 42, 50)
}

#  This code requires you to set the directory to the directory in which you've unpacked the supplement.

source('R/misc.functionsv1.4.R')

#   Get the us & canadian lines into the right projection
umw.domain <- spTransform(readOGR('../../data/input/shapes/usa/us.shp', 'us'), 
                          CRSobj = CRS(model.proj)) %>% 
  fortify(region = 'STATE_NAME')

can.domain <- spTransform(readOGR('../../data/input/shapes/canada', 'PROVINCE'), 
                          CRSobj = CRS(model.proj)) %>% fortify(region = 'NAME')

#  This file was processed in the first file in this chain.
#  replace with 'minn.wisc.mich.clean_v1_2'

used.data <- readOGR('../../data/output//aggregated_midwest//minn.wisc.mich.clean_v1_8.shp', 'minn.wisc.mich.clean_v1_8')

if (is.na(proj4string(used.data))) {
  # Just in case there happens to be no projection information associated with the data.
  proj4string(used.data) <- CRS('+proj=longlat +ellps=WGS84')
}

#########################################################################
#  Clean the tree data:

diams <- used.data@data[,4:7]
angles <- used.data@data[,16:19]
dists <- floor(used.data@data[,8:11])
species <- apply(used.data@data[,12:15], 2, as.character)
species[is.na(species)] <- 'No tree'

#  Points within a township are either sections or quartersections.  This
#  is the list of points that are sections.  All others are quarter-sections.
sections <- c(2, 5, 8, 11, 14, 18, 21, 24, 27, 30,
              34, 37, 40, 43, 46, 50, 53, 56, 59, 62,
              66, 70, 74, 78, 82,
              87, 89, 91, 93, 95, 98, 100, 102, 104, 106, 108,
              109, 111, 113, 115, 117, 119, 122, 123, 124, 125, 126)

#  These are the points on the outside of each township.
external <- c(109:120, 97:108, 87, 89, 91, 93, 95, 122:126)

#  These correction values are derived empirically and are described in the supplementary material.
#  One issue right now is the lack of an empirical theta value.

corr.vals <- read.csv('../../data/input/relation_tables/cogbill_corrections.csv')

correction <- data.frame(kappa = rep(NA, length(used.data)),
                         theta = rep(NA, length(used.data)),
                         zeta  = rep(NA, length(used.data)),
                         phi   = rep(NA, length(used.data)))

plot.trees <- rowSums(!(species == 'Water' | species == 'NonTree'), na.rm = TRUE)

point.no <- as.numeric(as.character(used.data$Point))

#  So there are a set of classes here, we can match them all up:

internal <- ifelse(!point.no %in% external, 'internal', 'external')
trees    <- ifelse(plot.trees == 2, 'P', '2NQ')
section  <- ifelse(point.no %in% sections, 'section', 'quartersection')
state    <- ifelse(substr(used.data@data$Township, 1, 2) == 'wi', 'Wisconsin',
                   ifelse(substr(used.data@data$Township, 1, 2) == 'mi', 'Michigan', 'Minnesota'))

corr.year     <- as.character(used.data@data$year)
corr.year[state == 'Michigan'] <- 'all'
corr.year[state == 'Wisconsin' & corr.year %in% c('1832-1834', '1834-1846')] <- 1845
corr.year[state == 'Wisconsin' & !(corr.year %in% c('1832-1834', '1834-1846'))] <- 1907
corr.year[state == 'Minnesota' & (corr.year %in% (as.character(1847:1855)))] <- 1855
corr.year[state == 'Minnesota' & !(corr.year %in% (as.character(1847:1855)))] <- 1907

match.vec <- apply(corr.vals[,1:4], 1, paste, collapse = '')
to.match <- apply(data.frame(state, corr.year, internal, section, stringsAsFactors = FALSE), 1, paste, collapse = '')

correction <- corr.vals[match(to.match, match.vec),]

#  Get rid of variables that aren't used anymore.
rm(external, sections, plot.trees, 
   corr.year, trees, section, state, 
   match.vec, to.match)

##########################################################
#  A set of helper functions for the paper:
ninefive <- function(x, to.rast = TRUE) {
  # Generate 95% CIs:
  
  if (class(x) == 'RasterLayer') {
    data <- getValues(x)
    
    data.nf <- quantile(data, c(0.025, 0.975), na.rm = TRUE)
    
    new.x <- x
    new.x[x == 0] <- NA
    new.x[x < data.nf[1]] <- NA
    new.x[x > data.nf[2]] <- data.nf[2]
    
    if (to.rast) {
      out <- new.x
    }
    else{
      out <- data.nf
    }
  }
  else {
    data.sort <- sort(x)
    out <- c(data.sort[length(data.sort) * 0.025], 
                 data.sort[length(data.sort) * 0.975])
  }
  return(out)
}

p <- function(x, fixed = FALSE) {
  
  # Pretty up the numbers:
  testthat::expect_true(class(x) %in% c('numeric', 'integer'))
  
  if (x > 100 & !fixed) {
    #  Round to three significant figures and add a comma:
    x <- round(x, -(floor(log(x, 10)) - 3))
    y <- formatC(x, format = 'd', digits = 2, big.mark = ',')
  }
  
  if (x > 100 & fixed) {
    y <- formatC(x, format = 'd', big.mark = ',')
  }
  
  if (x < 100) {
    y <- as.numeric(formatC(x, format = 'g', digits = 3))
  }
  
  return(y)
}

add_cells <- function(pts) {
  
  # This adds the `cell` column to a table with x & y coordinates:
  # Some assertions about `pts`
  testthat::expect_true((class(pts) %in% c("SpatialPoints", "SpatialPointsDataFrame") | all(c('x', 'y') %in% colnames(pts)) | ncol(pts) == 2))
  
  # A function to add cell numbers to a data.frame from the base raster:
  base.rast <- raster(xmn = -71000, xmx = 2297000, ncols = 296,
                      ymn = 58000,  ymx = 1498000, nrows = 180,
                      crs = '+init=epsg:3175')
  
  base.rast <- setValues(base.rast, 1:ncell(base.rast))
  
  points <- raster::extract(x = numbered.rast, y = pts)
  
  # Some assertions:
  # If this returns NAs then something is out of the domain.
  testthat::expect_true(all(!is.na(points)))
  
  return(points)

}

# A function - `check_vars` that basically looks at everything in memory and tells 
# you where it comes from and how often it's used.
source('R/check_var_use.R')
