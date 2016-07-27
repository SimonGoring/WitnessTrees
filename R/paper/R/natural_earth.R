#  Loading in the Natural Earth datasets:
#  This works in general, except that

library(rgeos)

umw.domain$paleon <- factor(umw.domain$id %in% c('Michigan', 'Minnesota', 'Wisconsin'))

#  Load & clip the base raster.  This is broader than the clipping extent for the 
#  domain, and it helps speed up the use of the natural earth layers later in the
#  analysis.
#  Natural earth data is stored in the Maps folder, and comes from: www.naturalearthdata.com
#  Here we are using the 10m data.

nat.earth <- crop(stack('../../data/input/NaturalEarth/BaseRaster/NE1_HR_LC_SR_W_DR.tif'), 
                  y = extent(c(-98 - 1, -83 + 1, 42 - 1, 50 + 1)))

if (model.proj == '+init=epsg:4326') {
  ext <- c(-98, -83, 42, 50)
}

if (model.proj == '+init=epsg:3175') {
  ext <- c(-100000, 1050000, 600000, 1600000)
}

quick.subset <- function(x, y, ext){
  # Subsets a shapefile (x = file, y = layer) to a given extent.
  shape <- readOGR(x, y)
  
  shape@data$id <- rownames(shape@data)
  
  x.f = fortify(shape, region = "id")
  x.join = join(x.f, shape@data, by = "id")
  
  aa <- x.join[,1:2]
  coordinates(aa) <- ~ long + lat
  proj4string(aa) <- CRS(proj4string(nat.earth))
  aa <- coordinates(spTransform(aa, CRS(model.proj)))
  
  x.join[,1:2] <- aa
  
  x.subset <- subset(x.join, x.join$long > xylimits[1] & x.join$long < xylimits[2] &
                           x.join$lat > xylimits[3] & x.join$lat < xylimits[4])
  
  x.subset
}

lakes.subset <- quick.subset('../../data/input/NaturalEarth/Lakes/ne_10m_lakes.shp',
                             'ne_10m_lakes', ext)
river.subset <- quick.subset('../../data/input/NaturalEarth/Rivers/ne_10m_rivers_lake_centerlines.shp',
                             'ne_10m_rivers_lake_centerlines', ext)
coast.subset <- quick.subset('../../data/input/NaturalEarth/Coasts/ne_10m_coastline.shp',
                             'ne_10m_coastline', ext)

nat.crop <- crop(projectRaster(nat.earth, crs = model.proj), y = extent(ext))

rast.table <- data.frame(xyFromCell(nat.crop, 1:ncell(nat.crop)),
                         getValues(nat.crop / 255))

base.map <- ggplot() + 
  geom_path(data = river.subset, 
            aes(x = long, y = lat, group = group), color = 'blue', alpha = 0.1) +
  geom_polygon(data = lakes.subset, 
               aes(x = long, y = lat, group = group), fill = '#ADD8E6') +
  coord_equal(xlim = ext[1:2], ylim = ext[3:4]) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  xlab('') + ylab('')

rm(rast.table)
