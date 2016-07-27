library(raster)
library(rgdal)

#  To do all the michigan stuff we need to first download the quads and
#  process them:

mich.quads <- readOGR('X:/students/Former_Students/BenS/Shapefiles/Michigan_Sections/statewide_plss_sections.shp', 'statewide_plss_sections')

all.quads <- list.files('X:/students/Former_Students/BenS/DigitizingMI/quadPts', full.names=TRUE, recursive=TRUE)

all.quads <- all.quads[regexpr('gml$',all.quads)>0]

for(i in 1:length(all.quads)){

  file.name <- substr(all.quads[i],  
                      max((gregexpr('/', all.quads[i], fixed = TRUE))[[1]]) + 1,
                      regexpr('.', all.quads[i], fixed = TRUE)-1)
  
  output.loc <- paste0('C:/MichOut/', file.name, '.shp')
  
  command <- paste('ogr2ogr -f \"ESRI Shapefile\" ', output.loc, ' ', all.quads[i])
  command <- gsub('/', '\\', command, fixed = TRUE)
  
  try(system(command))
}

shapes <- list.files('C:/MichOut/', full.names = TRUE)

shapes <- shapes[regexpr('.shp', shapes, fixed = TRUE) > 0]
start.stop <- data.frame(start = sapply(gregexpr('/', shapes, fixed = TRUE), max) + 1,
                         stop = sapply(gregexpr('.', shapes, fixed = TRUE), max) - 1)

shape.layer <- substr(shapes, start.stop[,1], start.stop[,2])

for(i in 1:length(shapes)){
  aa <- try(readOGR(shapes[i], shape.layer[i]))
  
  
  if(i == 1 & length(aa) > 1){
    aa <- aa[,regexpr('_', names(aa), fixed = TRUE)<0]
    
    if(!'fid' %in% names(aa)) aa@data$fid <- NA
    
    mich.layer <- aa
    mich.layer$quad <- shape.layer[i]
  }
  if(i > 1  & length(aa) > 1){
    aa <- aa[,regexpr('_', names(aa), fixed = TRUE)<0]
    
    if(!'fid' %in% names(aa)) aa@data$fid <- NA
    
    aa$quad <- shape.layer[i]
    
    mich.layer <- rbind(mich.layer, 
                        aa)
  }

}

mich.test <- mich.layer
proj4string(mich.test) <- CRS("+proj=omerc +lat_0=45.30916666666667 +lonc=-86 +alpha=337.25556 +k=0.9996 +x_0=2546731.496 +y_0=-4354009.816 +datum=NAD83 +units=m +no_uoff +no_defs +ellps=GRS80 +towgs84=0,0,0")

mich.sub <- mich.test[!is.na(mich.test$species1),]

mich.repro <- spTransform(mich.sub, CRSobj = CRS('+init=epsg:3175'))

#  Some points overlap between southern MI and the northern region:  
#  We need to add those points to northern MI and remove them from southern MI:

michigan.old <- readOGR('data/raw_data/mich/michigan.shp', 'michigan')

michigan.re <- spTransform(michigan.old, CRSobj = CRS('+init=epsg:3175'))

min.dist <- function(x, xy){
  dist <- sqrt(x[1] - xy[,1])^2 + (x[2] - xy[,2])^2/1000
  data.frame(close = min(dist, na.rm=TRUE),
             point = which.min(dist))
}

#  mich.repro is the new southern data (and northern samples).  This code returns
#  a two column data frame with the distance in meters (not kilometers)
closest <- apply(coordinates(mich.repro), 1, min.dist, xy = coordinates(michigan.re))
close <- do.call(rbind.data.frame, closest)

#  So, we want to replace things in Michigan that don't have close matches
#  and don't have identified trees with the newly digitized points.  So, 
#  species 1 is NA and it's one of the points close to the newly digitized data:

# It's got to be within at least a kilometer and have a missing SPP
chooser <- subset(close, close < 400)$point[is.na(michigan.re@data$SPP1)[subset(close, close < 400)$point]]
new.match <- which(close$close < 400 & is.na(michigan.re@data$SPP1)[close$point])

#  It's also got to be in the same township and range:
same.twp <- as.numeric(substr(mich.repro@data$Township[new.match], 1, 2)) == as.numeric(as.character(michigan.re@data$twp[chooser]))
same.rng <- as.numeric(substr(mich.repro@data$Range[new.match], 1, 2)) == as.numeric(as.character(michigan.re@data$rng[chooser]))

add.values <- function(x){
  if(regexpr('AZIM|SPP', x[1])>0)michigan.re@data[,x[1]] <<- 
    as.character(michigan.re@data[,x[1]])
  
  michigan.re@data[chooser[same.twp & same.rng],x[1]] <<- 
    as.character(mich.repro@data[new.match[same.twp & same.rng],x[2]])
}

alignment <-   cbind(c('SPP1', 'DBH1', 'AZIMUTH',  'DIST1',
                       'SPP2', 'DBH2', 'AZIMUTH2', 'DIST2',
                       'SPP3', 'DBH3', 'AZIMUTH3', 'DIST3',
                       'SPP4', 'DBH4', 'AZIMUTH4', 'DIST4'),
                     c('species1', 'diam1', 'az1', 'dist1',
                       'species2', 'diam2', 'az2', 'dist2',
                       'species3', 'diam3', 'az3', 'dist3',
                       'species4', 'diam4', 'az4', 'dist4'))

apply(alignment, 1, add.values)

writeOGR(michigan.re, 
         'data/output/aggregated_midwest/michigan_filled/michigan_filled.shp',
         'michigan_filled', 'ESRI Shapefile', overwrite = TRUE)

#  now clean out the points in southern michigan that are close but are repeats:
new.match <- which(close$close < 1 & !is.na(michigan.re@data$SPP1)[close$point])

drop.mich <- !(1:nrow(mich.repro)) %in% new.match

mich.repro <- mich.repro[drop.mich,]

writeOGR(mich.repro, 
         'data/output/southern_MI/so_michigan.shp',
         'so_michigan', 'ESRI Shapefile', overwrite = TRUE)
