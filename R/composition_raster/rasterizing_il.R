#  Just rasterizing for composition.
library(rgdal)
library(raster)

brugam <- read.csv('data/raw_data/illinois/brugam_il_pls_11-19-13_v1.csv', stringsAsFactor=FALSE)
arboretum <- read.csv('data/raw_data/illinois/arboretum_il_pls_11-19-13_v1.csv', stringsAsFactor=FALSE)
nd_il <- read.csv('data/raw_data/illinois/ndilpls_v1.2.csv', stringsAsFactor=FALSE)

#  Converting the ND data bit:
#  Clear leading and trailing spaces.
nd_il$L1.tree1 <- gsub('^[ \t]+|[ \t]+$', '', nd_il$L1.tree1)
nd_il$L1.tree1 <- gsub('[ ]{2}', ' ', nd_il$L1.tree1)
nd_il$L1.tree2 <- gsub('^[ \t]+|[ \t]+$', '', nd_il$L1.tree2)
nd_il$L1.tree2 <- gsub('[ ]{2}', ' ', nd_il$L1.tree2)

#  Conversion table load and convert.
taxon.conv <- read.csv('data/input/relation_tables/fullpaleon_conversion_v0.3_1.csv', header=TRUE, stringsAsFactors = FALSE)

#  Get rid of the wonky taxon conversions we know about:
taxon.conv$Level.3a[taxon.conv$Level.3a == '']       <- NA
taxon.conv$Level.3a[taxon.conv$Level.3a %in% '#N/A'] <- 'Other hardwood'

#  Create a long data frame by lat& long:
long.data <- data.frame(long = nd_il$PointX, #c(brugam$long, ),# arboretum$POINT_X),
                        lat  = nd_il$PointY, #c(brugam$Lat, ),# arboretum$POINT_Y),
                        SP1  = as.character(nd_il$L1.tree1),#c(brugam$Tree, )),# arboretum$Tree)),
                        SP2  = as.character(nd_il$L1.tree2),#c(rep(NA, nrow(brugam)), ),# rep(NA, nrow(arboretum)))),
                        stringsAsFactors=FALSE)

long.data$SP1[long.data$SP1 %in% c('No tree', 'No Tree', 'No Trees')] <- 'No tree'
long.data$SP2[long.data$SP2 %in% c('No tree', 'No Tree', 'No Trees')] <- 'No tree'

long.data$l3a.tree1 <-  taxon.conv$Level.3a[match(tolower(gsub('.', '',long.data$SP1, fixed = TRUE)), 
                                                  tolower(gsub('.', '',taxon.conv$Level.1, fixed = TRUE)))]

long.data$l3a.tree2 <-  taxon.conv$Level.3a[match(tolower(gsub('.', '',long.data$SP2, fixed = TRUE)), 
                                                  tolower(gsub('.', '',taxon.conv$Level.1, fixed = TRUE)))]


#  We then have to assign the points to appropriate raster cells:
coordinates(long.data) <- ~ long + lat

proj4string(long.data) <- CRS('+init=epsg:4326')

base.rast <- raster(xmn = -71000, xmx = 2297000, ncols=296,
                    ymn = 58000,  ymx = 1498000, nrows = 180,
                    crs = '+init=epsg:3175')

long.repro <- spTransform(long.data[,c('l3a.tree1', 'l3a.tree2')], CRS('+init=epsg:3175'))

rast.cell <- function(x){
  #  Assign the values for a single taxon to individual raster cells
  nc <- long.repro
  nc$test <- rowSums(long.repro@data == x, na.rm=TRUE)
  getValues(rasterize(nc, base.rast, field = 'test', fun = 'sum', na.rm=TRUE))
}

#  Get the unique list of taxa, excluding the NAs that get added from
#  the repetition of the Brugam and Arboretum data above:
taxa <- na.omit(as.character(unique(unlist(long.repro@data))))

rast.set <- sapply(taxa, rast.cell)
colnames(rast.set) <- taxa

rast.coord <- xyFromCell(base.rast, 1:ncell(base.rast))
rast.out <- cbind(rast.coord, rast.set)

write.csv(rast.out, 'data/output/gridded/IllinoisData_v0.2-2.csv')
