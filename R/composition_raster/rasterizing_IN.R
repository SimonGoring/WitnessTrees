#  Just rasterizing for composition.
library(rgdal)
library(raster)

indiana <- read.csv('data/raw_data/indiana/ndinpls_v1.3.csv')

#  We need to check to make sure the closest two trees are the first two.

#  Clear leading and trailing spaces.
indiana$L1.tree1 <- gsub('^[ \t]+|[ \t]+$', '', indiana$L1.tree1)
indiana$L1.tree1 <- gsub('[ ]{2}', ' ', indiana$L1.tree1)
indiana$L1.tree2 <- gsub('^[ \t]+|[ \t]+$', '', indiana$L1.tree2)
indiana$L1.tree2 <- gsub('[ ]{2}', ' ', indiana$L1.tree2)

#  Conversion table load and convert.
taxon.conv <- read.csv('data/input/relation_tables/fullpaleon_conversion_v0.3-2.csv', header=TRUE, stringsAsFactors = FALSE)

taxon.conv[taxon.conv[,6] == '', 6] <- NA
taxon.conv[taxon.conv[,6] %in% '#N/A', 6] <- 'Other hardwood'

indiana$l3a.tree1 <-  taxon.conv[match(tolower(indiana$L1.tree1), 
                                       tolower(taxon.conv$Level.1)),6]

indiana$l3a.tree2 <-  taxon.conv[match(tolower(indiana$L1.tree2), 
                                       tolower(taxon.conv$Level.1)),6]

coordinates(indiana) <- ~ POINT_X + POINT_Y

proj4string(indiana) <- CRS('+init=epsg:4326')

base.rast <- raster(xmn = -71000, xmx = 2297000, ncols=296,
                    ymn = 58000,  ymx = 1498000, nrows = 180,
                    crs = '+init=epsg:3175')

indiana.sp <- spTransform(indiana[,c('l3a.tree1', 'l3a.tree2')], CRS('+init=epsg:3175'))

rast.cell <- function(x){
  nc <- indiana.sp
  nc$test <- rowSums(nc@data == x, na.rm=TRUE)
  getValues(rasterize(nc, base.rast, field = 'test', fun = 'sum', na.rm=TRUE))
}

taxa <- as.character(unique(unlist(indiana.sp@data)))

rast.set <- sapply(taxa, rast.cell)
colnames(rast.set) <- taxa

rast.coord <- xyFromCell(base.rast, 1:ncell(base.rast))
rast.out <- cbind(rast.coord, rast.set)

write.csv(rast.out, 'data/output/gridded/IndianaData_v0.2.csv')
