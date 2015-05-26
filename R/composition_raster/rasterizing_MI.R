library(rgdal)

mich.shape <- readOGR('data/output/southern_MI/so_michigan.shp', 'so_michigan')

taxon.conv <- read.csv('data/input/relation_tables/fullpaleon_conversion_v0.3-3.csv')

mich.table <- data.frame(rbind(coordinates(mich.shape), coordinates(mich.shape)),
                         taxa = c(as.character(mich.shape@data$species1), 
                                  as.character(mich.shape@data$species2)))
                         
mich.table$trans.var <- taxon.conv$Level.3a[match(tolower(gsub('.', '',mich.table$taxa, fixed = TRUE)), 
                                                  tolower(gsub('.', '',taxon.conv$Level.1, fixed = TRUE)))]

check.table <- table(mich.table$taxa, mich.table$trans.var, useNA='always')

mich.table$trans.var[is.na(mich.table$trans.var)] <- 'No tree'
mich.table$trans.var[mich.table$trans.var==''] <- 'No tree'

coordinates(mich.table) <- ~ coords.x1 + coords.x2
proj4string(mich.table) <- CRS('+init=epsg:3175')

base.rast <- raster(xmn = -71000, xmx = 2297000, ncols=296,
                    ymn = 58000,  ymx = 1498000, nrows = 180,
                    crs = '+init=epsg:3175')

rast.cell <- function(x){
  #  Assign the values for a single taxon to individual raster cells
  nc <- mich.table
  nc$test <- rowSums(mich.table@data == x, na.rm=TRUE)
  getValues(rasterize(nc, base.rast, field = 'test', fun = 'sum', na.rm=TRUE))
}

taxa <- as.character(unique(unlist(mich.table$trans.var)))

rast.set <- sapply(taxa, rast.cell)
colnames(rast.set) <- taxa

rast.coord <- xyFromCell(base.rast, 1:ncell(base.rast))
mich.out <- cbind(rast.coord, rast.set)

write.csv(mich.out, 'data/output/gridded/so_Michigan_v0.3.csv')
