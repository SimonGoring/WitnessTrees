#  This script attempts to cluster the lost and novel forest data in the Upper Midwest.
#  We use k-medoid clustering to assign classes since it's more stable than k-means.

# We're trying to subset the datasets down to just the ones that are either lost or novel.
# Isolate the distances and then also isolate the compositional data.

# PLS types are the distances that are in the PLSS-FIA comparison & more different than the
# 95%ile difference.
pls.types      <- comp.grid[distances$dist[distances$class == 'PLSS-FIA'] > intern.quant[95],-1]
lost.distances <- subset(distances, dist > intern.quant[95] & class == 'PLSS-FIA')
pls.points     <- lost.distances[,1:2]

#  Matching the FIA points seems to be a bit more difficult than I anticipated.
#  It seems like the cell assignments are wrong. . .
novel.distances <- subset(distances, distances$class == 'FIA-PLSS' & 
                            distances$dist > intern.quant[95])

novel.distances$cell <- get_cells(novel.distances[,1:2])

# Get rid of the southern stuff. . .
novel.distances <- subset(novel.distances, cell %in% comp.grid$cell)

fia.types <- fia.aligned[match(novel.distances$cell, fia.aligned$cell),-1]
fia.points <- novel.distances[,1:2]

# There was a cell that was 100% willow.  It was weird, but doesn't seem to be there anymore.
fia.points <- fia.points[fia.types$Willow < 1,]
fia.types <- fia.types[fia.types$Willow < 1,]

#  About 5 clusters for each taxon:
#pls.types.clust <- kmeans((pls.types), centers = 5)
#fia.types.clust <- kmeans((fia.types), centers = 5)

pls.types.clust <- pam((pls.types), k = 5)
fia.types.clust <- pam((fia.types), k = 5)

pls.forest.types <- apply(pls.types.clust$medoid, 1, 
                          function(x){
                            paste(names(sort(x[x>.1],decreasing=TRUE)), collapse='-')
                          })

fia.forest.types <- apply(fia.types.clust$medoid, 1, 
                          function(x){
                            paste(names(sort(x[x>.1],decreasing=TRUE)), collapse='-')
                          })

pls.dominants <- data.frame(types = pls.forest.types, 
                            pct   = as.numeric(100 * (table(pls.types.clust$clustering) / length(pls.types.clust$clustering))),
                            tot   = as.numeric(100 * (table(pls.types.clust$clustering) / nrow(comp.grid))),
                            stringsAsFactors = FALSE)

#pls.dominants <- pls.dominants[order(pls.dominants$pct, decreasing = TRUE),]

fia.dominants <- data.frame(types = fia.forest.types, 
                            pct   = as.numeric(100 * (table(fia.types.clust$clustering) / length(fia.types.clust$clustering))),
                            tot   = as.numeric(100 * (table(fia.types.clust$clustering) / nrow(fia.aligned))),
                            stringsAsFactors = FALSE)

#fia.dominants <- fia.dominants[order(fia.dominants$pct, decreasing = TRUE),]

pls.lost <- SpatialPointsDataFrame(coords = pls.points, 
                                   data   = data.frame(type = pls.types.clust$clustering,
                                                       class = pls.dominants[pls.types.clust$clustering,1]),
                                   proj4string = CRS('+init=epsg:3175'))

fia.gain <- SpatialPointsDataFrame(coords = fia.points, 
                                   data   = data.frame(type = fia.types.clust$clustering,
                                                       class = fia.dominants[fia.types.clust$clustering,1]),
                                   proj4string = CRS('+init=epsg:3175'))

writeOGR(pls.lost, 
         dsn = paste0('../../data/output/wiki_outputs/mapped_lost_v',version,'.shp'),
         layer = paste0('mapped_lost_v',version),
         driver = 'ESRI Shapefile', overwrite=TRUE)

writeRaster(rasterize(pls.lost, numbered.rast, field = 'type'), 
            filename = 'fig_rasters/pls_lost.tif',
            overwrite = TRUE)

pls.lost.rast <- rasterize(pls.lost, numbered.rast, field = 'type')

writeOGR(fia.gain, 
         dsn = paste0('../../data/output/wiki_outputs/mapped_gain_v',version,'.shp'),
         layer = paste0('mapped_gain_v',version),
         driver = 'ESRI Shapefile', overwrite=TRUE)

writeRaster(rasterize(fia.gain, numbered.rast, field = 'type'), 
            filename = 'fig_rasters/fia_gainv0.9.9.tif',
            overwrite = TRUE)

fia.gain.rast <- rasterize(fia.gain, numbered.rast, field = 'type')

small_gainloss <- na.omit(data.frame(gain = getValues(fia.gain.rast), 
                                     loss = getValues(pls.lost.rast),
                                     cell = getValues(numbered.rast)))

mod_loss <- do.call(cbind.data.frame, lapply(1:5, function(x){
                colMeans(fia.aligned[match(small_gainloss$cell[small_gainloss$loss == x],
                                           fia.aligned$cell), -1], na.rm=TRUE)
                }))

mod_gain <- do.call(cbind.data.frame, lapply(1:5, function(x){
  colMeans(fia.aligned[match(small_gainloss$cell[small_gainloss$gain == x],
                             fia.aligned$cell), -1], na.rm=TRUE)
}))

pls_loss <- do.call(cbind.data.frame, lapply(1:5, function(x){
  colMeans(comp.grid[match(small_gainloss$cell[small_gainloss$loss == x],
                           comp.grid$cell), -1], na.rm=TRUE)
}))

pls_gain <- do.call(cbind.data.frame, lapply(1:5, function(x){
  colMeans(comp.grid[match(small_gainloss$cell[small_gainloss$gain == x],
                             comp.grid$cell), -1], na.rm=TRUE)
}))

colnames(pls_gain) <- 1:5
colnames(pls_loss) <- 1:5
colnames(mod_gain) <- 1:5
colnames(mod_loss) <- 1:5