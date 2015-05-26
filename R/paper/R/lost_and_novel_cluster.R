#  This script attempts to cluster the lost and novel forest data in the Upper Midwest.
#  We use k-medoid clustering to assign classes since it's more stable than k-means.

pls.types <- comp.grid[distances$dist[distances$class == 'PLSS-FIA'] > intern.quant[95],-1]

lost.distances <- distances[distances$dist > intern.quant[95] & distances$class == 'PLSS-FIA',]

pls.points <- distances[distances$dist > intern.quant[95] & distances$class == 'PLSS-FIA', 1:2]

#  Matching the FIA points seems to be a bit more difficult than I anticipated:
novel.distances <- distances[distances$class == 'FIA-PLSS' & distances$dist > intern.quant[95]  & distances$cell %in% count.table$cell, ]

fia.types <- fia.aligned[match(novel.distances$cell, agg.basa$cell),-1]
fia.points <- novel.distances[,1:2]

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
                            pct   = 100 * (table(pls.types.clust$clustering) / length(pls.types.clust$clustering)),
                            tot   = 100 * (table(pls.types.clust$clustering) / nrow(comp.grid)),
                            stringsAsFactors = FALSE)

fia.dominants <- data.frame(types = fia.forest.types, 
                            pct   = 100 * (table(fia.types.clust$clustering) / length(fia.types.clust$clustering)),
                            tot   = 100 * (table(fia.types.clust$clustering) / nrow(fia.aligned)),
                            stringsAsFactors = FALSE)


pls.lost <- SpatialPointsDataFrame(coords = pls.points, 
                                   data   = data.frame(type = pls.types.clust$clustering,
                                                       class = pls.dominants[pls.types.clust$clustering,1]),
                                   proj4string = CRS('+init=epsg:3175'))

fia.gain <- SpatialPointsDataFrame(coords      = fia.points, 
                                   data   = data.frame(type = fia.types.clust$clustering,
                                                       class = fia.dominants[fia.types.clust$clustering,1]),
                                   proj4string = CRS('+init=epsg:3175'))

writeOGR(pls.lost, 
         dsn = paste0('../../data/output/wiki_outputs/mapped_lost_v',version,'.shp'),
         layer = paste0('mapped_lost_v',version),
         driver = 'ESRI Shapefile', overwrite=TRUE)

writeOGR(fia.gain, 
         dsn = paste0('../../data/output/wiki_outputs/mapped_gain_v',version,'.shp'),
         layer = paste0('mapped_gain_v',version),
         driver = 'ESRI Shapefile', overwrite=TRUE)

#  The problem with k-means is that because there's a random start the order can change.
#  So we need to do something:

pls.dominants <- pls.dominants[order(pls.dominants[,3], decreasing = TRUE),]
fia.dominants <- fia.dominants[order(fia.dominants[,3], decreasing = TRUE),]
