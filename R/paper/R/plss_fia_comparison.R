#  PLSS and FIA comparisons.
#  Main outputs:
#  plss.dist - vector of minimum squared-chord distances between plss raster cells.
#  fia.dist  - vector of minimum squared-chord distances between fia raster cells.
#  plss.fia.dist - minimum squared-chord distances between plss cells and closest fia cells.

fia.trans <- read.csv('../../data/input/relation_tables/FIA_conversion_v0.2.csv', stringsAsFactors = FALSE)
fia.trans[fia.trans == ''] <- NA

#  Load FIA composition data:
get.fia.stack <- function(x, type){
  stack.files <- list.files(x, full = TRUE)
  stack.fia <- stack(stack.files)
  
  fia.vals <- getValues(stack.fia)
  name.end <- regexpr(paste('_', type, sep=''), colnames(fia.vals), fixed=TRUE)
  colnames(fia.vals) <- gsub("_$","", substr(colnames(fia.vals), 1, name.end-1))
  colnames(fia.vals) <- gsub('_', '.', colnames(fia.vals))
  
  fia.vals <- as.data.frame(fia.vals)
  
  fia.vals$x <- xyFromCell(stack.fia, 1:ncell(stack.fia))[,1]
  fia.vals$y <- xyFromCell(stack.fia, 1:ncell(stack.fia))[,2]
  
  fia.vals
}

fia.dens <- get.fia.stack('../../data/input/rasters/albers/density/', type = 'density')
fia.basa <- get.fia.stack('../../data/input/rasters/albers/basal_area/', type = 'basal')
fia.biom  <- get.fia.stack('../../data/input/rasters/albers/biomass/', type = 'biomass')
fia.diam  <- get.fia.stack('../../data/input/rasters/albers/diameter/', type = 'diameter')

re.agg <- function(x){
  # This takes the fia data (by species) and aggregates it to PalEON taxa.
  orig <- 1:nrow(x)
  
  x.test <- x[,!colnames(x) %in% c('x', 'y')]
  
  x.test$orig <- orig
  
  x.melt <- melt(x.test, na.rm=TRUE, value.name= 'value', id = c('orig'))
  
  x.melt$equiv <- fia.trans$PalEON[match(gsub('.', ' ', x.melt$variable, fixed = TRUE), 
                                         gsub('.', ' ', fia.trans$scientific, fixed = TRUE))]
  
  x.melt <- x.melt[x.melt$value > 0,]
  
  x.cast <- dcast(x.melt, orig ~ equiv, 
                  sum, drop = FALSE, 
                  value.var = 'value')
  
  x.cast$x <- x$x[x.cast$orig]
  x.cast$y <- x$y[x.cast$orig]
  
  x.cast$cell <- extract(setValues(dens, 1:ncell(dens)), x.cast[,c('x', 'y')])
  
  x.cast <- x.cast[!is.na(x.cast$cell),!is.na(colnames(x.cast))]
  x.melt_2 <- melt(x.cast[ ,!colnames(x.cast)%in% c('orig', 'x', 'y')], 
                   na.rm = TRUE, 
                   value.name= 'value', 
                   id = c('cell'))
  
  x.melt_2 <- x.melt_2[x.melt_2$value > 0,]
  
  x.cast <-  dcast(x.melt_2, 
                   cell ~ variable, 
                   mean, 
                   drop = FALSE, 
                   value.var = 'value')
  
  if(any(! unique(fia.trans$PalEON) %in% colnames(x.cast))){
    missing <- unique(fia.trans$PalEON)[!unique(fia.trans$PalEON) %in% colnames(x.cast) &
                                          !is.na(unique(fia.trans$PalEON))]
    new.cols <- length(missing)
    append.it <- matrix(ncol = new.cols, nrow = nrow(x.cast))
    colnames(append.it) <- missing
    x.cast <- cbind(x.cast, append.it)
  }
  
  x.cast[,c('cell', unique(fia.trans$PalEON)[!is.na(unique(fia.trans$PalEON))])]
}

agg.dens <- re.agg(fia.dens)
agg.basa <- re.agg(fia.basa)
agg.biom <- re.agg(fia.biom)

#  Diameters need to be done a little differently since they ought to be weighted by sample:

push.rast <- function(x){
  fia.rast <- setValues(dens, NA)
  rast.out <- rep(NA, ncell(dens))
  
  rast.out[x$cell] <- rowSums(x[,-1], na.rm=TRUE)
  
  rast.out[rast.out == 0] <- NA
  
  setValues(dens, rast.out)
}

fiadens <- push.rast(agg.dens)
fiabasa <- push.rast(agg.basa)
fiabiom <- push.rast(agg.biom)
fiadiam <- sqrt((fiabasa / fiadens)/pi) * 2 * 100

#  Get the taxa that are in common between the two datasets (e.g., exclude NonTrees):
name.set <- colnames(agg.dens)[colnames(agg.dens) %in% colnames(count.table)]

#  comp.grid is the PLSS data limited to cells with data.
comp.grid <- as.data.frame(composition.table[ , name.set])

good.points <- rowSums(is.na(composition.table)) < (ncol(composition.table)-3)

comp.pts <- composition.table[good.points, 1:2]

comp.grid <- comp.grid[good.points, ]

comp.grid[is.na(comp.grid)] <- 0

comp.grid[,2:ncol(comp.grid)] <- comp.grid[,-1] / rowSums(comp.grid[,-1], na.rm=TRUE)

#  This is the FIA data.
fia.agg <- agg.basa[!rowSums(is.na(agg.basa)) == ncol(agg.basa), name.set]
fia.agg[,2:ncol(fia.agg)] <- fia.agg[,-1] / rowSums(fia.agg[,-1], na.rm=TRUE)
fia.agg[is.na(fia.agg)] <- 0

sample.test <- data.frame(stat = rep(NA, 15),
                          est = rep(NA, 15),
                          df = rep(NA, 15),
                          p.val = rep(NA,15))

for(i in 2:ncol(fia.agg)){

  #  need to match up the FIA and PLSS rows:
  plss.vec <- match(comp.grid$cell, fia.agg$cell)
  
  diss.vec <- (fia.agg[na.omit(plss.vec),i] - comp.grid[!is.na(plss.vec),i])
  test.vec <- (fia.agg[na.omit(plss.vec),i] + comp.grid[!is.na(plss.vec),i]) > 0
  
  #  Test whether the difference in proportions between the PLSS and FIA
  #  are net-different than zero.
  
  tester <- try(t.test(diss.vec[test.vec], var.equal=TRUE))
  
  if(!class(tester) == 'try-error'){
    sample.test[i, ] <- c(tester$statistic, tester$estimate, 
                          tester$parameter, tester$p.value*ncol(fia.agg))
  }
}

sample.test <- sample.test[-1,]
rownames(sample.test) <- colnames(comp.grid)[-1]

biom.table <- read.csv('../../data/input/relation_tables/plss.pft.conversion_v0.1-1.csv')

biomass_fia <- exp(biom.table[match(rownames(sample.test), biom.table[,1]),2] + 
                   biom.table[match(rownames(sample.test), biom.table[,1]),3] * 
                   log(15 * 2.54))

#  Some code to compare the FIA and PLSS versus their differences.
#  We need to get a common table with FIA and PLSS data and then assign values
#  for them.

fia.aligned <- fia.agg[!rowSums(fia.agg[,-1]) == 0 & rowSums(is.na(fia.agg)) == 0,]
comp.grid  <- comp.grid[!rowSums(comp.grid[,-1]) == 0 & rowSums(is.na(comp.grid)) == 0,]

#rm(agg.basa, agg.biom)
if(paste0('distances_v', version, '.Rds') %in% list.files('../../data/output/aggregated_midwest/')){
  distances <- readRDS(paste0('../../data/output/aggregated_midwest/distances_v', version, '.Rds'))
} else {
  #  Check to see what the distances are like within each of the FIA and PLSS datasets
  #  and then between datasets.

  big.mat <- rbind(comp.grid[,-1], fia.aligned[,-1])
  
  big.dist <- vegdist(big.mat, 'bray')
  
  big.dist <- as.matrix(big.dist)
  
  plss.set <- big.dist[1:nrow(comp.grid), 1:nrow(comp.grid)]
  
  fia.set  <- big.dist[(nrow(comp.grid) + 1):nrow(big.dist), (nrow(comp.grid) + 1):nrow(big.dist)]
  
  plss.fia <- big.dist[1:nrow(comp.grid), (nrow(comp.grid) + 1):nrow(big.dist)]
  
  fia.plss <- big.dist[(nrow(comp.grid) + 1):nrow(big.dist), 1:nrow(comp.grid)]
  
  rm(big.dist)
  
  check.min <- function(x){
    if(nrow(x) == ncol(x)){
      diag(x) <- max(x)
    }
    apply(x, 1, min)
  }
  
  fia.pts <- xyFromCell(fiabasa, fia.aligned$cell)
  comp.pts <- xyFromCell(fiabasa, comp.grid$cell)
  
  fia.dist <- data.frame(fia.pts, 
                         dist = check.min(fia.set),
                         class = 'FIA', 
                         stringsAsFactors = FALSE)
  
  fia.plss.dist <- data.frame(fia.pts,
                              dist = check.min(fia.plss), 
                              class = 'FIA-PLSS', 
                              stringsAsFactors = FALSE)
  
  plss.dist <- data.frame(comp.pts, 
                          dist = check.min(plss.set), 
                          class = 'PLSS', 
                          stringsAsFactors = FALSE)
  
  plss.fia.dist <- data.frame(comp.pts, 
                              dist = check.min(plss.fia), 
                              class = 'PLSS-FIA', 
                              stringsAsFactors = FALSE)
  
  
  distances <- rbind(fia.dist, 
                     plss.dist, 
                     fia.plss.dist, 
                     plss.fia.dist)  
  
  distances$cell <- get_cells(distances[,1:2])
  
  saveRDS(distances, file = paste0('../../data/output/aggregated_midwest/distances_v', version, '.Rds'))
  rm(plss.dist, fia.dist, plss.fia, fia.plss)
}

#  Output the distances as raster layers:
#

########
  null.rast <- setValues(numbered.rast, NA)
  plss.rast <- null.rast
  dist_subset <- subset(distances, class == "PLSS")
  dist_subset$cell <- extract(numbered.rast, dist_subset[,1:2])
  plss.rast[dist_subset$cell] <- dist_subset$dist
########

for(i in unique(distances$class)){
  null.rast <- setValues(numbered.rast, NA)
  out.rast <- null.rast
  
  dist_subset <- subset(distances, class == i)
  
  dist_subset$cell <- extract(numbered.rast, dist_subset[,1:2])
  
  out.rast[dist_subset$cell] <- dist_subset$dist
  out.rast[is.na(plss.rast)] <- NA
  
  writeRaster(out.rast, 
              filename = paste0('fig_rasters/distRast_',i,'_v', version, '.tif'),
              overwrite = TRUE)
}

####

#  Test against the number of FIA plots:
fia.plot.num <- data.frame(x = subset(distances, class == 'FIA')[,1],
                           y = subset(distances, class == 'FIA')[,2],
                           dists = subset(distances, class == 'FIA')$dist,
                           plots = extract(fia.plots, subset(distances, class == 'FIA')[,1:2]),
                           rich  = rowSums(fia.aligned[,-1]>0))

library(statmod)

fia.highdist <- subset(fia.plot.num, dists > 0)
fia.num.model <- gam(dists ~ s(x, y, plots), data = fia.highdist, family = Gamma)

fia.rich.model <- gam(rich ~ s(x, y, by = plots), data = subset(fia.plot.num, rich > 0), family = poisson)

tester <- rbind(fia.plot.num, fia.plot.num, fia.plot.num, fia.plot.num, fia.plot.num)
tester$plots <- rep(c(1, 3, 5, 8, 10), each = nrow(fia.plot.num))

model.pred <- predict(fia.num.model, newdata = tester, type = 'response', se.fit = TRUE)
model.rich <- predict(fia.rich.model, newdata = tester, type = 'response', se.fit = TRUE)

tester$predicted <- model.pred[[1]]
tester$predrich  <- model.rich[[1]]

ggplot(tester, aes(x = x, y = predicted)) + geom_smooth(aes(color = factor(plots)))
ggplot(tester, aes(x = x, y = predrich)) + geom_smooth(aes(color = factor(plots)))

fia.cor<- try(t(cor(distances$dist[distances$class=='FIA'],fia.aligned)))
plss.cor<- try(t(cor(distances$dist[distances$class=='PLSS'],comp.grid)))

#  Moran Analysis:
library(spdep)

moran_test <- function(class.in){
  set.diss <- subset(distances, class %in% class.in)
  coordinates(set.diss) <- ~ x + y
  W.nb <- knn2nb(knearneigh(set.diss, k=4))
  # We choose a row standardized spatial weight matrix :
  W.listw <- nb2listw(W.nb,style="W")
  moran.test(x=set.diss$dist, listw=W.listw)
}

all.tests <- lapply(unique(distances$class), moran_test)
names(all.tests) <- unique(distances$class)

moran.stats <- sapply(all.tests, function(x)x$estimate[1])
