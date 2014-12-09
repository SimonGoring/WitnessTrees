##################################################################
#  morisita is a function found in misc.functions.r
estimates <- morisita(used.data, correction, veil = FALSE)

stem.density <- estimates[[1]]
basal.area <- estimates[[2]]

stem.density <- SpatialPointsDataFrame(coordinates(used.data), 
                                       data=data.frame(density = estimates[[1]],
                                                       basal   = estimates[[2]],
                                                       diams = rowMeans(diams[,1:2], na.rm=TRUE) * 2.54))

proj4string(stem.density) <- proj4string(used.data)

#  Based on our knowledge of the dataset, assuming NA values are density 0 is
#  justified in some regions.  Effectively this doesn't do anything, because we've already
#  checked it, but this is just a double check.
zero.trees <- is.na(stem.density$density) & (species[,2] %in% c('No tree', 'Water') | species[,1] %in% c('No tree', 'Water'))

stem.density$density[zero.trees] <- 0
stem.density$basal[zero.trees] <- 0

#  We use two different projection systems here.  This is the test to create the
#  base resolution.
if(model.proj == '+init=epsg:4326'){
  #  lat/long
  base.rast <- raster(xmn = -98.6, xmx = -66.1, ncols=391,
                      ymn = 36.5,  ymx = 49.75, nrows = 160,
                      crs = '+init=epsg:4326')
}

if(model.proj == '+init=epsg:3175'){
  base.rast <- raster(xmn = -71000, xmx = 2297000, ncols=296,
                          ymn = 58000,  ymx = 1498000, nrows = 180,
                          crs = '+init=epsg:3175')
}

#  This ought to be a better way of aggregating the cell densities, diameters, basal areas &cetera
#  The idea is to use cast/melt to put everything together.
#  This is a large & time consuming step.  Again, we place it within an if statement to help save time:

if('pointwise.ests_v0_91.Rdata' %in% list.files('data/output')){
  load('data/output/pointwise.ests_v0_91.Rdata')
}

if(!'pointwise.ests_v0_91.Rdata' %in% list.files('data/output')){
  numbered.rast <- setValues(base.rast, 1:ncell(base.rast))
  numbered.cell <- extract(numbered.rast, spTransform(stem.density,CRSobj=CRS(model.proj)))
  
  spec.table <- data.frame(xyFromCell(base.rast, numbered.cell),
                           cell = numbered.cell,
                           spec = c(as.character(species[,1]), as.character(species[,2])),
                           count = 1,
                           point = 1:nrow(species),
                           density = rep(stem.density$density/2, 2),
                           basal =  rep(stem.density$basal/2, 2),
                           diams = c(diams[,1], diams[,2]),
                           stringsAsFactors = FALSE)
  
  biom.table <- read.csv('data/input/plss.pft.conversion_v0.1-1.csv', 
                         stringsAsFactors=FALSE)
  
  form <- function(x){
    
    eqn <- match(x$spec, biom.table[,1])
    eqn[is.na(eqn)] <- 1  #  Sets it up for non-tree.
    
    b0 <- biom.table[eqn,2]
    b1 <- biom.table[eqn,3]
    
    biomass <- exp(b0 + b1 * log(x$diams * 2.54))
    biomass
    
  }
  
  #  This is the biomass of individual trees.  It needs to be converted into
  #  a stand level value, through the stem density estimate?  The values are
  #  in kg.
  
  biomass <- rep(NA, nrow(spec.table))
  
  for(i in 1:nrow(spec.table)){
    biomass[i] <- form(spec.table[i,])
    cat(i,'\n')
  }
  
  # convert to Mg.
  spec.table$biom <- biomass * spec.table$density / 1000
  spec.table <- spec.table[!is.na(spec.table$density), ]
  
  spec.table$spec[spec.table$spec == 'No Tree'] <- 'No tree'
  
  save(spec.table, file='data/output/pointwise.ests_v0_91.Rdata')
}

nine.nine.pct <- apply(spec.table[,6:9], 2, quantile, probs = 0.99, na.rm=TRUE)

spec.table$density[spec.table$density > nine.nine.pct['density']] <- nine.nine.pct['density']
spec.table$basal[spec.table$basal > nine.nine.pct['basal']] <- nine.nine.pct['basal']
#spec.table$biom[spec.table$biom > nine.nine.pct['biom']] <- nine.nine.pct['biom']

# These are not the full tables since the include only the cells with points in the database.
count.table <- dcast(spec.table, x + y + cell ~ spec, sum, na.rm=TRUE, value.var = 'count')

#61 3 62 49 43 53

unique.len <- function(x){length(unique(x))}

biomass.trees  <- dcast(spec.table, x + y + cell ~ spec, sum, na.rm=TRUE, value.var = 'count')
biomass.points <- dcast(spec.table, x + y + cell ~ spec, unique.len, value.var = 'point')

points.by.cell <- rowSums(count.table[,4:ncol(count.table)], na.rm=TRUE)
trees.by.cell  <- rowSums(count.table[,!colnames(count.table) %in% c('x', 'y', 'cell', 'No tree', 'Water')], na.rm=TRUE)

#  The problem with dcast is that it sums or averages by the count of each species, so the
#  mean per cell is not weighted properly.  We have to use the 'normalize' function defined
#  below to make sure that the densities &cetera sum properly.

density.table <- dcast(spec.table, x + y  + cell ~ spec, sum, na.rm=TRUE, value.var = 'density')
basal.table <- dcast(spec.table, x + y  + cell~ spec, sum, na.rm=TRUE, value.var = 'basal')
biomass.table <- dcast(spec.table, x + y  + cell~ spec, sum, na.rm=TRUE, value.var = 'biom')
diam.table <-  dcast(spec.table, x + y  + cell~ spec, sum, na.rm=TRUE, value.var = 'diams')

#  This function is to 
normalize <- function(x, mult = 2, value = points.by.cell) {x[,4:ncol(x)] <-  x[,4:ncol(x)] / value * mult; x}

density.table <- normalize(density.table)
basal.table <- normalize(basal.table)
biomass.table <- normalize(biomass.table)
diam.table <- normalize(diam.table, mult = 2.54, trees.by.cell)

#  We want rasterized versions of these tables with sums:
rast.fun <- function(x){
  
  to_grid <- data.frame(cell = x$cell, 
                        total = rowSums(x[,4:ncol(x)], na.rm=TRUE))
  
  empty <- rep(NA, ncell(base.rast))
  empty[to_grid$cell] <- to_grid$total
  setValues(base.rast, empty)

}

dens     <- rast.fun(density.table)
count.up <- rast.fun(count.table)
basal    <- rast.fun(basal.table)
biomass  <- rast.fun(biomass.table)
mdiam    <- rast.fun(diam.table); mdiam[mdiam==0] <- NA

#  To get standard deviations within cells we need to do things a bit differently:
rowset <- cbind(1:(nrow(spec.table)/2), (nrow(spec.table)/2+1):nrow(spec.table))

sd.table <- data.frame(cell = spec.table$cell[rowset[,1]],
                       density = spec.table$density[rowset[,1]] * 2,
                       basal   = spec.table$basal[rowset[,1]] + spec.table$basal[rowset[,2]],
                       biomass = spec.table$biom[rowset[,1]] + spec.table$biom[rowset[,2]])

dens.sd <- dcast(sd.table, cell ~ ., fun.aggregate=sd, na.rm=TRUE, value.var = 'density')
basal.sd <- dcast(sd.table, cell ~ ., sd, na.rm=TRUE, value.var = 'basal')
biomass.sd <- dcast(sd.table, cell ~ ., sd, na.rm=TRUE, value.var = 'biomass')

data.table <- data.frame(xyFromCell(dens, 1:ncell(dens)), 
                         stem.dens = getValues(dens),
                         basal.area = getValues(basal),
                         biomass = getValues(biomass),
                         diam = getValues(mdiam))

data.table <- data.table[!is.na(data.table[,3]), ]

if(model.proj == '+init=epsg:3175'){
  write.csv(data.table, 'data/output/density.basal.biomass_v1_9alb.csv')
}
if(model.proj == '+init=epsg:4326'){
  write.csv(data.table, 'data/output/density.basal.biomass_v1_9ll.csv')
}

#  Now we need to add zero cells to the dataframe:
comp.table <- data.frame(xyFromCell(base.rast, 1:ncell(base.rast)), 
                         cell = 1:ncell(base.rast),
                         matrix(ncol = ncol(density.table)-3, nrow = ncell(base.rast)))

colnames(comp.table)[4:ncol(comp.table)] <- colnames(count.table)[4:ncol(comp.table)]

reform <- function(x){
  comp.table[x$cell,] <- x
  comp.table
}

composition.table <- reform(basal.table)
composition.table[,4:ncol(composition.table)] <- composition.table[,4:ncol(composition.table)]/rowSums(composition.table[,4:ncol(composition.table)], na.rm=TRUE)

if(model.proj == '+init=epsg:4326'){
  write.csv(count.table, 'data/output/glo.forest.composition_v1_9ll.csv')
}

if(model.proj == '+init=epsg:3175'){
  write.csv(count.table, 'data/output/glo.forest.composition_v1_9alb.csv')  
}

pft.trans <- read.csv('data/input/pft.trans.table.csv', stringsAsFactor = FALSE)

to.pft <- function(x){
  pft.melt <- melt(x, id = c('x', 'y', 'cell'))
  pft.melt$pft <- pft.trans[match(gsub(' ', '.', pft.melt$variable), gsub(' ', '.', pft.trans$Taxon)),2]
  pft.pft <- dcast(pft.melt, x + y + cell ~ pft, fun.aggregate = sum, na.rm=TRUE)
  pft.pft <- pft.pft[,!colnames(pft.pft) %in% 'NA']

  colnames(pft.pft) <- c('x', 'y', 'cell', 'Grass', 'TBDT', 'TNDT', 'TNET')
  pft.pft
}

spec.table$pft <- pft.trans$PFT[match(spec.table$spec, pft.trans$Taxon)]

dens.pft <- to.pft(density.table)
basal.pft <- to.pft(basal.table)
biomass.pft <- to.pft(biomass.table)

unique.len <- function(x){length(unique(x))}

biomass.trees.pft  <- dcast(spec.table, x + y + cell ~ pft, sum, na.rm=TRUE, value.var = 'count')
biomass.points.pft <- dcast(spec.table, x + y + cell ~ pft, unique.len, value.var = 'point')

rm(spec.table, stem.density, diams, sd.table)
