##################################################################
#  morisita is a function found in misc.functions.r
estimates <- morisita(used.data, correction, veil = FALSE)

stem.density <- estimates[[1]]
basal.area <- estimates[[2]]

stem.density <- SpatialPointsDataFrame(coordinates(used.data), 
                                       data = data.frame(density = estimates[[1]],
                                                         basal   = estimates[[2]],
                                                         diams = rowMeans(diams[,1:2], 
                                                                          na.rm = TRUE) * 2.54))

proj4string(stem.density) <- proj4string(used.data)

#  Based on our knowledge of the dataset, assuming NA values are density 0 is
#  justified in some regions.  Effectively this doesn't do anything, because we've already
#  checked it, but this is just a double check.
zero.trees <- is.na(stem.density$density) & (species[,2] %in% c('No tree', 'Water') | species[,1] %in% c('No tree', 'Water'))

stem.density$density[zero.trees] <- 0
stem.density$basal[zero.trees] <- 0

#  We use two different projection systems here.  This is the test to create the
#  base resolution.
if (model.proj == '+init=epsg:4326') {
  #  lat/long
  base.rast <- raster(xmn = -98.6, xmx = -66.1, ncols = 391,
                      ymn = 36.5,  ymx = 49.75, nrows = 160,
                      crs = '+init=epsg:4326')
  numbered.rast <- setValues(base.rast, 1:ncell(base.rast))
  numbered.cell <- extract(numbered.rast, spTransform(stem.density, CRSobj = CRS(model.proj)))

}

if (model.proj == '+init=epsg:3175') {
  base.rast <- raster(xmn = -71000, xmx = 2297000, ncols = 296,
                          ymn = 58000,  ymx = 1498000, nrows = 180,
                          crs = '+init=epsg:3175')
  numbered.rast <- setValues(base.rast, 1:ncell(base.rast))
  numbered.cell <- extract(numbered.rast, spTransform(stem.density, CRSobj = CRS(model.proj)))
}

#  This ought to be a better way of aggregating the cell densities, diameters, basal areas &cetera
#  The idea is to use cast/melt to put everything together.
#  This is a large & time consuming step.  Again, we place it within an if statement to help save time:
#
#  This is where the biomass code is estimated:

if (paste0('pointwise.ests','_v',version, '.RDS') %in% list.files('../../data/output/aggregated_midwest/')) {

  spec.table <- readRDS(paste0('../../data/output/aggregated_midwest/pointwise.ests','_v',version, '.RDS'))

} else {
  
  spec.table <- data.frame(xyFromCell(base.rast, numbered.cell),
                           cell = numbered.cell,
                           spec = c(as.character(species[,1]), as.character(species[,2])),
                           count = 1,
                           point = 1:nrow(species),
                           density = rep(stem.density$density/2, 2),
                           basal =  rep(stem.density$basal/2, 2),
                           diams = c(diams[,1], diams[,2]),
                           stringsAsFactors = FALSE)
  
  biom.table <- read.csv('../../data/input/relation_tables/plss.pft.conversion_v0.1-1.csv', 
                         stringsAsFactors = FALSE)
  
  form <- function(x) {
    
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
  
  for (i in 1:nrow(spec.table)) {
    biomass[i] <- form(spec.table[i,])
    cat(i,'\n')
  }
  
  # convert to Mg.
  spec.table$biom <- biomass * spec.table$density / 1000
  spec.table      <- spec.table[!is.na(spec.table$density), ]
  
  spec.table$spec[spec.table$spec == 'No Tree'] <- 'No tree'
  
  saveRDS(spec.table, file = paste0('../../data/output/aggregated_midwest/pointwise.ests','_v',version, '.RDS'))
}

# Any plot with only one tree (and one No tree) should be removed.
# Start with plot numbers:
#  Any one tree plot will have the first tree listed, the second not:
one_tree <- c(which(!species[,1] %in% "No tree" & species[,2] %in% "No tree"),
              which(species[,1] %in% "No tree" & !species[,2] %in% "No tree"))

spec.table <- subset(spec.table, !point %in% one_tree)

nine.nine.pct <- apply(spec.table[,6:9], 2, quantile, probs = 0.99, na.rm = TRUE)

spec.table$density[spec.table$density > nine.nine.pct['density']] <- nine.nine.pct['density']
spec.table$basal[spec.table$basal > nine.nine.pct['basal']] <- nine.nine.pct['basal']
#spec.table$biom[spec.table$biom > nine.nine.pct['biom']] <- nine.nine.pct['biom']

#  Add the PFT classes to `spec.table`
pft.trans <- read.csv('../../data/input/relation_tables/pft.trans.table.csv', stringsAsFactor = FALSE)

to.pft <- function(x) {
  pft.melt     <- melt(x, id = c('x', 'y', 'cell'))
  pft.melt$pft <- pft.trans[match(gsub(' ', '.', pft.melt$variable), gsub(' ', '.', pft.trans$Taxon)),2]
  pft.pft      <- dcast(pft.melt, x + y + cell ~ pft, fun.aggregate = sum, na.rm = TRUE)
  pft.pft      <- pft.pft[,!colnames(pft.pft) %in% 'NA']
  
  colnames(pft.pft) <- c('x', 'y', 'cell', 'Grass', 'TBDT', 'TNDT', 'TNET')
  pft.pft
}

spec.table$pft <- pft.trans$PFT[match(spec.table$spec, pft.trans$Taxon)]

# These are not the full tables since they include only the cells with points in the database.
count.table <- dcast(spec.table, x + y + cell ~ spec, sum, na.rm = TRUE, value.var = 'count')

# This is to get count of plots with tree X by cell.
#  we expect to get '1' if the species are X & X or X & Y, but 0 if it's Y & Y (where X
#  is the focal species).

unique.len <- function(x) {length(unique(x))}

# This ultimately gets output to 
biomass.points <- dcast(spec.table, x + y + cell ~ spec, unique.len, value.var = 'point')
biomass.points     <- biomass.points[order(biomass.points$cell), ]

# There are a set of NAs in the PFTs from "Unknown Tree"  
biomass.points.pft <- dcast(spec.table, x + y + cell ~ pft, unique.len, value.var = 'point')
biomass.points.pft <- biomass.points.pft[order(biomass.points.pft$cell), ]

# This is to get the points with trees:
spec.adj <- spec.table
spec.adj$tree <- spec.table$spec %in% "No tree"

treed.points <- dcast(spec.adj, x + y + cell ~ tree, unique.len, value.var = 'point')
treed.points <- treed.points[order(treed.points$cell),]
colnames(treed.points) <- c('x', 'y', 'cell', 'tree', 'no tree')

# Now get the total number of plots per cell:
plots_per_cell <- dcast(spec.table, x + y + cell ~ ., unique.len, value.var = 'point')
plots_per_cell <- plots_per_cell[order(plots_per_cell$cell),]

points_per_cell_df <- data.frame(plots_per_cell,
                                 treed = treed.points$tree,
                                 untreed = treed.points$`no tree`)
colnames(points_per_cell_df)[4] <- "total"

#  The problem with dcast is that it sums or averages by the count of each species, so the
#  mean per cell is not weighted properly.  We have to use the 'normalize' function defined
#  below to make sure that the densities &cetera sum properly.

density.table <- dcast(spec.table, x + y  + cell ~ spec, sum, na.rm = TRUE, value.var = 'density')
basal.table   <- dcast(spec.table, x + y  + cell ~ spec, sum, na.rm = TRUE, value.var = 'basal')
biomass.table <- dcast(spec.table, x + y  + cell ~ spec, sum, na.rm = TRUE, value.var = 'biom')
diam.table    <- dcast(spec.table, x + y  + cell ~ spec, sum, na.rm = TRUE, value.var = 'diams')

## This is to normalize.  Diameters need to be averaged by the tree number, otherwise, by
##  the point number.
points.by.cell <- rowSums(count.table[,4:ncol(count.table)], na.rm = TRUE)
trees.by.cell  <- rowSums(count.table[,!colnames(count.table) %in% c('x', 'y', 'cell', 'No tree', 'Water')], na.rm = TRUE)

#  The function averages the estimates to a point level estimate from the aggregated sum.
#  Why is the multiplier * 2? Because there are two trees per point and we would underestimate otherwise.
normalize <- function(x, mult = 2, value = points.by.cell) {x[,4:ncol(x)] <-  x[,4:ncol(x)] / value * mult; x}

density.table <- normalize(density.table)
basal.table   <- normalize(basal.table)
biomass.table <- normalize(biomass.table)
diam.table    <- normalize(diam.table, mult = 2.54, trees.by.cell)

#  We want rasterized versions of these tables with sums:
rast.fun <- function(x) {
  
  to_grid <- data.frame(cell = x$cell, 
                        total = rowSums(x[,4:ncol(x)], na.rm = TRUE))
  
  empty <- rep(NA, ncell(base.rast))
  empty[to_grid$cell] <- to_grid$total
  setValues(base.rast, empty)

}

dens     <- rast.fun(density.table)
count.up <- rast.fun(count.table)
basal    <- rast.fun(basal.table)
biomass  <- rast.fun(biomass.table)
mdiam    <- rast.fun(diam.table); mdiam[mdiam == 0] <- NA

#  To get standard deviations within cells we need to do things a bit differently:
rowset <- cbind(1:(nrow(spec.table)/2), (nrow(spec.table)/2 + 1):nrow(spec.table))

sd.table <- data.frame(cell    = spec.table$cell[rowset[,1]],
                       density = spec.table$density[rowset[,1]] * 2,
                       basal   = spec.table$basal[rowset[,1]] + spec.table$basal[rowset[,2]],
                       biomass = spec.table$biom[rowset[,1]] + spec.table$biom[rowset[,2]])

dens.sd <- dcast(sd.table, cell ~ ., fun.aggregate = sd, na.rm = TRUE, value.var = 'density')
basal.sd <- dcast(sd.table, cell ~ ., sd, na.rm = TRUE, value.var = 'basal')
biomass.sd <- dcast(sd.table, cell ~ ., sd, na.rm = TRUE, value.var = 'biomass')

data.table <- data.frame(xyFromCell(dens, 1:ncell(dens)), 
                         stem.dens = getValues(dens),
                         basal.area = getValues(basal),
                         biomass = getValues(biomass),
                         diam = getValues(mdiam))

data.table <- data.table[!is.na(data.table[,3]), ]

if (model.proj == '+init=epsg:3175') {
  write.csv(data.table, paste0('../../data/output/aggregated_midwest/density.basal.biomass_alb','_v',version, '.csv'))
}
if (model.proj == '+init=epsg:4326') {
  write.csv(data.table, paste0('../../data/output/aggregated_midwest/density.basal.biomass_ll','_v',version, '.csv'))
}

#  Now we need to add zero cells to the dataframe:
comp.table <- data.frame(xyFromCell(base.rast, 1:ncell(base.rast)), 
                         cell = 1:ncell(base.rast),
                         matrix(ncol = ncol(density.table) - 3, nrow = ncell(base.rast)))

colnames(comp.table)[4:ncol(comp.table)] <- colnames(count.table)[4:ncol(comp.table)]

reform <- function(x) {
  comp.table[x$cell,] <- x
  comp.table
}

#  We use basal area throughout as our measure of composition:
composition.table <- reform(basal.table)
composition.table[,4:ncol(composition.table)] <- composition.table[,4:ncol(composition.table)]/rowSums(composition.table[,4:ncol(composition.table)], na.rm = TRUE)

## Load in the table to match the PFTs to PLSS taxa:

dens.pft <- to.pft(density.table)
basal.pft <- to.pft(basal.table)
biomass.pft <- to.pft(biomass.table)

unique.len <- function(x) {length(unique(x))}

#  Write.outputs:
add.v <- function(x, name) {
  
  if ("cell" %in% colnames(x)) {
    x <- x[order(x$cell), ]
  }
  
  #  Quick file name formatter:
  
  if (model.proj == '+init=epsg:4326') {
    p.ext <- '_ll'
  } else {
    p.ext <- '_alb'
  }
  
  write.csv(x, paste0('../../data/output/wiki_outputs/', name, p.ext, '_v',version, '.csv'),
            row.names = FALSE)
}

#  Write everything out:

## This is for the biomass model:
add.v(points_per_cell_df, 'plss_plots')
add.v(biomass.points, 'plss_plots_taxa')
add.v(biomass.points.pft, 'plss_plots_pft')

## This is everything else:
add.v(diam.table, 'plss_diam')
add.v(count.table, 'plss_composition')
add.v(density.table, 'plss_density')
add.v(basal.table, 'plss_basal')
add.v(biomass.table, 'plss_biomass')

add.v(dens.pft, 'plss_density_pft')
add.v(basal.pft, 'plss_basal_pft')
add.v(biomass.pft, 'plss_biomass_pft')
