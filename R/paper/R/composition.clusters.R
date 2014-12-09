library(igraph)

which.forest <- function(){
  #  This classifies forest composition based on the Anderson 1972 classes.
  cutoffs <- data.frame(prairie = 0,
                        savanna = 1,
                        open = 19,
                        forest = 40)
  
  condition <- getValues(dens)
  
  #  Set the forest class type.
  for.class <- findInterval(condition, t(cutoffs))
  
  setValues(base.rast, for.class)

}

#  Get the Anderson forest classes from our analysis.
anderson.class <- which.forest()

#  Rule based classification:
#  Prairie: density class is 1, no other real classification.
#  Savanna: density class is 2, no other classification.
#  Temperate Deciduous:  density class is 3 or 4 and TD > 70%
#  Temperate Evergreen:  density class is 3 or 4 and TE > 70%
#  Mixed Wood:           density class is 3 or 4 and each is < 70%

test.class <- function(pfts){
  #  Assigns forest classes based on the Prentice and Haxeltine model with the same numbering as
  #  in the potential vegetation model.

  x <- pfts[,4:7]
  cell <- pfts$cell
  class <- getValues(anderson.class)[cell]
  
  rf.classes <- rep(NA, ncell(anderson.class))
  
  #  Grassland
  rf.classes[cell][class == 1] <- 10
  
  # Savanna
  rf.classes[cell][class == 2 & apply(x, 1, which.max) %in% c(1, 2)] <- 9
  
  # Temperate Deciduous Forest
  rf.classes[cell][class > 2 & (x/rowSums(x, na.rm=TRUE))$TBDT > .7] <- 5
  
  # Evergreen/Deciduous Mixed Forest/Woodland
  rf.classes[cell][class > 2 & ((x/rowSums(x, na.rm=TRUE))$TNET + (x/rowSums(x, na.rm=TRUE))$TNDT)  > .7] <- 8
  
  # Temperate Needleleaf Evergreen Forest/Woodland
  rf.classes[cell][is.na(rf.classes[cell])] <- 4

  setValues(dens, rf.classes)
}

dens.class <- test.class(dens.pft)
basal.class <- test.class(basal.pft)
biom.class <- test.class(biomass.pft)

patch.size <- function(table, class){
  
  rast <- crop(table, pot.veg)
  
  clumped <- getValues(clump(rast == class, directions=4))
  
  if(model.proj == '+init=epsg:3175'){
    output <- table(clumped) * 64
  }
  
  output
  
}

base <- resample(pot.veg, dens, method = 'ngb')
base[is.na(basal.class)] <- NA

#  This gives me the mean patch size for each zone.
patch.ests <- lapply(list(base, basal.class, dens.class, basal.class, biom.class),
                     function(y) unlist(sapply(c(4, 5, 8, 9, 10), function(x) patch.size(y, x))))
       
#  iGraph has a function edge that works differently than the raster package's edge.
detach('package:igraph')

edge.size <- function(rast){
  if(!model.proj == '+init=epsg:3175'){
    rast <- projectRaster(rast, crs='+init=epsg:3175')
  }
  edge.rast <- boundaries(rast, directions = 4, classes = TRUE) - boundaries(rast, directions = 4)
  
  sum(getValues(edge.rast), na.rm=TRUE)
}

base.edge <- sum(getValues(boundaries(base, directions=4)) == 0, na.rm=TRUE)

edge.est <-  sapply(list(base, basal.class, dens.class, basal.class, biom.class),
                   function(x) edge.size(x))

edge.est <- round(edge.est / base.edge * 100, 1)
