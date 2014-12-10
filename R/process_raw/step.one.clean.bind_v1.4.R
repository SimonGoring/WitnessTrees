#  Binding and cleaning the Wisconsin and Minnesota Public Lands Surveys, data is
#  sourced from the Mladenoff Lab at the University of Wisconsin.  The lab has
#  granted us permission to use the data, but not permission to distribute the
#  original datasets.  A version of the Minnesota data can be obtained from 
#  http://deli.dnr.state.mn.us/metadata/pveg_btreept3.html
#  The wisconsin data may be obtained by contacting David Mladenoff at:
#  mladenoff@wisc.edu
#
#  This file opens the Wisconsin and Minnesota datasets, renames the columns of
#  the Minnesota shapefile to match those of the Wisconsin dataset, and then
#  binds a number of columns from both datasets together (but not the complete
#  set of columns, since there are some columns unique to each dataset.  The
#  The ultimate dataset has the following columns:
#  Point:  PLS point number
#  Township:  Township line
#  Range:  Range line
#  diam (1 through 4):  Bearing tree diameter
#  dist  (1 through 4):  Distance to bearing tree
#  species  (1 through 4):  Species of the bearing tree based on the lumping 
#                           described in 'witnesstreecodesv4.csv'
#  az (1 through 4): Azimuth to the bearing trees.

#  Data:  The following files are distributed with this code:
#  Maps/glo_corn_ex.shp : This dataset represents a subset of the full Wisconsin data,
#  and uses only 1% of the original dataset.
#  Maps/Minnesota_ex.shp : 5% of the original Minnesota dataset.
#  Maps/glpotveg_5min.asc : The potential vegetation map of Ramankutty and Foley.
#  witnesstreecodesv4.csv : The lumping codes used with the witness tree dataset.

#  Simon Goring, January 28, 2012.
#  Updated {1.3): Nov 18, 2014
#  Updated (1.4): Dec  8, 2014

library(sp)
library(spdep)
library(rgdal)
library(raster)

wisc <- readOGR('data/raw_data/wisc/glo_corn.shp', 'glo_corn')
minn <- readOGR('data/raw_data/minn/Minnesota.shp', 'Minnesota')
mich <- readOGR('data/output/aggregated_midwest/michigan_filled/michigan_filled.shp', 'michigan_filled')

#  The files are in unique projections, this standardizes the projections to
#  a long-lat projection:

wisc <- spTransform(wisc, CRS('+proj=longlat +ellps=WGS84'))
minn <- spTransform(minn, CRS('+proj=longlat +ellps=WGS84'))
mich <- spTransform(mich, CRS('+proj=longlat +ellps=WGS84'))

#  The wisconsin Range is set as a single value, the 'E' and 'W' codes are in
#  RANGDIR.  Looking at the data it also looks like there are a few ranges
#  that are miscoded.
wisc$RANGDIR[wisc$RANGDIR %in% '2'] <- 'W'
wisc$RANGDIR[wisc$RANGDIR %in% '4'] <- 'E'

wisc$RANGDIR[wisc$RANGDIR %in% 'W' & coordinates(wisc)[,1] > 5e+05] <- 'E'

#  Character vectors are read into R as factors, to merge them they need to
#  be converted to character strings first and then bound together.  To ensure
#  township and ranges are unique we add a state abbreviation to the front of
#  the twonship name.
twp <- c(paste('mn', as.character(minn$TWP)), 
         paste('wi', as.character(wisc$TOWNSHIP)), 
         paste('mi', as.character(mich$twp)))
rng <- c(as.character(minn$RNG), 
         paste(as.character(wisc$RANGE),as.character(wisc$RANGDIR), sep=''), 
               as.character(mich$rng))

#  The index of column names in Minnesota becomes the same as the Wisconsin.
names(minn)[c(8, 10:25)] <- names(wisc)[c(5, 13:28)]
names(mich)[c(36, 13:28)] <- names(wisc)[c(5, 13:28)]

#  This step throws up four warnings, the warnings seem to come from 
#  the distance fields, which are stored as a factor for some reason.
minn$DIST1 <- as.numeric(levels(minn$DIST1)[minn$DIST1])
minn$DIST2 <- as.numeric(levels(minn$DIST2)[minn$DIST2])
minn$DIST3 <- as.numeric(levels(minn$DIST3)[minn$DIST3])
minn$DIST4 <- as.numeric(levels(minn$DIST4)[minn$DIST4])

#  We have made a choice to say that all taxa labelled 'Beech' in Minnesota are likely
#  Bluebeech, or, in our dataset, Ironwood.
minn@data[minn@data == 'BE'] <- 'IR'

#  The merged dataset is called nwmw, Minnesota comes first, then Wisconsin.
nwmw <- rbind(minn[,c(8, 10:25)], wisc[,c(5, 13:28)], mich[,c(36, 13:28)])

 #  There are a set of 9999 values for distances which I assume are meant to be NAs.  Also, there are a set of points where
#  the distance to the tree is 1 or 2 feet.  They cause really big density estimates!
nwmw@data [ nwmw@data == '9999'] <- NA
nwmw@data [ nwmw@data == '8888'] <- NA
nwmw@data [ nwmw@data == '_'] <- NA
nwmw@data [ nwmw@data == '99999'] <- NA
nwmw@data [ nwmw@data == '999999'] <- NA
nwmw@data [ nwmw@data == '6666'] <- NA
nwmw@data [ nwmw@data == '999'] <- NA

# There is some cleaning to do.  A bit frustrating.  We can't confirm the diameters of
#  a number of points, although we hope to at some point in the future:
#  No stem density removals, none of the plots look like they have 'weird' points.
#  Basal area removals:
nwmw@data[456454,] <- rep(NA, ncol(nwmw))  #  in Michigan, second tree diameter 120", too big!
nwmw@data[458684,] <- rep(NA, ncol(nwmw))  #  in Michigan, second tree diameter 114", too big!
nwmw@data[459133,] <- rep(NA, ncol(nwmw))  #  in Michigan, second tree diameter 1014", too big!
nwmw@data[512571,] <- rep(NA, ncol(nwmw))  #  in Michigan, second tree diameter 61", but trees are listed as NA
nwmw@data[516526,] <- rep(NA, ncol(nwmw))  #  second tree diameter 112". In WI?
nwmw@data[516526,] <- rep(NA, ncol(nwmw))  #  second tree diameter 112".  In WI?

diams <-  cbind(as.numeric(nwmw$DIAM1), 
                as.numeric(nwmw$DIAM2), 
                as.numeric(nwmw$DIAM3), 
                as.numeric(nwmw$DIAM4))

dists <-  cbind(as.numeric(nwmw$DIST1), 
                as.numeric(nwmw$DIST2), 
                as.numeric(nwmw$DIST3), 
                as.numeric(nwmw$DIST4))

azimuths <- cbind(as.character(nwmw$AZ1), 
                  as.character(nwmw$AZ2),
                  as.character(nwmw$AZ3),
                  as.character(nwmw$AZ4))

#  getAngle converts the four character azimuth (e.g. N43E) to a numeric, 360
#  degree angle.  It also has to deal with a number of special cases.
#  The code for getAngles is a bit scuzzy, but it leaves only 231 azimuths 
#  untranslated, this is a manageable number.
source('R/process_raw/get_angle.R')
azimuths <- apply(azimuths, 2, get_angle)

#####  Cleaning Trees:  
#      Changing tree codes to lumped names:
spec.codes <- read.csv('data/input/relation_tables/fullpaleon_conversion_v0.3_1.csv', stringsAsFactor = FALSE)
spec.codes <- subset(spec.codes, Domain %in% 'Upper Midwest')

lumped <- data.frame(abbr = as.character(spec.codes$Level.1),
                     lump = as.character(spec.codes$Level.3a))

species.old <- data.frame(as.character(nwmw$SP1), 
                          as.character(nwmw$SP2), 
                          as.character(nwmw$SP3), 
                          as.character(nwmw$SP4), stringsAsFactors = FALSE)

species <- t(apply(species.old, 1, 
                   function(x) lumped[match(x, lumped[,1]), 2]))

#  Here there needs to be a check, comparing species.old against species.
test.table <- table(unlist(species.old), unlist(species), useNA='always')
write.csv(test.table, 'data/output/tests/clean.bind.test.csv')

#  There are a set of dead taxa (DA, DB & cetera) that we exclude.  Only AM is
#  unknown at this point.  This excludes 213 trees.
species[species %in% ''] <- 'No tree'

#  Now we assign species that don't fit to the 'No tree' category.
species[is.na(species)] <- 'No tree'

######
#  Some annoying things that need to be done:
#  First, there are some points where the first tree has a distance of zero
#  since it is the plot center.  
#  In these cases, the first azimuth is sometimes listed in a strange way, either
#  it's an NA (obviously) or it's a strange value.  In any case, we need to
#  ensure the azimuth is something recognized, I set it to 0.  It doesn't really
#  matter though.

treed.center <- (dists[,1] == 0 & !is.na(azimuths[,1]) & diams[,1] > 0)
treed.center[is.na(treed.center)] <- FALSE

azimuths[treed.center,1] <- 0

#  Another special case, two trees of distance 1.  What's up with that?!
dists[rowSums(dists == 1, na.rm=T) > 1, ] <- rep(NA, 4)

#  When the object is NA, or the species is not a tree (NonTree or Water), set
#  the distance to NA.
dists[is.na(species) | species == 'No tree' | species == 'Water'] <- NA

#  Now we're looking at 36843 points with no usable data,
#  5022 points with only one tree
#  524 points with only two trees
#  59 points with three trees
#  380645 points with four trees

#  At this point we need to make sure that the species are ordered by distance
#  so that trees one and two are actually the closest two trees.

#  There's an annoying problem that has to do with having a character string in
#  the subsequent sort/order function, in that it converts everything to a
#  character string.  To fix it I change the string 'species' into a numeric
#  where each number is associated with a factor level.

sp.levels <- levels(factor(unlist(species)))

species.num <- t(apply(species, 1, function(x) match(x, sp.levels)))

usable.data <- data.frame(diams,
                          dists,
                          species.num,
                          azimuths,
                          stringsAsFactors = FALSE)

rank.fun <- function(x){
  #  This is the function we use to re-order the points so that the closest is first, based on distances.
  
  test.dists <- as.vector(x[c(5:8)])
  
  ranker <- order(test.dists, na.last=TRUE)
  
  return(x[c(ranker, ranker+4, ranker + 8, ranker + 12)])
}

colnames(usable.data) <- c(paste('diam', 1:4, sep =''),
                         paste('dist', 1:4, sep = ''), 
                         paste('species', 1:4, sep = ''),
                         paste('az', 1:4, sep = ''))

ranked.data <- matrix(nrow = nrow(usable.data),
                      ncol = ncol(usable.data))

for(i in 1:nrow(ranked.data)){
  if( sum(is.na(ranked.data[i,5:8]))<2 ){
    ranked.data[i,] <- unlist(usable.data[i,])
  } else{
    ranked.data[i,] <- unlist(rank.fun(usable.data[i,]))
  }
  if(i%%6500 == 0)cat('.')
}

#ranked.data <- t(apply(usable.data, 1, rank.fun)) # need to drop 'id'

species <- data.frame(species1 = sp.levels[ranked.data[, 9]],
                      species2 = sp.levels[ranked.data[,10]],
                      species3 = sp.levels[ranked.data[,11]],
                      species4 = sp.levels[ranked.data[,12]])

#  We need to bin the year information so that we can use it to calculate
#  appropriate Cottam Correction factors.  The survey instructions for the PLS
#  change at a number of points during the sirveys in Wisconsin, but are
#  considered to be fixed by the time.
#  Some wisconsin samples don't have a year.  Look this up and figure out why.
#  It causes a problem with the Cottam correction factor.

mn_survey <- read.csv('data/raw_data/minn/MN_Surveys.csv')
mn_survey$TOWN <- paste('T', formatC(mn_survey$TOWN, width=3, flag='0'), 'N', sep ='')
mn_survey$RANG <- paste('R', formatC(mn_survey$RANG, width=2, flag='0'), mn_survey$RDIR, sep ='')

get.minn.year <- function(x)which(paste(mn_survey$TOWN, mn_survey$RANG) == paste(x$TWP, x$RNG))

minn.year <- rep(NA, nrow(minn@data))
for(i in 1:nrow(minn@data)){
  if(is.na(minn.year[i])){
    minn.test <- get.minn.year(minn@data[i,])
    if(length(minn.test) == 1)minn.year[i] <- mn_survey$YEAR[minn.test]
  } 
}

wisc.year <- ifelse(wisc@data$YEAR_ > 1851, '1851+',
                    ifelse(wisc@data$YEAR_ > 1846, '1846-1851',
                           ifelse(wisc@data$YEAR_ > 1834, '1834-1846',
                                  ifelse(wisc@data$YEAR_ > 1832, '1832-1834','None'))))

#  Michigan has 47 instances where the year is '2', and 12 where the year is '9999'
#  The 9999s don't clean up because we're not using the nwmw data here:

mich.year <- as.numeric(as.character(mich@data$SURVYR))
mich.year[mich.year == '9999'] <- 2

mich.year <- ifelse(mich.year > 1851, '1851+',
                    ifelse(mich.year > 1846, '1846-1851',
                           ifelse(mich.year > 1834, '1834-1846',
                                  ifelse(mich.year > 1832, '1832-1834','None'))))
mich.year[is.na(mich.year)] <- 'None'

survey.year <- factor(c(minn.year, wisc.year, mich.year))

#  These are the columns for the final dataset.

final.data <- data.frame(nwmw$POINT,
                        twp,
                        rng,
                        ranked.data[,1:8],
                        species,
                        ranked.data[,13:16],
                        survey.year,
                        stringsAsFactors = FALSE)

colnames(final.data) <- c('Point', 'Township', 'Range',
                          paste('diam',    1:4, sep =''),
                          paste('dist',    1:4, sep = ''), 
                          paste('species', 1:4, sep = ''),
                          paste('az',      1:4, sep = ''), 'year')
           
#  Turn it into a SpatialPointsDataFrame:
coordinates(final.data) <- coordinates(nwmw)

#  Write the data out as a shapefile.
writeOGR(final.data, 
         'data/output/aggregated_midwest/minn.wisc.mich.clean_v1_5.shp', 
         'minn.wisc.mich.clean_v1_5', 'ESRI Shapefile',
         overwrite_layer = TRUE, check_exists = TRUE)
