#  Read in each of the tables, read in the transformation table, 
#  load the PalEON composition table and then sum everything together.

source('R/composition_raster/rasterizing_MI.R')
source('R/composition_raster/rasterizing_IN.R')
source('R/composition_raster/rasterizing_il.R')

strip.index <- function(x){
  # Sometimes the write.csv gets used without removing the rownames.
  #  This is a short function to make sure we're not accidentally
  #  including them.
  if(all(diff(x[,1])==1)){
    x <- x[,-1]
  }
  x
}

#  Load the data:
indiana <- strip.index(read.csv('data/output/gridded//IllinoisData_v0.2.csv'))
michigan <- strip.index(read.csv('data/output/gridded/so_Michigan_v0.3.csv'))
illinois <- strip.index(read.csv('data/output/gridded/IndianaData_v0.2.csv'))
paleon <- strip.index(read.csv('data/output/aggregated_midwest/glo.forest.composition_v1_91alb.csv'))

paleon <- paleon[,!colnames(paleon) %in% 'cell']

length(unique(c(nrow(indiana), nrow(michigan), nrow(illinois))))==1
final.taxa <- unique(c(colnames(paleon), colnames(indiana), colnames(michigan), colnames(illinois)))

output <- matrix(nrow = nrow(paleon), ncol = length(final.taxa))
colnames(output) <- final.taxa

library(reshape2)

add.to <- function(x, input){
  #  function to take input data and add it to 
  
  out.melt <- melt(input, id=c('x', 'y'))
  inp.melt <- melt(x, id=c('x', 'y'))
  
  if(all(is.na(input))){
    out.melt <- inp.melt
  } else {
    out.melt <- rbind(out.melt, inp.melt)
  }
  
  out.melt$variable <- tolower(gsub('[ ]|[[:punct:]]', '.', out.melt$variable))
    
  out.cast <- dcast(out.melt, x + y ~ variable, fun.aggregate = sum, value.var = 'value', na.rm=TRUE)
  
  return(out.cast)

}

output <- add.to(indiana, output)
output <- add.to(michigan, output)
output <- add.to(illinois, output)
output <- add.to(paleon, output)

taxon.conv <- read.csv('data/input/relation_tables/fullpaleon_conversion_v0.3_1.csv', header=TRUE, stringsAsFactors = FALSE)
restore.cols <- data.frame(input = tolower(gsub('[ ]|[[:punct:]]', '.', unique(taxon.conv$Level.3a))),
                           output = unique(taxon.conv$Level.3a), stringsAsFactors = FALSE)

#  plot(y~x, data = subset(output, output[,i]>0), 
#       main = colnames(output)[i],
#       xlim=range(output$x), ylim=range(output$y), pch=19, cex = 0.2); i <- i+1

colnames(output)[3:ncol(output)] <- restore.cols[match(colnames(output)[3:ncol(output)],
                                                       restore.cols[,1]),2]

write.csv(output, 'data/output/gridded/western_comp_v0.4-2.csv', row.names=FALSE)
