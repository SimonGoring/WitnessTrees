good.taxa <- c('Tamarack', 'Cedar/juniper', 'Fir', 'Hemlock', 'Pine', 
               'Spruce', 'Ash', 'Basswood', 'Beech', 'Birch', 'Elm', 
               'Ironwood', 'Maple', 'Oak', 'Poplar')

aa <- cooccur(mat = comp.grid[ ,good.taxa], type = 'site_spp', spp_names = TRUE, thresh = FALSE)
bb <- cooccur(mat = fia.aligned[ ,good.taxa], type = 'site_spp', spp_names = TRUE, thresh = FALSE)

aa$results$era <- 'plss'
bb$results$era <- 'fia'

cc <- do.call(rbind.data.frame, list(aa$results, bb$results))

pls.mat <- matrix(nrow = length(good.taxa), ncol = length(good.taxa))

