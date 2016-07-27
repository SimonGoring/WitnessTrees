#  I stripped the ninefive out of here.

val.plot <- data.frame(xyFromCell(dens, 1:ncell(dens)),
                       dens = c(getValues(dens), 
                                getValues(basal), 
                                getValues(biomass)),
                       class = factor(rep(c('Density', 'Basal', 'Biomass'), each = ncell(dens))))

dens.facet <- base.map +
  geom_tile(data = val.plot[val.plot$class == 'Density',], 
            aes(x = x, y = y, fill = dens), alpha = 0.8) +
  scale_fill_gradient(low='white', high='black', na.value = NA, trans = 'sqrt') +
  geom_path(data = umw.domain, aes(x = long, y = lat, group = group, linesize = paleon)) +
  geom_path(data = can.domain, aes(x = long, y = lat, group = group)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = 'top',
        legend.title = element_blank())
  
basa.facet <-  base.map +
  geom_tile(data = val.plot[val.plot$class == 'Basal',], 
            aes(x = x, y = y, fill = dens), 
            alpha = 0.8) +
  scale_fill_gradient(low='white', high='black', trans='sqrt', na.value = NA) +
  geom_path(data = umw.domain, aes(x = long, y = lat, group = group, linesize = paleon)) +
  geom_path(data = can.domain, aes(x = long, y = lat, group = group)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = 'top',
        legend.title = element_blank())

biom.facet <-   base.map +
  geom_tile(data = val.plot[val.plot$class == 'Biomass',], 
            aes(x = x, y = y, fill = dens), 
            alpha = 0.8) +
  scale_fill_gradient(low='white', high='black', trans='sqrt', na.value = NA) +
  geom_path(data = umw.domain, aes(x = long, y = lat, group = group, linesize = paleon)) +
  geom_path(data = can.domain, aes(x = long, y = lat, group = group)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = 'top',
        legend.title = element_blank())

anders.data <- val.plot[val.plot$class == 'Density',]
anders.data$dclass <- NA
anders.data$dclass[anders.data$dens < 47] <- 'Savanna'
anders.data$dclass[anders.data$dens < 0.5] <- 'Prairie'
anders.data$dclass[anders.data$dens >47] <- 'Forest'


anders.facet <- base.map +
  geom_tile(data = anders.data, 
            aes(x = x, y = y, fill = dclass), 
            alpha = 0.8) +
  scale_fill_brewer(type='qual', na.value = NA) +
  geom_path(data = umw.domain, aes(x = long, y = lat, group = group, linesize = paleon)) +
  geom_path(data = can.domain, aes(x = long, y = lat, group = group)) +
  geom_path(data=river.subset, aes(x = long, y = lat, group = group), color = 'blue', alpha=0.6) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = 'top',
        legend.title = element_blank())

writeRaster(dens, paste0('fig_rasters/dens_v', version, '.tif'), overwrite=TRUE)
writeRaster(basal, paste0('fig_rasters/basal_v', version, '.tif'), overwrite=TRUE)
writeRaster(biomass, paste0('fig_rasters/biomass_v', version, '.tif'), overwrite=TRUE)
