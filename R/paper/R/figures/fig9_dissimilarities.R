histogram <- ggplot(distances, aes(x = dist, linetype = class, size = class)) + 
  geom_density(alpha = 0.5) +
  xlab('Composition Dissimilarity') +
  ylab('Density') +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_sqrt(expand = c(0,0)) +
  scale_linetype_manual(values = c(1, 1, 2, 2)) +
  scale_size_manual(values = c(1, 2, 1, 2)) +
  theme_bw() +
  theme(axis.text = element_text(family='serif', size = 20),
        axis.title = element_text(family='serif', size = 24, face = 'bold'),
        legend.position = 'none')

ggsave(histogram, filename = paste0('figures/fig9_densities_v', version, '.tiff'), width = 9, height = 4, dpi = 300)

plss.dissim <- distances

intern.dissim <- ninefive(plss.dissim$dist[plss.dissim$class %in% c('FIA', 'PLSS')])

plss.ninefive <- distances

#  Anything really big gets assigned a value of 100:
plss.ninefive$dist[plss.ninefive$class %in% c('FIA-PLSS', 'PLSS-FIA') & 
                     plss.ninefive$dist > intern.dissim[2] & 
                     !is.na(plss.ninefive$dist)] <- 100

plss.ninefive$dist[!(plss.ninefive$class %in% c('FIA-PLSS', 'PLSS-FIA') & 
                       plss.ninefive$dist > intern.dissim[2] & 
                       !is.na(plss.ninefive$dist))] <- NA

plss.ninefive <- plss.ninefive[!is.na(plss.ninefive$dist),]

plss.ninefive$dist <- factor(plss.ninefive$dist)

spatial.internal <- base.map + # base.map is defined in 'R/naturalearth.R'
  geom_tile(data = plss.dissim[plss.dissim$class %in% c('FIA', 'PLSS'),], 
            aes(x = x, y = y, fill = dist)) +
  scale_fill_gradient2(low='white', high='darkred', na.value = NA) +
  geom_path(data = umw.domain, aes(x = long, y = lat, group = group, linesize = paleon)) +
  geom_path(data = can.domain, aes(x = long, y = lat, group = group)) +
  facet_wrap(~class, nrow = 1)

spatial.external <- base.map + # base.map is defined in 'R/naturalearth.R'
  geom_tile(data = plss.dissim[plss.dissim$class %in% c('FIA-PLSS', 'PLSS-FIA'),], aes(x = x, y = y, fill = dist)) +
  facet_wrap(~class, nrow = 1) +
  scale_fill_gradient(low='white', high='red', na.value = NA) +
  geom_point(data = plss.ninefive[plss.ninefive$class %in% c('FIA-PLSS', 'PLSS-FIA'),], aes(x = x, y = y), color = 'blue', 
             alpha = 0.5, shape = 15) +
  geom_path(data = umw.domain, aes(x = long, y = lat, group = group, linesize = paleon)) +
  geom_path(data = can.domain, aes(x = long, y = lat, group = group)) +
  theme(legend.position = 'none')
