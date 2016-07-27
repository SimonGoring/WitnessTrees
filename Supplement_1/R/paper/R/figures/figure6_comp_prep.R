#  Setting up the dataframe for the Figure X plot.

pct <- function(x) x / rowSums(x, na.rm=TRUE)

pftplots <- data.frame(x = rep(dens.pft[,1], 3 * 3),
                       y = rep(dens.pft[,2], 3 * 3),
                       pcts = c(unlist(pct(dens.pft[,5:7])),
                                unlist(pct(basal.pft[,5:7])),
                                unlist(pct(biomass.pft[,5:7]))),
                       pfts = rep(rep(c('BDT', 'NDT', 'NET'), 
                                      each = nrow(dens.pft)),3),
                       measure = factor(rep(c('Density', 'Basal Area', 'Biomass'), each = nrow(dens.pft) * 3),
                                        levels=c('Density', 'Basal Area', 'Biomass')))

pftplots$pfts <- factor(pftplots$pfts, levels = c('BDT', 'NET', 'NDT'))
#pftplots$measure <- factor(pftplots$measure, levels = c('Biomass', 'Basal Area', 'Density'))

out.map <- base.map +
  geom_tile(data = pftplots, aes(x = x, y = y, fill=pcts), alpha = 0.8) +
  scale_fill_gradient(low='white', high='black', trans='sqrt', 
                      limits=c(0, 1), na.value = NA) +
  geom_path(data = umw.domain, aes(x = long, y = lat, group = group, linesize = paleon)) +
  geom_path(data = can.domain, aes(x = long, y = lat, group = group)) +
  facet_wrap(~measure + pfts, nrow = 3) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        strip.text = element_text(family='serif', size = 18, face = 'bold'),
        legend.position = 'none')

ggsave(out.map, filename = paste0('figures/fig6_pftComp_v', version, '.tiff'), dpi = 300, width = 7, height = 7)
