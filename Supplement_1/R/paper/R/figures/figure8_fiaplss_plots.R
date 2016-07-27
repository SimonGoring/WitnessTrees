
plss.diam <- getValues(mdiam)
plss.diam[plss.diam < 20] <- NA

allometry.plot <- data.frame(plss = c(plss.diam,
                                      getValues(dens),
                                      getValues(basal),
                                      getValues(biomass)),
                             fia =  c(getValues(fiadiam),
                                      getValues(fiadens),
                                      getValues(fiabasa),
                                      getValues(fiabiom/1000)),
                             vars = rep(c('Diameter', 'Density', 'Basal Area', 'Biomass'), 
                                        each = ncell(dens)),
                             lowval = rep(getValues(dens), 4))

allometry.plot <- allometry.plot[!(rowSums(is.na(allometry.plot) | allometry.plot == 0)>0), ]

all.pl <- function(x, label.x, label.y){
  
  range.vals <- ninefive(c(x$plss, x$fia))
  
  ggplot(x, aes(x = plss, y = fia)) + 
  geom_point(alpha = 0.3) +
  geom_abline(intercept = 0, slope = 1, size = 1.5, color = 'gray', alpha = 0.5) +
  theme_bw() +
  xlab(label.x) + ylab(label.y) +
  scale_x_continuous(expand=c(0,5)) +
  scale_y_continuous(expand=c(0,5)) +
  theme(axis.text = element_text(family='serif', size = 14),
        axis.title = element_text(family='serif', size = 14, face = 'bold')) +
  coord_equal(xlim=range.vals, ylim=range.vals)
}

dens.fp <- all.pl(allometry.plot[allometry.plot$vars == 'Density',],
                  expression('PLSS Stem Density stems'~ha^-1), 
                  expression('FIA Stem Density stems'~ha^-1))
basa.fp <- all.pl(allometry.plot[allometry.plot$vars == "Basal Area",],
                  expression('PLSS Basal Area'~m^2~ha^-1), 
                  expression('FIA Basal Area'~m^2~ha^-1))
biom.fp <- all.pl(allometry.plot[allometry.plot$vars == "Biomass" ,],
                  expression('PLSS Biomass'~Mg~ha^-1), 
                  expression('FIA Biomass'~Mg~ha^-1))

biomdiff <- data.frame(xyFromCell(biomass, 1:ncell(biomass)),
                       higher = factor(getValues(biomass-fiabiom/1000)>0),
                       value  = getValues(biomass-fiabiom/1000))

writeRaster(setValues(dens, biomdiff$value), 
            filename = paste0('fig_rasters/biomdiff_v', version,'.tiff'),
            overwrite = TRUE)

biomdiff <- biomdiff[!is.na(biomdiff[,3]),]

biomass.comp <- base.map +
                  geom_tile(data = biomdiff, aes(x = x, y = y, fill = higher)) +
                  scale_fill_brewer(type = 'seq', palette='Set3') +
                  geom_path(data = umw.domain, aes(x = long, y = lat, group = group, linesize = paleon)) +
                  geom_path(data = can.domain, aes(x = long, y = lat, group = group)) +
                  theme_bw() +
                  theme(axis.text.x = element_blank(),
                        axis.text.y = element_blank()) +
                  xlab('') + ylab('')

