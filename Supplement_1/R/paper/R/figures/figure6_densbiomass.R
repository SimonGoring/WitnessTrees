
allometry.plot <- data.frame(x.vals = c(getValues(ninefive(dens)),
                                        getValues(ninefive(dens)),
                                        getValues(ninefive(basal))),
                             y.vals = c(getValues(ninefive(basal)),
                                        getValues(ninefive(biomass)),
                                        getValues(ninefive(biomass))),
                             vars = rep(c('Density vs Basal Area', 
                                          'Density vs Biomass', 
                                          'Basal Area vs Biomass'), 
                                        each=ncell(dens)))

allometry.plot <- allometry.plot[rowSums(is.na(allometry.plot)) == 0,]

all.pl <- function(x, label.x, label.y){
  
  ggplot(data=x, aes(x = x.vals, y = y.vals)) +
    geom_point(alpha=0.2) +
    stat_smooth(size = 2, color = 'black') + 
    stat_smooth(size = .5, linetype=2, color = 'gray70') + 
    theme_bw() +
    xlab(label.x) + ylab(label.y) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    theme(axis.text = element_text(family='serif', size = 14),
          axis.title = element_text(family='serif', size = 14, face = 'bold')) 
}

ba.biom <- all.pl(allometry.plot[allometry.plot$vars == 'Basal Area vs Biomass',],
                  expression('Basal area'~m^2~ha^-1), 
                  expression('Biomass'~Mg~ha^-1))
dens.ba <- all.pl(allometry.plot[allometry.plot$vars == "Density vs Basal Area",],
                  expression('Stem Density stems'~ha^-1), 
                  expression('Basal area'~m^2~ha^-1))
dens.biom <- all.pl(allometry.plot[allometry.plot$vars == "Density vs Biomass" ,],
                    expression('Stem Density stems '~ha^-1), 
                    expression('Biomass'~Mg~ha^-1))

# grid.arrange(dens.biom, ba.biom, dens.ba, nrow = 1)