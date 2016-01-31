## ----setup-knitr, echo = FALSE, message = FALSE--------------------------
knitr::opts_chunk$set(
  comment = " ",
  error = FALSE,
  eval = TRUE,
  cache = FALSE,
  tidy = TRUE,
  fig.path = "figure/",
  cache.path = "cache/"
)



library("pander")
library("plyr")
panderOptions('table.style', 'rmarkdown')
panderOptions('table.split.table', 300)
ps.options(fonts='serif')


## ----dataLoad, echo = FALSE, message=FALSE, results = 'hide', warning=FALSE----

version <- '0.9-7'

#  Choose one of the two below (3175 is the Albers, 4326 is lat/long):
model.proj <- '+init=epsg:3175'
#model.proj <- '+init=epsg:4326'

if(model.proj == '+init=epsg:3175'){
  xylimits <- c(-100000, 1050000, 600000, 1600000)
}
if(model.proj == '+init=epsg:4326'){
  xylimits <- c(-98, -83, 42, 50)
}

source('R/load_data_v0.2.R')
source('R/natural_earth.R')


## ----fig2_angles, echo = FALSE, warning=FALSE, message=FALSE-------------

source('R/figures/figure2_angles.R')
source('R/figures/figureX_angles.R')


## ----loadFIAplotsperCell, echo = FALSE, message=FALSE, results = 'hide', warning=FALSE----

fia.plots <- raster('../../data/input/rasters/numplots_alb.tif')


## ----CountHist, echo = FALSE, warning=FALSE, message=FALSE---------------
#  Note that this will not produce the same results as presented in the paper since the dataset attached here represents only a 30% subset of the total dataset used in analysis.
#  This code does not produce any plots.

quadrant <- floor(apply(used.data@data[,16:19], 2, function(x)as.numeric(as.character(x)))/90)

two.quads <- !(quadrant[,1] - quadrant[,2] == 0)

good.trees <- with(used.data@data, (!(diam1 > 254 | diam2 > 254)) | (species1 %in% c('No Tree', 'Water') & species2 %in% c('No Tree', 'Water')))

used.data$good <- two.quads | good.trees

base.rast <- raster(xmn = -71000, xmx = 2297000, ncols=296,
                        ymn = 58000,  ymx = 1498000, nrows = 180,
                        crs = '+init=epsg:3175')

count.plot <- rasterize(spTransform(used.data, CRSobj = CRS(model.proj)), 
                        base.rast, 'good', sum)
  
count.df <- data.frame(xyFromCell(count.plot, 1:ncell(count.plot)),
                       cell = get_cells(xyFromCell(count.plot, 1:ncell(count.plot))),
                       points = getValues(count.plot))

write.csv(na.omit(count.df),
          paste0('../../data/output/wiki_outputs/plss_counts_alb_',version,'.csv'))

source('R/base_calculationsv_1.4.R')

rm(quadrant)


## ----loadtaxon, echo = FALSE, warning=FALSE, message=FALSE---------------
taxon_table <- read.csv('../../data//output//tests//clean.bind.test.csv')


## ----fig3_calculations, echo = FALSE, warning=FALSE, message=FALSE-------

source('R/figures/figure3_panel_means.R')


## ----fig4_calc, echo=FALSE, message=FALSE, warning=FALSE-----------------

resids <- data.frame(xyFromCell(dens, 1:ncell(dens)),
                     dens = getValues(ninefive(dens)),
                     biom = getValues(ninefive(biomass)))

mean.mod <- lm(biom ~ dens - 1, data = resids)

resids$resid <- NA
resids$resid[!is.na(resids$dens)] <- resid(mean.mod, type = 'response')

biom_sd <- base.map +
            geom_tile(data = resids, aes(x = x, y = y, fill = resid)) +
            scale_fill_gradient2(low='blue', mid = 'white', high='red', na.value = NA) +
            geom_path(data = umw.domain, aes(x = long, y = lat, group = group, linesize = paleon)) +
            geom_path(data = can.domain, aes(x = long, y = lat, group = group)) +
            theme_bw() +
            theme(axis.text.x = element_blank(),
                  axis.text.y = element_blank(),
                  legend.position = 'none') +
            xlab('') + ylab('')

rast.output <- writeRaster(setValues(dens, resids$resid), 
                           paste0('fig_rasters/fig4resid_v', version,'.tif'),
                           overwrite = TRUE)



## ----StemDensityCors, echo=FALSE, message=FALSE, warning=FALSE-----------

taxon.cor <- data.frame(taxon = colnames(density.table)[4:ncol(density.table)],
                        cor = NA,
                        pval = NA,
                        df = NA,
                        stat = NA)

for(i in 4:ncol(density.table)){
  is.good <- basal.table[,i] > 0 & rowSums(basal.table[, 4:ncol(density.table)], na.rm=TRUE) > 0
  
  sample.prop <- basal.table[is.good,i] / rowSums(basal.table[is.good, 4:ncol(density.table)], na.rm=TRUE)
  
  test <- try(cor.test(sample.prop, rowSums(basal.table[is.good, 4:ncol(density.table)], na.rm=TRUE)))
  if(!class(test) == 'try-error'){
    taxon.cor$cor[i-3]  <- test$estimate
    taxon.cor$pval[i-3] <- test$p.value
    taxon.cor$stat[i-3] <- test$statistic
    taxon.cor$df[i-3]   <- sum(is.good) - 2
  }
}

taxon.cor <- taxon.cor[!is.na(taxon.cor[,2]),]
rownames(taxon.cor) <- taxon.cor[,1]
taxon.cor[,2] <- round(taxon.cor[,2], 2)

taxon.pct <- data.frame(all = colMeans(composition.table[,4:ncol(composition.table)]*100,
                                       na.rm=TRUE),
                        forested = colMeans(composition.table[getValues(dens)>47,4:ncol(composition.table)]*100,
                                       na.rm=TRUE))

taxon.pct <- round(taxon.pct, 0)


## ----fig5_calc, echo = FALSE, warning=FALSE, message=FALSE---------------

#  Composition.table should actually be the basal table, and we've set that already in base_calculations.
excl <- !colnames(composition.table) %in% c('x', 'y', 'cell', 'Water')

comp.intoplot <- composition.table[,excl]
comp.intoplot[comp.intoplot == 0 ] <- NA

# Top taxa:
top.set <- colSums(comp.intoplot, na.rm=TRUE) > 
  sort(colSums(comp.intoplot, na.rm=TRUE))[ncol(comp.intoplot) - 15]

comp.intoplot <- comp.intoplot[,top.set]
comp.intoplot <- comp.intoplot / rowSums(comp.intoplot, na.rm=TRUE)

comp.intoplot <- data.frame(composition.table[,1:2], comp.intoplot)

compplots <- melt(comp.intoplot, id = c('x', 'y'))

compplots$species <- factor(compplots$variable, 
                            levels = c('Tamarack', 
                                       'Cedar.juniper', 
                                       'Fir', 'Hemlock', 'Pine', 'Spruce',
                                       'Ash', 'Basswood', 'Beech', 'Birch', 'Elm',
                                       'Ironwood', 'Maple', 'Oak', 'Poplar'))

composition.plots <- base.map +
      geom_tile(data = compplots, aes(x = x, y = y, fill=value), alpha = 0.8) +
      scale_fill_gradient(low='white', high='black', trans='sqrt', na.value = NA) +
      geom_path(data = umw.domain, aes(x = long, y = lat, group = group, linesize = paleon)) +
      geom_path(data = can.domain, aes(x = long, y = lat, group = group)) +
      theme_bw() +
      facet_wrap(~species, nrow=3) +
      theme(axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            strip.text = element_text(family='serif', size = 20, face = 'bold'),
            legend.position = 'none')
            
#  This needs to get fixed so that the strip.background = element_rect(fill = c('#ffff99', '#80b1d3', rep('#bebada', 5), rep('#8dd3c7', 8))))
ggsave(composition.plots, 
       filename = paste0('figures/fig5_compplots_v', version, '.tiff'), 
                         dpi = 300, width = 15, height = 9)
rm(compplots)


## ----ForestClusterPlots, echo = FALSE, warning=FALSE, message=FALSE------

source('R/load_fia_data.R')
source('R/figures/Figure6_clusterPlots.R')


## ----tTests, echo = FALSE, warning=FALSE, message=FALSE------------------
comp.pft <- basal.pft[,4:7] / rowSums(basal.pft[,4:7], na.rm=TRUE)

bn.test <- t.test(comp.pft$TBDT[comp.pft$TNDT > 0 & basal.pft$TBDT > 0], 
                  comp.pft$TNDT[comp.pft$TNDT > 0 & basal.pft$TBDT > 0])

comp.means <- round(apply(comp.pft, 2, function(x) mean(x[x>0], na.rm=TRUE)) * 100, 1)
comp.sd <-  round(apply(comp.pft, 2, function(x) sd(x[x>0], na.rm=TRUE)) * 100, 1)


## ----fig6_calc, dpi=150, echo = FALSE, message = FALSE, warning=FALSE, results='hide'----

source('R/figures/figure6_comp_prep.R')


## ----fig7_calc, echo = FALSE, warning=FALSE, message=FALSE, fig.width=7, fig.height = 5, dpi=150----
#  Plotting out the qualitative composition data.
#  Note that I think the Basal Class assignment is wrong for this figure and that MW and TN are reversed.

source('R/composition.clusters.R')

factor.vals <- factor(c(getValues(base), 
                        getValues(dens.class), 
                        getValues(basal.class), 
                        getValues(biom.class)))

levels(factor.vals) <- c('TN', 'TD', 'MW', 'Sa', 'Pr')
factor.vals <- factor(factor.vals, levels = c('Sa', 'Pr', 'TN', 'TD', 'MW'))


plot.frame <- data.frame(xyFromCell(base, 1:ncell(base)),
                         values = factor.vals,
                         model = factor(rep(c('RK1999', 'Density', 'Basal Area', 'Biomass'), each = ncell(base)),
                                        levels = c('Density', 'RK1999', 'Basal Area', 'Biomass')))

plot.frame <- plot.frame[plot.frame$model %in% c('Density'),]

forest.class <- base.map +
        geom_tile(data = plot.frame, aes(x = x, y = y, fill = values)) +
        scale_fill_brewer(type = 'qual', palette='Set3') +
        geom_path(data = umw.domain, aes(x = long, y = lat, group = group, linesize = paleon)) +
        geom_path(data = can.domain, aes(x = long, y = lat, group = group)) +
        theme_bw() + facet_wrap(~model, nrow=1) +
        theme(axis.text.x = element_blank(),
              axis.text.y = element_blank()) +
        xlab('') + ylab('')


## ----fig8_calc, echo = FALSE, message = FALSE, warning=FALSE-------------
source('R/plss_fia_comparison.R')

source('R/figures/figure8_fiaplss_plots.R')

diff.ts <- list(density = t.test(getValues((fiadens)) - getValues((dens))),
                basal   = t.test(getValues((fiabasa)) - getValues((basal))),
                biomass = t.test(getValues((fiabiom))/1000 - getValues((biomass))),
                diam    = t.test(getValues((fiadiam))[is.finite(getValues(mdiam))] - 
                                 getValues((mdiam))[is.finite(getValues(mdiam))]))
tot_biom_diff <- sum(rowSums(cbind(-getValues(fiabiom)/1000, getValues(biomass))), na.rm=TRUE)


## ----diff_by_zone, echo = FALSE, message = FALSE, warning=FALSE----------
#comparing stuff:
test_clusts <- clust_plot

num.rast <- setValues(fiadens, 1:ncell(fiadens))

test_clusts$cell <- extract(num.rast, test_clusts[,1:2])

test_clusts$cells <- (!is.na(getValues((fiadens))) &
                        !is.na(getValues((dens))))[test_clusts$cell]
test_clusts$dens <- (getValues((fiadens)) - getValues((dens)))[test_clusts$cell]
test_clusts$basa <- (getValues((fiabasa)) - getValues((basal)))[test_clusts$cell]
test_clusts$biom <- (getValues((fiabiom))/1000 - getValues((biomass)))[test_clusts$cell]
test_clusts$plss_dens <- getValues((dens))[test_clusts$cell]
test_clusts$plss_basa <- getValues((basal))[test_clusts$cell]
test_clusts$plss_biom <- getValues((biomass))[test_clusts$cell]

zone_melt <- melt(test_clusts, id.vars = c('x', 'y', 'cluster', 'cell'), variable.name = 'estimate')

zone_cast <- dcast(zone_melt, cluster ~ estimate, value.var = 'value', mean, na.rm=TRUE)
zone_cast$cells <- dcast(zone_melt, 
                         cluster ~ estimate, 
                         value.var = 'value', sum, na.rm=TRUE)$cells

zone_cast <- zone_cast[order(zone_cast$cells, decreasing = TRUE),]

zone_cast[,3:8] <- round(zone_cast[,3:8], 1)

zone_cast$dens <- paste0(zone_cast$dens, ' (',zone_cast$plss_dens,')')
zone_cast$basa <- paste0(zone_cast$basa, ' (',zone_cast$plss_basa,')')
zone_cast$biom <- paste0(zone_cast$biom, ' (',zone_cast$plss_biom,')')

zone_cast <- zone_cast[,1:5]

colnames(zone_cast) <- c('Forest Type', 'Number of Cells', 
                         'Stem Density stems ha^-1^',
                         'Basal Area m^2^ ha^-1^', 'Biomass Mg ha^-1^')

zone_cast <- zone_cast[,1:5]

## ----diff_table, echo = FALSE, message = FALSE, warning=FALSE, results='asis'----

pander::pandoc.table(zone_cast, style = "rmarkdown")


## ----diamdiff, echo = FALSE, message = FALSE, warning=FALSE--------------
  diam.diff <- t.test(getValues(mdiam) - getValues(fiadiam))

  diam.diff.fia <- t.test((getValues(mdiam) - getValues(fiadiam))[getValues(mdiam) < getValues(fiadiam)])
                          
  diam.diff.plss<- t.test((getValues(mdiam) - getValues(fiadiam))[getValues(mdiam) > getValues(fiadiam)])

  ranger <- ninefive(getValues(mdiam) - getValues(fiadiam))

  cor.tests <- list(to.PLSS= cor.test((getValues(mdiam) - getValues(fiadiam)),
                                       getValues(mdiam), 
                                       use = 'pairwise.complete.obs'),
                    to.fia = cor.test((getValues(mdiam) - getValues(fiadiam)),
                                       getValues(fiadiam), 
                                       use = 'pairwise.complete.obs'))

  biom.test <- cor.test((getValues(mdiam) - getValues(fiadiam)),
                                       getValues(biomass), 
                                       use = 'pairwise.complete.obs')

## ----fig9_calc, echo = FALSE, warning=FALSE, message=FALSE, results='hide'----

source('R/figures/fig9_dissimilarities.R')

intern.quant <- quantile(distances$dist[distances$class %in% c('PLSS','FIA')], 
                        seq(0, 1, by = 0.01), na.rm=TRUE)

PLS.fia.quant <- quantile(distances$dist[distances$class %in% c('PLSS-FIA')], 
                        seq(0, 1, by = 0.01), na.rm=TRUE)

fia.PLS.quant <- quantile(distances$dist[distances$class %in% c('FIA-PLSS')], 
                        seq(0, 1, by = 0.01), na.rm=TRUE)

diffs <- (1 - (c(findInterval(intern.quant[96], PLS.fia.quant) - 1, findInterval(intern.quant[96], fia.PLS.quant) - 1)/100)) * 100

PLS.na <- cor(comp.grid, distances$dist[distances$class == 'PLSS-FIA'] > intern.quant[95])

#  We need to align the FIA and PLSS data with the distance table.  This seems more problematic than it ought to be.

good.cells <- data.frame(cell = unique(c(fia.agg$cell, 
                                         subset(distances, class == 'FIA-PLSS')$cell)))
good.cells$fia <- match(good.cells$cell, fia.agg$cell)
good.cells$dist <- match(good.cells$cell, subset(distances, class == 'FIA-PLSS')$cell)

good.cells <- na.omit(good.cells)
 
fia.data <- fia.agg[good.cells$fia,-1]
fia.na <-  subset(distances, class == 'FIA-PLSS')$dist[good.cells$dist] > intern.quant[95]

fia.na <- cor(fia.data, fia.na)


## ----quantifyingLost_n_found, echo = FALSE, warning=FALSE, message=FALSE, results='hide'----

source('R/lost_and_novel_cluster.R')


## ----figure9, echo = FALSE, message = FALSE, warning=FALSE---------------

source('R/figures/Fig9_noveltyDistance.R')

fig9_output <- figure_9()


## ----fig11_calc, echo = FALSE, message = FALSE, warning=FALSE------------

source('R/figures/fig10_transect_plots.R')


## ----fig1_output,echo=FALSE,warning=FALSE,message=FALSE,fig.width=4,fig.height = 3,dpi=300, dev=c('png','tiff')----

val.plot <- data.frame(xyFromCell(dens, 1:ncell(dens)),
                       dens = !is.na(getValues(dens)))

val.plot$dens[!val.plot$dens] <- NA

fig1plot <- base.map +
  geom_tile(data = na.omit(val.plot), 
            aes(x = x, y = y), alpha = 0.8, fill = 'lightgray') +
  geom_path(data=river.subset, 
            aes(x = long, y = lat, group = group), color = 'blue', alpha = 0.1) +
    geom_polygon(data=lakes.subset, 
                 aes(x = long, y = lat, group = group), fill = '#ADD8E6') +
  geom_path(data = umw.domain, aes(x = long, y = lat, group = group, linesize = paleon)) +
  geom_path(data = can.domain, aes(x = long, y = lat, group = group)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = 'none',
        legend.title = element_blank())

fig1plot

## ----fig2_output, echo = FALSE, warning=FALSE, message=FALSE, fig.width=5, fig.height = 5, dpi=300, dev=c('png','tiff')----

grid.arrange(dens.facet, anders.facet, basa.facet, biom.facet, nrow = 2)


## ----fig3_output, echo=FALSE, message=FALSE, warning=FALSE, width=5, fig.height = 5, dpi=300, dev=c('png','tiff')----

  biom_sd


## ----fig4_output, echo = FALSE, warning=FALSE, message=FALSE, fig.width=7, fig.height = 5, dpi=300, dev=c('png','tiff')----
composition.plots

## ----fig5_clusters, echo = FALSE, warning=FALSE, message=FALSE, fig.width=7, fig.height = 5, dpi=300, dev=c('png','tiff')----

map_plots


## ----fig6_output, echo = FALSE, message = FALSE, warning=FALSE, results='hide', fig.width=6, fig.height=6, dpi=300, dev=c('png','tiff')----

out.map


## ----fig7_output, echo = FALSE, message = FALSE, warning=FALSE, fig.width=7, fig.height = 7, dpi=300, dev=c('png','tiff')----

# ggsave(arrangeGrob(dens.fp, basa.fp, biom.fp, biomass.comp, nrow = 2), 
#        filename = paste0('figures/fig8_fourpanelbiomass', version, '.tiff'), 
#        width=8, height = 8, dpi=300)

plot(arrangeGrob(dens.fp, basa.fp, biom.fp, biomass.comp, nrow = 2))
       

## ----fig8_output, echo = FALSE, warning=FALSE, message=FALSE, results='hide', fig.width=7, fig.height = 5, dpi=300, dev=c('png','tiff')----

grid.arrange(histogram, spatial.internal, spatial.external, ncol = 1)


## ----fig9_output, echo = FALSE, message = FALSE, warning=FALSE, fig.width=7, fig.height = 6, dpi=300, dev=c('png','tiff')----

plot(fig9_output$plot)


## ----fig10_output, echo = FALSE, message = FALSE, warning=FALSE, fig.width=7, fig.height = 6, dpi=300, dev=c('png','tiff')----
trans.plot

