#  This is the code to plot the Transect data.
#  Transect axes:
transect.lims <- list(one = data.frame(x = c(-18093, 631094),
                                       y = c(1470422, 668889)),
                      two = data.frame(x = c(-48780, 675121),
                                       y = c(832084, 1151701)))

#  We want to get the x & y coordinates of all points along the transects and then
#  pull the data.

get_points <- function(x){
  
  #  Input comes from the x & y limits of the transect (`transect`)
  
  cells <- extract(numbered.rast, 
                   data.frame(x = seq(x$x[1], x$x[2], length.out = 500),
                              y = seq(x$y[1], x$y[2], length.out = 500)))
  
  cells <- unique(cells)
  
  transect.plss <- count.table[match(cells, count.table$cell),3:ncol(count.table)]
  
  #  Remove non-tree and empty plots:
  transect.plss <- transect.plss[!(rowSums(is.na(transect.plss))>0),]
  
  transect.out  <- transect.plss[,!colnames(transect.plss) %in% c('No tree', 'cell')]
  transect.out  <- transect.out / rowSums(transect.out)
  transect.out$cell <- transect.plss$cell
  transect.plss <- transect.out[!is.na(rowSums(transect.out[,!colnames(transect.out) == 'cell'])),]
  
  transect.fia <- agg.dens[match(cells, agg.dens$cell),]
  transect.fia <- transect.fia[!is.na(transect.fia$cell),]
  transect.fia[is.na(transect.fia)] <- 0
  
  #  We are just using the "x" axis of the data as the parameter on the x axis.
  transect.plss$cell <- xyFromCell(numbered.rast, transect.plss$cell)[,1]
  transect.fia$cell  <- xyFromCell(numbered.rast, transect.fia$cell)[,1]
  
  transect.plss$class <- 'PLSS'
  transect.fia$class <- 'FIA'
  
  transect.fia[,2:(ncol(transect.fia)-1)] <- transect.fia[,2:(ncol(transect.fia)-1)] / rowSums(transect.fia[,2:(ncol(transect.fia)-1)], na.rm=TRUE)
  
  good.taxa <- c('cell', 'class', 'Tamarack', 'Pine', 'Birch', 'Elm', 
                 'Maple', 'Oak', 'Poplar', 'Spruce', 'Hemlock')
  
  sample.out <- rbind(melt(transect.plss[,good.taxa],
                           id.vars = c('cell', 'class')),
                      melt(transect.fia[,good.taxa],
                           id.vars = c('cell', 'class')))
  
  sample.out <- sample.out[!is.na(sample.out$value),]
  
  sample.out$PFT <- biom.table$PFT[match(sample.out$variable, biom.table[,1])]
  sample.out$PFT <- gsub(pattern = 'Temperate', replacement = '', x = sample.out$PFT)
  sample.out$PFT <- gsub(pattern = 'Tree', replacement = '', x = sample.out$PFT)
  
  sample.out          
}

transect.taxa <- lapply(transect.lims, get_points)
transect.taxa[[1]]$transect <- 'One'
transect.taxa[[2]]$transect <- 'Two'

transects <- do.call(rbind.data.frame, transect.taxa)

transects$class <- factor(transects$class, levels = c('PLSS', 'FIA'))
transects$PFT2  <- gsub('Evergreen|Deciduous', '', transects$PFT)
transects$PFT2  <- gsub(' ', '', transects$PFT2, fixed = TRUE)

transects$full <- apply(transects, 1, function(x){
      paste(x[c('transect','PFT2', 'class')], collapse = ', ')})

getPalette <- colorRampPalette(brewer.pal(9, "Paired"))

transects <- na.omit(transects)

trans.plot <- ggplot(transects, aes(x = cell, y = value, linetype = PFT)) +
  geom_line(alpha = 0.4, aes(group = variable, color = variable)) +
  scale_x_continuous(expand = c(0,0), 
                     breaks = c(seq(0, 6e+05, by = 5e+04)),
                     labels = c('', '', '1e+05', '', '', '', '3e+05',
                                '', '', '', '5e+05', '', '')) +
  scale_y_sqrt() +
  coord_cartesian(ylim=c(0,1)) +
  geom_smooth(method = 'gam', method.args=list(family="betar"), se = FALSE,
              formula = y ~ te(x), size = 1.2, aes(color = variable)) + 
  scale_color_manual(values = getPalette(9)) +
  facet_wrap(~full, ncol = 2) +
  theme_bw() +
  xlab('Meters East') + ylab('Percent Composition') +
  theme(axis.text = element_text(family='serif', size = 16),
        axis.title = element_text(family='serif', size = 18, face = 'bold'),
        legend.text = element_text(family='serif', size = 14),
        strip.text =  element_text(family='serif', size = 18, face = 'bold'))

ggsave(trans.plot, filename = paste0('figures/fig10_transectplot_v', version, '.tiff'), dpi = 300, width = 10, height = 8)

colnames(transects)[3] <- 'taxon'
transects$taxon <- as.character(transects$taxon)

gam.test <- function(i){
  tax <- taxon.spec$taxon[i]
  tran <- taxon.spec$transect[i]
  
  input <- na.omit(subset(transects, transect == tran & taxon == tax))
  
  model1 <- gam(value ~ s(cell, by = class),
                data = input,
                family = betar)
  model2 <- gam(value ~ s(cell),
                data = input,
                family = betar)
  
  diff(AIC(model1, model2)[,2])
  
}

taxon.spec <- data.frame(transect = c('One', 'Two'),
                         taxon    = rep(unique(transects$taxon), each = 2),
                         AIC = NA)

for(i in 1:nrow(taxon.spec)){
  output <- gam.test(i)
  taxon.spec$AIC[i] <- output[1]
}

#  Beta diversity along the gradient:
comp.min <- comp.grid[comp.grid$cell %in% agg.dens$cell,]
rownames(comp.min) <- comp.grid$cell[comp.grid$cell %in% agg.dens$cell]

fia.min <- fia.aligned[agg.dens$cell %in% comp.grid$cell,]
rownames(fia.min) <- agg.dens$cell[agg.dens$cell %in% comp.grid$cell]

pls.beta <- betadiver(comp.min, 'sor')
fia.beta <- betadiver(fia.min, 'sor')


get_cells <- function(x){
  cells <- extract(numbered.rast, 
                   data.frame(x = seq(x$x[1], x$x[2], length.out = 500),
                              y = seq(x$y[1], x$y[2], length.out = 500)))
  
  cells <- cells[!duplicated(cells)]
  cells
}

trans.cells <- lapply(transect.lims, get_cells)


get_match <- function(cells, x, table){
  matches   <- match(cells[[x]], rownames(table))
  t(rbind(matches[-length(matches)], matches[-1]))
}


get_trans_t <- function(x){
  pls.matches <- get_match(trans.cells, x, comp.min)
  pls.trans.a <- as.matrix(pls.beta)[t(pls.matches)]
  
  fia.matches <- get_match(trans.cells, x, fia.min)
  fia.trans.a <- as.matrix(fia.beta)[t(fia.matches)]
  
  t.test(fia.trans.a, pls.trans.a, paired=TRUE)
}

trans.t.tests <- lapply(1:2, get_trans_t)
