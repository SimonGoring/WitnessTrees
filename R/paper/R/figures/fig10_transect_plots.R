library(RColorBrewer)

#  Transect axes:
transects <- list(one = data.frame(x = c(-18093, 631094),
                                   y = c(1470422, 668889)),
                  two = data.frame(x = c(-48780, 675121),
                                   y = c(832084, 1151701)))

get_cells <- function(x){
  cells <- extract(num.rast, 
                   data.frame(x = seq(x$x[1], x$x[2], length.out = 500),
                              y = seq(x$y[1], x$y[2], length.out = 500)))
  
  cells <- cells[!duplicated(cells)]
  cells
}

trans.cells <- lapply(transects, get_cells)

get_points <- function(x){
  cells <- extract(num.rast, 
                   data.frame(x = seq(x$x[1], x$x[2], length.out = 500),
                              y = seq(x$y[1], x$y[2], length.out = 500)))
  
  cells <- cells[!duplicated(cells)]
  
  transect.plss <- count.table[match(cells, count.table$cell),3:ncol(count.table)]
  
  #  Remove non-tree and empty plots:
  transect.plss <- transect.plss[!(rowSums(is.na(transect.plss))>0),]
  
  transect.out  <- transect.plss[,!colnames(transect.plss) %in% c('No tree', 'cell')]
  transect.out  <- transect.out / rowSums(transect.out)
  transect.out$cell <- transect.plss$cell
  transect.plss <- transect.out[!is.na(rowSums(transect.out[,!colnames(transect.out) == 'cell'])),]
  
  transect.fia <- agg.dens[match(cells, agg.dens$cell),]
  
  transect.plss$cell <- xyFromCell(num.rast, transect.plss$cell)[,1]
  transect.fia$cell <- xyFromCell(num.rast, transect.fia$cell)[,1]

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

transect.taxa <- lapply(transects, get_points)
transect.taxa[[1]]$transect <- 'One'
transect.taxa[[2]]$transect <- 'Two'

transects <- do.call(rbind.data.frame, transect.taxa)

transects$class <- factor(transects$class, levels = c('PLSS', 'FIA'))
transects$PFT2  <- gsub('Evergreen|Deciduous', '', transects$PFT)

getPalette <- colorRampPalette(brewer.pal(9, "Paired"))

trans.plot <- ggplot(transects, aes(x = cell, 
                               y = value, 
                               linetype = PFT)) +
  geom_line(alpha=0.2, aes(group = variable, color = variable)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_sqrt(expand = c(0,0)) +
  #  stat_smooth(se = FALSE, method = 'gam', 
  #              formula = y ~ te(x), family = binomial, size = 2) + 
  geom_smooth(method = 'gam', family = betar, se = FALSE,
              formula = y ~ s(x, k = 10), size = 1.2, aes(color = variable)) + 
  scale_color_manual(values = getPalette(9)) +
  facet_wrap(~transect+PFT2+class, ncol = 2) +
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
comp.min <- comp.grid[plss.cells %in% fia.rows,]
rownames(comp.min) <- plss.cells[plss.cells %in% fia.rows]

fia.min <- fia.aligned[fia.rows%in% plss.cells,]
rownames(fia.min) <- fia.rows[fia.rows %in% plss.cells]

pls.beta <- betadiver(comp.min, 'sor')
fia.beta <- betadiver(fia.min, 'sor')

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
