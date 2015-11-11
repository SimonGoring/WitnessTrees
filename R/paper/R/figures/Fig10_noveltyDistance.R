#  Distance, novelty & dispersal

# Are "novel" cells further from remnant forest than expected?
figure_10a <- function(){
  
  breaks <- quantile(plss.dissim$dist[plss.dissim$class %in% c('PLSS')], c(seq(0, 1, by = 0.01)))
  breaks2 <- quantile(plss.dissim$dist[plss.dissim$class %in% c('PLSS')], c(seq(0, 1, by = 0.1)))
  
  fia.distances <- distances[distances$class == 'FIA-PLSS',]
  
  point.distances <- as.matrix(dist(fia.distances[,1:2]))
  
  #diag(point.distances) <- max(point.distances)
  
  bb <- apply(point.distances[fia.distances$dist < breaks[25],], 2, min)
  
  fia_by_dist <- data.frame(x = fia.distances$x,
                            y = fia.distances$y,
                            sp_dist = bb/1000,
                            eco_dist = fia.distances$dist,
                            eco_quant = findInterval(fia.distances$dist,breaks,all.inside=TRUE)/100)
  
  #fia_by_dist <- fia_by_dist[fia_by_dist$sp_dist>0,]
  
  old_classes <- pam(comp.grid[,-1], k = 4)
  rem_class <- factor(old_classes$clustering[match(fia.aligned[,1],comp.grid[,1])],
                      labels=c('Maple/Cedar/Hemlock/Birch',
                               'Oak/Poplar/Maple',
                               'Pine/Spruce',
                               'Oak Savanna'))
    
  log_mod <- glm(eco_quant ~ sp_dist*rem_class, family=binomial(logit), data = fia_by_dist)
  
  log_pred <- predict(log_mod, newdata=data.frame(sp_dist = rep(seq(0, 100),4), 
                                                  rem_class = factor(rep(levels(rem_class), each=101))), 
                      se.fit=TRUE, type= 'response')
  
  log_output <- data.frame(dist = rep(seq(0, 100),4),
                           rem_class =  factor(rep(levels(rem_class), each=101)),
                           values = log_pred[[1]],
                           val_min = log_pred[[1]] - log_pred[[2]],
                           val_max = log_pred[[1]] + log_pred[[2]])
  
  dist_plot <- ggplot(data = fia_by_dist, aes(x = sp_dist, y = eco_quant)) +
    geom_jitter(position = position_jitter(width = .5)) +
    geom_hline(yintercept = 0.95, linetype = 2, color = 'red') +
    geom_ribbon(data=log_output, aes(x = dist, ymin=val_min, ymax=val_max, color = rem_class), 
                inherit.aes=FALSE, alpha = 0.5) +
    geom_path(data=log_output, aes(x = dist, y=values,group = rem_class, color = rem_class), 
              inherit.aes=FALSE, size = 1, name = 'Forest Type') +
    theme_bw() +
    coord_cartesian(xlim=c(0, 80), ylim=c(0, 1)) +
    theme(axis.title.y = element_blank(),
          axis.text.y  = element_blank(),
          axis.title.x = element_text(family='serif', size = 16, face='bold'),
          axis.text.x  = element_text(family='serif', size = 14),
          legend.position = "none") +
    xlab('Distance from Remnant (25%ile) Forest')
  
  dist_dist <- ggplot(fia_by_dist) +
    geom_bar(aes(x = eco_quant, y = (..count..)/sum(..count..)), 
             breaks = seq(0,1,0.05), color = 'red', fill = 'lightgray', alpha = 0.5) +
    coord_flip(xlim=c(-0.0001, 1)) +
    theme_bw() +
    scale_y_reverse() +
    geom_vline(xintercept = 0.95, linetype = 2, color = 'red') +
    geom_vline(xintercept = 0.25, linetype = 2, size = 2) +
    theme(axis.title = element_text(family='serif', size = 18, face='bold'),
          axis.text  = element_text(family='serif', size = 14)) +
    xlab('Dissimilarity Quantile') +
    ylab('Proportion of Cells')
  
  #Okay, so what species dominate these cells?
  
  return(list(plot = grid.arrange(dist_dist, dist_plot, ncol=2, widths = c(.8, 1.6)),
         fun = log_mod,
         quants = fia_by_dist,
         paths = log_output))
}