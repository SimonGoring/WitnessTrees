#  Distance, novelty & dispersal

# Are "novel" cells further from remnant forest than expected?
figure_10a <- function(){
  
  breaks <- quantile(plss.dissim$dist[plss.dissim$class %in% c('PLSS')], c(seq(0, 1, by = 0.01)))
  breaks2 <- quantile(plss.dissim$dist[plss.dissim$class %in% c('PLSS')], c(seq(0, 1, by = 0.1)))
  
  #  Ecological distances in fia.distances$dist
  fia.distances <- distances[distances$class == 'FIA-PLSS',]
  
  #  Spatial distances in the first two columns
  point.distances <- as.matrix(dist(fia.distances[,1:2]))
  
  old_classes <- pam(comp.grid[,-1], k = 5)
  
  
  rem_class <- factor(old_classes$clustering[match(fia.aligned[,1],comp.grid[,1])],
                      labels=c('Tamarack/Pine/Spruce/Poplar',
                               'Oak/Poplar/Basswood/Maple',
                               'Pine',
                               'Hemlock/Cedar/Birch/Maple',
                               'Oak Savanna'))
  
  get_glmmodel <- function(selection){
    #  What is the closest 'remnant' point?
    
    fia.dist <- fia.distances
    
    if(!selection == 'breaks'){
      # This randomizes the distances if we're doing the null model:
      fia.dist$dist <- sample(fia.distances$dist)
    } 

    #  Then find the 
    sp_wmin <- apply(point.distances[fia.dist$dist < breaks[25],], 2, which.min)
    sp_dist <- apply(point.distances[fia.dist$dist < breaks[25],], 2, min)
    
    fia_by_dist <- data.frame(x = fia.dist$x,
                              y = fia.dist$y,
                              rem_class = rem_class,
                              sp_dist = sp_dist/1000,
                              eco_dist = fia.dist$dist,
                              eco_quant = findInterval(fia.dist$dist,breaks,all.inside=TRUE)/100)
    
    list(glm(eco_quant ~ sp_dist*rem_class, family=quasibinomial, data = na.omit(fia_by_dist)),
         fia_by_dist)
  }
  
  log_mod <- get_glmmodel('breaks')
    
  #  Here we want to test whether these coefficients are different than random:
  mods <- lapply(1:100, function(x)get_glmmodel('random'))
  
  get_thresh <- function(x){
    #
    #  Get the cutoff thresholds for reaching 95%ile dissimilarity:
    #
    model <- x[[1]]
    
    data_in <- data.frame(sp_dist = rep(seq(0, 100),5), 
                          rem_class = rep(levels(rem_class), each=101))
    
    data_model <- predict(model, newdata= data_in, type= 'response', se.fit=TRUE)
                  
    data_in$log_pred <- data_model[[1]]
    data_in$err      <- data_model[[2]] * 1.96
    
    do.call(rbind.data.frame, lapply(levels(rem_class),
                                     function(x){
                                       data.frame(low = which.min(abs(0.95-subset(data_in, rem_class == x)$log_pred - 
                                                                       subset(data_in, rem_class == x)$err)),
                                                  high = which.min(abs(0.95- 
                                                                      subset(data_in, rem_class == x)$log_pred + 
                                                                        subset(data_in, rem_class == x)$err)),
                                                  zone = x)
                                                  
    }))
    
  }
  
  null_thresh  <- do.call(rbind.data.frame,lapply(mods,get_thresh))
  model_thresh <- get_thresh(log_mod)
  
  thresh_bounds <- aggregate(x = null_thresh, by = list(null_thresh$zone), FUN = mean)
  
  signif <- sapply(mods, function(x){
      low  <- coef(log_mod[[1]]) < confint.default(x[[1]])[,1]
      high <- coef(log_mod[[1]]) > confint.default(x[[1]])[,2]
      
      output <- rep(NA, length(low))
      output[low]  <- 'low'
      output[high] <- 'high'
      output[is.na(output)] <- 'ns'
      output
    })
  
  plot_output <- function(input){
    log_pred <- predict(input[[1]], newdata=data.frame(sp_dist = rep(seq(0, 100),5), 
                                                       rem_class = rep(levels(rem_class), each=101)),
                                    type= 'response', se.fit=TRUE)
    
    log_output <- data.frame(dist = rep(seq(0, 100), 5),
                             rem_class =  factor(rep(levels(rem_class), each=101)),
                             values = log_pred,
                             val_min = log_pred[[1]] - log_pred[[2]]*1.96,
                             val_max = log_pred[[1]] + log_pred[[2]]*1.96)
    
    log_output$rem_class <- factor(log_output$rem_class,levels(log_output$rem_class)[c(5,3,4,1,2)])

    ggplot(data = input[[2]], aes(x = sp_dist, y = eco_quant)) +
      geom_jitter(position = position_jitter(width = 1), alpha=0.2) +
      geom_hline(yintercept = 0.95, linetype = 2, color = 'red') +
      geom_ribbon(data=log_output, aes(x = dist, ymin=val_min, ymax=val_max, color = rem_class),
                  inherit.aes=FALSE, alpha = 0.5) +
      geom_path(data=log_output, aes(x = dist, y=values.fit, group = rem_class, color = rem_class),
                inherit.aes=FALSE, size = 1, name = 'Forest Type') +
      scale_color_brewer(type='qual') +
      theme_bw() +
      coord_cartesian(xlim=c(0, 80), ylim=c(0, 1)) +
      theme(axis.title.y = element_blank(),
            axis.text.y  = element_blank(),
            axis.title.x = element_text(family='serif', size = 12, face='bold'),
            axis.text.x  = element_text(family='serif', size = 10)) +
      xlab('Distance from Remnant (25%ile) Forest')
  }
  
  
  dist_dist <- ggplot(fia_by_dist) +
    geom_bar(aes(x = eco_quant, y = (..count..)/sum(..count..)), 
             breaks = seq(0,1,0.05), color = 'red', fill = 'lightgray', alpha = 0.5) +
    coord_flip(xlim=c(-0.0001, 1)) +
    theme_bw() +
    scale_y_reverse() +
    geom_vline(xintercept = 0.95, linetype = 2, color = 'red') +
    geom_vline(xintercept = 0.25, linetype = 2, size = 2) +
    theme(axis.title = element_text(family='serif', size = 12, face='bold'),
          axis.text  = element_text(family='serif', size = 10)) +
    xlab('Dissimilarity Quantile') +
    ylab('Proportion of Cells')
  
  clust_map <- na.omit(data.frame(fia_by_dist[,1:2], cluster = rem_class))

  map_plots <- base.map + 
    geom_tile(data = clust_map, aes(x = x, y = y, fill=cluster)) +
    scale_fill_brewer(type='qual') +
    geom_path(data=river.subset, aes(x = long, y = lat, group = group), color = 'blue', alpha = 0.1) +
    geom_polygon(data=lakes.subset, aes(x = long, y = lat, group = group), fill = '#ADD8E6') +
    geom_path(data = umw.domain, aes(x = long, y = lat, group = group, linesize = paleon)) +
    geom_path(data = can.domain, aes(x = long, y = lat, group = group)) +
    theme_bw() +
    theme(axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        legend.position = "none")
  
  #Okay, so what species dominate these cells?
  
  return(list(plot = grid.arrange(dist_dist, dist_plot, ncol=2, widths = c(.8, 1.6)),
         fun = log_mod,
         quants = fia_by_dist,
         paths = log_output))
}