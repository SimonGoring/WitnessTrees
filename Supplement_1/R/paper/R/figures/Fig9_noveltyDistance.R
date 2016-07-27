#  Distance, novelty & dispersal

# Are "novel" cells further from remnant forest than expected?
figure_9 <- function(){

  #  Get the values for the quantile breaks from the PLSS internal dissimilarities:
  breaks  <- quantile(plss.dissim$dist[plss.dissim$class %in% c('PLSS')], c(seq(0, 1, by = 0.01)))
  breaks2 <- quantile(plss.dissim$dist[plss.dissim$class %in% c('PLSS')], c(seq(0, 1, by = 0.1)))

  #  Ecological distances in fia.distances$dist for the FIA to PLSS minimum dissimilarities
  #  Everything above the 95%ile is considered novel.
  fia.distances <- distances[distances$class == 'FIA-PLSS',]
  
  #  Spatial distances are botained from dissimilarity matrix of the first two columns of
  #  fia.distances:
  point.distances <- as.matrix(dist(fia.distances[,1:2]))
  
  #  old_classes is calculated higher in the paper using the same call.
  #  We have to modify the call because we're using a smaller subset of the points.
  #  Just want to make sure we're doing things right.
  
  rem_class_fun <- factor(old_classes$clustering[match(fia.aligned[,1],comp.grid[,1])],
                          labels=c('Tamarack/Pine/Spruce/Poplar',
                                   'Oak/Poplar/Basswood/Maple',
                                   'Pine',
                                   'Hemlock/Cedar/Birch/Maple',
                                   'Oak Savanna'))
  
  #  What is the model that defines the link between distance from a remnant plot
  #  and novelty?
  get_glmmodel <- function(selection, thresh){
    #  This function has two options, one is whether or not we're calculating a
    #  null model.  We use "thresh" to define the threshold, currently at 0.25.
    
    fia.dist <- fia.distances
    
    if(!selection == 'breaks'){
      # This randomizes the distances if we're doing the null model:
      fia.dist$dist <- sample(fia.distances$dist)
    } 

    #  Then find the spatial distance to all points below the "thresh" percentile.
    sp_wmin <- apply(point.distances[fia.dist$dist < breaks[thresh],], 2, which.min)
    sp_dist <- apply(point.distances[fia.dist$dist < breaks[thresh],], 2, min)
  
    fia_by_dist <- data.frame(x = fia.dist$x,
                              y = fia.dist$y,
                              rem_class = rem_class_fun,
                              sp_dist = sp_dist/1000,
                              eco_dist = fia.dist$dist,
                              eco_quant = findInterval(fia.dist$dist,breaks,all.inside=TRUE)/100 > 0.95)
    
    fia_by_dist <- na.omit(fia_by_dist)
    
    eco_mod <- glm(eco_quant ~ sp_dist * rem_class,
                      family=binomial, 
                      data = fia_by_dist)
    
    list(eco_mod,
         fia_by_dist)
  }
  
  log_mod <- get_glmmodel('breaks', thresh=25)
    
  #  Here we want to test whether these coefficients are different than random:
  mods <- lapply(1:100, function(x)get_glmmodel('random'))
  
  get_thresh <- function(model_in, na_cut){
    #
    #  Get the cutoff thresholds for reaching 95%ile dissimilarity:
    #
    model <- model_in[[1]]
    
    data_model <- predict(model, type= 'response', se.fit=TRUE)
                  
    data_in          <- model_in[[2]]
    data_in$log_pred <- data_model[[1]]
    data_in$err      <- data_model[[2]] * 1.96  # convert from SE to CI.
    
    do.call(rbind.data.frame, lapply(levels(rem_class_fun),
                                     function(x){
                                       
                                       #  Approximate the upper and lower bounds of the confidence interval:
                                       low_range <- approx(x = data_in$sp_dist[data_in$rem_class %in% x],
                                                           y = data_in$log_pred[data_in$rem_class %in% x] - 
                                                                               data_in$err[data_in$rem_class %in% x], 
                                                           xout = seq(0,100, by = 1))
                                       
                                       high_range <- approx(x = data_in$sp_dist[data_in$rem_class %in% x],
                                                           y = data_in$log_pred[data_in$rem_class %in% x] + 
                                                             data_in$err[data_in$rem_class %in% x], 
                                                           xout = seq(0,100, by = 1))
                                      
                                       #  Return the lowest and highest distances for the model output:
                                       data.frame(low  = ifelse(all(high_range$y>0.5, na.rm=TRUE),0,
                                                                high_range$x[which.min(abs(na_cut - high_range$y))]),
                                                  high =  ifelse(all(low_range$y>0.5, na.rm=TRUE),0,
                                                                 low_range$x[which.min(abs(na_cut -  low_range$y))]),
                                                  zone = x)
    }))
    
  }
  
  null_thresh  <- do.call(rbind.data.frame,lapply(mods,get_thresh, na_cut=0.5))
  model_thresh <- get_thresh(log_mod, 0.5)
  
  # This is the null model.
  thresh_bounds <- aggregate(x = null_thresh, by = list(null_thresh$zone), FUN = mean)
  
  plot_output <- function(input){
    #  Just to get the ggplot plotted nicely:
    
    log_pred <- predict(input[[1]], type= 'response', se.fit=TRUE)
    
    log_output <- data.frame(dist = na.omit(input[[2]])$sp_dist,
                             rem_class =  na.omit(input[[2]])$rem_class,
                             values = log_pred$fit,
                             val_min = log_pred$fit - log_pred$se.fit*1.96,
                             val_max = log_pred$fit + log_pred$se.fit*1.96)
    
    log_output <- log_output[order(log_output$values),]
    
    ggplot(data = input[[2]], aes(x = sp_dist, y = as.numeric(eco_quant))) +
      geom_jitter(position = position_jitter(width = 1, height = 0.2), alpha=0.5, aes(color = rem_class)) +
      geom_ribbon(data=log_output, aes(x = dist, ymin=val_min, ymax=val_max, color = rem_class),
                  inherit.aes=FALSE, alpha = 0.5) +
      geom_path(data=log_output, aes(x = dist, y = values, 
                                     group = rem_class, color = rem_class),
                inherit.aes=FALSE, size = 1) +
      scale_color_brewer(type='qual', name = 'Forest Type') +
      geom_hline(yintercept = 0.5, linetype = 2, color = 'red') +
      theme_bw() +
      coord_cartesian(xlim=c(-1, 80), ylim=c(-0.1, 1), expand = FALSE) +
      theme(axis.title.y = element_blank(),
            axis.text.y  = element_blank(),
            axis.title.x = element_text(family='serif', size = 12, face='bold'),
            axis.text.x  = element_text(family='serif', size = 10),
            legend.position = "none") +
      xlab('Distance from Remnant (25%ile) Forest')
  }
  
  model_quantile <- ecdf(plss.dissim$dist[plss.dissim$class %in% c('PLSS')])
  
  log_mod[[2]]$eco_quantile <- model_quantile(log_mod[[2]]$eco_dist)
    
  
  dist_plot <- plot_output(log_mod)
  dist_dist <- ggplot(log_mod[[2]]) +
    geom_histogram(aes(x = eco_quantile, y = (..count..)/sum(..count..)), 
             binwidth = 0.05, color = 'red', fill = 'lightgray', alpha = 0.5) +
    coord_flip(xlim=c(-0.0001, 1), ylim = c(0, .4), expand = FALSE) +
    theme_bw() +
    scale_y_reverse() +
    geom_vline(xintercept = 0.95, linetype = 2, color = 'red') +
    geom_vline(xintercept = 0.25, linetype = 2, size = 2) +
    theme(axis.title = element_text(family='serif', size = 12, face='bold'),
          axis.text  = element_text(family='serif', size = 10)) +
    xlab('Dissimilarity Quantile') +
    ylab('Proportion of Cells')
  
  return(list(plot = dist_plot,
         fun = log_mod[[1]],
         quants = log_mod[[2]],
         distances = list(observed = model_thresh, null = thresh_bounds)))
}