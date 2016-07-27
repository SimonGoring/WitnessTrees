min.azm <- function(az1, az2){
  acos(cos(((az1 - az2)/360) * 2*pi)) * 360 / (2*pi)
}

angle.dat <- data.frame(azm.diff = min.azm(angles$az1, angles$az2),
                        correction)

angle.dat <- angle.dat[!(angles$az1 %in% c(45, 135, 225, 315) | angles$az2 %in%  c(45, 135, 225, 315)),]

angle.dat <- na.omit(angle.dat)

angle.dat$theta_name <- paste0('theta = ', angle.dat$theta)

angle.diff.plot <- ggplot(angle.dat, aes(x = azm.diff)) + 
      geom_bar(aes(y = (..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]),
               binwidth = 10) +
      facet_wrap(~theta_name) +
      theme_bw() +
      ylab('Proportion of Trees at Azimuth') +
      xlab('Azimuthal Offset') +    
      scale_x_continuous(expand=c(0,0), breaks=c(0, 90, 180, 270, 360), 
                         labels = c('', '90', '', '270', ''), limits=c(0,180)) +
      scale_y_continuous(expand=c(0,0), breaks = seq(0, .1, by = 0.025), 
                         labels = c(0, '', 0.05,'', .1), limits = c(0, 0.2)) +
      theme(axis.text = element_text(family='serif', size = 12),
            axis.title = element_text(family='serif', size = 18, face = 'bold'),
            strip.text = element_text(family='serif', size = 14, face = 'bold'))
