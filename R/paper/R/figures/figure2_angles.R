angle.dat <- data.frame(azm = c(angles$az1, angles$az2), 
                        rbind(correction, correction))

angle.dat <- na.omit(angle.dat[angle.dat$azm > 0, ])
angle.dat$zeta <- paste0('zeta = ', angle.dat$zeta)

angle.plot <- ggplot(angle.dat, aes(x = azm)) + 
  geom_bar(aes(y = (..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]),
           binwidth = 10) +
  facet_wrap(~zeta) +
  theme_bw() +
  ylab('Proportion of Trees at Azimuth') +
  xlab('Azimuth') +    
  scale_x_continuous(expand=c(0,0), breaks=c(0, 90, 180, 270, 360), 
                     labels = c('', '90', '', '270', '')) +
  scale_y_continuous(expand=c(0,0), breaks = seq(0, .075, by = 0.025), 
                     labels = c('0', '', '0.05', '')) +
  coord_cartesian(ylim=c(0, 0.075)) +
  theme(axis.text = element_text(family='serif', size = 12),
        axis.title = element_text(family='serif', size = 24, face = 'bold'),
        strip.text = element_text(family='serif', size = 20, face = 'bold'))

ggsave(plot = angle.plot, filename = 'figures/Fig2_zeta.tiff', dpi = 300, width = 8, height = 6)
