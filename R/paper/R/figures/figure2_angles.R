angle.dat <- data.frame(azm = c(angles$az1, angles$az2), 
                        rbind(correction, correction))

angle.dat <- na.omit(angle.dat[angle.dat$azm > 0, ])

angle.plot <- ggplot(angle.dat, aes(x = azm)) + 
  geom_bar(aes(y = (..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]),
           binwidth = 10) +
  facet_wrap(~zeta) +
  theme_bw() +
  ylab('Proportion of Trees at Azimuth') +
  xlab('Azimuth') +    
  scale_x_continuous(expand=c(0,0), breaks=c(0, 90, 180, 270, 360), 
                     labels = c('', '90', '', '270', '')) +
  scale_y_continuous(expand=c(0,0), breaks = seq(0, .1, by = 0.025), 
                     labels = c(0, '', 0.05,'', .1)) +
  theme(axis.text = element_text(family='serif', size = 12),
        axis.title = element_text(family='serif', size = 18, face = 'bold'),
        strip.text = element_text(family='serif', size = 14, face = 'bold'))
