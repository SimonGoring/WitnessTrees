#  Testing the angles:
angle.dat <- data.frame(azm = c(angles$az1, angles$az2), 
                        rbind(correction, correction))

angle.dat <- na.omit(angle.dat)

ggplot(angle.dat, aes(x = azm)) + 
  geom_bar(aes(y = (..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..])) +
  facet_wrap(~zeta) +
  theme_bw() +
  ylab('Proportion of Trees at Azimuth') +
  xlab('Azimuth')

angle.dat2 <- data.frame(azm = abs(angles$az1 - angles$az2), 
                         correction)

angle.dat2 <- na.omit(angle.dat2)

ggplot(angle.dat2, aes(x = azm)) + 
  geom_bar(aes(y = (..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..])) +
  facet_wrap(~theta) +
  theme_bw() +
  ylab('Proportion of Trees at Azimuth') +
  xlab('Azimuth')

angle.dat <- data.frame(azm = c(diam), 
                        rbind(correction, correction))

angle.dat <- na.omit(angle.dat)

ggplot(angle.dat, aes(x = azm)) + 
  geom_bar(aes(y = (..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..])) +
  facet_wrap(~zeta) +
  theme_bw() +
  ylab('Proportion of Trees at Azimuth') +
  xlab('Azimuth')
