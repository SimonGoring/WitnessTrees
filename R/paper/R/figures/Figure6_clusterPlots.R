# Figure 6 Forest clusters plot:

old_classes <- pam(comp.grid[,-1], k = 5) # Based on a distance of 0.32

rem_class <- factor(old_classes$clustering,
                    labels=c('Tamarack/Pine/Spruce/Poplar',
                             'Oak/Poplar/Basswood/Maple',
                             'Pine',
                             'Hemlock/Cedar/Birch/Maple',
                             'Oak Savanna'))

clust_plot <- data.frame(comp.pts, 
                         cluster = rem_class,
                         clustNum = as.numeric(rem_class))

cluster_out <- clust_plot
coordinates(cluster_out) <- ~ x + y
proj4string(cluster_out) <- CRS('+init=epsg:3175')

clust_rast <- rasterize(cluster_out, base.rast, field = "clustNum")
writeRaster(clust_rast, 'fig_rasters/plss_clusters.tif', overwrite = TRUE)

map_plots <- base.map + 
  geom_tile(data = clust_plot, aes(x = x, y = y, fill=cluster)) +
  scale_fill_brewer(type = 'qual') +
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

ggsave(filename = "figures/fig_5_ClusterPlots.tiff", plot = map_plots, dpi=300, width=4, height=4)
