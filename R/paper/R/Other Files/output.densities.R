data.out <- stack(dens, dens.sd, basal, basal.sd)
names(data.out) <- c('Stem Density st/ha', 'Density sd st/ha', 'Basal area m2/ha', 'Basal sd m2/ha')
writeRaster(data.out, '../dens_ba_ll_01.tif')

##  Set up the netCDF output.
dim_lon = ncdim_def('longitude', 'degrees_east', seq(-98.6, -66.1, by = 5/60))
dim_lat = ncdim_def('latitude', 'degrees_north', seq(36.5,  49.75, by = 5/60))

var_sd = ncvar_def('cell stem density', 'stems/ha', list(dim_lon, dim_lat), -999)
var_sdsd = ncvar_def('cell stem density sd', 'stems/ha', list(dim_lon, dim_lat), -999)
var_ba = ncvar_def('cell basal area', 'm2/ha', list(dim_lon, dim_lat), -999)
var_basd = ncvar_def('cell basal area sd', 'm2/ha', list(dim_lon, dim_lat), -999)

nc <- nc_create('../dens_ba_ll_01.nc', list(var_sd, var_sdsd, var_ba, var_basd))

for.nc <- function(x){
  x.out <- t(as.matrix(x))
  x.out[, ncol(x.out):1]
}

ncvar_put(nc, var_sd, for.nc(dens))
ncvar_put(nc, var_sdsd, for.nc(dens.sd))
ncvar_put(nc, var_ba, for.nc(basal))
ncvar_put(nc, var_basd, for.nc(basal.sd))

nc_close(nc)

