###  
# Generate composition stack and save as netcdf and tiff:
compstack.df <- composition.table
coordinates(compstack.df) <- ~ x + y
proj4string(compstack.df) <- CRS(proj4string(dens))

composition.table[,3:ncol(composition.table)] <- composition.table[,3:ncol(composition.table)] /
  rowSums(composition.table[,3:ncol(composition.table)], na.rm=TRUE)

for(i in 1:(ncol(composition.table) - 2)){
  if(i == 1){
    comp.stack <- setValues(base.rast, composition.table[,i+2])
  }
  else{
    comp.stack <- stack(comp.stack, setValues(base.rast, composition.table[,i+2]))
  }
}

names(comp.stack) <- colnames(composition.table)[3:ncol(composition.table)]

writeRaster(comp.stack, '../fullcomposition_ll_01.tif')

dim_lon = ncdim_def('longitude', 'degrees_east', seq(-98.6, -66.1, by = 5/60))
dim_lat = ncdim_def('latitude', 'degrees_north', seq(36.5,  49.75, by = 5/60))

var.list <- list()

for(i in 1:nlayers(comp.stack)){
  var.list[[i]] <- ncvar_def(names(comp.stack)[i], '% composition', list(dim_lon, dim_lat), -999)
}

nc <- nc_create('../fullcomposition_ll_01.nc', var.list)

for.nc <- function(x){
  x.out <- t(as.matrix(x))
  x.out[, ncol(x.out):1]
}

for(i in 1:nlayers(comp.stack)){
  ncvar_put(nc, var.list[[i]], for.nc(comp.stack[[i]]))
}

nc_close(nc)
