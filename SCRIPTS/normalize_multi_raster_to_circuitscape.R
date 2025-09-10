################################################################################
######################### Normalizing multiple rasters ########################
################################################################################

##### normalize +1 raster, by LCM

## INSTALL PACKAGES
#install.packages("raster", dep = TRUE)

## LOAD PACKAGES
library(raster)


###WORLDCLIM

## SET YOUR WORKING DIRECTORY (adjust the path to your folder)
#setwd("D:/Unifesp/Mestrado/analises/rasters_ib/earthenv/heterogeneity/cs_gdm_earth_het_andes/")
#setwd("D:/Unifesp/Mestrado/analises/rasters_ib/earthenv/land_cover/cs_gdm_earth_lc_andes")
#setwd("D:/Unifesp/Mestrado/analises/rasters_ib/earthenv/topography/cs_gdm_earth_tp_andes")
#setwd("D:/Unifesp/Mestrado/analises/rasters_ib/envirem/cs_gdm_env_andes")
#setwd("D:/Unifesp/Mestrado/analises/rasters_ib/hwsd/cs_gdm_hwsd_andes")
setwd("D:/Unifesp/Mestrado/analises/rasters_ib/worldclim/cs_gdm_wc_andes")
##SOLARGIS FAZER URGENTE
setwd("D:/Unifesp/Mestrado/analises/rasters_ib/solargis/normalizados")
#setwd(D:\Unifesp\Mestrado\analises\rasters_ib\solargis\DIF)
#setwd(D:\Unifesp\Mestrado\analises\rasters_ib\solargis\DNI)
#setwd(D:\Unifesp\Mestrado\analises\rasters_ib\solargis\GHI)
#setwd(D:\Unifesp\Mestrado\analises\rasters_ib\solargis\GTI)
#setwd(D:\Unifesp\Mestrado\analises\rasters_ib\solargis\OPTA)
#setwd(D:\Unifesp\Mestrado\analises\rasters_ib\solargis\PVOUT)
#setwd(D:\Unifesp\Mestrado\analises\rasters_ib\solargis\TEMP)



####PITCAIRNIA LANUGINOSA (ANDES E ESPINHACO)####
setwd("D:/Unifesp/Mestrado/analises/rasters_ib/PIT_LANU_all_rasters")

## Function to rescale a raster
raster01 = function(r){
  
  # get the min max values
  minmax_r = range(values(r), na.rm=TRUE)
  
  # rescale
  return( (r - minmax_r[1]) / (diff(minmax_r)))
}

## List all raster files in the directory (for example, all .asc files)
raster_files <- list.files(pattern = ".asc")

## Loop through each raster file, normalize it, and save it
for (file in raster_files) {
  
  # Open the raster
  r <- raster(file)
  
  # Plot the original raster (optional)
  plot(r)
  
  # Rescale the raster
  rn <- raster01(r)
  
  # Replace 0 values with 0.0001
  values(rn)[values(rn) == 0] <- 0.0001
  
  # Define the output filename (adding "ed_" as a prefix)
  output_file <- paste0("pit_normalizado_", file)
  
  # Save the normalized raster
  raster::writeRaster(rn, output_file, options = c("COMPRESS=DEFLATE"), format = "ascii", overwrite = TRUE)
  
  # Plot the rescaled raster (optional)
  plot(rn)
}



###################################DAQUI PRA BAIXO ADEQUAR DE ACORDO COM OS BANCOS DE DADOS QUE VC TIVER#######################

################LEMBRA QUE TEM QUE SER NO FORMATO ASCII!!!!!!!!!!!!######################


#normalizar rasters de variaveis para dora
################################################################################
setwd("D:/Unifesp/Mestrado/analises/GDM/DORA/raster/corte_rasters_amz_dora")

## Function to rescale a raster
raster01 = function(r){
  
  # get the min max values
  minmax_r = range(values(r), na.rm=TRUE)
  
  # rescale
  return( (r - minmax_r[1]) / (diff(minmax_r)))
}

## List all raster files in the directory (for example, all .asc files)
raster_files <- list.files(pattern = "\\.asc$")

## Loop through each raster file, normalize it, and save it
for (file in raster_files) {
  
  # Open the raster
  r <- raster(file)
  
  # Plot the original raster (optional)
  plot(r)
  
  # Rescale the raster
  rn <- raster01(r)
  
  # Replace 0 values with 0.0001
  values(rn)[values(rn) == 0] <- 0.0001
  
  # Define the output filename (adding "ed_" as a prefix)
  output_file <- paste0("norm_", file)
  
  # Save the normalized raster
  raster::writeRaster(rn, output_file, options = c("COMPRESS=DEFLATE"), format = "ascii", overwrite = TRUE)
  
  # Plot the rescaled raster (optional)
  plot(rn)
}
####################################
