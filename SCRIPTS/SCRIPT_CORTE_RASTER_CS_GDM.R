#################SCRIPT CORTE RASTERS CIRCUITSCAPE E GDM###############

library(raster)

######WORLDCLIM##### 
##lr:demorou 3h

# Definir o diretório de trabalho e a pasta de destino
setwd("D:/Unifesp/Mestrado/analises/rasters_ib/worldclim")
dir.create("cs_gdm_wc")
pasta_destino <- "cs_gdm_wc"

# Definir a extensão da região neotropical
extent_neotropical <- extent(-120, -25, -60, 30)

# Listar todos os arquivos .tif na pasta
arquivos <- list.files(pattern = ".tif")

# Iterar sobre cada arquivo e cortar
for (arquivo in arquivos) {
  # Carregar o raster
  raster_original <- raster(arquivo)
  
  # Verificar a projeção do raster
  projection_raster <- projection(raster_original)
  
  # Criar um RasterLayer a partir da extensão
  raster_extensao <- raster(extent_neotropical)
  crs(raster_extensao) <- projection_raster
  
  # Cortar o raster
  raster_cortado <- crop(raster_original, raster_extensao)
  
  # Salvar o raster cortado na nova pasta
  writeRaster(raster_cortado,
              filename = file.path(pasta_destino, arquivo),
              format = "ascii",
              overwrite = TRUE)
  print(paste("Raster", arquivo, "cortado com sucesso!"))
}


######EARTHENV#####
# Definir o diretório de trabalho e a pasta de destino
setwd("D:/Unifesp/Mestrado/analises/rasters_ib/earthenv")
dir.create("cs_gdm_earth")
pasta_destino <- "cs_gdm_earth"

# Definir a extensão da região neotropical
extent_neotropical <- extent(-120, -25, -60, 30)

# Listar todos os arquivos .tif na pasta
arquivos <- list.files(pattern = ".tif")

# Iterar sobre cada arquivo e cortar
for (arquivo in arquivos) {
  # Carregar o raster
  raster_original <- raster(arquivo)
  
  # Verificar a projeção do raster
  projection_raster <- projection(raster_original)
  
  # Criar um RasterLayer a partir da extensão
  raster_extensao <- raster(extent_neotropical)
  crs(raster_extensao) <- projection_raster
  
  # Cortar o raster
  raster_cortado <- crop(raster_original, raster_extensao)
  
  # Salvar o raster cortado na nova pasta
  writeRaster(raster_cortado,
              filename = file.path(pasta_destino, arquivo),
              format = "ascii",
              overwrite = TRUE)
  print(paste("Raster", arquivo, "cortado com sucesso!"))
}


######ENVIREM#####
# Definir o diretório de trabalho e a pasta de destino
setwd("D:/Unifesp/Mestrado/analises/rasters_ib/envirem")
dir.create("cs_gdm_env")
pasta_destino <- "cs_gdm_env"

# Definir a extensão da região neotropical
extent_neotropical <- extent(-120, -25, -60, 30)

# Listar todos os arquivos .tif na pasta
arquivos <- list.files(pattern = ".tif")

# Iterar sobre cada arquivo e cortar
for (arquivo in arquivos) {
  # Carregar o raster
  raster_original <- raster(arquivo)
  
  # Verificar a projeção do raster
  projection_raster <- projection(raster_original)
  
  # Criar um RasterLayer a partir da extensão
  raster_extensao <- raster(extent_neotropical)
  crs(raster_extensao) <- projection_raster
  
  # Cortar o raster
  raster_cortado <- crop(raster_original, raster_extensao)
  
  # Salvar o raster cortado na nova pasta
  writeRaster(raster_cortado,
              filename = file.path(pasta_destino, arquivo),
              format = "ascii",
              overwrite = TRUE)
  print(paste("Raster", arquivo, "cortado com sucesso!"))
}


######HARMONIZED WORLD DATABASE#####

#ELE JÁ É .ASC - LOOP ESTÁ DIFERENTE JA

# Definir o diretório de trabalho e a pasta de destino
setwd("D:/Unifesp/Mestrado/analises/rasters_ib/hwsd")
dir.create("cs_gdm_hwsd")
pasta_destino <- "cs_gdm_hwsd"

# Definir a extensão da região neotropical
extent_neotropical <- extent(-120, -25, -60, 30)

# Listar todos os arquivos .tif na pasta
arquivos <- list.files(pattern = ".asc")

# Iterar sobre cada arquivo e cortar
for (arquivo in arquivos) {
  # Carregar o raster
  raster_original <- raster(arquivo)
  
  # Verificar a projeção do raster
  projection_raster <- projection(raster_original)
  
  # Criar um RasterLayer a partir da extensão
  raster_extensao <- raster(extent_neotropical)
  crs(raster_extensao) <- projection_raster
  
  # Cortar o raster
  raster_cortado <- crop(raster_original, raster_extensao)
  
  # Salvar o raster cortado na nova pasta
  writeRaster(raster_cortado,
              filename = file.path(pasta_destino, arquivo),
              overwrite = TRUE)
  print(paste("Raster", arquivo, "cortado com sucesso!"))
}


######SOLARGIS#####

##REVER##

# Definir o diretório de trabalho e a pasta de destino
setwd("D:/Unifesp/Mestrado/analises/rasters_ib/worldclim")
dir.create("cs_gdm_wc")
pasta_destino <- "cs_gdm_wc"

# Definir a extensão da região neotropical
extent_neotropical <- extent(-120, -25, -60, 30)

# Listar todos os arquivos .tif na pasta
arquivos <- list.files(pattern = ".tif")

# Iterar sobre cada arquivo e cortar
for (arquivo in arquivos) {
  # Carregar o raster
  raster_original <- raster(arquivo)
  
  # Verificar a projeção do raster
  projection_raster <- projection(raster_original)
  
  # Criar um RasterLayer a partir da extensão
  raster_extensao <- raster(extent_neotropical)
  crs(raster_extensao) <- projection_raster
  
  # Cortar o raster
  raster_cortado <- crop(raster_original, raster_extensao)
  
  # Salvar o raster cortado na nova pasta
  writeRaster(raster_cortado,
              filename = file.path(pasta_destino, arquivo),
              format = "ascii",
              overwrite = TRUE)
  print(paste("Raster", arquivo, "cortado com sucesso!"))
}


