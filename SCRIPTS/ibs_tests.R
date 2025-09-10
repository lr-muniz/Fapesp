######SCRIPT DANTAS 2 - VRIESEA OLIGANTHA

#remotes::install_github("dinmatias/reconproGS")
library(adegenet)
library(reconproGS)
library(poppr)
library(vcfR)
library(hierfstat)
library(pegas)
library(dplyr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(Cairo)
library(radiator)
library(ape)
library(fossil)
library(BiodiversityR) # also loads vegan 
library(ecodist)
library(raster)
library(geobr)
library(reshape2)
library(genepop)
library(graph4lg)
library(apex)
library(pegas)
library(mmod)
library(poppr)
library(ade4)
library(sp)
library(tseries)
library(maptools)

#data(wrld_simpl)
#plot(wrld_simpl)

setwd("D:/Unifesp/Mestrado/analises/art_especies/dantas2/IB")
dir()
#input genetic distance
data1 <- read.table("input_fst_ib_vriesea_oligantha.txt", sep=';', header = F)

data1 <- data1[-1, ]

colnames(data1) = c("JAC" , "MKA" , "MCH" , "DIA" , "MUC" , "RCO" , "DIM" , "SGO" , "ABA" , "LIC" , "GMO" , "DIB" , "OUR" , "CIP")

rownames(data1) = c("JAC" , "MKA" , "MCH" , "DIA" , "MUC" , "RCO" , "DIM" , "SGO" , "ABA" , "LIC" , "GMO" , "DIB" , "OUR" , "CIP")

rgen <- as.dist(data1)

rgen

#input geographical distance (order: long/lat)
data <- read.table("input_coord_ib_vriesea_oligantha.txt", sep=';', header = T)

dist <- earth.dist(data, dist = T)

rgeo <- as.dist(dist)
rgeo

wrgeo <- melt(as.matrix(rgeo), varnames = c("row", "col"))

write.table(wrgeo, file="geographic_distance.txt")

#IBD - Mantel

gen_geo <- mantel.rtest(rgen, dist, nrepet = 9999)

gen_geo

A <- data.frame(gen_geo$obs, gen_geo$pvalue)

write.table(A, file="results_IBD_vriesea_oligantha.txt")

########################################################################
#IBE
# run IBD first

#IBE_envirem
data2 <-data.frame(data$longitude,data$latitude)

# import environmental database
setwd("D:/Unifesp/Mestrado/analises/rasters_ib/envirem")
dir()

td10 <- dir(pattern = "current")

r <- raster::stack(td10)

for(i in 1:length(td10)){ 
  # Extrair valores dos dados raster
  values <- extract(r[[i]], data2)
  
  # Opcional: tratar valores NA (se necessário)
  values[is.na(values)] <- 0
  
  # Criar um data frame com as coordenadas e valores extraídos
  df <- cbind.data.frame(coordinates(data2), values)
  
  # Calcular a distância ecológica
  eco_dist <- bcdist(values, rmzero = FALSE)
  
  # Converter o objeto dist em uma matriz antes de salvar
  eco_dist_matrix <- as.matrix(eco_dist)
  
  # Salvar o resultado de eco_dist em um arquivo .txt
  eco_dist_file <- file.path(paste0("ecodist_envirem_", i, ".txt"))
  write.table(eco_dist_matrix, file = eco_dist_file, row.names = FALSE, col.names = FALSE)
  
  # Realizar o teste de Mantel
  gen_env_total <- mantel.rtest(rgen, eco_dist, nrepet = 10000)
  
  # Criar um data frame com o resultado do Mantel
  A <- data.frame(gen_env_total$obs, gen_env_total$pvalue, td10[[i]])
  
  # Salvar o resultado do Mantel em um arquivo .txt
  mantel_file <- file.path(paste0("Mantel_IBE_envirem_", i, ".txt"))
  write.table(A, file = mantel_file, row.names = FALSE, col.names = TRUE)
}


#IBE_earthenv_freswater
setwd("D:/Unifesp/Mestrado/analises/rasters_ib/earthenv/freshwater/")
dir()

td10 <- dir(pattern = ".nc")

r <- raster::stack(td10)

for(i in 1:length(td10)){ 
  # Extrair valores dos dados raster
  values <- extract(r[[i]], data2)
  
  # Opcional: tratar valores NA (se necessário)
  values[is.na(values)] <- 0
  
  # Criar um data frame com as coordenadas e valores extraídos
  df <- cbind.data.frame(coordinates(data2), values)
  
  # Calcular a distância ecológica
  eco_dist <- bcdist(values, rmzero = FALSE)
  
  # Converter o objeto dist em uma matriz antes de salvar
  eco_dist_matrix <- as.matrix(eco_dist)
  
  # Salvar o resultado de eco_dist em um arquivo .txt
  eco_dist_file <- file.path(paste0("ecodist_earth_fw_", i, ".txt"))
  write.table(eco_dist_matrix, file = eco_dist_file, row.names = FALSE, col.names = FALSE)
  
  # Realizar o teste de Mantel
  gen_env_total <- mantel.rtest(rgen, eco_dist, nrepet = 10000)
  
  # Criar um data frame com o resultado do Mantel
  A <- data.frame(gen_env_total$obs, gen_env_total$pvalue, td10[[i]])
  
  # Salvar o resultado do Mantel em um arquivo .txt
  mantel_file <- file.path(paste0("Mantel_IBE_earth_fw_", i, ".txt"))
  write.table(A, file = mantel_file, row.names = FALSE, col.names = TRUE)
}

#IBE_earthenv_heterogeneity
setwd("D:/Unifesp/Mestrado/analises/rasters_ib/earthenv/heterogeneity/")
dir()
td10 <- dir(pattern = ".tif")

r <- raster::stack(td10)

for(i in 1:length(td10)){ 
  # Extrair valores dos dados raster
  values <- extract(r[[i]], data2)
  
  # Opcional: tratar valores NA (se necessário)
  # values[is.na(values)] <- 0
  
  # Criar um data frame com as coordenadas e valores extraídos
  df <- cbind.data.frame(coordinates(data2), values)
  
  # Calcular a distância ecológica
  eco_dist <- bcdist(values, rmzero = FALSE)
  
  # Converter o objeto dist em uma matriz antes de salvar
  eco_dist_matrix <- as.matrix(eco_dist)
  
  # Salvar o resultado de eco_dist em um arquivo .txt
  eco_dist_file <- file.path(paste0("ecodist_earth_het_", i, ".txt"))
  write.table(eco_dist_matrix, file = eco_dist_file, row.names = FALSE, col.names = FALSE)
  
  # Realizar o teste de Mantel
  gen_env_total <- mantel.rtest(rgen, eco_dist, nrepet = 10000)
  
  # Criar um data frame com o resultado do Mantel
  A <- data.frame(gen_env_total$obs, gen_env_total$pvalue, td10[[i]])
  
  # Salvar o resultado do Mantel em um arquivo .txt
  mantel_file <- file.path(paste0("Mantel_IBE_earth_het_", i, ".txt"))
  write.table(A, file = mantel_file, row.names = FALSE, col.names = TRUE)
}

#IBE_earthenv_land_cover
setwd("D:/Unifesp/Mestrado/analises/rasters_ib/earthenv/land_cover/")
dir()

td10 <- dir(pattern = ".tif")

r <- raster::stack(td10)

for(i in 1:length(td10)){ 
  # Extrair valores dos dados raster
  values <- extract(r[[i]], data2)
  
  # Opcional: tratar valores NA (se necessário)
  # values[is.na(values)] <- 0
  
  # Criar um data frame com as coordenadas e valores extraídos
  df <- cbind.data.frame(coordinates(data2), values)
  
  # Calcular a distância ecológica
  eco_dist <- bcdist(values, rmzero = FALSE)
  
  # Converter o objeto dist em uma matriz antes de salvar
  eco_dist_matrix <- as.matrix(eco_dist)
  
  # Salvar o resultado de eco_dist em um arquivo .txt
  eco_dist_file <- file.path(paste0("ecodist_earth_lc_", i, ".txt"))
  write.table(eco_dist_matrix, file = eco_dist_file, row.names = FALSE, col.names = FALSE)
  
  # Realizar o teste de Mantel
  gen_env_total <- mantel.rtest(rgen, eco_dist, nrepet = 10000)
  
  # Criar um data frame com o resultado do Mantel
  A <- data.frame(gen_env_total$obs, gen_env_total$pvalue, td10[[i]])
  
  # Salvar o resultado do Mantel em um arquivo .txt
  mantel_file <- file.path(paste0("Mantel_IBE_earth_lc_", i, ".txt"))
  write.table(A, file = mantel_file, row.names = FALSE, col.names = TRUE)
}

#IBE_earthenv_topography
setwd("D:/Unifesp/Mestrado/analises/rasters_ib/earthenv/topography/")
dir()

td10 <- dir(pattern = ".tif")

r <- raster::stack(td10)

for(i in 1:length(td10)){ 
  # Extrair valores dos dados raster
  values <- extract(r[[i]], data2)
  
  # Opcional: tratar valores NA (se necessário)
  # values[is.na(values)] <- 0
  
  # Criar um data frame com as coordenadas e valores extraídos
  df <- cbind.data.frame(coordinates(data2), values)
  
  # Calcular a distância ecológica
  eco_dist <- bcdist(values, rmzero = FALSE)
  
  # Converter o objeto dist em uma matriz antes de salvar
  eco_dist_matrix <- as.matrix(eco_dist)
  
  # Salvar o resultado de eco_dist em um arquivo .txt
  eco_dist_file <- file.path(paste0("ecodist_earth_topo_", i, ".txt"))
  write.table(eco_dist_matrix, file = eco_dist_file, row.names = FALSE, col.names = FALSE)
  
  # Realizar o teste de Mantel
  gen_env_total <- mantel.rtest(rgen, eco_dist, nrepet = 10000)
  
  # Criar um data frame com o resultado do Mantel
  A <- data.frame(gen_env_total$obs, gen_env_total$pvalue, td10[[i]])
  
  # Salvar o resultado do Mantel em um arquivo .txt
  mantel_file <- file.path(paste0("Mantel_IBE_earth_topo_", i, ".txt"))
  write.table(A, file = mantel_file, row.names = FALSE, col.names = TRUE)
}

#IBE_hwsd
setwd("D:/Unifesp/Mestrado/analises/rasters_ib/hwsd/")
dir()

td10 <- dir(pattern = "sq")

r <- raster::stack(td10)

for(i in 1:length(td10)){ 
  # Extrair valores dos dados raster
  values <- extract(r[[i]], data2)
  
  # Opcional: tratar valores NA (se necessário)
  # values[is.na(values)] <- 0
  
  # Criar um data frame com as coordenadas e valores extraídos
  df <- cbind.data.frame(coordinates(data2), values)
  
  # Calcular a distância ecológica
  eco_dist <- bcdist(values, rmzero = FALSE)
  
  # Converter o objeto dist em uma matriz antes de salvar
  eco_dist_matrix <- as.matrix(eco_dist)
  
  # Salvar o resultado de eco_dist em um arquivo .txt
  eco_dist_file <- file.path(paste0("ecodist_hwsd_", i, ".txt"))
  write.table(eco_dist_matrix, file = eco_dist_file, row.names = FALSE, col.names = FALSE)
  
  # Realizar o teste de Mantel
  gen_env_total <- mantel.rtest(rgen, eco_dist, nrepet = 10000)
  
  # Criar um data frame com o resultado do Mantel
  A <- data.frame(gen_env_total$obs, gen_env_total$pvalue, td10[[i]])
  
  # Salvar o resultado do Mantel em um arquivo .txt
  mantel_file <- file.path(paste0("Mantel_IBE_hwsd_", i, ".txt"))
  write.table(A, file = mantel_file, row.names = FALSE, col.names = TRUE)
}

#IBE_solargis_DIF
setwd("D:/Unifesp/Mestrado/analises/rasters_ib/solargis/DIF")
dir()

td10 <- dir(pattern = ".tif")

r <- raster::stack(td10)

for(i in 1:length(td10)){ 
  # Extrair valores dos dados raster
  values <- extract(r[[i]], data2)
  
  # Opcional: tratar valores NA (se necessário)
  values[is.na(values)] <- 0
  
  # Criar um data frame com as coordenadas e valores extraídos
  df <- cbind.data.frame(coordinates(data2), values)
  
  # Calcular a distância ecológica
  eco_dist <- bcdist(values, rmzero = FALSE)
  
  # Converter o objeto dist em uma matriz antes de salvar
  eco_dist_matrix <- as.matrix(eco_dist)
  
  # Salvar o resultado de eco_dist em um arquivo .txt
  eco_dist_file <- file.path(paste0("ecodist_solargis_dif_", i, ".txt"))
  write.table(eco_dist_matrix, file = eco_dist_file, row.names = FALSE, col.names = FALSE)
  
  # Realizar o teste de Mantel
  gen_env_total <- mantel.rtest(rgen, eco_dist, nrepet = 10000)
  
  # Criar um data frame com o resultado do Mantel
  A <- data.frame(gen_env_total$obs, gen_env_total$pvalue, td10[[i]])
  
  # Salvar o resultado do Mantel em um arquivo .txt
  mantel_file <- file.path(paste0("Mantel_IBE_solargis_dif_", i, ".txt"))
  write.table(A, file = mantel_file, row.names = FALSE, col.names = TRUE)
}


#IBE_solargis_DNI
setwd("D:/Unifesp/Mestrado/analises/rasters_ib/solargis/DNI")
dir()

td10 <- dir(pattern = ".tif")

r <- raster::stack(td10)

for(i in 1:length(td10)){ 
  # Extrair valores dos dados raster
  values <- extract(r[[i]], data2)
  
  # Opcional: tratar valores NA (se necessário)
  values[is.na(values)] <- 0
  
  # Criar um data frame com as coordenadas e valores extraídos
  df <- cbind.data.frame(coordinates(data2), values)
  
  # Calcular a distância ecológica
  eco_dist <- bcdist(values, rmzero = FALSE)
  
  # Converter o objeto dist em uma matriz antes de salvar
  eco_dist_matrix <- as.matrix(eco_dist)
  
  # Salvar o resultado de eco_dist em um arquivo .txt
  eco_dist_file <- file.path(paste0("ecodist_solargis_dni_", i, ".txt"))
  write.table(eco_dist_matrix, file = eco_dist_file, row.names = FALSE, col.names = FALSE)
  
  # Realizar o teste de Mantel
  gen_env_total <- mantel.rtest(rgen, eco_dist, nrepet = 10000)
  
  # Criar um data frame com o resultado do Mantel
  A <- data.frame(gen_env_total$obs, gen_env_total$pvalue, td10[[i]])
  
  # Salvar o resultado do Mantel em um arquivo .txt
  mantel_file <- file.path(paste0("Mantel_IBE_solargis_dni_", i, ".txt"))
  write.table(A, file = mantel_file, row.names = FALSE, col.names = TRUE)
}

#IBE_solargis_GHI
setwd("D:/Unifesp/Mestrado/analises/rasters_ib/solargis/GHI")
dir()

td10 <- dir(pattern = ".tif")

r <- raster::stack(td10)

for(i in 1:length(td10)){ 
  # Extrair valores dos dados raster
  values <- extract(r[[i]], data2)
  
  # Opcional: tratar valores NA (se necessário)
  values[is.na(values)] <- 0
  
  # Criar um data frame com as coordenadas e valores extraídos
  df <- cbind.data.frame(coordinates(data2), values)
  
  # Calcular a distância ecológica
  eco_dist <- bcdist(values, rmzero = FALSE)
  
  # Converter o objeto dist em uma matriz antes de salvar
  eco_dist_matrix <- as.matrix(eco_dist)
  
  # Salvar o resultado de eco_dist em um arquivo .txt
  eco_dist_file <- file.path(paste0("ecodist_solargis_ghi_", i, ".txt"))
  write.table(eco_dist_matrix, file = eco_dist_file, row.names = FALSE, col.names = FALSE)
  
  # Realizar o teste de Mantel
  gen_env_total <- mantel.rtest(rgen, eco_dist, nrepet = 10000)
  
  # Criar um data frame com o resultado do Mantel
  A <- data.frame(gen_env_total$obs, gen_env_total$pvalue, td10[[i]])
  
  # Salvar o resultado do Mantel em um arquivo .txt
  mantel_file <- file.path(paste0("Mantel_IBE_solargis_ghi_", i, ".txt"))
  write.table(A, file = mantel_file, row.names = FALSE, col.names = TRUE)
}

#IBE_solargis_GTI
setwd("D:/Unifesp/Mestrado/analises/rasters_ib/solargis/GTI")
dir()

td10 <- dir(pattern = ".tif")

r <- raster::stack(td10)

for(i in 1:length(td10)){ 
  # Extrair valores dos dados raster
  values <- extract(r[[i]], data2)
  
  # Opcional: tratar valores NA (se necessário)
  values[is.na(values)] <- 0
  
  # Criar um data frame com as coordenadas e valores extraídos
  df <- cbind.data.frame(coordinates(data2), values)
  
  # Calcular a distância ecológica
  eco_dist <- bcdist(values, rmzero = FALSE)
  
  # Converter o objeto dist em uma matriz antes de salvar
  eco_dist_matrix <- as.matrix(eco_dist)
  
  # Salvar o resultado de eco_dist em um arquivo .txt
  eco_dist_file <- file.path(paste0("ecodist_solargis_gti_", i, ".txt"))
  write.table(eco_dist_matrix, file = eco_dist_file, row.names = FALSE, col.names = FALSE)
  
  # Realizar o teste de Mantel
  gen_env_total <- mantel.rtest(rgen, eco_dist, nrepet = 10000)
  
  # Criar um data frame com o resultado do Mantel
  A <- data.frame(gen_env_total$obs, gen_env_total$pvalue, td10[[i]])
  
  # Salvar o resultado do Mantel em um arquivo .txt
  mantel_file <- file.path(paste0("Mantel_IBE_solargis_gti_", i, ".txt"))
  write.table(A, file = mantel_file, row.names = FALSE, col.names = TRUE)
}

#IBE_solargis_OPTA
setwd("D:/Unifesp/Mestrado/analises/rasters_ib/solargis/OPTA")
dir()

td10 <- dir(pattern = ".tif")

r <- raster::stack(td10)

for(i in 1:length(td10)){ 
  # Extrair valores dos dados raster
  values <- extract(r[[i]], data2)
  
  # Opcional: tratar valores NA (se necessário)
  values[is.na(values)] <- 0
  
  # Criar um data frame com as coordenadas e valores extraídos
  df <- cbind.data.frame(coordinates(data2), values)
  
  # Calcular a distância ecológica
  eco_dist <- bcdist(values, rmzero = FALSE)
  
  # Converter o objeto dist em uma matriz antes de salvar
  eco_dist_matrix <- as.matrix(eco_dist)
  
  # Salvar o resultado de eco_dist em um arquivo .txt
  eco_dist_file <- file.path(paste0("ecodist_solargis_opta_", i, ".txt"))
  write.table(eco_dist_matrix, file = eco_dist_file, row.names = FALSE, col.names = FALSE)
  
  # Realizar o teste de Mantel
  gen_env_total <- mantel.rtest(rgen, eco_dist, nrepet = 10000)
  
  # Criar um data frame com o resultado do Mantel
  A <- data.frame(gen_env_total$obs, gen_env_total$pvalue, td10[[i]])
  
  # Salvar o resultado do Mantel em um arquivo .txt
  mantel_file <- file.path(paste0("Mantel_IBE_solargis_opta_", i, ".txt"))
  write.table(A, file = mantel_file, row.names = FALSE, col.names = TRUE)
}

#IBE_solargis_PVOUT
setwd("D:/Unifesp/Mestrado/analises/rasters_ib/solargis/PVOUT")
dir()

td10 <- dir(pattern = ".tif")

r <- raster::stack(td10)

for(i in 1:length(td10)){ 
  # Extrair valores dos dados raster
  values <- extract(r[[i]], data2)
  
  # Opcional: tratar valores NA (se necessário)
  values[is.na(values)] <- 0
  
  # Criar um data frame com as coordenadas e valores extraídos
  df <- cbind.data.frame(coordinates(data2), values)
  
  # Calcular a distância ecológica
  eco_dist <- bcdist(values, rmzero = FALSE)
  
  # Converter o objeto dist em uma matriz antes de salvar
  eco_dist_matrix <- as.matrix(eco_dist)
  
  # Salvar o resultado de eco_dist em um arquivo .txt
  eco_dist_file <- file.path(paste0("ecodist_solargis_pvout_", i, ".txt"))
  write.table(eco_dist_matrix, file = eco_dist_file, row.names = FALSE, col.names = FALSE)
  
  # Realizar o teste de Mantel
  gen_env_total <- mantel.rtest(rgen, eco_dist, nrepet = 10000)
  
  # Criar um data frame com o resultado do Mantel
  A <- data.frame(gen_env_total$obs, gen_env_total$pvalue, td10[[i]])
  
  # Salvar o resultado do Mantel em um arquivo .txt
  mantel_file <- file.path(paste0("Mantel_IBE_solargis_pvout_", i, ".txt"))
  write.table(A, file = mantel_file, row.names = FALSE, col.names = TRUE)
}

#IBE_solargis_TEMP
setwd("D:/Unifesp/Mestrado/analises/rasters_ib/solargis/TEMP")
dir()

td10 <- dir(pattern = ".tif")

r <- raster::stack(td10)

for(i in 1:length(td10)){ 
  # Extrair valores dos dados raster
  values <- extract(r[[i]], data2)
  
  # Opcional: tratar valores NA (se necessário)
  # values[is.na(values)] <- 0
  
  # Criar um data frame com as coordenadas e valores extraídos
  df <- cbind.data.frame(coordinates(data2), values)
  
  # Calcular a distância ecológica
  eco_dist <- bcdist(values, rmzero = FALSE)
  
  # Converter o objeto dist em uma matriz antes de salvar
  eco_dist_matrix <- as.matrix(eco_dist)
  
  # Salvar o resultado de eco_dist em um arquivo .txt
  eco_dist_file <- file.path(paste0("ecodist_solargis_temp_", i, ".txt"))
  write.table(eco_dist_matrix, file = eco_dist_file, row.names = FALSE, col.names = FALSE)
  
  # Realizar o teste de Mantel
  gen_env_total <- mantel.rtest(rgen, eco_dist, nrepet = 10000)
  
  # Criar um data frame com o resultado do Mantel
  A <- data.frame(gen_env_total$obs, gen_env_total$pvalue, td10[[i]])
  
  # Salvar o resultado do Mantel em um arquivo .txt
  mantel_file <- file.path(paste0("Mantel_IBE_solargis_temp_", i, ".txt"))
  write.table(A, file = mantel_file, row.names = FALSE, col.names = TRUE)
}

#IBE_worldclim
setwd("D:/Unifesp/Mestrado/analises/rasters_ib/worldclim/")
dir()

td10 <- dir(pattern = ".tif")

r <- raster::stack(td10)

for(i in 1:length(td10)){ 
  # Extrair valores dos dados raster
  values <- extract(r[[i]], data2)
  
  # Opcional: tratar valores NA (se necessário)
  # values[is.na(values)] <- 0
  
  # Criar um data frame com as coordenadas e valores extraídos
  df <- cbind.data.frame(coordinates(data2), values)
  
  # Calcular a distância ecológica
  eco_dist <- bcdist(values, rmzero = FALSE)
  
  # Converter o objeto dist em uma matriz antes de salvar
  eco_dist_matrix <- as.matrix(eco_dist)
  
  # Salvar o resultado de eco_dist em um arquivo .txt
  eco_dist_file <- file.path(paste0("ecodist_wc_", i, ".txt"))
  write.table(eco_dist_matrix, file = eco_dist_file, row.names = FALSE, col.names = FALSE)
  
  # Realizar o teste de Mantel
  gen_env_total <- mantel.rtest(rgen, eco_dist, nrepet = 10000)
  
  # Criar um data frame com o resultado do Mantel
  A <- data.frame(gen_env_total$obs, gen_env_total$pvalue, td10[[i]])
  
  # Salvar o resultado do Mantel em um arquivo .txt
  mantel_file <- file.path(paste0("Mantel_IBE_wc_", i, ".txt"))
  write.table(A, file = mantel_file, row.names = FALSE, col.names = TRUE)
}




##############################################################################

########### IBR #####

#nesse caso, dá pra rodar todas as matrizes de resistência juntas, uma vez que, no meu caso, todas as matrizes de resistência foram colocadas na mesma pasta.
##ele faz a correlação com a distancia genética e a matriz de resistencia das variáveis que foram geradas no circuitscape

library(ade4)

# 1. Definir diretório e arquivos
setwd("D:/Unifesp/Mestrado/analises/GDM/vriesea_oligantha_e/")

# Matrizes de resistência (Circuitscape)
resistance_files <- list.files(pattern = "*resistances.out", full.names = TRUE)

# Matriz genética (FST)
genetic_file <- list.files(pattern = "input_fst_cs_vriesea_oligantha.out")  # Substitua pelo seu arquivo

# 2. Ler a matriz genética
genetic <- as.dist(read.table(genetic_file, sep = " ", header = TRUE, row.names = 1))

# 3. Criar arquivo de resultados
output_file <- "resultados_IBR_vrie.txt"
write("Resultados de IBR (Mantel entre FST e matrizes de resistência)\n", file = output_file)

# 4. Loop para testar cada matriz de resistência
for (r_file in resistance_files) {
  # Ler matriz de resistência
  rdist <- as.dist(read.csv(r_file, sep = " ", row.names = 1, header = TRUE))
  
  # Teste de Mantel
  mantel_result <- mantel.rtest(genetic, rdist, nrepet = 9999)
  
  # Salvar resultados
  write(paste("\n--- Arquivo de resistência:", r_file, "---"), file = output_file, append = TRUE)
  write(paste("Correlação (R):", mantel_result$obs), file = output_file, append = TRUE)
  write(paste("Valor-p:", mantel_result$pvalue), file = output_file, append = TRUE)
}

# Opcional: Visualizar resultados no console
print(readLines(output_file))

