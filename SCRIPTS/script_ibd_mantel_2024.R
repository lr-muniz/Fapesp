#remotes::install_github("dinmatias/reconproGS")
library(adegenet)
library(reconproGS)
library("poppr")
library("vcfR")
library("hierfstat")
library("pegas")
library(dplyr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(Cairo)
library(radiator)
library(ape)
library("fossil")
library(BiodiversityR) # also loads vegan 
library(ecodist)
library(raster)
library(geobr)
library(reshape2)
library(genepop)
library("graph4lg")
library("apex")
library("pegas")
library("mmod")
library("poppr")
library(ade4)
library(sp)
library(tseries)
library(maptools)
#install.packages("maptools")
data(wrld_simpl)
plot(wrld_simpl)
setwd("D:/Unifesp/Mestrado/analises/jauru")

data1 <- read.table("input_fst_jauru.txt", sep='\t', header = F)

data1 <- data1[-1,-1 ]

colnames(data1) = c("ALC", "RVE", "PET", "MIN", "PGO") 
rownames(data1) = c("ALC", "RVE", "PET", "MIN", "PGO") 

rgen <- as.dist(data1)

rgen
#entrada dado distancia genetica

#geographical distance

data <- read.table("input_coord_jauru.txt", header = T)

dist <- earth.dist(data, dist = T)

rgeo <- as.dist(dist)
rgeo

wrgeo <- melt(as.matrix(rgeo), varnames = c("row", "col"))

write.table(wrgeo, file="distancia_geografica.txt")

## MANTEL

#IBD

gen_geo <- mantel.rtest(rgen, dist, nrepet = 9999)

gen_geo

A <- data.frame(gen_geo$obs, gen_geo$pvalue)

write.table(A, file="Mantel_IBD.txt")
A
##IBR
# Make a place holder for the cs_run.exe path
CS_exe <- 'C:/Users/pc/Desktop/imputs/imputs/Circuitscape/cs_run.exe' # Don't forget the "Program Files" problem

#IBR_envirem

setwd("C:/Users/pc/Desktop/imputs/jauri/envirem/")

td3 <- dir(pattern = ".asc")
td3

data
var.t <- raster::stack(td3)

# Plot it
# Rasterize points using the cost extent

sites <- rasterize(x = data,y = var.t)

values(sites)[values(sites) == 0] <- 0.0001
values(var.t)[values(var.t) == 0] <- 0.0001

# Write rasters to your working directory
writeRaster(sites,"sites_rast.asc",overwrite=TRUE)

plot(var.t[[1]])
points(data$longitude, data$latitude, col="red", pch=20, cex=0.9)

for(i in 1:length(td3)){ 
writeRaster(var.t[[i]],file.path(paste0("cost_rast_",i,".asc")),overwrite=TRUE)

# Make an .ini file
CS_ini <- c("[circuitscape options]",            
            "data_type = raster",
            "scenario = pairwise",
            paste(c("point_file =",
                    "habitat_file =",
                    "output_file ="),
                  paste(getwd(),c("sites_rast.asc",
                                  file.path(paste0("cost_rast_",i,".asc")),
                                  file.path(paste0("CS_",i,".out"))),
                        sep="/")))

# Write it to your working directory
writeLines(CS_ini,file.path(paste0("myini_",i,".ini")))

# Make the CS run cmd
CS_run <- paste(CS_exe, paste(getwd(),file.path(paste0("myini_",i,".ini")),sep="/")) # Make the cmd

CS_run
# Run the command
system(CS_run)

# Import the effective resistance
rdist <- as.dist(read.csv(file.path(paste0("CS_",i,"_resistances.out")),sep=" ",row.names=1,header=1))

wrdist <- melt(as.matrix(rdist), varnames = c("row", "col"))

write.table(wrdist, file=file.path(paste0("distancia_envirem_",i,"_resistencia.txt")))

gen_ibr_clima <- mantel.rtest(rgen, rdist, nrepet = 9999)

A <- data.frame(gen_ibr_clima$obs, gen_ibr_clima$pvalue)
  
write.table(A, file=file.path(paste0("Mantel_envirem_",i,".txt")))
}

#############################################################################
#IBR_solarGis

setwd("C:/Users/pc/Desktop/imputs/jauri/solargis/")

td4 <- dir(pattern = ".asc")

var.t2 <- raster::stack(td4)

# Rasterize points using the cost extent
sites <- rasterize(x = data,y = var.t2)

values(sites)[values(sites) == 0] <- 0.0001
values(var.t2)[values(var.t2) == 0] <- 0.0001

# Write rasters to your working directory
writeRaster(sites,"sites_rast.asc",overwrite=TRUE)

for(i in 1:length(td4)){ 
  writeRaster(var.t2[[i]],file.path(paste0("cost_rast_",i,".asc")),overwrite=TRUE)
  
  # Make an .ini file
  CS_ini <- c("[circuitscape options]",            
              "data_type = raster",
              "scenario = pairwise",
              paste(c("point_file =",
                      "habitat_file =",
                      "output_file ="),
                    paste(getwd(),c("sites_rast.asc",
                                    file.path(paste0("cost_rast_",i,".asc")),
                                    file.path(paste0("CS_",i,".out"))),
                          sep="/")))
  
  # Write it to your working directory
  writeLines(CS_ini,file.path(paste0("myini_",i,".ini")))
  
  # Make the CS run cmd
  CS_run <- paste(CS_exe, paste(getwd(),file.path(paste0("myini_",i,".ini")),sep="/")) # Make the cmd
  
  # Run the command
  system(CS_run)
  
  # Import the effective resistance
  rdist <- as.dist(read.csv(file.path(paste0("CS_",i,"_resistances.out")),sep=" ",row.names=1,header=1))
  
  wrdist <- melt(as.matrix(rdist), varnames = c("row", "col"))
  
  write.table(wrdist, file=file.path(paste0("distancia_solargis_",i,"_resistencia.txt")))
  
  gen_ibr_solargis1 <- mantel.rtest(rgen, rdist, nrepet = 9999)
  
  A <- data.frame(gen_ibr_solargis1$obs, gen_ibr_solargis1$pvalue)
  
  write.table(A, file=file.path(paste0("Mantel_solargis_",i,".txt")))
}

#############################################################################

#IBR_solarGis2

setwd("C:/Users/pc/Desktop/imputs/jauri/solargis2/")

td8 <- dir(pattern = ".asc")
td8
var.t8 <- raster::stack(td8)

# Rasterize points using the cost extent
sites <- rasterize(x = data,y = var.t8)

values(sites)[values(sites) == 0] <- 0.0001
values(var.t8)[values(var.t8) == 0] <- 0.0001

# Write rasters to your working directory
writeRaster(sites,"sites_rast.asc",overwrite=TRUE)

for(i in 1:length(td8)){ 
  writeRaster(var.t8[[i]],file.path(paste0("cost_rast_",i,".asc")),overwrite=TRUE)
  
  # Make an .ini file
  CS_ini <- c("[circuitscape options]",            
              "data_type = raster",
              "scenario = pairwise",
              paste(c("point_file =",
                      "habitat_file =",
                      "output_file ="),
                    paste(getwd(),c("sites_rast.asc",
                                    file.path(paste0("cost_rast_",i,".asc")),
                                    file.path(paste0("CS_",i,".out"))),
                          sep="/")))
  
  # Write it to your working directory
  writeLines(CS_ini,file.path(paste0("myini_",i,".ini")))
  
  # Make the CS run cmd
  CS_run <- paste(CS_exe, paste(getwd(),file.path(paste0("myini_",i,".ini")),sep="/")) # Make the cmd
  
  # Run the command
  system(CS_run)
  
  # Import the effective resistance
  rdist <- as.dist(read.csv(file.path(paste0("CS_",i,"_resistances.out")),sep=" ",row.names=1,header=1))
  
  wrdist <- melt(as.matrix(rdist), varnames = c("row", "col"))
  
  write.table(wrdist, file=file.path(paste0("distancia_solargis2_",i,"_resistencia.txt")))
  
  gen_ibr_solargis2 <- mantel.rtest(rgen, rdist, nrepet = 9999)
  
  A <- data.frame(gen_ibr_solargis2$obs, gen_ibr_solargis2$pvalue)
  
  write.table(A, file=file.path(paste0("Mantel_solargis2_",i,".txt")))
}

#############################################################################

#IBR_Solo

setwd("C:/Users/pc/Desktop/imputs/jauri/hwsd/")

td5 <- dir(pattern = ".asc")
td5
var.t3 <- raster::stack(td5)

# Rasterize points using the cost extent
sites <- rasterize(x = data,y = var.t3)

values(sites)[values(sites) == 0] <- 0.0001
values(var.t3)[values(var.t3) == 0] <- 0.0001

# Write rasters to your working directory
writeRaster(sites,"sites_rast.asc",overwrite=TRUE)

for(i in 1:length(td5)){ 
  writeRaster(var.t3[[i]],file.path(paste0("cost_rast_",i,".asc")),overwrite=TRUE)
  
  # Make an .ini file
  CS_ini <- c("[circuitscape options]",            
              "data_type = raster",
              "scenario = pairwise",
              paste(c("point_file =",
                      "habitat_file =",
                      "output_file ="),
                    paste(getwd(),c("sites_rast.asc",
                                    file.path(paste0("cost_rast_",i,".asc")),
                                    file.path(paste0("CS_",i,".out"))),
                          sep="/")))
  
  # Write it to your working directory
  writeLines(CS_ini,file.path(paste0("myini_",i,".ini")))
  
  # Make the CS run cmd
  CS_run <- paste(CS_exe, paste(getwd(),file.path(paste0("myini_",i,".ini")),sep="/")) # Make the cmd
  
  # Run the command
  system(CS_run)
  
  # Import the effective resistance
  rdist <- as.dist(read.csv(file.path(paste0("CS_",i,"_resistances.out")),sep=" ",row.names=1,header=1))
  
  wrdist <- melt(as.matrix(rdist), varnames = c("row", "col"))
  
  write.table(wrdist, file=file.path(paste0("distancia_solo_",i,"_resistencia.txt")))
  
  gen_ibr_radiacao <- mantel.rtest(rgen, rdist, nrepet = 9999)
  
  A <- data.frame(gen_ibr_radiacao$obs, gen_ibr_radiacao$pvalue)
  
  write.table(A, file=file.path(paste0("Mantel_solo_",i,".txt")))
}
#############################################################################
#IBR_Clima

setwd("C:/Users/pc/Desktop/imputs/jauri/clima/")

td6 <- dir(pattern = ".asc")
td6

var.t6 <- raster::stack(td6)

# Rasterize points using the cost extent
sites <- rasterize(x = data,y = var.t6)

plot(sites)

values(sites)[values(sites) == 0] <- 0.0001
values(var.t6)[values(var.t6) == 0] <- 0.0001

# Write rasters to your working directory
writeRaster(sites,"sites_rast.asc",overwrite=TRUE)

for(i in 1:length(td6)){ 
  writeRaster(var.t6[[i]],file.path(paste0("cost_rast_",i,".asc")),overwrite=TRUE)
  
  # Make an .ini file
  CS_ini <- c("[circuitscape options]",            
              "data_type = raster",
              "scenario = pairwise",
              paste(c("point_file =",
                      "habitat_file =",
                      "output_file ="),
                    paste(getwd(),c("sites_rast.asc",
                                    file.path(paste0("cost_rast_",i,".asc")),
                                    file.path(paste0("CS_",i,".out"))),
                          sep="/")))
  
  # Write it to your working directory
  writeLines(CS_ini,file.path(paste0("myini_",i,".ini")))
  
  # Make the CS run cmd
  CS_run <- paste(CS_exe, paste(getwd(),file.path(paste0("myini_",i,".ini")),sep="/")) # Make the cmd
  
  # Run the command
  system(CS_run)
  
  # Import the effective resistance
  rdist <- as.dist(read.csv(file.path(paste0("CS_",i,"_resistances.out")),sep=" ",row.names=1,header=1))
  
  wrdist <- melt(as.matrix(rdist), varnames = c("row", "col"))
  
  write.table(wrdist, file=file.path(paste0("distancia_clima_",i,"_resistencia.txt")))
  
  gen_ibr_solo <- mantel.rtest(rgen, rdist, nrepet = 9999)
  
  A <- data.frame(gen_ibr_solo$obs, gen_ibr_solo$pvalue)
  
  write.table(A, file=file.path(paste0("Mantel_clima_",i,".txt")))
}
