### script download occurrences - spocc ###

## Abrir os dados baixados e acertar as colunas da espécie####

## packages
library(colorRamps)#colore mapas
library(mapr)#faz mapas
library(raster)#GIS
library(sf)#GIs
library(tidyverse)#tabelas
library(readr)
# check loaded packages
search()

###---------------------------------------------------------------------------------------###
## data
## occ
# directory
setwd("D:/Unifesp/Mestrado/analises/modelagem/gbif")#direcionar para a sua pasta

# import data
list.files()

#"Euphorbia attastoma.csv"    
[11] "Lychnophora ericoides.csv"   "Oxalis laciniata.csv"       
[13] "Pitcairnia lanuginosa.csv"   "Richterago discoidea.csv"   
[15] "Tillandsia capillaris.csv"   "Tillandsia virescens.csv"   
[17] "Vellozia auriculata.csv"     "Vriesea oligantha.csv"  


occ <- read.delim("Tillandsia virescens.csv", header = TRUE, stringsAsFactors = FALSE)
occ<- data.frame(occ$species, occ$decimalLongitude, occ$decimalLatitude)
colnames(occ)<-c("species","decimalLongitude","decimalLatitude")
occ
## variaveis
# directory

# import data
#var <- raster::raster("wc2.1_2.5m_bio_10.tif")
#var

###---------------------------------------------------------------------------------------###

## clear data

# 1. registers and date
occ.cl <- occ %>% 
  dplyr::distinct(decimalLongitude, decimalLatitude, .keep_all = TRUE) %>% # remove duplicates
  dplyr::filter(!is.na(decimalLongitude)) %>% # remove decimalLongitude with NAs
  dplyr::filter(!is.na(decimalLatitude)) # remove decimalLatitude with NAs
  # filter by date
occ.cl

###---------------------------------------------------------------------------------------###

## 2. points inside limit
# plot
ggplot() +
  geom_point(data = occ.cl %>% dplyr::select(decimalLongitude, decimalLatitude), aes(x = decimalLongitude, y = decimalLatitude), 
             color = "black", alpha = .3) +
  theme_minimal()

readr::write_csv(occ.cl, "Tillandsia_clean.csv")
getwd()
dir()
###---------------------------------------------------------------------------------------###


