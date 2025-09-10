### script frequency ensemble ###

# packages
library(colorRamps)
library(raster)
library(tidyverse)
library(Cairo)

# list packages
search()

###----------------------------------------------------------------------------###

# import data
# directory
setwd("D:/Unifesp/Mestrado/analises/modelagem/raster_novo/03_cedrela_50/02_modelos_multiplos/")
getwd()

#ja foi: 

# enms
tif <- dir(patt = ".tif$")
tif

pres = grep("_pres_",tif, value=T)
pres

# import evaluates
csv <- purrr::map_dfr(dir(patt = ".csv"), readr::read_csv)
csv

enm <- var <- raster::stack(pres[1])

###----------------------------------------------------------------------------###

## frequency ensemble 
# species
sp <- csv %>% 
  dplyr::select(species) %>% 
  dplyr::distinct() %>% 
  dplyr::pull()
sp <- sp[!is.na(sp)]
sp
sp2 <- gsub("\\s*\\(.*\\)", "", sp)
sp2
# algorithms
al <- csv %>% 
  dplyr::select(algorithm) %>% 
  dplyr::distinct() %>%
  dplyr::pull() %>% 
  stringr::str_replace("_", "")
al <- al[!is.na(al)]
al

# ensembles
ens <- enm
ens[] <- 0
ens
plot(ens)
# directory
dir.create("00_ensemble_freq")

# for presente
for(i in sp2){
  # select model by species
  tif.sp <- grep(i, pres, value = TRUE)
  eva.sp <- csv %>% dplyr::filter(species == sp)
  tif.sp
  eva.sp
  # information
  print(paste0("The ensemble for ", i, " started, relax, take a coffee, it may take awhile..."))
  
  for(j in al){

    # select model by algorithms
    tif.al <- grep(j, tif.sp, value = TRUE)
    eva.al <- eva.sp %>% dplyr::filter(algorithm == j)
    tif.al
    eva.al
    # information
    print(paste0("The ensemble for '", i, "', algorithm '", j, "' are going!"))
    
    # import raster
    enm.al <- stack(tif.al)
    enm.al
    for(k in seq(length(tif.al))){
   
        # sum
      ens <- sum(ens, enm.al[[k]] >= eva.al[k, 5] %>% dplyr::pull())
      eva.al[k, 5]
    }
    
  }
  
  # export
  setwd("00_ensemble_freq")
  writeRaster(ens / (length(tif.sp)), paste0("ensemble_freq_", i, "_present.tif"), 
              options = c("COMPRESS=DEFLATE"), format = "GTiff", overwrite = TRUE)
  setwd("..")
  
  # information
  print(paste0("Nice! The ensemble of ", i, " it's done!"))
  
  
  ens[] <- 0
  
  print("Yeh! It's over!!!")
  
}

 ###----------------------------------------------------------------------------###

# directory
setwd("00_ensemble_freq")
getwd()
# import
tif <- dir(patt = ".tif$")
tif

pres = grep("_present",tif, value=T)
pres

mo <- raster::stack(pres)
mo

# map
# occurrences
cores <- colorRampPalette(c("blue", "red"))(20)
plot(mo, col = cores)

CairoPDF(width = 36, height = 34, file = "Present.pdf", pointsize=40, bg = "white", title = "", paper = "special", pagecentre=TRUE) #
plot(mo, col = cores)
dev.off()
###----------------------------------------------------------------------------###

#HS
setwd("../")
tif <- dir(patt = ".tif$")
tif

hs = grep("_hs_",tif, value=T)
hs

# for presente
for(i in sp2){
  # select model by species
  tif.sp <- grep(i, hs, value = TRUE)
  eva.sp <- csv %>% dplyr::filter(species == sp)

  # information
  print(paste0("The ensemble for ", i, " started, relax, take a coffee, it may take awhile..."))
  
  for(j in al){
    
    # select model by algorithms
    tif.al <- grep(j, tif.sp, value = TRUE)
    eva.al <- eva.sp %>% dplyr::filter(algorithm == j)
    tif.al

    # information
    print(paste0("The ensemble for '", i, "', algorithm '", j, "' are going!"))
    
    # import raster
    enm.al <- stack(tif.al)
    tif.al
    for(k in seq(length(tif.al))){
      
      # sum
      ens <- sum(ens, enm.al[[k]] >= eva.al[k, 5] %>% dplyr::pull())
      eva.al[k, 5]
    }
    
  }
  
  # export
  setwd("00_ensemble_freq")
  writeRaster(ens / (length(tif.sp)), paste0("ensemble_freq_", i, "_hs.tif"), 
              options = c("COMPRESS=DEFLATE"), format = "GTiff", overwrite = TRUE)
  setwd("..")
  
  # information
  print(paste0("Nice! The ensemble of ", i, " it's done!"))
  
  
  ens[] <- 0
  
  print("Yeh! It's over!!!")
  
}

###----------------------------------------------------------------------------###

# directory
setwd("00_ensemble_freq")
getwd()
# import

tif <- dir(patt = ".tif$")
tif

hsl = grep("_hs",tif, value=T)
hsl
mhs <- stack(hsl)
mhs

# map
# occurrences

plot(mhs, col = cores)

CairoPDF(width = 36, height = 34, file = "HS.pdf", pointsize=40, bg = "white", title = "", paper = "special", pagecentre=TRUE) #
plot(mhs, col = cores)
dev.off()

##########################
#LGM
setwd("../")
tif <- dir(patt = ".tif$")
tif

lgm = grep("_lgm_",tif, value=T)
lgm

# for presente
for(i in sp2){
  # select model by species
  tif.sp <- grep(i, lgm, value = TRUE)
  eva.sp <- csv %>% dplyr::filter(species == sp)
  
  # information
  print(paste0("The ensemble for ", i, " started, relax, take a coffee, it may take awhile..."))
  
  for(j in al){
    
    # select model by algorithms
    tif.al <- grep(j, tif.sp, value = TRUE)
    eva.al <- eva.sp %>% dplyr::filter(algorithm == j)
    tif.al
    
    # information
    print(paste0("The ensemble for '", i, "', algorithm '", j, "' are going!"))
    
    # import raster
    enm.al <- stack(tif.al)
    tif.al
    for(k in seq(length(tif.al))){
      
      # sum
      ens <- sum(ens, enm.al[[k]] >= eva.al[k, 5] %>% dplyr::pull())
      eva.al[k, 5]
    }
    
  }
  
  # export
  setwd("00_ensemble_freq")
  writeRaster(ens / (length(tif.sp)), paste0("ensemble_freq_", i, "_lgm.tif"), 
              options = c("COMPRESS=DEFLATE"), format = "GTiff", overwrite = TRUE)
  setwd("..")
  
  # information
  print(paste0("Nice! The ensemble of ", i, " it's done!"))
  
  
  ens[] <- 0
  
  print("Yeh! It's over!!!")
  
}

###----------------------------------------------------------------------------###

# directory
setwd("00_ensemble_freq")
getwd()
# import

tif <- dir(patt = ".tif$")
tif

lgml = grep("_lgm",tif, value=T)
lgml
mlgm <- stack(lgml)
mlgm

# map
# occurrences

plot(mlgm, col = cores)

CairoPDF(width = 36, height = 34, file = "LGM.pdf", pointsize=40, bg = "white", title = "", paper = "special", pagecentre=TRUE) #
plot(mlgm, col = cores)
dev.off()

##########################
#LIG
setwd("../")
tif <- dir(patt = ".tif$")
tif

lig = grep("_lig_",tif, value=T)
lig

# for presente
for(i in sp2){
  # select model by species
  tif.sp <- grep(i, lig, value = TRUE)
  eva.sp <- csv %>% dplyr::filter(species == sp)
  
  # information
  print(paste0("The ensemble for ", i, " started, relax, take a coffee, it may take awhile..."))
  
  for(j in al){
    
    # select model by algorithms
    tif.al <- grep(j, tif.sp, value = TRUE)
    eva.al <- eva.sp %>% dplyr::filter(algorithm == j)
    tif.al
    
    # information
    print(paste0("The ensemble for '", i, "', algorithm '", j, "' are going!"))
    
    # import raster
    enm.al <- stack(tif.al)
    tif.al
    for(k in seq(length(tif.al))){
      
      # sum
      ens <- sum(ens, enm.al[[k]] >= eva.al[k, 5] %>% dplyr::pull())
      eva.al[k, 5]
    }
    
  }
  
  # export
  setwd("00_ensemble_freq")
  writeRaster(ens / (length(tif.sp)), paste0("ensemble_freq_", i, "_lig.tif"), 
              options = c("COMPRESS=DEFLATE"), format = "GTiff", overwrite = TRUE)
  setwd("..")
  
  # information
  print(paste0("Nice! The ensemble of ", i, " it's done!"))
  
  
  ens[] <- 0
  
  print("Yeh! It's over!!!")
  
}

###----------------------------------------------------------------------------###

# directory
setwd("00_ensemble_freq")
getwd()
# import

tif <- dir(patt = ".tif$")
tif

ligl = grep("_lig",tif, value=T)
ligl
mlig <- stack(ligl)
mlig

# map
# occurrences

plot(mlig, col = cores)

CairoPDF(width = 36, height = 34, file = "LIG.pdf", pointsize=40, bg = "white", title = "", paper = "special", pagecentre=TRUE) #
plot(mlig, col = cores)
dev.off()

##########################
#MH
setwd("../")
tif <- dir(patt = ".tif$")
tif

mh = grep("_mh_",tif, value=T)
mh

# for presente
for(i in sp2){
  # select model by species
  tif.sp <- grep(i, mh, value = TRUE)
  eva.sp <- csv %>% dplyr::filter(species == sp)
  
  # information
  print(paste0("The ensemble for ", i, " started, relax, take a coffee, it may take awhile..."))
  
  for(j in al){
    
    # select model by algorithms
    tif.al <- grep(j, tif.sp, value = TRUE)
    eva.al <- eva.sp %>% dplyr::filter(algorithm == j)
    tif.al
    
    # information
    print(paste0("The ensemble for '", i, "', algorithm '", j, "' are going!"))
    
    # import raster
    enm.al <- stack(tif.al)
    tif.al
    for(k in seq(length(tif.al))){
      
      # sum
      ens <- sum(ens, enm.al[[k]] >= eva.al[k, 5] %>% dplyr::pull())
      eva.al[k, 5]
    }
    
  }
  
  # export
  setwd("00_ensemble_freq")
  writeRaster(ens / (length(tif.sp)), paste0("ensemble_freq_", i, "_mh.tif"), 
              options = c("COMPRESS=DEFLATE"), format = "GTiff", overwrite = TRUE)
  setwd("..")
  
  # information
  print(paste0("Nice! The ensemble of ", i, " it's done!"))
  
  
  ens[] <- 0
  
  print("Yeh! It's over!!!")
  
}


###----------------------------------------------------------------------------###

# directory
setwd("00_ensemble_freq")
getwd()
# import

tif <- dir(patt = ".tif$")
tif

mhl = grep("_mh",tif, value=T)
mhl
mmh <- stack(mhl)
mmh

# map
# occurrences

plot(mmh, col = cores)

CairoPDF(width = 36, height = 34, file = "MH.pdf", pointsize=40, bg = "white", title = "", paper = "special", pagecentre=TRUE) #
plot(mmh, col = cores)
dev.off()

##########################
#YD
setwd("../")
tif <- dir(patt = ".tif$")
tif

yd = grep("_yd_",tif, value=T)
yd

# for presente
for(i in sp2){
  # select model by species
  tif.sp <- grep(i, yd, value = TRUE)
  eva.sp <- csv %>% dplyr::filter(species == sp)
  
  # information
  print(paste0("The ensemble for ", i, " started, relax, take a coffee, it may take awhile..."))
  
  for(j in al){
    
    # select model by algorithms
    tif.al <- grep(j, tif.sp, value = TRUE)
    eva.al <- eva.sp %>% dplyr::filter(algorithm == j)
    tif.al
    
    # information
    print(paste0("The ensemble for '", i, "', algorithm '", j, "' are going!"))
    
    # import raster
    enm.al <- stack(tif.al)
    tif.al
    for(k in seq(length(tif.al))){
      
      # sum
      ens <- sum(ens, enm.al[[k]] >= eva.al[k, 5] %>% dplyr::pull())
      eva.al[k, 5]
    }
    
  }
  
  # export
  setwd("00_ensemble_freq")
  writeRaster(ens / (length(tif.sp)), paste0("ensemble_freq_", i, "_yd.tif"), 
              options = c("COMPRESS=DEFLATE"), format = "GTiff", overwrite = TRUE)
  setwd("..")
  
  # information
  print(paste0("Nice! The ensemble of ", i, " it's done!"))
  
  
  ens[] <- 0
  
  print("Yeh! It's over!!!")
  
}

  ###----------------------------------------------------------------------------###

# directory
setwd("00_ensemble_freq")
getwd()
# import

tif <- dir(patt = ".tif$")
tif

ydl = grep("_yd",tif, value=T)
ydl
myd <- stack(ydl)
myd

# map
# occurrences

plot(myd, col = cores)

CairoPDF(width = 36, height = 34, file = "YD.pdf", pointsize=40, bg = "white", title = "", paper = "special", pagecentre=TRUE) #
plot(myd, col = cores)
dev.off()
  
library(beepr)
beep(2)
beep(2)

    