### script enm dismo - multiple algorithms ###

# packages
library(colorRamps)
library(dismo)
library(kernlab)
library(randomForest)
library(raster)
library(rJava)
library(tidyverse)
library(viridis)
library(ggplot2)
library(dplyr)
#install.packages("rJava")
###---------------------------------------------------------------------------###
## import data
# 1. occ
# directory
setwd("D:/Unifesp/Mestrado/analises/modelagem/gbif/gbif+artigos/modelagem/")
dir()
getwd()
# occurrences
occ <- readr::read_csv("Vriesea_modelagem_limpa.csv")
occ


#ja foi: 

#"cedrela_modelagem_limpa.csv"         
#"Euphorbia_modelagem_limpa.csv"      
#"Lychnophora_modelagem_limpa.csv"                               
#"Oxalis_modelagem_limpa.csv"         
#"Pitcairnia_modelagem_limpa.csv"      
#"Richterago_modelagem_limpa.csv"
#"Tillandsia_capi_modelagem_limpa.csv" 
#"Tillandsia_vire_modelagem_limpa.csv"
#"Vellozia_modelagem_limpa.csv"        
#"Vriesea_modelagem_limpa.csv"        

# create shapefile
occ.sh <-  occ %>% 
  sf::st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
occ.sh

# plot
ggplot() + 
  geom_sf(data = occ.sh, col = "black", size = 2.5, alpha = .4) +
  theme_minimal()

# export shapefile
sf::st_write(occ.sh, "occ.shp", append = FALSE)

###---------------------------------------------------------------------------###

##  2. variables

setwd("D:/Unifesp/Mestrado/analises/modelagem/raster_novo/CHELSA_cur_V1_2B_r2_5m/")
# list files
pres <- dir(patt = ".tif$")
pres

# import rasters
var <- raster::stack(pres)
var

# plot
plot(var, col = viridis::viridis(100))

#LGM
setwd("D:/Unifesp/Mestrado/analises/modelagem/raster_novo/chelsa_LGM_v1_2B_r2_5m/")
getwd()
lgm <- dir(patt = ".tif$")
lgm

# import rasters
var.lgm <- raster::stack(lgm)
var.lgm

# plot
plot(var.lgm, col = viridis::viridis(100))

#HS
setwd("D:/Unifesp/Mestrado/analises/modelagem/raster_novo/HS1_v1_2_5m/")

hs <- dir(patt = ".tif$")
hs

# import rasters
var.hs <- raster::stack(hs)
var.hs

# plot
plot(var.hs, col = viridis::viridis(100))

#MH
setwd("D:/Unifesp/Mestrado/analises/modelagem/raster_novo/MH_v1_2_5m/")

mh <- dir(patt = ".tif$")
mh

# import rasters
var.mh <- raster::stack(mh)
var.mh

# plot
plot(var.mh, col = viridis::viridis(100))


#YD
setwd("D:/Unifesp/Mestrado/analises/modelagem/raster_novo/YDS_v1_2_5m/")

yd <- dir(patt = ".tif$")
yd

# import rasters
var.yd <- raster::stack(yd)
var.yd

# plot
plot(var.yd, col = viridis::viridis(100))

#LIG
setwd("D:/Unifesp/Mestrado/analises/modelagem/raster_novo/LIG_v1_2_5m/")

lig <- dir(patt = ".tif$")
lig

# import rasters
var.lig <- raster::stack(lig)
var.lig

# plot
plot(var.lig, col = viridis::viridis(100))

names(var) <- c("layer.1", "layer.2","layer.3","layer.4","layer.5","layer.6","layer.7","layer.8")
names(var.lgm) <- c("layer.1", "layer.2","layer.3","layer.4","layer.5","layer.6","layer.7","layer.8")
names(var.hs) <- c("layer.1", "layer.2","layer.3","layer.4","layer.5","layer.6","layer.7","layer.8")
names(var.mh) <- c("layer.1", "layer.2","layer.3","layer.4","layer.5","layer.6","layer.7","layer.8")
names(var.yd) <- c("layer.1", "layer.2","layer.3","layer.4","layer.5","layer.6","layer.7","layer.8")
names(var.lig) <- c("layer.1", "layer.2","layer.3","layer.4","layer.5","layer.6","layer.7","layer.8")

###---------------------------------------------------------------------------###

## extract coordinates for background
# coordinates
## background coordinates
bc <- tibble::as.tibble(raster::rasterToPoints(var)[, 1:2])
bc
colnames(bc) <- c("lon", "lat")
bc

###---------------------------------------------------------------------------###

# verify maxent
# copy maxent.jar in "C:\Users\seu_nome\Documents\R\win-library\3.5.1\dismo\java"
file.exists(paste0(system.file(package = "dismo"), "/java/maxent.jar"))

###---------------------------------------------------------------------------###

### enms ###
# diretory
getwd()
setwd("../")

dir.create("03_vrisea_50")
setwd("03_vrisea_50")


dir.create("02_modelos_multiplos")
setwd("02_modelos_multiplos")
getwd()

# enms

for(i in 1:length(unique(occ[, 1]))){ # for to each specie
  
  # graphics
  dir.create(paste0("graphics_",i))
  
  # variables for evaluate
  eval.Bioclim <- NULL
  #eval.Gower <- NULL
  eval.GLM <- NULL
  eval.RandomForest <- NULL
  eval.Maxent <- NULL
  eval.SVM <- NULL
  eval.names <- NULL
  maxent.results <- matrix()
  unique(occ[1, 1])[[1]]
  # selecting presence and absence
  id.specie <- as.character(unique(occ[1, 1]))[[i]]
  pr.specie <- occ[2:3]
  id.background <- sample(nrow(bc), nrow(pr.specie))
  bc.specie <- bc[id.background, ]
  
  # for
  for(r in 1:50){	# numero de replicas --> o ideal é acima de 10
    
    ## preparing the models
    # train and test data	
    pr.sample.train <- sample(nrow(pr.specie), round(0.7 * nrow(pr.specie)))
    bc.sample.train <- sample(nrow(bc.specie), round(0.7 * nrow(bc.specie)))
    train <- na.omit(dismo::prepareData(x = var, p = pr.specie[pr.sample.train, ], b = bc.specie[bc.sample.train, ]))
    test <- na.omit(dismo::prepareData(x = var, p = pr.specie[-pr.sample.train, ], b = bc.specie[-bc.sample.train, ]))
    
    
    ### algorithms ###
    
    ## 1. bioclim
    print(paste(id.specie, "Bioclim", ifelse(r < 10, paste0("0", r), r)))
    Bioclim <- dismo::bioclim(train[which(train[, 1] == 1), -1])
    raster::writeRaster(dismo::predict(var, Bioclim, progress = "text"), paste0("bioclim_pres_", id.specie, ifelse(r < 10, paste0("0", r), r), ".tif"), options = c("COMPRESS=DEFLATE"), format = "GTiff")
    raster::writeRaster(dismo::predict(var.lgm, Bioclim, progress = "text"), paste0("bioclim_lgm_", id.specie, ifelse(r < 10, paste0("0", r), r), ".tif"), options = c("COMPRESS=DEFLATE"), format = "GTiff")
    raster::writeRaster(dismo::predict(var.hs, Bioclim, progress = "text"), paste0("bioclim_hs_", id.specie, ifelse(r < 10, paste0("0", r), r), ".tif"), options = c("COMPRESS=DEFLATE"), format = "GTiff")
    raster::writeRaster(dismo::predict(var.mh, Bioclim, progress = "text"), paste0("bioclim_mh_", id.specie, ifelse(r < 10, paste0("0", r), r), ".tif"), options = c("COMPRESS=DEFLATE"), format = "GTiff")
    raster::writeRaster(dismo::predict(var.yd, Bioclim, progress = "text"), paste0("bioclim_yd_", id.specie, ifelse(r < 10, paste0("0", r), r), ".tif"), options = c("COMPRESS=DEFLATE"), format = "GTiff")
    raster::writeRaster(dismo::predict(var.lig, Bioclim, progress = "text"), paste0("bioclim_lig_", id.specie, ifelse(r < 10, paste0("0", r), r), ".tif"), options = c("COMPRESS=DEFLATE"), format = "GTiff")
    eBioclim <- dismo::evaluate(p = test[test[, 1] == 1, -1], a = test[test[, 1] == 0, -1], model = Bioclim)
    idBioclim <- which(eBioclim@t == as.numeric(threshold(eBioclim, "spec_sens")))
    eval.Bioclim.sp <- c(id.specie, ifelse(r < 10, paste0("0", r), r), "bioclim", eBioclim@t[idBioclim], eBioclim@auc, (eBioclim@TPR[idBioclim] + eBioclim@TNR[idBioclim] - 1))
    eval.Bioclim <- rbind(eval.Bioclim, eval.Bioclim.sp)
    
    setwd("graphics_1")
    tiff(paste0("bioclim_response_", id.specie, ifelse(r < 10, paste0("0", r), r), ".tif")); dismo::response(Bioclim); dev.off()
    tiff(paste0("bioclim_auc_", id.specie, ifelse(r < 10, paste0("0", r), r), ".tif")); plot(eBioclim, "ROC"); dev.off()
    setwd("..")
    
    
    ## 2. gower
    #print(paste(id.specie, "Gower", ifelse(r < 10, paste0("0", r), r)))
    #Gower <- dismo::domain(train[which(train[, 1] == 1), -1])	
    #raster::writeRaster(dismo::predict(var, Gower, progress = "text"), paste0("gower_pres_", id.specie, ifelse(r < 10, paste0("0", r), r), ".tif"), options = c("COMPRESS=DEFLATE"), format = "GTiff") 
    #raster::writeRaster(dismo::predict(var.f, Gower, progress = "text"), paste0("gower_future_", id.specie, ifelse(r < 10, paste0("0", r), r), ".tif"), options = c("COMPRESS=DEFLATE"), format = "GTiff") 
    #eGower <- dismo::evaluate(p = test[test[, 1] == 1, -1], a = test[test[, 1] == 0, -1], model = Gower)
    #idGower <- which(eGower@t == as.numeric(threshold(eGower, "spec_sens")))
    #eval.Gower.sp <- c(id.specie, ifelse(r < 10, paste0("0", r), r), "gower", eGower@t[idGower], eGower@auc, (eGower@TPR[idGower] + eGower@TNR[idGower] - 1))
    #eval.Gower <- rbind(eval.Gower, eval.Gower.sp)
    
    #setwd("graphics")
    #tiff(paste0("gower_response_", id.specie, ifelse(r < 10, paste0("0", r), r), ".tif")); dismo::response(Gower); dev.off()
    #tiff(paste0("gower_auc_", id.specie, ifelse(r < 10, paste0("0", r), r), ".tif")); plot(eGower, "ROC"); dev.off()
    #setwd("..")
    
  
    ## 3. glm
    print(paste(id.specie, "GLM", ifelse(r < 10, paste0("0", r), r)))
    GLM <- glm(pb ~ ., data = train)	
    raster::writeRaster(dismo::predict(var, GLM, progress = "text"), paste0("glm_pres_", id.specie, ifelse(r < 10, paste0("0", r), r), ".tif"), options = c("COMPRESS=DEFLATE"), format = "GTiff") 
    raster::writeRaster(dismo::predict(var.lgm, GLM, progress = "text"), paste0("glm_lgm_", id.specie, ifelse(r < 10, paste0("0", r), r), ".tif"), options = c("COMPRESS=DEFLATE"), format = "GTiff")
    raster::writeRaster(dismo::predict(var.hs, GLM, progress = "text"), paste0("glm_hs_", id.specie, ifelse(r < 10, paste0("0", r), r), ".tif"), options = c("COMPRESS=DEFLATE"), format = "GTiff")
    raster::writeRaster(dismo::predict(var.mh, GLM, progress = "text"), paste0("glm_mh_", id.specie, ifelse(r < 10, paste0("0", r), r), ".tif"), options = c("COMPRESS=DEFLATE"), format = "GTiff")
    raster::writeRaster(dismo::predict(var.yd, GLM, progress = "text"), paste0("glm_yd_", id.specie, ifelse(r < 10, paste0("0", r), r), ".tif"), options = c("COMPRESS=DEFLATE"), format = "GTiff")
    raster::writeRaster(dismo::predict(var.lig, GLM, progress = "text"), paste0("glm_lig_", id.specie, ifelse(r < 10, paste0("0", r), r), ".tif"), options = c("COMPRESS=DEFLATE"), format = "GTiff")
    eGLM <- dismo::evaluate(p = test[test[, 1] == 1, -1], a = test[test[, 1] == 0, -1], model = GLM)
    idGLM <- which(eGLM@t == as.numeric(threshold(eGLM, "spec_sens")))
    eval.GLM.sp <- c(id.specie, ifelse(r < 10, paste0("0", r), r), "glm", eGLM@t[idGLM], eGLM@auc, (eGLM@TPR[idGLM] + eGLM@TNR[idGLM] - 1))
    eval.GLM <- rbind(eval.GLM, eval.GLM.sp)
    
    setwd("graphics_1")
    tiff(paste0("glm_response_", id.specie, ifelse(r < 10, paste0("0", r), r), ".tif")); dismo::response(GLM); dev.off()
    tiff(paste0("glm_auc_", id.specie, ifelse(r < 10, paste0("0", r), r), ".tif")); plot(eGLM, "ROC"); dev.off()
    setwd("..")
    
    
    ## 4. random forest
    print(paste(id.specie, "Random Forest", ifelse(r < 10, paste0("0", r), r)))
    RandomForest <- randomForest::randomForest(pb ~ ., data = train)
    raster::writeRaster(dismo::predict(var, RandomForest, progress = "text"), paste0("randomforest_pres_", id.specie, ifelse(r < 10, paste0("0", r), r), ".tif"), options = c("COMPRESS=DEFLATE"), format = "GTiff") 
    raster::writeRaster(dismo::predict(var.lgm, RandomForest, progress = "text"), paste0("randomforest_lgm_", id.specie, ifelse(r < 10, paste0("0", r), r), ".tif"), options = c("COMPRESS=DEFLATE"), format = "GTiff")
    raster::writeRaster(dismo::predict(var.hs, RandomForest, progress = "text"), paste0("randomforest_hs_", id.specie, ifelse(r < 10, paste0("0", r), r), ".tif"), options = c("COMPRESS=DEFLATE"), format = "GTiff")
    raster::writeRaster(dismo::predict(var.mh, RandomForest, progress = "text"), paste0("randomforest_mh_", id.specie, ifelse(r < 10, paste0("0", r), r), ".tif"), options = c("COMPRESS=DEFLATE"), format = "GTiff")
    raster::writeRaster(dismo::predict(var.yd, RandomForest, progress = "text"), paste0("randomforest_yd_", id.specie, ifelse(r < 10, paste0("0", r), r), ".tif"), options = c("COMPRESS=DEFLATE"), format = "GTiff")
    raster::writeRaster(dismo::predict(var.lig, RandomForest, progress = "text"), paste0("randomforest_lig_", id.specie, ifelse(r < 10, paste0("0", r), r), ".tif"), options = c("COMPRESS=DEFLATE"), format = "GTiff")
    eRandomForest <- dismo::evaluate(p = test[test[, 1] == 1, -1], a = test[test[, 1] == 0, -1], model = RandomForest)
    idRandomForest <- which(eRandomForest@t == as.numeric(threshold(eRandomForest, "spec_sens")))
    eval.RandomForest.sp <- c(id.specie, ifelse(r < 10, paste0("0", r), r), "randomforest", eRandomForest@t[idRandomForest], eRandomForest@auc, (eRandomForest@TPR[idRandomForest] + eRandomForest@TNR[idRandomForest] - 1))
    eval.RandomForest <- rbind(eval.RandomForest, eval.RandomForest.sp)
    
    setwd("graphics_1")
    tiff(paste0("randomforest_auc_", id.specie, ifelse(r < 10, paste0("0", r), r), ".tif")); plot(eRandomForest, "ROC"); dev.off()
    setwd("..")
    
    ## 5. maxent	
    print(paste(id.specie, "Maxent", ifelse(r < 10, paste0("0", r), r)))
    Maxent <- dismo::maxent(train[, -1], train[, 1])	
    raster::writeRaster(dismo::predict(var, Maxent, progress = "text"), paste0("maxent_pres_", id.specie, ifelse(r < 10, paste0("0", r), r), ".tif"), options = c("COMPRESS=DEFLATE"), format = "GTiff") 
    raster::writeRaster(dismo::predict(var.lgm, Maxent, progress = "text"), paste0("maxent_lgm_", id.specie, ifelse(r < 10, paste0("0", r), r), ".tif"), options = c("COMPRESS=DEFLATE"), format = "GTiff")
    raster::writeRaster(dismo::predict(var.hs, Maxent, progress = "text"), paste0("maxent_hs_", id.specie, ifelse(r < 10, paste0("0", r), r), ".tif"), options = c("COMPRESS=DEFLATE"), format = "GTiff")
    raster::writeRaster(dismo::predict(var.mh, Maxent, progress = "text"), paste0("maxent_mh_", id.specie, ifelse(r < 10, paste0("0", r), r), ".tif"), options = c("COMPRESS=DEFLATE"), format = "GTiff")
    raster::writeRaster(dismo::predict(var.yd, Maxent, progress = "text"), paste0("maxent_yd_", id.specie, ifelse(r < 10, paste0("0", r), r), ".tif"), options = c("COMPRESS=DEFLATE"), format = "GTiff")
    raster::writeRaster(dismo::predict(var.lig, Maxent, progress = "text"), paste0("maxent_lig_", id.specie, ifelse(r < 10, paste0("0", r), r), ".tif"), options = c("COMPRESS=DEFLATE"), format = "GTiff")
    eMaxent <- dismo::evaluate(p = test[test[, 1] == 1, -1], a = test[test[, 1] == 0, -1], model = Maxent)
    idMaxent <- which(eMaxent@t == as.numeric(threshold(eMaxent, "spec_sens")))
    eval.Maxent.sp <- c(id.specie, ifelse(r < 10, paste0("0", r), r), "maxent", eMaxent@t[idMaxent], eMaxent@auc, (eMaxent@TPR[idMaxent] + eMaxent@TNR[idMaxent] - 1))
    eval.Maxent <- rbind(eval.Maxent, eval.Maxent.sp)
    
    setwd("graphics_1")
    tiff(paste0("maxent_response_",  id.specie, ifelse(r < 10, paste0("0", r), r), ".tif")); dismo::response(Maxent); dev.off()
    tiff(paste0("maxent_contribution_",  id.specie, ifelse(r < 10, paste0("0", r), r), ".tif")); plot(Maxent); dev.off()
    tiff(paste0("maxent_auc_",  id.specie, ifelse(r < 10, paste0("0", r), r), ".tif")); plot(eMaxent, "ROC"); dev.off()
    maxent.results <- tibble::as.tibble(data.frame(maxent.results, as.matrix(Maxent@results)))
    setwd("..")
    
    ## 6. svm	
    print(paste(id.specie, "SVM", ifelse(r < 10, paste0("0", r), r)))
    SVM <- kernlab::ksvm(pb ~ ., data = train)
    raster::writeRaster(dismo::predict(var, SVM, progress = "text"), paste0("svm_pres_", id.specie, ifelse(r < 10, paste0("0", r), r), ".tif"), options = c("COMPRESS=DEFLATE"), format = "GTiff") 
    raster::writeRaster(dismo::predict(var.lgm, SVM, progress = "text"), paste0("svm_lgm_", id.specie, ifelse(r < 10, paste0("0", r), r), ".tif"), options = c("COMPRESS=DEFLATE"), format = "GTiff")
    raster::writeRaster(dismo::predict(var.hs, SVM, progress = "text"), paste0("svm_hs_", id.specie, ifelse(r < 10, paste0("0", r), r), ".tif"), options = c("COMPRESS=DEFLATE"), format = "GTiff")
    raster::writeRaster(dismo::predict(var.mh, SVM, progress = "text"), paste0("svm_mh_", id.specie, ifelse(r < 10, paste0("0", r), r), ".tif"), options = c("COMPRESS=DEFLATE"), format = "GTiff")
    raster::writeRaster(dismo::predict(var.yd, SVM, progress = "text"), paste0("svm_yd_", id.specie, ifelse(r < 10, paste0("0", r), r), ".tif"), options = c("COMPRESS=DEFLATE"), format = "GTiff")
    raster::writeRaster(dismo::predict(var.lig, SVM, progress = "text"), paste0("svm_lig_", id.specie, ifelse(r < 10, paste0("0", r), r), ".tif"), options = c("COMPRESS=DEFLATE"), format = "GTiff")
    eSVM <- dismo::evaluate(p = test[test[, 1] == 1, -1], a = test[test[, 1] == 0, -1], model = SVM)
    idSVM <- which(eSVM@t == as.numeric(threshold(eSVM, "spec_sens")))
    eval.SVM.sp <- c(id.specie, ifelse(r < 10, paste0("0", r), r), "svm", eSVM@t[idSVM], eSVM@auc, (eSVM@TPR[idSVM] + eSVM@TNR[idSVM] - 1))
    eval.SVM <- rbind(eval.SVM, eval.SVM.sp)
    
    setwd("graphics_1")
    tiff(paste0("svm_auc_", id.specie, ifelse(r < 10, paste0("0", r), r), ".tif")); plot(eSVM, "ROC"); dev.off()
    setwd("..")
    
    eval.names <- c(eval.names, paste0(id.specie, ifelse(r < 10, paste0("0", r), r)))	
    
  } # ends for "r"
  
  # maxent results
  setwd("graphics_1")
  na <- attributes(Maxent@results)[[2]][[1]]
  maxent.results <- tibble::as.tibble(data.frame(na, maxent.results[, -1]))
  colnames(maxent.results) <- c("names", paste0("rep", 1:r))
  readr::write_csv(maxent.results, paste0("_maxent_results", id.specie, ".csv", overwrite=TRUE))
  setwd("..")
  
  # evaluations
  dimnames(eval.Bioclim) <- list(eval.names, c("species", "replica", "algorithm", "thrs", "AUC", "TSS"))
  #dimnames(eval.Gower) <- list(eval.names, c("species", "replica", "algorithm", "thrs", "AUC", "TSS"))  
  dimnames(eval.GLM) <- list(eval.names, c("species", "replica", "algorithm", "thrs", "AUC", "TSS"))  
  dimnames(eval.RandomForest) <- list(eval.names, c("species", "replica", "algorithm", "thrs", "AUC", "TSS"))
  dimnames(eval.Maxent) <- list(eval.names, c("species", "replica", "algorithm", "thrs", "AUC", "TSS"))
  dimnames(eval.SVM) <- list(eval.names, c("species", "replica", "algorithm", "thrs", "AUC", "TSS"))
  
  write.csv(eval.Bioclim, paste0("zEval_", "bioclim_", id.specie, ".csv"))
  #write.csv(eval.Gower, paste0("zEval_", "gower_", id.specie, ".csv"))
  write.csv(eval.GLM, paste0("zEval_", "glm_", id.specie, ".csv"))
  write.csv(eval.RandomForest, paste0("zEval_", "randomforest_", id.specie, ".csv"))
  write.csv(eval.Maxent, paste0("zEval_", "maxent_", id.specie, ".csv"))
  write.csv(eval.SVM, paste0("zEval_", "svm_", id.specie, ".csv"))
  
} # ends for"i"

###----------------------------------------------------------------------------###

##### RODAR ENSEMBLE #####

