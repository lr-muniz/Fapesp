# Carregar pacotes
library(usdm)
library(raster)


# directory
setwd("D:/Unifesp/Mestrado/analises/rasters_ib/SEM/temp_precip/sa_temp_extent/wc/")

dir()
###---------------------------------------------------------------------------------------###

# list files
ti1 <- dir(patt = ".tif")

ti1
#present

var1 <- raster::stack(ti1)

# vif 10- fizemos para th=10 e th=5. LCM - optamos pelo 5 pois ele é mais adequado para objetivos de publicacao (MUDAR O NOME NA HORA DE SALVAR E NAO CONFUNDIR VALOR DE TH) 
vi.10 <- usdm::vifstep(var1, th = 5)
vi.10
vi.10@results


#save TH = 5
write.table(
  vi.10@results, 
  file = "resultados_5_vif.txt", 
  sep = "\t", 
  row.names = FALSE, 
  quote = FALSE
)


vi.10 <- usdm::vifstep(var1, th = 10)
vi.10
vi.10@results


#save TH = 10
#write.table(
#  vi.10@results, 
#  file = "resultados_10_vif.txt", 
#  sep = "\t", 
#  row.names = FALSE, 
#  quote = FALSE
#  )
