############ PCA OFICIAL - VERSÃO ATUALIZADA ############

library(terra)
library(sf)
library(proxy)     
library(factoextra) 

# ========================================
# 1. Configurar diretório e verificar arquivos - mudar conforme categoria
# ========================================
setwd("D:/Unifesp/Mestrado/analises/SEM/vellozia/")

dir.create("land")
dir.create("precip")
dir.create("solo")
dir.create("topo")
#dir.create("radia")


dir()  # Verificar se os arquivos esperados estão presentes

# ========================================
# 2. Carregar dados de rasters - 5 categorias (LCM) - land, topo, solo, precip, radia
# ========================================

wc_files <- list.files(
  path = "D:/Unifesp/Mestrado/analises/rasters_ib/SEM/solo/sa_solo",
  pattern = ".tif$",
  full.names = TRUE
)
wc <- rast(wc_files)
print(wc)  # Verificar metadados

###caminhos
#D:/Unifesp/Mestrado/analises/rasters_ib/SEM/landcover/sa_land
#D:/Unifesp/Mestrado/analises/rasters_ib/SEM/temp_precip/sa_temp/wc
#D:/Unifesp/Mestrado/analises/rasters_ib/SEM/solo/sa_solo/
#D:/Unifesp/Mestrado/analises/rasters_ib/SEM/topogra/sa_topo

#D:/Unifesp/Mestrado/analises/rasters_ib/SEM/radia/sa_radia

# ========================================
# 3. Carregar coordenadas das amostras
# ========================================
coords <- read.table(
  "input_coord_sem_vellozia_auriculata_new.txt",
  header = FALSE,
  sep = "",
  col.names = c("ID", "longitude", "latitude")
)

# Converter para objeto espacial (sf)
points_sf <- st_as_sf(coords, coords = c("longitude", "latitude"), crs = 4326)


# Plot de verificação
plot(wc[[1]], main = "Verificação de sobreposição")
plot(st_geometry(points_sf), add = TRUE, col = "red", pch = 19)

# ========================================
# 4. Extrair e padronizar variáveis ambientais
# ========================================
env_values <- terra::extract(wc, vect(points_sf))[, -1]  # Remove coluna ID

# Verificar e remover colunas com variância zero antes de padronizar
vars_removidas <- colnames(env_values)[apply(env_values, 2, var, na.rm = TRUE) == 0]
if (length(vars_removidas) > 0) {
  message("Variáveis removidas (variância zero): ", paste(vars_removidas, collapse = ", "))
  env_values <- env_values[, !colnames(env_values) %in% vars_removidas]
}

env_scaled <- scale(env_values)  # Padronização (média = 0, DP = 1)

# Checar NAs
if (any(is.na(env_scaled))) {
  warning("Dados contêm NAs! Verificar antes de calcular distâncias.")
  env_scaled <- na.omit(env_scaled)  # Remove NAs ou trate conforme necessário
}

# ========================================
# 5. Matriz de distância ambiental (Euclidiana)
# ========================================
dist_env <- dist(env_scaled, method = "euclidean")
dist_env_matrix <- as.matrix(dist_env)

# Visualização parcial
print(round(dist_env_matrix[1:5, 1:5], 2))

# Salvar matriz de distância
write.csv(dist_env_matrix, "matriz_distancia_ambiental_.csv")

# ========================================
# 6. PCA e distância nos componentes principais
# ========================================
pca_env <- prcomp(env_scaled, scale. = FALSE)

# Visualização dos resultados da PCA
# 6.1. Variância explicada
fviz_eig(pca_env, addlabels = TRUE, 
         main = "Variância explicada pelos PCs",
         barfill = "#4E84C4") +
  theme_minimal()

ggsave("variancia_explicada_PCA_.png", device = "png", width = 8, height = 6)

# 6.2. Biplot para interpretação das variáveis
fviz_pca_biplot(pca_env, 
                repel = TRUE,
                col.var = "red",
                col.ind = "#696969",
                title = "Biplot - Variáveis Ambientais",
                addEllipses = FALSE)


# 6.3. Scores dos indivíduos nos PCs
fviz_pca_ind(pca_env,
             col.ind = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)

# Extrair e salvar os scores da PCA
env_pca_scores <- pca_env$x[, 1:2]  # Primeiros dois componentes
write.csv(env_pca_scores, "sem_PCA_scores_env_.csv")

# 6.4. Matriz de distância baseada nos PCs
dist_pca <- dist(env_pca_scores, method = "euclidean")
dist_pca_matrix <- as.matrix(dist_pca)

# Salvar matriz de distância da PCA
write.csv(dist_pca_matrix, "matriz_distancia_PCA_.csv")

# ========================================
# 7. Informações adicionais da PCA (opcional)
# ========================================
# 7.1. Autovalores
eigenvals <- pca_env$sdev^2
write.csv(eigenvals, "autovalores_PCA_.csv")

# 7.2. Cargas fatoriais (loadings)
str(pca_env)

loadings <- pca_env$rotation[, 1, drop = FALSE]

#loadings <- pca_env$rotation[, 1:2]  # Loadings dos dois primeiros PCs
#loadings <- pca_env$rotation[, 1:min(2, ncol(pca_env$rotation))] #pega de todos os pcs disponiveis
write.csv(loadings, "loadings_PCA_.csv")

contrib_PC1 <- (loadings[, 1]^2) * 100
contrib_PC1_df <- data.frame(Variavel = rownames(loadings), Contribuicao_PC1 = contrib_PC1)
write.csv(contrib_PC1_df, "contribuicao_variaveis_PC1_.csv")

fviz_pca_var(pca_env,
             col.var = "contrib", # Colorir pela contribuição
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,
             title = "Contribuição das variáveis para os PCs")


# 7.3. Resumo completo da PCA
summary(pca_env)
# Salvar o sumário da PCA em um arquivo de texto
capture.output(
  summary(pca_env), 
  file = "resumo_completo_PCA_.txt"
)

##################### vai ser alterado no script sem, pra pairwise##########

### repetir para cada agrupamento de variaveis definidas, e colocar na pasta criada anteriormente lan/topo/precip/radia/solo #####



