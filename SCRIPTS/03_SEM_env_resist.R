############# SEM com scores de PCA por categoria ################

# 1. Pacotes

#if(!require(lavaan)) install.packages("lavaan")
#if(!require(tidySEM)) install.packages("tidySEM")
#if(!require(dplyr)) install.packages("dplyr")
#if(!require(lavaanPlot)) install.packages("lavaanPlot")

library(lavaan)
library(tidySEM)
library(dplyr)
library(lavaanPlot)

# 2. Diretório
# Este deve ser o diretório base onde os arquivos de scores (e outros inputs) estão.
setwd("D:/Unifesp/Mestrado/analises/SEM/pitcairnia/") 


# 3. Dados genéticos - PCoA (para dist_fst_sem)
fst_matrix <- read.table("input_fst_sem_pitcairnia_lanuginosa.out", header = TRUE)[,-1]
dist_fst <- as.dist(as.matrix(fst_matrix))
dist_fst_sem <- as.numeric(dist_fst) # <-- Esta é a variável de 91 linhas para o SEM

# 4. Coordenadas geográficas (para geo_dist)
coords <- read.table("input_5_coord_cs_pitcairnia_lanuginosa.txt", 
                     header = FALSE, 
                     col.names = c("ID", "longitude", "latitude"))
dist_coords <- dist(coords[, c("longitude", "latitude")]) 
geo_dist <- as.numeric(dist_coords) # <-- Esta é a variável de 91 linhas para o SEM

# 5. Carregar scores PCA ambientais (sem radiação) E CONVERTER PARA DISTÂNCIAS ENTRE PARES

### LCM - ATENCAO!! 
### se por um acaso foi retirada alguma categoria, adequar os grupos_env ou grupos_resist

env_PCAs_paired <- list()
grupos_env <- c("land", "topo", "solo", "clima")

for (g in grupos_env) {
  file_path <- paste0("sem_PCA_scores_env_", g, ".csv") 
  
  # Lendo os scores da PCA ambiental 
  # Confirme se a coluna é 'PC1' ou 'Dim.1'
  env_scores_linhas <- read.csv(file_path, row.names = 1)$PC1 # Seus scores de PCA ambiental
  
  # CONVERTER PARA DISTÂNCIA EUCLIDIANA ENTRE OS SCORES DOS PARES 
  # A distância é calculada entre os valores de PC1 de cada par de populações
  dist_env_pc1 <- dist(as.matrix(env_scores_linhas), method = "euclidean")
  env_PCAs_paired[[g]] <- scale(as.numeric(dist_env_pc1)) # Escalar a distância
}

# 6. Carregar scores PCA de resistência (sem radiação) 
grupos_resist <- c("land", "topo", "solo", "clima")
resist_PCAs_paired <- list()
for (g in grupos_resist) { 
  file_path <- paste0("sem_PCA_scores_resist_", g, ".csv") 
  
  # Estes arquivos já contêm os scores da PCA para os 91 pares (como Dim.1 ou PC1)
  # Com base em image_3b085a.png, o nome da coluna é "Dim.1"
  resist_PCAs_paired[[g]] <- scale(read.csv(file_path, row.names = 1)$Dim.1) 
}

# 7. Construir o data.frame para o SEM
# AGORA TODAS AS VARIÁVEIS DEVEM TER 91 LINHAS
sem_data <- data.frame(
  genetic_dist = dist_fst_sem, # <-- CORRIGIDO: use o nome da variável correta
  geo_dist = geo_dist,         
  env_land = env_PCAs_paired$land, 
  env_topo = env_PCAs_paired$topo, 
  env_solo = env_PCAs_paired$solo, 
  env_clima = env_PCAs_paired$clima, 
  res_land = resist_PCAs_paired$land,
  res_topo = resist_PCAs_paired$topo,
  res_solo = resist_PCAs_paired$solo,
  res_clima = resist_PCAs_paired$clima
)


# 8. Definir modelo SEM em formato relacional
model <- '
  # Efeitos sobre distância genética
  genetic_dist ~ env_land + 
                 env_topo + 
                 env_solo + 
                 env_clima +
                 res_land + 
                 res_topo + 
                 res_solo + 
                 res_clima +
                 geo_dist # Distância geográfica como preditor

  # Efeitos da distância geográfica sobre ambiente (PCs de pares de pontos)
  env_land ~ geo_dist
  env_topo ~ geo_dist
  env_solo ~ geo_dist
  env_clima ~ geo_dist

  # Efeitos da distância geográfica sobre resistência (PCs de pares de pontos)
  res_land ~ geo_dist
  res_topo ~ geo_dist
  res_solo ~ geo_dist
  res_clima ~ geo_dist
'
#####################################################
# --- CHECAR E REMOVER VARIÁVEIS ALTAMENTE CORRELACIONADAS ---
# Defina o limiar de correlação (ex: 0.95 ou 0.99)
limiar_corr <- 0.97

# Calcular matriz de correlação
cor_matrix <- cor(sem_data, use = "pairwise.complete.obs")

# Função para identificar colunas redundantes
find_redundant_vars <- function(cor_mat, limiar) {
  remove_vars <- c()
  for (i in 1:(ncol(cor_mat) - 1)) {
    for (j in (i + 1):ncol(cor_mat)) {
      if (abs(cor_mat[i, j]) > limiar) {
        var_remove <- colnames(cor_mat)[j]  # remove a segunda para manter consistência
        remove_vars <- c(remove_vars, var_remove)
      }
    }
  }
  unique(remove_vars)
}

# Identificar variáveis redundantes
redundantes <- find_redundant_vars(cor_matrix, limiar_corr)

# Remover se existir
if (length(redundantes) > 0) {
  message("Variáveis removidas por alta correlação (>", limiar_corr, "): ", 
          paste(redundantes, collapse = ", "))
  sem_data <- sem_data[, !(colnames(sem_data) %in% redundantes)]
} else {
  message("Nenhuma variável com correlação acima de ", limiar_corr)
}

# --- Fim da checagem ---

#####################################################
# 9. Ajustar modelo
fit <- sem(model, data = sem_data)

# 10. Resultados
summary(fit, fit.measures = TRUE, standardized = TRUE)

# 11. Salvar coeficientes e sumário (opcional)
write.csv(standardizedSolution(fit), "SEM_coefficients_results.csv")
capture.output(summary(fit, fit.measures = TRUE, standardized = TRUE), 
               file = "SEM_full_results.txt")

# 12. Visualização com lavaanPlot
pdf("SEM_sem_radiacao.pdf", width = 12, height = 10)
lavaanPlot(
  model = fit,
  coef = TRUE,
  stand = TRUE,
  graph_options = list(rankdir = "TB", fontsize = 10),
  edge_options = list(color = "darkblue"),
  node_options = list(shape = "box", fontname = "Helvetica")
)
dev.off()


####SALVAR DIRETO DO VIEWER AO LADO####
## EXPORT > SAVE AS IMAGE > NOMEAR > ESCOLHER DIRETORIO CORRETO


#library(ggplot2)
# pdf("SEM_sem_radiacao.pdf", width = 20, height = 15) # Dimensões em polegadas para PDF
# tryCatch({
#   lavaanPlot(
#     model = fit,
#     coef = TRUE,
#     stand = TRUE,
#     graph_options = list(rankdir = "TB", fontsize = 10),
#     edge_options = list(color = "darkblue"),
#     node_options = list(shape = "box", fontname = "Helvetica")
#   )
#   Sys.sleep(1) 
# }, finally = {
#   dev.off()
# })

