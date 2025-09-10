####### PCA RESIST GEMINI #######

# --- 0. Instalar e Carregar Pacotes (se necessário) ---
# Se você ainda não tem esses pacotes, remova o '#' da linha de instalação e execute.
# install.packages("tidyverse")
# install.packages("factoextra")
# install.packages("FactoMineR")

library(tidyverse)
library(factoextra)
library(FactoMineR)

# --- 1. Carregar Suas Matrizes de Resistência do Circuitscape Automaticamente ---

# Define o diretório onde seus arquivos de resistência estão localizados.
# SUBSTITUA ESTE CAMINHO PELO DIRETÓRIO REAL ONDE OS ARQUIVOS DA FOTO ESTÃO.
diretorio_arquivos <- "D:/Unifesp/Mestrado/analises/SEM/euphorbia/RESIST/topo/"


## LCM- mudar a pasta de acordo com o grupo de variaveis que fiz

# OU, se os arquivos estiverem na mesma pasta do seu script R:
# diretorio_arquivos <- getwd()

# Define o diretório onde os resultados da PCA serão salvos.
# É uma boa prática criar uma pasta separada para os outputs.
output_dir <- file.path(diretorio_arquivos, "PCA_Outputs")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  message(paste0("Diretório de saída criado: ", output_dir))
}


# Lista todos os arquivos no diretório que terminam com ".out"
caminhos_arquivos <- list.files(path = diretorio_arquivos, 
                                pattern = "\\.out$", 
                                full.names = TRUE)

# Verifica se foram encontrados arquivos
if (length(caminhos_arquivos) == 0) {
  stop("Nenhum arquivo .out encontrado no diretório especificado. Verifique o caminho e o padrão.")
}

# Gera os nomes para as suas colunas de forma automática a partir dos nomes dos arquivos.
nome_das_matrizes <- basename(caminhos_arquivos) %>%
  str_remove_all("\\.out$") # Remove a extensão ".out"

cat("Arquivos encontrados e nomes das colunas gerados:\n")
print(caminhos_arquivos)
print(nome_das_matrizes)
cat("\n")

# Lista para armazenar as matrizes carregadas do Circuitscape
todas_as_matrizes_carregadas <- list()

# O número de populações será inferido do primeiro arquivo.
num_populacoes <- NULL 
pop_names <- NULL

# Loop para carregar cada arquivo e pré-processá-lo
for (i in seq_along(caminhos_arquivos)) {
  filepath <- caminhos_arquivos[i]
  
  matriz_df <- read.table(filepath, 
                          header = TRUE, 
                          sep = "",      
                          row.names = 1) 
  
  matriz_num <- as.matrix(matriz_df)
  
  # Verificar e definir num_populacoes e pop_names na primeira iteração
  if (is.null(num_populacoes)) {
    num_populacoes <- nrow(matriz_num)
    pop_names <- rownames(matriz_num)
    cat(paste0("Número de Populações/Pontos detectado: ", num_populacoes, "\n"))
    cat("Nomes dos Pontos/Populações:\n")
    print(pop_names)
    cat("\n")
  } else {
    if (nrow(matriz_num) != num_populacoes || !all(rownames(matriz_num) == pop_names)) {
      stop(paste0("Inconsistência nas dimensões ou nomes de pontos entre o arquivo: ", 
                  basename(filepath), " e o primeiro arquivo lido. Todas as matrizes devem ser idênticas em estrutura."))
    }
  }
  
  todas_as_matrizes_carregadas[[i]] <- matriz_num
  
  if (i == 1) {
    cat("Exemplo da Primeira Matriz de Resistência Carregada (primeiras linhas/colunas):\n")
    print(head(matriz_num))
    cat("\n")
  }
}

# --- 2. Extrair Valores Únicos e Concatenar para Matriz Unificada para PCA ---

lista_df_valores_resistencia <- list()

for (i in seq_along(todas_as_matrizes_carregadas)) {
  matriz_atual <- todas_as_matrizes_carregadas[[i]]
  
  valores_matriz_atual <- as.numeric(as.dist(matriz_atual))
  
  pares_pontos <- combn(pop_names, 2, FUN = function(x) paste(x[1], x[2], sep = "-"))
  
  df_temp <- tibble(
    `Par_de_Pontos` = pares_pontos,
    !!sym(nome_das_matrizes[i]) := valores_matriz_atual 
  )
  lista_df_valores_resistencia[[i]] <- df_temp
}

matriz_pca_input <- reduce(lista_df_valores_resistencia, full_join, by = "Par_de_Pontos") %>%
  column_to_rownames(var = "Par_de_Pontos") 

cat("\nMatriz de Entrada para a PCA (antes da remoção de variáveis com variância zero - primeiras 5 linhas):\n")
print(head(matriz_pca_input))
cat("\nDimensões da Matriz de Entrada para PCA (antes da remoção):", dim(matriz_pca_input), "\n")

# --- NOVA ETAPA: Remover colunas com variância zero ---
# A variância zero indica que a variável não tem nenhuma variação e, portanto, não contribuirá para a PCA.
vars_com_variancia_zero <- colnames(matriz_pca_input)[apply(matriz_pca_input, 2, var, na.rm = TRUE) == 0]

if (length(vars_com_variancia_zero) > 0) {
  message("Variáveis de resistência removidas (variância zero): ", paste(vars_com_variancia_zero, collapse = ", "))
  matriz_pca_input <- matriz_pca_input[, !colnames(matriz_pca_input) %in% vars_com_variancia_zero]
  
  # Atualiza os nomes das matrizes, caso algum tenha sido removido
  nome_das_matrizes <- colnames(matriz_pca_input) 
  
  cat("\nMatriz de Entrada para a PCA (APÓS remoção de variáveis com variância zero - primeiras 5 linhas):\n")
  print(head(matriz_pca_input))
  cat("\nDimensões da Matriz de Entrada para PCA (APÓS remoção):", dim(matriz_pca_input), "\n")
  
  if (ncol(matriz_pca_input) == 0) {
    stop("Todas as variáveis de resistência foram removidas devido a variância zero. PCA não pode ser realizada.")
  }
} else {
  message("Nenhuma variável de resistência com variância zero encontrada. Nenhuma variável foi removida.")
}


# --- 3. Realizando a PCA ---
# A função PCA() do FactoMineR já faz a padronização (scaling) por padrão.
res.pca <- PCA(matriz_pca_input, scale.unit = TRUE, ncp = ncol(matriz_pca_input), graph = FALSE)

# --- 4. Análise dos Resultados da PCA ---

# 4.1. Variância Explicada
cat("\nVariância Explicada por cada Componente Principal:\n")
print(res.pca$eig)

# Salvar a tabela de variância explicada por cada PC
# Isso inclui autovalores, porcentagem de variância e porcentagem cumulativa.
write.csv(res.pca$eig, file.path(output_dir, "pca_variancia_explicada.csv"), row.names = TRUE)
message(paste0("Tabela de variância explicada salva em: ", file.path(output_dir, "pca_variancia_explicada.csv")))


# Scree Plot
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50)) +
  labs(title = "Scree Plot: Variância Explicada por Componente Principal") +
  theme_minimal()

# Salvar o Scree Plot
ggsave(file.path(output_dir, "scree_plot.png"), device = "png", width = 8, height = 6)
message(paste0("Scree Plot salvo em: ", file.path(output_dir, "scree_plot.png")))


# 4.2. Contribuição de Cada Matriz de Resistência (Variável) para os PCs
cat("\nContribuição de cada Matriz de Resistência para os Componentes Principais (%):\n")
contribuicoes_var <- get_pca_var(res.pca)$contrib 
print(round(contribuicoes_var, 2))

# Salvar a tabela de contribuições das variáveis para os PCs
write.csv(contribuicoes_var, file.path(output_dir, "pca_contribuicoes_variaveis.csv"), row.names = TRUE)
message(paste0("Tabela de contribuições das variáveis salva em: ", file.path(output_dir, "pca_contribuicoes_variaveis.csv")))


# Visualizar a contribuição das variáveis para o PC1
fviz_contrib(res.pca, choice = "var", axes = 1, top = ncol(matriz_pca_input),
             title = "Contribuição das Matrizes de Resistência para o PC1") +
  theme_minimal()
ggsave(file.path(output_dir, "contribuicao_pc1_plot.png"), device = "png", width = 10, height = 7)


# Visualizar a contribuição das variáveis para o PC2
fviz_contrib(res.pca, choice = "var", axes = 2, top = ncol(matriz_pca_input),
             title = "Contribuição das Matrizes de Resistência para o PC2") +
  theme_minimal()
ggsave(file.path(output_dir, "contribuicao_pc2_plot.png"), device = "png", width = 10, height = 7)


# 4.3. Cargas (Loadings) das Variáveis
cat("\nCargas (Loadings) das Variáveis em cada Componente Principal:\n")
cargas_var <- get_pca_var(res.pca)$coord 
print(round(cargas_var, 3))

# Salvar a tabela de cargas (loadings) das variáveis
write.csv(cargas_var, file.path(output_dir, "pca_cargas_variaveis.csv"), row.names = TRUE)
message(paste0("Tabela de cargas das variáveis salva em: ", file.path(output_dir, "pca_cargas_variaveis.csv")))


# Visualizar as variáveis em um biplot (PCs 1 e 2)
fviz_pca_var(res.pca, col.var = "black", repel = TRUE,
             title = "Biplot das Variáveis (Matrizes de Resistência)") +
  theme_minimal()
ggsave(file.path(output_dir, "biplot_variaveis.png"), device = "png", width = 8, height = 8)


# 4.4. Coordenadas dos Indivíduos (Pares de Pontos) nos PCs (Scores para SEM)
cat("\nCoordenadas dos Pares de Pontos (Indivíduos) nos Componentes Principais (primeiras 5 linhas):\n")
scores_individuos <- get_pca_ind(res.pca)$coord
print(head(round(scores_individuos, 3)))

# Salvar os scores da PCA (input para o SEM)
# Estes são os valores que você usará como variáveis no seu modelo de Equações Estruturais.
write.csv(scores_individuos, file.path(output_dir, "sem_PCA_scores_resist_.csv"), row.names = TRUE)
message(paste0("Scores da PCA para SEM salvos em: ", file.path(output_dir, "sem_PCA_scores_resist_.csv")))


# Visualizar os indivíduos (pares de pontos) no plano dos PCs 1 e 2
fviz_pca_ind(res.pca, col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,
             title = "Gráfico de Indivíduos (Pares de Pontos) nos PCs") +
  theme_minimal()
ggsave(file.path(output_dir, "individuos_pc_plot.png"), device = "png", width = 8, height = 8)

############################# esta saindo lista adjacente ######################
#####sai por par, não 'matriz'. vai ser usado assim, pairwise

