library(tidyverse)
library(onemap)

dados_mapmaker <- read_table("originais/dados/m_feb06.raw", col_names = F)

# Separando o arquivo do mapmaker em dados novos e antigos e colocando
# no formato do onemap

# Arrumando os cabecalhos
cabecalho <- cabecalho_novos <- cabecalho_antigos <- dados_mapmaker %>%
  slice(1, 2)

info <- as.numeric(strsplit(pull(cabecalho[2, 1]), " ")[[1]])
cabecalho[2, 1] <- paste(info[1], info[2], 0, 0, info[3])

# Separando os marcadores MgSTS (novos)
# Deixando todas as informacoes fenotipicas no conjunto de dados antigo
# Senao o combine_onemap do CRAN nao funciona
cabecalho_novos[2, 1] <- paste(info[1], info[2] - 261, 0, 0, 0)
cabecalho_antigos[2, 1] <- paste(info[1], 261, 0, 0, info[3])

# Os marcadores CA174 e MgSTS16 estao formatados diferente dos outros
# (tudo em uma mesma linha)
CA174 <- dados_mapmaker %>%
  slice(195) %>%
  mutate_all(str_replace, " +", " ")

MgSTS16 <- dados_mapmaker %>%
  slice(522) %>%
  mutate_all(str_replace, " +", " ")

# Os outros marcadores estao com o nome em uma linha e as informacoes na de baixo
dados <- dados_mapmaker %>%
  slice(-195, -522) %>% # tirando o CA174 e o MgSTS16
  slice(seq(3, 836)) %>% # tirando o cabecalho e as caracteristicas fenotipicas
  mutate_all(str_remove_all, " ") %>%
  mutate_all(~paste(lag(.), .)) %>%  # trazendo as informacoes pra linha de cima
  slice(seq(2, 832, 2)) %>%  # apagando as linhas pares
  add_row(CA174, .before=195) %>%
  add_row(MgSTS16, .after=417) %>% 
  separate(X1, c("chrom", "gen"), sep=" ") %>%
  mutate_at(vars(gen), str_to_lower) %>% # formato do onemap
  # adicionando espaco entre os individuos
  pmap_dfr(function(chrom, gen) {
    tibble(X1 = paste(chrom, paste(strsplit(gen, "")[[1]], collapse=" ")))
  })

# Comparando o numero de marcadores com o artigo
dados %>%
  mutate_at(vars(X1), str_replace, " ", "_") %>%
  separate(X1, c("chrom", "gen"), sep = "_") %>%
  mutate_at(vars(chrom), str_remove, "\\*") %>%
  mutate_at(vars(chrom), str_to_upper) %>%
  mutate(tipo = case_when(chrom %in% c("CYCA", "CYCB", "LFY", "AP3") ~ "gen",
                          str_starts(chrom, "AAT") | str_starts(chrom, "AG") ~ "microssat",
                          str_starts(chrom, "MGSTS") ~ "novo",
                          T ~ "AFLP")) %>%
  mutate_at(vars(gen), function(x) {
    sapply(x, function(y) {
      paste(sort(unique(strsplit(y, " ")[[1]])), collapse=", ")
    })
  }) %>%
  group_by(tipo, gen) %>%
  count()

# Paper:
#  27 AFLP codominantes
# 197 AFLP   dominantes
#   4 gene-based
#  27 microssatélites
#----------------------
# 255 total

# Dados:
#  27 AFLP codominantes
# 204 AFLP   dominantes
#   4 gene-based
#  26 microssatélites
#----------------------
# 261 total +
# 157 MgSTS

# Vendo os que restaram após a filtragem do paper
source("analise/programas/03_mapa_antigo.R")
dados %>%
  mutate_at(vars(X1), str_replace, " ", "_") %>%
  separate(X1, c("chrom", "gen"), sep = "_") %>%
  mutate_at(vars(chrom), str_remove, "\\*") %>%
  mutate(tipo = case_when(chrom %in% c("CYCA", "CYCB", "LFY", "AP3") ~ "gen",
                          str_starts(chrom, "AAT") | str_starts(chrom, "AG") ~ "microssat",
                          str_starts(chrom, "MgSTS") ~ "novo",
                          T ~ "AFLP")) %>%
  mutate_at(vars(chrom), str_to_upper) %>%
  mutate_at(vars(gen), function(x) {
    sapply(x, function(y) {
      paste(sort(unique(strsplit(y, " ")[[1]])), collapse=", ")
    })
  }) %>%
  filter(chrom %in% mapa_antigo$marcador) %>%
  group_by(gen, tipo) %>%
  count()

# Total: 174, igual o paper

# Adicionando informacao do tipo de segregacao (formato do onemap)
mktype <- dados %>% 
  mutate_all(str_replace, " ", "_") %>%
  separate(X1, c("chrom", "gen"), sep="_") %>%
  # pegando os caracteres unicos que aparecem no gen (bd, ac ou abh)
  pmap_dfr(function(chrom, gen) {
    tibble(X1 = paste(chrom, paste(sort(unique(strsplit(gen, " ")[[1]])), collapse="")))
  }) %>%
  mutate_all(str_replace, " ", "_") %>%
  separate(X1, c("chrom", "gen"), sep="_") %>%
  mutate_at(vars(gen), str_remove, "^-") %>% # ignorando os dados perdidos
  mutate(mktype = case_when(gen=="bd" ~ "D.B",
                            gen=="ac" ~ "C.A",
                            gen=="abh" ~ "A.H.B")) %>%
  select(chrom, mktype)

# Adicionando o mktype aos dados
dados <- dados %>%
  mutate_all(str_replace, " ", "_") %>%
  separate(X1, c("chrom", "gen"), sep="_") %>%
  left_join(mktype, by=c("chrom")) %>%
  mutate_at(vars(gen), str_replace_all, "h", "ab") %>% # formato do onemap
  # deixando os nomes dos marcadores em caixa alta pra padronizar
  # (tem algum aat com o nome todo em minusculas)
  mutate_at(vars(chrom), str_to_upper) %>%
  mutate_at(vars(chrom), str_replace, "MGSTS", "MgSTS") %>%
  unite(X1, chrom, mktype, gen, sep=" ")

# Criando um data frame com as caracteristicas fenotipicas
linhas_variaveis <- which(str_detect(pull(dados_mapmaker[,1]), "^\\*"))[419:434]
variaveis <- tibble(inicio = linhas_variaveis) %>%
  mutate(fim = lead(inicio) - 1) %>% 
  mutate(fim = ifelse(is.na(fim), 1434, fim)) %>%
  pmap_dfr(function(inicio, fim) {
    dados_mapmaker %>%
      slice(seq(inicio, fim)) %>%
      mutate(linha=row_number()) %>%
      pivot_wider(names_from=linha, values_from=X1) %>%
      unite(X1, sep="\t") %>%
      mutate_all(str_replace_all, " +", " ") %>%
      mutate_all(str_replace_all, "\\t", " ")
  })

# linha com os individuos (formato do onemap)
individuos <- paste(sapply(seq(1, info[1]), function(x) paste0("I", x)), collapse=" ")

# juntando tudo pra fazer os conjuntos de dados
dados_onemap <- rbind(cabecalho, individuos, dados, variaveis)
dados_onemap_novos <- rbind(cabecalho_novos, individuos, filter(dados, str_starts(X1, "*MGSTS")))
dados_onemap_antigos <- rbind(cabecalho_antigos, individuos, filter(dados, !str_starts(X1, "*MGSTS")), variaveis)

# escrevendo os arquivos
conexao <- file("analise/dados/m_feb06_onemap.raw")
writeLines(pull(dados_onemap), conexao)
close(conexao)

conexao <- file("analise/dados/m_feb06_novos.raw")
writeLines(pull(dados_onemap_novos), conexao)
close(conexao)

conexao <- file("analise/dados/m_feb06_antigos.raw")
writeLines(pull(dados_onemap_antigos), conexao)
close(conexao)
