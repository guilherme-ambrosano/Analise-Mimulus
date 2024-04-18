library(tidyverse)
library(magrittr, include.only = "%T>%")
library(onemap)
library(gt)
library(patchwork)
library(ggrepel)

source("analise/programas/03_mapa_antigo.R")


# Carregando os dados ----

dados <- read_onemap(input="analise/dados/m_feb06_onemap.raw") # 418 markers
bins <- create_data_bins(dados, find_bins(dados, exact = F))   # 416

# uma funcao pra converter o nome do marcador no numero que o onemap entende
get_number <- function(marker) {
  which(colnames(bins$geno)==marker)
}
get_name <- function(number) {
  colnames(bins$geno)[number]
}

# filtrando marcadores para garantir que vao restar 14 grupos de ligacao
filtr <- c(402, 236, 220, 374, 258, 167, 393)
keep <- setdiff(seq(1, 416), filtr)

dois_pts <- rf_2pts(bins)
final <- make_seq(dois_pts, keep) # 407

# Fazendo dois mapas separados para dominantes e codominantes ----

# table(sapply(lapply(apply(bins$geno[,keep], 2, unique), sort), paste, collapse="_"))

# dominantes   <- which(unlist(lapply(lapply(apply(bins$geno, 2, unique), sort), paste, collapse="_"))!="0_1_2_3")
# codominantes <- which(unlist(lapply(lapply(apply(bins$geno, 2, unique), sort), paste, collapse="_"))=="0_1_2_3")
# 
# final_codominantes <- drop_marker(final, dominantes)
# final_dominantes   <- drop_marker(final, codominantes)
# 
# LGs_codominantes <- group(final_codominantes, LOD = 5)
# LGs_dominantes   <- group(final_dominantes, LOD = 4)

# Comparando diferentes valores de lod ----

# lista_lods <- lapply(c(5.41, 5.5, 6, 6.5, 7), function(lod) {
#   LGs <- group(final, LOD = lod)
#   rcd(make_seq(LGs, 10))
#   })
# draw_map(lista_lods, names=T)


# Fazendo o mapa com as heuristicas ----

suggest_lod(bins)
qchisq(0.05/ncol(combn(407, 2)), 1, lower=F)*0.2172
LGs <- group(final, LOD=5.41)

lista_mapa <- lapply(seq(1, LGs$n.groups), function(x) rcd(make_seq(LGs, x)))
draw_map(lista_mapa, names=T)

# Funcoes para obter o grupo de ligacao de um marcador
get_LG_onemap <- function(marker=NULL, number=NULL) {
  if(is.null(number)) {
    number <- get_number(marker)
  }
  which(unlist(lapply(seq(1, length(lista_mapa)), function(x) number %in% lista_mapa[[x]]$seq.num)))
}

get_LG_mapa_antigo <- function(marker) {
  mapa_antigo[mapa_antigo$marcador==marker, "LG"]
}
get_LG_onemap("BA153")
get_LG_mapa_antigo("BA153")

# Garantindo que os marcadores do paper permanecem nos mesmos grupos de ligacao
lapply(seq(1, 14), function(x) {
  mapa_antigo %>% 
    filter(LG==x) %>%
    pull(marcador) %>%
    sapply(get_LG_onemap) %>%
    unlist()
})

# Deixando na ordem do artigo
lista_mapa <- lista_mapa[c(1:9, 11:14, 10)]

draw_map(lista_mapa, names=T)

(Reduce(`+`, lapply(lista_mapa, rf_graph_table)) &
  scale_fill_gradientn(colours=rainbow(4), limits=c(0, 0.5), na.value = "white")) +
  plot_layout(guides="collect")


# Comparando diferentes heuristicas

# funcoes <- list(rcd, record, ug, mds_onemap)
# ordens <- lapply(funcoes, function(x) lapply(lapply(seq(1, LGs$n.groups), function(y) make_seq(LGs, y)), x))
# graficos_algoritmos <- lapply(seq(1, length(funcoes)), function(x) rf_graph_table(ordens[[x]][[1]]))
# Reduce(`+`, graficos_algoritmos)

cbind(mds_onemap(make_seq(LGs, 2))[[1]],
      seriation( make_seq(LGs, 2))[[1]],
      rcd(       make_seq(LGs, 2))[[1]],
      record(    make_seq(LGs, 2))[[1]],
      ug(        make_seq(LGs, 2))[[1]])


# Arrumando o mapa manualmente ----

# Função pra fazer o mapa com o ggplot
ggdraw_map_comp <- function(...) {
  grupos <- list(...)
  lapply(grupos, function(x) { 
    tibble(pos=c(0, cumsum(kosambi(x[[3]]))),
           marker=sapply(x[[1]], get_name))
    }) %>%
    bind_rows(.id="grupo") %>%
    mutate_at(vars(grupo), str_replace, "^", "Alternativa ") %>%
    mutate_at(vars(grupo), factor, levels=paste("Alternativa", seq(1, length(grupos)))) %>% 
    ggplot(aes(x=grupo, y=pos)) +
    geom_text_repel(aes(label=marker), hjust=0, nudge_x = 0.1) +
    geom_point(shape=3, stroke=1.5, size=0.1) +
    geom_line(size=1) +
    labs(x="") +
    scale_y_reverse(name="Distance (cM)") +
    theme_minimal() +
    theme()
}

# Função pra plotar o mapa junto com o heatmap
# Também é possível comparar com um mapa alternativo, lado a lado
plotar_mapa <- function(x, print=T, alternativa=NULL) {
  if(length(alternativa)!=0) {
    plot <- ((rf_graph_table(x) +
                rf_graph_table(alternativa) &
                scale_fill_gradientn(colours = rainbow(4), limits=c(0, 0.5),
                                     na.value = "white")) +
               plot_layout(guides="collect")) /
      ggdraw_map_comp(x, alternativa)
  } else {
    plot <- rf_graph_table(x) + 
      (ggdraw_map_comp(x) +
         theme(axis.text.x = element_blank())) +
      plot_layout(design = "112
                            112")
  }
  
  if(print) {
    print(plot)
  }
  return(x)
}

# Função para remover marcadores, rodar o try_seq do onemap e 
# escolher a opção com o LOD = 0
# Se houverem posições alternativas, ele vai avisar
auto_try_seq <- function(dados, ...) {
  marcadores <- list(0, ...)
  reduce(marcadores, .f=function(dados_atualizados, marcador) {
    if(is(dados_atualizados, "numeric")) {
      dados_atualizados <- dados
    }
    dados_atualizados %>%
      drop_marker(marcador) %>%
      map() %>%
      try_seq(marcador) %T>%
      {
        lods <- .[[2]]
        alternativas <- lods != 0 & lods > -3
        if(sum(alternativas) > 0) {
          cat(paste(sum(alternativas), "alternativa (s):\t"))
          cat(paste(which(alternativas), collapse = ", "))
          cat("\n")
        }
      } %>%
      make_seq(which(.[[2]]==0))
  })
}

# marcadores que foram filtrados no paper
filtr_paper <- which(!(colnames(bins$geno) %in% mapa_antigo$marcador |
                         str_starts(colnames(bins$geno), "MgSTS")))

# marcadores distorcidos (nao seguem a segregacao prevista pela lei de mendel)
nums_distorted <-  select_segreg(test_segregation(bins), 
                                 distorted = T, numbers = T)

# Linkage group 1
lista_mapa[[1]]$seq.num[lista_mapa[[1]]$seq.num %in% filtr_paper]
lista_mapa[[1]]$seq.num[lista_mapa[[1]]$seq.num %in% nums_distorted]

LG1 <- lista_mapa[[1]] %>%
  drop_marker(c(8, 1, 225, 219, 130, 33, 325, 49)) %>%
  map() %>%
  drop_marker(c(257, 28)) %>%
  map() %>%
  auto_try_seq(257) %>%
  try_seq(28) %>%
  make_seq(2) %>%
  auto_try_seq(47, 325, 165, 257) %>%
  plotar_mapa(alt=lista_mapa[[1]])

# Linkage group 2
lista_mapa[[2]]$seq.num[lista_mapa[[2]]$seq.num %in% filtr_paper]
lista_mapa[[2]]$seq.num[lista_mapa[[2]]$seq.num %in% nums_distorted]

LG2 <- lista_mapa[[2]] %>%
  drop_marker(c(38, 211, 342, 
                2, 45, 243)) %>%
  map() %>%
  auto_try_seq(251, 293) %>%
  drop_marker(107) %>%
  map() %>%
  try_seq(107) %>%
  make_seq(9) %>%
  plotar_mapa()

# Linkage group 3
lista_mapa[[3]]$seq.num[lista_mapa[[3]]$seq.num %in% filtr_paper]
lista_mapa[[3]]$seq.num[lista_mapa[[3]]$seq.num %in% nums_distorted]

LG3 <- lista_mapa[[3]] %>%
  drop_marker(c(254, 113, 116, 228, 128, 182, 388,
                32)) %>%
  map() %>%
  auto_try_seq(319, 53, 3, 217, 100,
               361, 410, 129, 118, 216, 91, 123) %>%
  plotar_mapa()
  
# Linkage group 4
lista_mapa[[4]]$seq.num[lista_mapa[[4]]$seq.num %in% filtr_paper]
lista_mapa[[4]]$seq.num[lista_mapa[[4]]$seq.num %in% nums_distorted]

LG4 <- lista_mapa[[4]] %>%
  drop_marker(c(124, 222, 207, 80, 205, 62, 143, 153, 24, 42, 189, 133, 144,
                76, 135, 307, 333)) %>%
  map() %>%
  auto_try_seq(288, 105, 73, 157, 238, 109, 392, 344, 318) %>%
  map() %>%
  plotar_mapa()

# Linkage group 5
lista_mapa[[5]]$seq.num[lista_mapa[[5]]$seq.num %in% filtr_paper]
lista_mapa[[5]]$seq.num[lista_mapa[[5]]$seq.num %in% nums_distorted]

LG5 <- lista_mapa[[5]] %>%
  drop_marker(c(139, 111, 85, 90, 233, 150, 5)) %>%
  map() %>%
  auto_try_seq(54, 155, 127, 111, 61) %>%
  plotar_mapa()

# Linkage group 6
lista_mapa[[6]]$seq.num[lista_mapa[[6]]$seq.num %in% filtr_paper]
lista_mapa[[6]]$seq.num[lista_mapa[[6]]$seq.num %in% nums_distorted]

LG6 <- lista_mapa[[6]] %>%
  drop_marker(c(223, 9, 51, 12, 104, 39, 239, 97, 203, 212, 148,
                20, 11, 141, 245, 265, 309)) %>%
  map() %>%
  auto_try_seq(179, 253, 6, 82, 229, 253, 184, 415) %>%
  plotar_mapa()

# Linkage group 7
lista_mapa[[7]]$seq.num[lista_mapa[[7]]$seq.num %in% filtr_paper]
lista_mapa[[7]]$seq.num[lista_mapa[[7]]$seq.num %in% nums_distorted]

LG7 <- lista_mapa[[7]] %>%
  drop_marker(c(147, 117, 110, 237, 40, 70, 199, 234,
                220, 10, 136, 163)) %>%
  map() %>%
  auto_try_seq(117, 353, 338, 114, 160, 362) %>%
  plotar_mapa()

# Linkage group 8
lista_mapa[[8]]$seq.num[lista_mapa[[8]]$seq.num %in% filtr_paper]
lista_mapa[[8]]$seq.num[lista_mapa[[8]]$seq.num %in% nums_distorted]

LG8 <- lista_mapa[[8]] %>%
  drop_marker(c(14, 149, 178, 186,
                36, 252, 171, 406, 321, 300)) %>%
  map() %>%
  auto_try_seq(68, 46, 164, 101, 232,
               349, 371, 313, 315) %>%
  plotar_mapa()

# Linkage group 9
lista_mapa[[9]]$seq.num[lista_mapa[[9]]$seq.num %in% filtr_paper]
lista_mapa[[9]]$seq.num[lista_mapa[[9]]$seq.num %in% nums_distorted]

LG9 <- lista_mapa[[9]] %>%
  drop_marker(c(226, 27, 221, 29, 17, 67)) %>%
  map() %>%
  auto_try_seq(59, 86, 120, 214, 172,
               29, 67,
               308, 37, 96) %>%
  plotar_mapa()

# Linkage group 10
lista_mapa[[10]]$seq.num[lista_mapa[[10]]$seq.num %in% filtr_paper]
lista_mapa[[10]]$seq.num[lista_mapa[[10]]$seq.num %in% nums_distorted]

LG10 <- lista_mapa[[10]] %>%
  drop_marker(c(196, 191, 210, 131, 169, 94, 
                227, 151, 102, 357, 337, 298)) %>%
  map() %>%
  auto_try_seq(81, 69, 81, 56, 25, 177, 58) %>%
  auto_try_seq(191, 255, 22, 122, 373, 320) %>% 
  plotar_mapa()

# Linkage group 11
lista_mapa[[11]]$seq.num[lista_mapa[[11]]$seq.num %in% filtr_paper]
lista_mapa[[11]]$seq.num[lista_mapa[[11]]$seq.num %in% nums_distorted]

LG11 <- lista_mapa[[11]] %>%
  map() %>% 
  auto_try_seq(250, 43, 284, 26, 218) %>%
  plotar_mapa()

# Linkage group 12
lista_mapa[[12]]$seq.num[lista_mapa[[12]]$seq.num %in% filtr_paper]
lista_mapa[[12]]$seq.num[lista_mapa[[12]]$seq.num %in% nums_distorted]

LG12 <- lista_mapa[[12]] %>%
  drop_marker(c(187, 197, 202, 192, 204, 241, 195, 134, 174, 180, 379, 328)) %>%
  map() %>%
  auto_try_seq(181, 183, 77) %>%
  plotar_mapa()

# Linkage group 13
lista_mapa[[13]]$seq.num[lista_mapa[[13]]$seq.num %in% filtr_paper]
lista_mapa[[13]]$seq.num[lista_mapa[[13]]$seq.num %in% nums_distorted]

LG13 <- lista_mapa[[13]] %>%
  drop_marker(c(240, 279, 327)) %>%
  map() %>%
  auto_try_seq(78, 206, 63, 380, 301) %>%
  plotar_mapa()

# Linkage group 14
lista_mapa[[14]]$seq.num[lista_mapa[[14]]$seq.num %in% filtr_paper]
lista_mapa[[14]]$seq.num[lista_mapa[[14]]$seq.num %in% nums_distorted]

LG14 <- lista_mapa[[14]] %>%
  drop_marker(c(98, 175, 
                274, 343, 332)) %>%
  map() %>% 
  auto_try_seq(258, 60, 273, 322, 350, 416, 376, 304,
               391, 305, 374, 358) %>%
  plotar_mapa()

# Mapa final ----

lista_mapa_final <- list(LG1, LG2, LG3, LG4, LG5, LG6, LG7,
                         LG8, LG9, LG10, LG11, LG12, LG13, LG14)
inverter_mapa <- function(x) {
  y <- x
  y$seq.num <- rev(x$seq.num)
  y$seq.rf <- rev(x$seq.rf)
  return(y)
}

# Invertendo alguns grupos que ficaram de cabeça pra baixo
lista_mapa_final[c(1, 3, 4, 8, 9, 10, 11, 12, 13, 14)] <-
  lapply(lista_mapa_final[c(1, 3, 4, 8, 9, 10, 11, 12, 13, 14)], inverter_mapa)

# Salvando o mapa
save(bins, dois_pts, lista_mapa_final, file="analise/saida/LGs.RData")
