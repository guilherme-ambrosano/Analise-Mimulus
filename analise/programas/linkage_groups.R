library(tidyverse)
library(onemap)
library(gt)
library(patchwork)

source("analise/programas/tabelas.R")
source("analise/programas/mapa_antigo.R")

dados <- read_onemap(input="analise/dados/m_feb06_onemap.raw") # 418 markers
bins <- create_data_bins(dados, find_bins(dados, exact = F))   # 416

get_number <- function(marker) {
  which(colnames(bins$geno)==marker)
}

nums_distorted <-  select_segreg(test_segregation(bins), 
                                 distorted = T, numbers = T)

# Estranhos:
# Generally, such markers placed in a single
# interval with high probability (LOD >> 2), but their
# placement resulted in substantial (< 7 cM) increases in
# map length, nonadditive two-point distances within the
# group, and many apparent genotyping errors. 

nums_estranhos <- c(which(!(colnames(bins$geno) %in% mapa_antigo$marcador |
                              str_starts(colnames(bins$geno), "MGSTS"))),
                    318, 402, 36, 4, 6, 135, 76, 245, 20,
                    167, 274)

keep <- setdiff(seq(1, 416), nums_estranhos)

dois_pts <- rf_2pts(bins)
final <- make_seq(dois_pts, keep) # 356

# Fazendo dois mapas separados para dominantes e codominantes

# table(sapply(lapply(apply(bins$geno[,keep], 2, unique), sort), paste, collapse="_"))

# dominantes   <- which(unlist(lapply(lapply(apply(bins$geno, 2, unique), sort), paste, collapse="_"))!="0_1_2_3")
# codominantes <- which(unlist(lapply(lapply(apply(bins$geno, 2, unique), sort), paste, collapse="_"))=="0_1_2_3")
# 
# final_codominantes <- drop_marker(final, dominantes)
# final_dominantes   <- drop_marker(final, codominantes)
# 
# LGs_codominantes <- group(final_codominantes, LOD = 5)
# LGs_dominantes   <- group(final_dominantes, LOD = 4)

# Comparando diferentes valores de lod para um LG especifico

# lista_lods <- lapply(c(4.5, 5, 5.5, 6, 6.5, 7), function(lod) {
#   LGs <- group(final, LOD = lod)
#   print(LGs$n.groups)
#   rcd(make_seq(LGs, 3))
#   })
# draw_map(lista_lods, names=T)

LGs <- group(final, LOD=7)

lista_mapa <- lapply(seq(1, LGs$n.groups), function(x) rcd(make_seq(LGs, x)))
draw_map(lista_mapa[c(6, 1, 2, 7, 3, 5, 4, 13, 10, 9, 11, 12, 14, 8)], names=T)

(Reduce(`+`, lapply(lista_mapa[c(6, 1, 2, 7, 3, 5, 4, 13, 10, 9, 11, 12, 14, 8)],
                    rf_graph_table)) &
  scale_fill_gradientn(colours=rainbow(4), limits=c(0, 0.5), na.value = "white")) +
  plot_layout(guides="collect")

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

mapa_antigo %>% 
  filter(LG==10) %>%
  pull(marcador) %>%
  sapply(get_LG_onemap) %>%
  unlist()

# Arrumando as posições de alguns marcadores, de acordo com o artigo

# Linkage group 6
LG6 <- lista_mapa[[6]] %>%
  drop_marker(c(33, 165, 257, 16)) %>%
  map() %>%
  try_seq(33) %>%
  make_seq(1) %>%
  try_seq(165) %>%
  make_seq(1) %>%
  try_seq(257) %>%
  make_seq(4) %>%
  try_seq(16) %>%
  make_seq(9)

# draw_map(list(LG6, lista_mapa[[6]]), names=T)
# Reduce(`+`, lapply(list(LG6, lista_mapa[[6]]), rf_graph_table))

lista_mapa[[6]] <- LG6

# Linkage group 1
LG1 <- lista_mapa[[1]] %>%
  drop_marker(c(251, 45, 342, 107, 314, 243)) %>%
  map() %>%
  try_seq(251) %>%
  make_seq(1) %>%
  try_seq(107) %>%
  make_seq(8) %>%
  try_seq(314) %>%
  make_seq(12) %>%
  try_seq(243) %>%
  make_seq(14) %>%
  try_seq(342) %>%
  make_seq(11) %>%
  try_seq(45) %>%
  make_seq(8)

lista_mapa[[1]] <- LG1

# Linkage group 2
LG2 <- lista_mapa[[2]] %>%
  drop_marker(100) %>%
  map() %>%
  try_seq(100) %>%
  make_seq(5) %>%
  drop_marker(c(32, 118, 129, 388)) %>%
  map() %>%
  try_seq(32) %>%
  make_seq(5) %>%
  try_seq(118) %>%
  make_seq(8) %>%
  try_seq(129) %>%
  make_seq(11) %>%
  try_seq(388) %>%
  make_seq(15)

lista_mapa[[2]] <- LG2

# Linkage group 7
LG7 <- lista_mapa[[7]] %>%
  drop_marker(c(73, 238)) %>%
  map() %>%
  try_seq(73) %>%
  make_seq(6) %>%
  try_seq(238) %>%
  make_seq(20)

lista_mapa[[7]] <- LG7

# Linkage group 3
LG3 <- lista_mapa[[3]] %>%
  drop_marker(c(5, 127, 155)) %>%
  map() %>%
  try_seq(5) %>%
  make_seq(5) %>%
  try_seq(127) %>%
  make_seq(7) %>%
  try_seq(155) %>%
  make_seq(12) %>%
  drop_marker(5) %>%
  map() %>%
  try_seq(5) %>%
  make_seq(8) %>%
  drop_marker(61) %>%
  map() %>%
  try_seq(61) %>%
  make_seq(6)
  
lista_mapa[[3]] <- LG3

# Linkage group 5
LG5 <- lista_mapa[[5]] %>%
  drop_marker(c(229, 82, 141, 253, 137, 179, 11)) %>%
  map() %>%
  try_seq(229) %>%
  make_seq(13) %>%
  try_seq(82) %>% 
  make_seq(13) %>%
  try_seq(141) %>% 
  make_seq(17) %>%
  try_seq(253) %>%
  make_seq(17) %>%
  try_seq(137) %>%
  make_seq(27) %>% # 29
  try_seq(179) %>%
  make_seq(31) %>%
  try_seq(11) %>%
  make_seq(25)

lista_mapa[[5]] <- LG5

# Linkage group 4
LG4 <- lista_mapa[[4]] %>%
  drop_marker(c(10, 114)) %>%
  map() %>%
  try_seq(10) %>%
  make_seq(4) %>%
  try_seq(114) %>%
  make_seq(6)

lista_mapa[[4]] <- LG4

# Linkage group 13
LG13 <- lista_mapa[[13]] %>%
  drop_marker(c(101, 349)) %>%
  map() %>%
  try_seq(101) %>%
  make_seq(16) %>%
  try_seq(349) %>% 
  make_seq(19)

lista_mapa[[13]] <- LG13

# Linkage group 9
LG9 <- lista_mapa[[9]] %>%
  drop_marker(c(227, 56, 69, 173, 106)) %>%
  map() %>%
  try_seq(227) %>%
  make_seq(8) %>%
  try_seq(56) %>%
  make_seq(12) %>%
  try_seq(69) %>%
  make_seq(15) %>%
  try_seq(173) %>%
  make_seq(13) %>% 
  try_seq(106) %>%
  make_seq(32) %>%
  drop_marker(c(373, 357, 292)) %>%
  map() %>%
  try_seq(373) %>%
  make_seq(18) %>%
  try_seq(357) %>%
  make_seq(5) %>%
  try_seq(292) %>%
  make_seq(4)

lista_mapa[[9]] <- LG9

# Linkage group 11
LG11 <- lista_mapa[[11]] %>%
  drop_marker(c(284, 250, 249, 218)) %>%
  map() %>%
  try_seq(284) %>%
  make_seq(7) %>%
  try_seq(250) %>%
  make_seq(10) %>%
  try_seq(218) %>%
  make_seq(17) %>%
  try_seq(249) %>% 
  make_seq(16)

lista_mapa[[11]] <- LG11

# Linkage group 12
LG12 <- lista_mapa[[12]] %>%
  drop_marker(c(181, 183, 328, 44)) %>% 
  map() %>% 
  try_seq(181) %>%
  make_seq(6) %>%
  try_seq(183) %>%
  make_seq(9) %>%
  try_seq(328) %>% 
  make_seq(9) %>%
  try_seq(44) %>%
  make_seq(6)

lista_mapa[[12]] <- LG12
  
# Linkage group 14
LG14 <- lista_mapa[[14]] %>%
  drop_marker(c(48, 206)) %>%
  map() %>% 
  try_seq(48) %>% 
  make_seq(15) %>% 
  try_seq(206) %>% 
  make_seq(16)

lista_mapa[[14]] <- LG14

# Linkage group 8
LG8 <- lista_mapa[[8]] %>%
  drop_marker(c(258, 235, 215, 125, 248)) %>%
  map() %>%
  try_seq(258) %>%
  make_seq(1) %>%
  try_seq(235) %>%
  make_seq(5) %>%
  try_seq(215) %>%
  make_seq(6) %>%
  try_seq(125) %>%
  make_seq(13) %>%
  try_seq(248) %>%
  make_seq(13)

lista_mapa[[8]] <- LG8

# Mapa final

lista_mapa_final <- lista_mapa[c(6, 1, 2, 7, 3, 5, 4, 13, 10, 9, 11, 12, 14, 8)]
lista_mapa_final[[1]]$seq.num <- rev(lista_mapa_final[[1]]$seq.num)
lista_mapa_final[[3]]$seq.num <- rev(lista_mapa_final[[3]]$seq.num)
lista_mapa_final[[6]]$seq.num <- rev(lista_mapa_final[[6]]$seq.num)
lista_mapa_final[[9]]$seq.num <- rev(lista_mapa_final[[9]]$seq.num)
lista_mapa_final[[10]]$seq.num <- rev(lista_mapa_final[[10]]$seq.num)
lista_mapa_final[[11]]$seq.num <- rev(lista_mapa_final[[11]]$seq.num)
lista_mapa_final[[13]]$seq.num <- rev(lista_mapa_final[[13]]$seq.num)
lista_mapa_final[[14]]$seq.num <- rev(lista_mapa_final[[14]]$seq.num)

draw_map(lista_mapa_final, names=T)
do.call(draw_map2, list(lista_mapa_final, output="analise/saida/map_ref.pdf"))

(Reduce(`+`, lapply(lista_mapa_final, rf_graph_table)) &
  scale_fill_gradientn(colours=rainbow(4), limits=c(0, 0.5), na.value = "white")) +
  plot_layout(guides="collect") +
  plot_annotation(tag_prefix = "LG", tag_levels = "1") +
  ggsave("analise/saida/heatmaps.pdf", width = 13, height = 8, dpi = 300)

# Comparando diferentes heuristicas

# funcoes <- list(rcd, record, ug, seriation, mds_onemap)
# ordens <- lapply(funcoes, function(x) lapply(lapply(seq(1, LGs$n.groups), function(y) make_seq(LGs, y)), x))
# 
# graficos_algoritmos <- lapply(seq(1, 3), function(x) rf_graph_table(ordens[[x]][[1]]))
# Reduce(`+`, graficos_algoritmos)
# 
# cbind(mds_onemap(make_seq(LGs, 2))[[1]],
#       seriation( make_seq(LGs, 2))[[1]],
#       rcd(       make_seq(LGs, 2))[[1]],
#       record(    make_seq(LGs, 2))[[1]],
#       ug(        make_seq(LGs, 2))[[1]])
# 
