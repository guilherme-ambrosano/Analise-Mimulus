library(tidyverse)
library(onemap)
library(gt)

library(ggvenn)
library(ggrepel)
library(patchwork)

source("analise/programas/02_tabelas.R")
source("analise/programas/03_mapa_antigo.R")
load("analise/saida/LGs.RData")

dados <- read_onemap(input="analise/dados/m_feb06_onemap.raw") # 418 markers

# A function to get the marker's number from its code
get_number <- function(marker) {
  which(colnames(bins$geno)==marker)
}
# A function to get the marker's code from its number
get_name <- function(number) {
  colnames(bins$geno)[number]
}

mapa_novo <- lapply(lista_mapa_final, function(x) tibble(marcador=x$seq.num)) %>%
  bind_rows(.id="LG") %>%
  mutate_at(vars(LG), as.numeric) %>%
  mutate_at(vars(marcador), ~sapply(., get_name))

length(unique(mapa_novo$marcador)) # 288

174/255 # 68,2% filtered in the paper
288/418 # 68,9% filtered

sum(sapply(lista_mapa_final, function(x) {
  cumsum(kosambi(x$seq.rf))[length(x$seq.rf)]
}))  # 2374.856 cM

# Map ----

# Plot the map
ggdraw_map <- function(...) {
  grupos <- list(...)
  lapply(grupos, function(x) { 
    tibble(pos=c(0, cumsum(kosambi(x[[3]]))),
           marker=sapply(x[[1]], get_name))
    }) %>%
    bind_rows(.id="grupo") %>%
    mutate_at(vars(grupo), str_replace, "^", "LG") %>%
    add_row(tibble(pos=c(0, 10, 20), marker=c("", "20 cM", ""), grupo=rep("", 3))) %>%
    mutate_at(vars(grupo), factor, levels=c(paste0("LG", seq(1, length(grupos))), "")) %>% 
    ggplot(aes(x=1, y=pos)) +
    geom_text_repel(aes(label=marker, size=grupo==""),
                    point.padding=1, box.padding=0.35, min.segment.length = 0.1,
                    max.time = 2, max.iter = 200000, seed=123) +
    geom_point(aes(color=grupo==""), shape=3, stroke=1.5, size=0.1) +
    scale_color_manual(values=c("black", "transparent")) +
    scale_size_manual(values = c(2, 4)) +
    guides(color=F, size=F) +
    geom_line(size=1) +
    facet_wrap(~grupo, ncol=5) +
    labs(x="") +
    scale_x_continuous(limits=c(0, 2)) +
    scale_y_reverse() +
    theme_void() +
    theme(strip.text = element_text(size=15, family="serif", hjust=0.5),
          panel.spacing = unit(1, "pt"),
          plot.margin = unit(c(0.75, 0.25, 0.15, 0.5), "in"))
}

# Plotting...
do.call(ggdraw_map, lista_mapa_final) +
  ggsave("analise/saida/mapa_final.pdf", width = 8, height = 11, dpi=300)
# do.call(draw_map2, list(lista_mapa_final, output="analise/saida/map_ref.pdf"))

# Heatmaps ----
((Reduce(`+`, lapply(lista_mapa_final, rf_graph_table)) &
    theme(axis.text = element_text(size=4)) &
    scale_fill_gradientn(colours=rainbow(4), limits=c(0, 0.5), na.value = "white")) +
    plot_layout(guides="collect", ncol=3) +
    plot_annotation(tag_prefix = "LG", tag_levels = "1", 
                    theme=theme(plot.margin = unit(c(0.9, 0.9, 0.9, 0.9), "in"))) &
    guides(fill = guide_colorbar(barheight = unit(0.1, "in"), barwidth = unit(5, "in"),
                                 title.vjust = 1)) &
    labs(fill="Recombination fraction") &
    theme(plot.tag.position = "top",
          plot.tag = element_text(family="serif"),
          legend.position = "bottom")) +
  ggsave("analise/saida/heatmaps.pdf", width = 8, height = 11, dpi = 300)

# Transmission rate distortion ----

# Table with segregation information
freqs <- tabela_segregacao(dados)

# Get locus position in the map
get_position <- function(LGs, numeros) {
  lista_mapa_invertido <- lapply(lista_mapa_final, inverter_mapa)
  posicoes <- sapply(seq(1, length(numeros)), function(x) {
    loco <- which(lista_mapa_invertido[[LGs[x]]][[1]] == numeros[x])
    c(0, cumsum(kosambi(lista_mapa_invertido[[LGs[x]]][[3]])))[loco]
  })
}

inverter_mapa <- function(x) {
  y <- x
  y$seq.num <- rev(x$seq.num)
  y$seq.rf <- rev(x$seq.rf)
  return(y)
}

# Figure 3
freqs %>%
  mutate(numero = sapply(marcador, get_number),
         LG=case_when(numero %in% lista_mapa_final[[1]]$seq.num  ~ 1,
                      numero %in% lista_mapa_final[[2]]$seq.num  ~ 2,
                      numero %in% lista_mapa_final[[3]]$seq.num  ~ 3,
                      numero %in% lista_mapa_final[[4]]$seq.num  ~ 4,
                      numero %in% lista_mapa_final[[5]]$seq.num  ~ 5,
                      numero %in% lista_mapa_final[[6]]$seq.num  ~ 6,
                      numero %in% lista_mapa_final[[7]]$seq.num  ~ 7,
                      numero %in% lista_mapa_final[[8]]$seq.num  ~ 8,
                      numero %in% lista_mapa_final[[9]]$seq.num  ~ 9,
                      numero %in% lista_mapa_final[[10]]$seq.num ~ 10,
                      numero %in% lista_mapa_final[[11]]$seq.num ~ 11,
                      numero %in% lista_mapa_final[[12]]$seq.num ~ 12,
                      numero %in% lista_mapa_final[[13]]$seq.num ~ 13,
                      numero %in% lista_mapa_final[[14]]$seq.num ~ 14,
                      T ~ 0)) %>%
  filter(LG!=0) %>%
  mutate(pos = get_position(LG, numero)) %>%
  mutate_at(vars(A, H, B, D, C), str_remove, " \\([:digit:]+,[:digit:]%\\)$") %>%
  mutate_at(vars(A, H, B, D, C), as.numeric) %>%
  # B = GG /A = NN
  mutate(GG = case_when(segregacao=="A.H.B" ~ 0.25 - B/(A+H+B),
                        segregacao=="D.B" ~ 0.25 - B/(B+D)),
         NN = case_when(segregacao=="A.H.B" ~ A/(A+H+B)-0.25,
                        segregacao=="C.A" ~ A/(A+C) -0.25),
         cor = ifelse(p.val.adj.signif=="ns", "Nonsignificant", "Significant")) %>%
  pivot_longer(cols=c(GG, NN)) %>%
  select(LG, segregacao, A, H, B, D, C, pos, name, value, cor) %>%
  mutate_at(vars(LG), ~factor(paste0("LG", .), levels=paste0("LG", seq(1, 14)))) %>% 
  group_by(LG) %>%
  mutate(max_pos=max(pos)) %>%
  ungroup() %>%
  ggplot(aes(x=pos, y=value)) +
  facet_wrap(~LG, ncol=2, scales="free_y", dir = "v") +
  geom_point(aes(shape=name, color=cor), stroke=1, size=1) +
  geom_segment(aes(x=0, y=0, xend=max_pos, yend=0), linetype="dotted") +
  scale_shape_manual(values=c(4, 3)) +
  scale_color_manual(values=c("lightgrey", "black")) +
  scale_y_continuous(limits=c(-0.2, 0.2), breaks=c(-.3, -.2, -.1, 0, .1, .2, .3)) +
  geom_text(data=tibble(x=c(180, 180),
                        y=c(0.13, -0.13),
                        LG=factor(c("LG1", "LG1"), levels=paste0("LG", seq(1, 14))),
                        label=c("NN excess/\nGG deficit",
                                "GG excess/\nNN deficit")),aes(x=x, y=y, label=label),
            family="serif", fontface="bold", size=3) +
  labs(x="Marker position (cM)", y="Deviation of homozygous genotype frequency from Mendelian expectation",
       shape="", color="") +
  theme_minimal(base_family="serif") +
  theme(panel.grid = element_blank(),
        axis.ticks = element_line(),
        axis.line = element_line(),
        strip.text = element_text(family="serif", hjust=0),
        plot.margin = unit(c(0.9, 1, 1, 0.9), "in"),
        legend.position="bottom") +
  ggsave("analise/saida/TRD.pdf", width = 8, height = 11, dpi=300)


# How many markers of each type ----

tibble(chrom=colnames(bins$geno),
       segr=bins$segr.type) %>%
  mutate(tipo = case_when(chrom %in% c("CYCA", "CYCB", "LFY", "AP3") ~ "gen",
                          str_starts(chrom, "AAT") | str_starts(chrom, "AG") ~ "microssat",
                          str_starts(chrom, "MgSTS") ~ "novo",
                          T ~ "AFLP")) %>%
  group_by(tipo, segr) %>%
  count()

# non-distorted
tibble(chrom=colnames(bins$geno),
       segr=bins$segr.type) %>%
  filter(chrom %in% select_segreg(test_segregation(bins), distorted = F, numbers = F)) %>%
  mutate(tipo = case_when(chrom %in% c("CYCA", "CYCB", "LFY", "AP3") ~ "gen",
                          str_starts(chrom, "AAT") | str_starts(chrom, "AG") ~ "microssat",
                          str_starts(chrom, "MgSTS") ~ "novo",
                          T ~ "AFLP")) %>%
  group_by(tipo, segr) %>%
  count()

# previous map
tibble(chrom=colnames(bins$geno),
       segr=bins$segr.type) %>%
  filter(chrom %in% select_segreg(test_segregation(bins), distorted = F, numbers = F)) %>%
  filter(chrom %in% mapa_antigo$marcador) %>%
  mutate(tipo = case_when(chrom %in% c("CYCA", "CYCB", "LFY", "AP3") ~ "gen",
                          str_starts(chrom, "AAT") | str_starts(chrom, "AG") ~ "microssat",
                          str_starts(chrom, "MgSTS") ~ "novo",
                          T ~ "AFLP")) %>%
  group_by(tipo, segr) %>%
  count()

# new map
tibble(chrom=colnames(bins$geno),
       segr=bins$segr.type) %>%
  filter(chrom %in% select_segreg(test_segregation(bins), distorted = F, numbers = F)) %>%
  filter(chrom %in% mapa_novo$marcador) %>%
  mutate(tipo = case_when(chrom %in% c("CYCA", "CYCB", "LFY", "AP3") ~ "gen",
                          str_starts(chrom, "AAT") | str_starts(chrom, "AG") ~ "microssat",
                          str_starts(chrom, "MgSTS") ~ "novo",
                          T ~ "AFLP")) %>%
  group_by(tipo, segr) %>%
  count()

# TRD ---
freqs %>%
  gt()

freqs %>%
  filter(marcador %in% mapa_novo$marcador) %>%
  group_by(segregacao) %>%
  count() %>%
  ungroup() %>%
  mutate(porc = n/sum(n))

# Total:
# segreg         n  porc
# A.H.B        213 0.510
# C.A           92 0.221
# D.B          113 0.270

# Mapa:
# segreg         n  porc
# A.H.B        172 0.597
# C.A           51 0.177
# D.B           65 0.226

freqs %>% 
  filter(marcador %in% mapa_novo$marcador) %>%
  mutate_at(vars(p.val.adj.signif), str_replace, "\\*+", "*") %>%
  group_by(p.val.adj.signif, segregacao) %>%
  count() %>%
  pivot_wider(names_from=segregacao, values_from=n)

freqs %>% 
  filter(marcador %in% mapa_novo$marcador) %>%
  mutate_at(vars(p.val.adj.signif), str_replace, "\\*+", "*") %>%
  group_by(p.val.adj.signif) %>%
  count() %>%
  ungroup() %>%
  mutate(porc = n/sum(n))

# Total: 62/418 (14.8%) w/ TRD - p < 0.05
# Map:   41/288 (14.2%) w/ TRD

# No paper, 49%
# Não corrigiram o p-valor?

freqs %>%
  filter(marcador %in% mapa_novo$marcador) %>%
  mutate_at(vars(p.val.adj.signif), str_replace, "\\*+", "*") %>%
  mutate_at(vars(segregacao), ~case_when(.=="A.H.B" ~ "Codominante",
                                         T ~ "Dominante")) %>%
  group_by(segregacao, p.val.adj.signif) %>%
  count() %>%
  group_by(segregacao) %>%
  mutate(porc = n/sum(n))

# Total: 45 (21.1%) of codominants e 17 (8.3%) of dominants
# Map:   33 (19.2%) of codominants e  8 (6.9%) of dominants

freqs %>%
  filter(marcador %in% mapa_novo$marcador) %>%
  group_by(p.val.adj.signif) %>%
  count() %>%
  ungroup() %>%
  mutate(porc = n/sum(n))

# Total: 39 (9.3%) w/ p < 0,001
# Map:   28 (9.7%) w/ p < 0,001

freqs %>%
  filter(marcador %in% mapa_novo$marcador) %>%
  mutate_at(vars(segregacao), ~case_when(.=="A.H.B" ~ "Codominante",
                                         T ~ "Dominante")) %>%
  group_by(segregacao, p.val.adj.signif) %>%
  count() %>%
  group_by(segregacao) %>%
  mutate(porc = n/sum(n))

# Total, codominants: 29/213 (13.6%) w/ p < 0,001
# Map,  codominants:  22/172 (12.8%) w/ p < 0,001

# Total, dominants:  10/205  (4.9%) w/ p < 0,001
# Map,  dominants:    6/116  (5.2%) w/ p < 0,001

mapa_antigo %>% 
  left_join(freqs) %>%
  mutate_at(vars(p.val.adj.signif), str_replace, "\\*+", "*") %>%
  group_by(LG, segregacao, p.val.adj.signif) %>%
  count() %>%
  pivot_wider(names_from=segregacao, values_from=n)

mapa_novo %>% 
  left_join(freqs) %>%
  mutate_at(vars(p.val.adj.signif), str_replace, "\\*+", "*") %>%
  group_by(LG, segregacao, p.val.adj.signif) %>%
  count() %>%
  pivot_wider(names_from=segregacao, values_from=n)

mapa_novo %>% 
  left_join(freqs) %>%
  mutate_at(vars(A, H, B, D, C), str_remove, "^[:digit:]+ \\(") %>%
  mutate_at(vars(A, H, B, D, C), str_replace, ",", ".") %>%
  mutate_at(vars(A, H, B, D, C), str_remove, "%\\)$") %>%
  mutate_at(vars(A, H, B, D, C), as.numeric) %>%
  filter(p.val.adj.signif != "ns") %>%
  mutate(fav = case_when(segregacao=="A.H.B" & B > 25 & A < 25 ~ "GG",
                         segregacao=="A.H.B" & A > 25 & B < 25 ~ "NN",
                         segregacao=="C.A" & A > 1/3 ~ "NN",
                         segregacao=="D.B" & B > 1/3 ~ "GG")) %>%
  mutate_at(vars(segregacao), ~case_when(.=="A.H.B" ~ "Codominante",
                                         T ~ "Dominante")) %>%
  group_by(segregacao, fav) %>%
  count() %>%
  group_by(segregacao) %>%
  mutate(porc = n/sum(n))

# Codominants: 24/33 (72,7%) favouring GG (16 w/ p < 0,001); 
#               7/33 (21,2%) favouring NN  (5 w/ p < 0,001)

# Dominants:    6 (75%) GG (4 w/ p < 0,001);
#               2 (25%) NN (2 w/ p < 0,001)

mapa_novo %>% 
  left_join(freqs) %>%
  mutate_at(vars(A, H, B, D, C), str_remove, "^[:digit:]+ \\(") %>%
  mutate_at(vars(A, H, B, D, C), str_replace, ",", ".") %>%
  mutate_at(vars(A, H, B, D, C), str_remove, "%\\)$") %>%
  mutate_at(vars(A, H, B, D, C), as.numeric) %>%
  filter(p.val.adj.signif != "ns") %>%
  mutate(fav = case_when(segregacao=="A.H.B" & B > 25 & A < 25 ~ "GG",
                         segregacao=="A.H.B" & A > 25 & B < 25 ~ "NN",
                         segregacao=="C.A" & A > 1/3 ~ "NN",
                         segregacao=="D.B" & B > 1/3 ~ "GG")) %>%
  group_by(segregacao, fav, p.val.adj.signif) %>%
  count() %>%
  pivot_wider(names_from=p.val.adj.signif, values_from=n)
