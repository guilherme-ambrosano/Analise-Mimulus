library(tidyverse)
library(rstatix)
library(janitor)
library(scales)

# Algumas funcoes para facilitar a analise descritiva

tabela_frequencias <- function(onemap.obj) {
  sapply(apply(onemap.obj$geno, 2, table), . %>%
           data.frame() %>%
           mutate_at(vars(Var1), factor, levels=seq(0, 5)) %>%
           complete(Var1, fill=list(Freq=NA)) %>%
           pull(Freq)) %>%
    data.frame() %>%
    t() %>%
    data.frame() %>%
    rownames_to_column("marcador") %>%
    as_tibble() %>%
    rename(Perdidos=X1, A=X2, H=X3, B=X4, D=X5, C=X6) %>%
    mutate(segregacao=onemap.obj$segr.type) %>%
    group_by(segregacao, marcador) %>%
    adorn_percentages("row") %>%
    adorn_pct_formatting() %>%
    adorn_ns("front") %>%
    ungroup() %>%
    mutate_at(vars(-marcador, -segregacao), str_replace_all, " +", " ") %>%
    mutate_at(vars(-marcador, -segregacao), str_replace, "\\.", ",") %>%
    mutate_at(vars(-marcador, -segregacao), str_remove, "^ ") %>%
    mutate_at(vars(-marcador, -segregacao), ~ifelse(.=="NA (-)", "\u2012", .))
}

tabela_segregacao <- function(onemap.obj) {
  
  testes_segregacao <- test_segregation(onemap.obj)
  t(testes_segregacao$Results.of.tests) %>%
    as_tibble() %>%
    mutate_all(unlist) %>%
    mutate(marcador = testes_segregacao$Marker) %>%
    left_join(tabela_frequencias(onemap.obj),
              by="marcador") %>%
    arrange(segregacao) %>%
    adjust_pvalue(p.col="p.val", method="bonferroni") %>%
    add_significance(p.col="p.val.adj") %>%
    mutate_at(vars(p.val.adj), pvalue, acc=0.0001, dec=",") %>%
    mutate_at(vars(qui.quad), number, acc=0.0001, dec=",") %>%
    mutate_at(vars(perc.genot), ~./100) %>%
    mutate_at(vars(perc.genot), percent, acc=0.1, dec=",") %>%
    select(marcador, segregacao, Hypothesis, perc.genot,
           Perdidos, A, H, B, D, C,
           qui.quad, p.val.adj, p.val.adj.signif)
}
