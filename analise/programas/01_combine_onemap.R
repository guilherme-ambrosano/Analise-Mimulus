# Arrumando a funcao combine_onemap que esta no CRAN

# Se os dois conjuntos de dados apresentam os mesmos dados fenotipicos,
# a funcao combine_onemap do CRAN dá pau

combine_onemap <- function (...) 
{
  onemap.objs <- list(...)
  n.objs <- length(onemap.objs)
  if (!n.objs) {
    stop("You must provide a list of OneMap objects as input.")
  }
  for (i in 1:n.objs) {
    if (!is(onemap.objs[[i]], "onemap")) 
      stop("All objects must be of class 'onemap'.")
  }
  if (n.objs == 1) {
    stop("Nothing to merge.")
  }
  crosstype <- class(onemap.objs[[1]])[2]
  for (i in 2:n.objs) {
    if (!is(onemap.objs[[i]], crosstype)) 
      stop("All objects must be of the same cross type.")
  }
  n.mar <- 0
  n.ind <- 0
  sampleIDs <- NULL
  labs.phe <- NULL
  sampleID.flag <- FALSE
  for (i in 1:n.objs) {
    n.mar <- n.mar + onemap.objs[[i]]$n.mar
    # Nao somar o numero de caracteristicas fenotipicas dos conjuntos
    # (podem ser as mesmas em conjuntos diferentes)
    labs.phe <- c(labs.phe, colnames(onemap.objs[[i]]$pheno))
    cur.sampleIDs <- rownames(onemap.objs[[i]]$geno)
    sampleIDs <- unique(c(sampleIDs, cur.sampleIDs))
    if (is.null(cur.sampleIDs)) {
      sampleID.flag <- TRUE
    }
  }
  # Pegar o numero de caracteristicas fenotipicas total de todos os conjuntos
  n.phe <- length(unique(labs.phe))
  if (sampleID.flag) {
    n.ind <- onemap.objs[[1]]$n.ind
    for (i in 2:n.objs) {
      if (onemap.objs[[i]]$n.ind != n.ind) 
        stop("Sample IDs are missing in at least one dataset. All objects must contain the same number of individuals and in the same order.")
    }
  }
  else {
    n.ind <- length(sampleIDs)
  }
  geno <- matrix(0, nrow = n.ind, ncol = n.mar)
  colnames(geno) <- rep(NA, n.mar)
  if (!sampleID.flag) {
    rownames(geno) <- sampleIDs
  }
  segr.type <- rep(NA, n.mar)
  segr.type.num <- rep(NA, n.mar)
  CHROM <- rep(NA, n.mar)
  POS <- rep(NA, n.mar)
  if (n.phe) {
    pheno <- matrix(NA, nrow = n.ind, ncol = n.phe)
    colnames(pheno) <- seq(1, n.phe)
  }
  else {
    pheno <- NULL
  }
  mrk.start <- 1
  phe.start <- 1
  for (i in 1:n.objs) {
    cur.n.mar <- onemap.objs[[i]]$n.mar
    mrk.end <- mrk.start + cur.n.mar - 1
    if (sampleID.flag) {
      ind.matches <- 1:n.ind
    }
    else {
      ind.matches <- match(rownames(onemap.objs[[i]]$geno), 
                           rownames(geno))
    }
    geno[ind.matches, mrk.start:mrk.end] <- onemap.objs[[i]]$geno
    colnames(geno)[mrk.start:mrk.end] <- colnames(onemap.objs[[i]]$geno)
    segr.type[mrk.start:mrk.end] <- onemap.objs[[i]]$segr.type
    segr.type.num[mrk.start:mrk.end] <- onemap.objs[[i]]$segr.type.num
    if (!is.null(onemap.objs[[i]]$CHROM)) {
      CHROM[mrk.start:mrk.end] <- onemap.objs[[i]]$CHROM
    }
    if (!is.null(onemap.objs[[i]]$POS)) {
      POS[mrk.start:mrk.end] <- onemap.objs[[i]]$POS
    }
    
    labs.pheno <- colnames(onemap.objs[[i]]$pheno)
    # cur.pheno = matriz com as caracteristicas fenotipicas desse conjunto
    # que ainda nao foram adicionadas a matriz pheno
    cur.pheno <- onemap.objs[[i]]$pheno[,!(labs.pheno %in% colnames(pheno)), drop=F]
    cur.n.phe <- ncol(cur.pheno)
    phe.end <- phe.start + cur.n.phe - 1
    
    if (cur.n.phe) {
      pheno[ind.matches, phe.start:phe.end] <- cur.pheno
      colnames(pheno)[phe.start:phe.end] <- colnames(cur.pheno)
    }
    mrk.start <- mrk.start + cur.n.mar
    phe.start <- phe.start + cur.n.phe
  }
  if (anyDuplicated(colnames(geno))) {
    warning("Duplicate marker names found. Please check.")
  }
  if (all(is.na(CHROM))) {
    CHROM <- NULL
  }
  if (all(is.na(POS))) {
    POS <- NULL
  }
  input <- "combined"
  structure(list(geno = geno, n.ind = n.ind, n.mar = n.mar, 
                 segr.type = segr.type, segr.type.num = segr.type.num, 
                 n.phe = n.phe, pheno = pheno, CHROM = CHROM, POS = POS, 
                 input = input), class = c("onemap", crosstype))
}
