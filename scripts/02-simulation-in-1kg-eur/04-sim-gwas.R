#!/usr/bin/env Rscript
library(data.table)

## get random seeds
get.seeds <- function(seed, size, maxs=999999){
  set.seed(seed)
  sample(1:maxs, size=size)
}

## set causal SNPs
sample.causal <- function(seed, M, ratio=0.1){
  size <- round(M*min(ratio, 1-ratio))
  set.seed(seed)
  idx.causal <- sort(sample(1:M, size=size))
  if(ratio > 0.5) idx.causal <- setdiff(1:M, idx.causal)
  return(idx.causal)
}

## subset SNPs
subset.snps <- function(X0, seed, size=1e4){
  breaks <- c(0, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
  labels <- c('(0,0.005]', '(0.005,0.01]', '(0.01,0.05]', '(0.05,0.1]',
              '(0.1,0.2]', '(0.2,0.3]', '(0.3,0.4]', '(0.4,0.5]')

  N <- dim(X0)[1]
  set.seed(seed)
  seeds <- get.seeds(seed, length(breaks)-1)
  df <- data.frame(maf=colSums(X0) / N / 2) %>% 
    mutate(maf=ifelse(maf > 0.5, 1-maf, maf)) %>% 
    mutate(maf_bin=cut(maf, breaks=breaks, labels=labels,
                       include.lowest = F, right = T)) %>% 
    mutate(isnp=1:n())
  
  sel.snps <- do.call('rbind', lapply(1:length(labels), function(i, seeds, labels, df, size){
      df.sub <- filter(df, maf_bin==labels[i])
      if(nrow(df.sub) <= size) return(df.sub)
      set.seed(seeds[i])
      return(sample_n(df.sub, size=size, replace = FALSE))
    }, seeds=seeds, labels=labels, df=df, size=size)) %>%
    arrange(isnp) %>%
    relocate(isnp)
  return(sel.snps)
}

## simulate phenotype
sim.phe <- function(i, setups, X0, out.dir){
  N <- dim(X0)[1]
  M <- dim(X0)[2]

  setup <- setups[i, ]
  h2 <- sprintf('%.2f', setup$h2)
  ratio <- sprintf('%.2f', setup$ratio)
  out.file <- paste0(out.dir, '/h2_', h2, '_ratio_', ratio, '.RData')
  idx <- sample.causal(setup[['causal.seed']], M, setup[['ratio']])
  Mc <- length(idx)
  # trueA <- rep(0, M)
  # trueA[idx] <- M / length(idx)
  seed <- as.integer(setup[['seed']])
  nphes <- as.integer(setup[['nphes']])

  h2 <- setup[['h2']]
  # sds <- sqrt(h2 / M * trueA)
  sds <- rep(sqrt(h2 / Mc), Mc)
  set.seed(seed)
  b <- t(mapply(rnorm, n=nphes, sd=sds)) ## beta
  y0 <- X0[, idx] %*% b

  sel.snps <- subset.snps(X0, seed, size=1e4)
  X0 <- X0[, sel.snps$isnp]
  X <- scale(X0, scale=FALSE)

  seeds <- get.seeds(seed, nphes)
  # e.sd <- sqrt(1 - h2)
  if(h2==0){
    e.sd <- rep(1, nphes)
  }else{
    e.sd <- sqrt((1-h2)/h2 * apply(y0, 2, var))
  }
  err <- do.call('cbind', lapply(1:nphes,
          function(i, seeds, N, e.sd){
            set.seed(seeds[i]);
            rnorm(N, sd=e.sd)
          }, seeds=seeds, N=N, e.sd=e.sd)) ## error term
  
  y <- scale(y0 + err) ## phenotype

  
  x2 <- colSums(X**2)
  bhat <- (t(X) %*% y) / x2  ## bhat

  # SE of bhat
  se <- 0 * bhat
  for(i in 1:nphes){
    se[, i] <- sqrt(colSums((t(t(X)*bhat[, i]) - y[, i])**2) / (N-2) / x2)
  }
  
  y2 <- crossprod(y) * diag(rep(1, nphes))
  b <- (bhat / ((N-2)*se**2 + bhat**2)) %*% y2 ## b = t(X) %*% y
  b0 <- t(X0) %*% y

  ## save results
  save(X0, y, b, b0, sel.snps, file=out.file)
}

## setup simulation params
setup.sim <- function(N=489, nreps=5, seed=437151, RN=c(1/N, 0.01, 0.05, 0.1, 0.5, 0.85, 1.0)){ 
  nt <- sum(round(RN * N))
  setups <- data.frame(
    h2 = c(0, 0.01, 0.2, 0.5, 0.8),
    ratio = c(0.1, 0.1, 0.1, 0.1, 0.1),
    # nphes = nt * nreps
    nphes = N * 2
  )

  setups$seed <- get.seeds(seed, nrow(setups))
  setups$causal.seed <- seed
  return(setups)
}

load('data/X0.RData')
out.dir <- './simulation'

if(!dir.exists(out.dir)) dir.create(out.dir, recursive=TRUE)
setups <- setup.sim(N=dim(X0)[1], nreps=5)
for(i in 1:nrow(setups)){
  sim.phe(i, setups=setups, X0=X0, out.dir=out.dir)
}
