#!/usr/bin/env Rscript
library(ggsankey)
library(data.table)
library(tidyverse)
library(patchwork)
library(scales)

mycols <- c('#FF6F00FF', '#C71000FF', '#008EA0FF', '#8A4198FF', '#5A9599FF', '#FF6348FF',
            '#84D7E1FF', '#FF95A8FF', '#3D3B25FF', '#ADE2D0FF', '#1A5354FF', '#3F4041FF')
mytheme <- theme_bw() + theme(text=element_text(size=18))

save.png <- function(p, out.png='png/test.png', w=8, h=6, res=100, units='in'){
  out.dir <- dirname(out.png)
  out.dir <- ifelse(out.dir=='', '.', out.dir)
  if(!dir.exists(out.dir)) dir.create(out.dir, recursive=TRUE)
  png(out.png, width=w, height=h, res=res, units=units)
  print(p)
  tmp <- dev.off()
}


sample.causal <- function(seed, M, ratio=0.1){
  size <- round(M*min(ratio, 1-ratio))
  set.seed(seed)
  idx.causal <- sort(sample(1:M, size=size))
  if(ratio > 0.5) idx.causal <- setdiff(1:M, idx.causal)
  return(idx.causal)
}

read.acc <- function(data.file, nreps=5, M=NULL, N=NULL, RN=NULL, sim.dir='./simulation', causal.seed=437151){
  data.file <- paste0(sim.dir, '/', data.file)
  out.dir <- gsub('\\.RData', '', data.file)
  h2 <- as.numeric(gsub('.*h2_(.*)_ratio.*', '\\1', data.file))
  ratio <- as.numeric(gsub('.*ratio_(.*)\\.RData', '\\1', data.file))
  idx.causal <- sample.causal(causal.seed, M, ratio)
  
  breaks <- c(-Inf, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 0.99999999999, 1)
  labels <- c('[0, 0.05]', '(0.05, 0.1]', '(0.1, 0.2]', '(0.2, 0.3]',
            '(0.3, 0.4]', '(0.4, 0.5]', '(0.5, 0.6]', '(0.6, 0.7]', '(0.7, 0.8]',
            '(0.8, 0.9]', '(0.9, 0.95]', '(0.95, 0.99]', '(0.99, 1)', '[1]')

  df <- list()
  RN <- sort(RN)
  idx.file <- paste0(out.dir, '/idx.RData') # idx, sel.snps
  load(idx.file)
  df <- sel.snps %>% 
    mutate(is_causal=isnp %in% idx.causal) %>%
    select(maf_bin, is_causal)
  for(i in 1:length(RN)){
    rn <- ifelse(RN[i]==1/N, '1/N', RN[i])
    acc.dir <- sprintf('%s/nt-%d', out.dir, round(RN[i]*N))
    acc <- do.call('cbind', lapply(1:nreps, function(irep){
      fread(paste0(acc.dir, '/', irep, '.tsv'), select = 3, header=F)
    })) %>%
      rowMeans(., na.rm=TRUE) %>% 
      as.data.frame() %>%
      set_names('acc') %>% 
      mutate(acc=cut(acc, breaks=breaks, labels=labels)) %>% 
      set_names(rn)
  
    df <- cbind(df, acc)
  }
  df$h2 <- h2
  # df$rn <- factor(df$rn, levels=unique(df$rn))
  return(df)
}

nreps <- 5
N <- fread('data/bed/data.fam', header=F) %>% nrow
M <- fread('data/bed/data.bim', header=F) %>% nrow
RN <- c(1/N, 0.01, 0.05, 0.1, 0.5, 0.85, 1.0)
sim.dir <- './simulation'

df <- do.call('rbind', lapply(list.files(sim.dir, pattern='*.RData'), read.acc, M=M, N=N, RN=RN, nreps=nreps))