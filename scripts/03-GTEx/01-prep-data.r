#!/usr/bin/env Rscript
setwd('/share/home/lanao/projects/06-solveX/GTEx/v7/04-solve')

library(data.table)
library(tidyverse)

out.dir <- './data'
if(!dir.exists(out.dir)) dir.create(out.dir, recursive = TRUE)
load('data/Y-norm.rda') # Y.norm
gwas.df <- fread('data/gwas.tsv') %>%
  select(variant_id, gene_id, slope, slope_se) %>%
  set_names(c('SNP', 'gene', 'BETA', 'SE')) %>%
  add_count(SNP, name='n.genes') %>%
  mutate(RN=n.genes/nrow(Y.norm))



Y2.df <- colSums(Y.norm**2) %>%
  as.data.frame() %>%
  rownames_to_column('gene') %>%
  set_names('gene', 'Y2')

# BETA = BETA / ((N-2)*SE**2 + BETA**2) * Y2

N <- nrow(Y.norm)

chunk.size <- 5e4
snps <- unique(gwas.df$SNP)
i.chunk <- 1
for(i in seq(1, length(snps), by=chunk.size)){
  snp.chunk <- snps[i:min(i+chunk.size-1, length(snps))]
  
  gwas.b <- gwas.df %>%
    filter(SNP %in% snp.chunk) %>%
    left_join(Y2.df, by='gene') %>%
    mutate(b = BETA / ((N-2)*SE**2 + BETA**2) * Y2) %>%
    select(SNP, gene, b) %>%
    spread(gene, b)
  out.file <- file.path(out.dir, paste0('b-', i.chunk, '.tsv.gz'))
  fwrite(gwas.b, out.file, na = 'NA',
    sep = '\t', quote = FALSE, col.names = TRUE, compress='gzip')
  i.chunk <- i.chunk + 1
}



## save normalized Y
Y.norm %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column('gene') %>%
  # filter(gene %in% colnames(gwas.b)) %>%
  fwrite(file.path(out.dir, 'Y.tsv.gz'),
    sep = '\t', quote = FALSE, col.names = TRUE, compress='gzip')
