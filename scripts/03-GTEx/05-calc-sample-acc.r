#!/usr/bin/env Rscript

setwd('/share/home/lanao/projects/06-solveX/GTEx/v7/04-solve')
library(BEDMatrix)

calc.acc <- function(X.file, geno, bim, samples){
  X <- fread(X.file, header = F, drop = 2)
  colnames(X) <- c('SNP', samples)
  X <- column_to_rownames(X, 'SNP')
  idx <- bim |> filter(SNP %in% rownames(X)) |> pull(i)
  
  if(length(idx) == 0){
    return(NULL)
  }
  Xtrue <- t(geno[, idx]) |>
    as.data.frame()
  kept.snps <- intersect(rownames(X), rownames(Xtrue))
  kept.samples <- intersect(colnames(X), colnames(Xtrue))
  X <- X[kept.snps, kept.samples]
  Xtrue <- Xtrue[kept.snps, kept.samples]
  return(as.data.frame(X==Xtrue))
}

N <- 369
df <- readRDS('data/acc.rds') |>
  rename('MAF.X'='MAF') |>
  left_join(fread('../data/snp-cnt.tsv'), by='SNP') |>
  left_join(fread('../data/snp-maf.tsv'), by='SNP') |>
  mutate(RN=n.genes/N)
kept.snps <- df |>
  filter(MAF <= 0.1) |>
  filter(RN > 0.2) |>
  pull(SNP)

samples <- fread('data/Y.tsv.gz', header=T, drop=1) |> colnames()
geno <- BEDMatrix::BEDMatrix('../genotype/Whole_Blood_Genotype.bed', simple_names = TRUE)
bim <- fread('../genotype/Whole_Blood_Genotype.bim', header=F, select=2, col.names='SNP') |>
  mutate(i=1:n()) |>
  filter(SNP %in% kept.snps)

X.dir <- 'data/solved-X'
X.files <- list.files(X.dir, full.names = T, pattern='*.tsv.gz')

X.comp <- parallel::mclapply(
    X.files, calc.acc, geno=geno, bim=bim,
    samples=samples, mc.cores=10) |> 
  do.call(what=rbind)

sample.acc <- colSums(X.comp, na.rm=T) / nrow(X.comp)
sample.acc |>
  as.data.frame() |>
  set_names('acc') |>
  mutate(N=n()) |>
  arrange(desc(acc)) |>
  mutate(sample.ratio=(1:n()) / N) |>
  filter(acc > 0.9) |>
  tail()
