#!/usr/bin/env Rscript

setwd('/share/home/lanao/projects/06-solveX/GTEx/v7/04-solve')
library(BEDMatrix)

calc.acc <- function(X.file, geno, bim, samples){
  X <- fread(X.file, header = F, drop = 2)
  colnames(X) <- c('SNP', samples)
  X <- column_to_rownames(X, 'SNP')
  idx <- bim |> filter(SNP %in% rownames(X)) |> pull(i)
  Xtrue <- t(geno[, idx]) |>
    as.data.frame()
  
  kept.snps <- intersect(rownames(X), rownames(Xtrue))
  kept.samples <- intersect(colnames(X), colnames(Xtrue))
  X <- X[kept.snps, kept.samples]
  Xtrue <- Xtrue[kept.snps, kept.samples]
  Xnull <- X * 0

  row.maf <- rowMeans(Xtrue, na.rm=TRUE) |> 
    as.data.frame() |> 
    set_names('MAF') |> 
    mutate(MAF=MAF/2) |>
    rownames_to_column('SNP')
  null.acc <- rowMeans(X == Xnull, na.rm = T) |>
    as.data.frame() |>
    set_names('acc.null') |>
    rownames_to_column('SNP')
  ## row-wise accuracy
  row.acc <- rowMeans(X == Xtrue, na.rm = T) |>
    as.data.frame() |>
    set_names('acc') |>
    rownames_to_column('SNP') |>
    left_join(null.acc, by='SNP') |>
    left_join(row.maf, by='SNP')
  return(row.acc)
}

samples <- fread('data/Y.tsv.gz', header=T, drop=1) |> colnames()
geno <- BEDMatrix::BEDMatrix('../genotype/Whole_Blood_Genotype.bed', simple_names = TRUE)
bim <- fread('../genotype/Whole_Blood_Genotype.bim', header=F, select=2, col.names='SNP') |>
  mutate(i=1:n())

X.dir <- 'data/solved-X'
X.files <- list.files(X.dir, full.names = T, pattern='*.tsv.gz')
df.acc <- parallel::mclapply(
    X.files, calc.acc, geno=geno, bim=bim,
    samples=samples, mc.cores=10) |> 
  do.call(what=rbind)

saveRDS(df.acc, file='data/acc.rds')
