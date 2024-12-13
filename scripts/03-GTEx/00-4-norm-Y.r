#!/usr/bin/env Rscript

setwd('/share/home/lanao/projects/06-solveX/GTEx/v7/04-solve')

library(data.table)
library(tidyverse)

norm.pheno <- function(x, inverse.gaussian=FALSE){
  # replace NA with mean
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  if(inverse.gaussian){
    x <- qnorm((rank(x, na.last = 'keep') - 0.5)/sum(!is.na(x)))
  }else{
    x <- (x - mean(x))/sd(x) # normalize to N(0, 1)
  }
  return(x)
}

Y <- fread('../GTEx_Analysis_v7_eQTL_expression_matrices/Whole_Blood.v7.normalized_expression.bed.gz')
covar <- fread('../GTEx_Analysis_v7_eQTL_covariates/Whole_Blood.v7.covariates.txt')


y <- Y %>%
  column_to_rownames('gene_id') %>%
  select(-c(1:3)) %>%
  t() %>%
  as.data.frame() %>%
  lapply(FUN=norm.pheno, inverse.gaussian=TRUE) %>%
  do.call(what=cbind)
row.names(y) <- colnames(Y)[-c(1:4)]

covs <- covar %>%
  column_to_rownames('ID') %>%
  t() %>%
  as.data.frame() %>%
  lapply(FUN=norm.pheno, inverse.gaussian=FALSE) %>%
  do.call(what=cbind)
row.names(covs) <- colnames(covar)[-1]
covs <- covs[row.names(y), ]

# residualize
res.y <- function(i.gene, y, covs){
  model <- lm(y[, i.gene] ~ covs)
  return(residuals(model))
}
Y.norm <- lapply(seq_len(ncol(y)), res.y, y=y, covs=covs) %>%
  do.call(what=cbind) %>%
  as.data.frame() %>%
  set_names(colnames(y))

save(Y.norm, file = 'data/Y-norm.rda')
