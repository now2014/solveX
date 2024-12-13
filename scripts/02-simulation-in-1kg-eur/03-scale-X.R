#!/usr/bin/env Rscript
require(BEDMatrix)

scale.bed <- function(bed.file, out.dir, overwrite=FALSE){
  X0.file <- paste0(out.dir, '/X0.RData')
  X.file <- paste0(out.dir, '/X.RData')
  if(file.exists(X0.file) & file.exists(X.file) & !overwrite) return(NULL)
  X0 <- as.matrix(BEDMatrix(bed.file))
  X0[is.na(X0)] <- 0
  X <- scale(X0)
  save(X0, file=X0.file)
  save(X, file=X.file)
}

bed.file <- 'data/bed/data.bed'
out.dir <- 'data'
scale.bed(bed.file, out.dir, FALSE)
