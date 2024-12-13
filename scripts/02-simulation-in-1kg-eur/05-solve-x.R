#!/usr/bin/env Rscript
# library(lpSolve)
# library(linprog)
library(data.table)
library(parallel)
library(BEDMatrix)

## get random seeds
get.seeds <- function(seed, size, maxs=999999){
  set.seed(seed)
  sample(1:maxs, size=size)
}

compareX <- function(Xhat, X){
  M <- ncol(Xhat)
  acc <- matrix(NA, nrow=M, ncol=3)
  acc[, 1] <- as.vector(colSums(Xhat == X, na.rm=TRUE))
  acc[, 2] <- as.vector(colSums(Xhat != X, na.rm=TRUE))
  acc[, 3] <- acc[,1]/(acc[, 1] + acc[, 2])
  return(acc)
}

solveX <- function(beta, Y=NULL, err=0.001, boundX=TRUE){
  NT <- nrow(Y)
  N <- ncol(Y)
  beta <- matrix(beta, ncol=1)
  if(sum(is.na(beta)) > 0 | sum(!is.finite(beta)) > 0 | sum(is.nan(beta)) > 0){
    return(matrix(rep(NA, N), ncol=1))
  }

  cv <- rep(1, N)
  if(length(err)==1){
    bv <- c((1+err) * beta, (1-err)*beta)
  }else{
    bv <- c(beta-err, beta+err)
  }
  cd <- c(bv[1:NT] > beta, bv[(NT+1):(2*NT)] > beta)
  cd <- ifelse(cd, '<=', '>=')

  Amat <- rbind(Y, Y)

  ## 0 <= x <= 2
  if(boundX){
    bv <- c(bv, rep(0, N), rep(2, N))
    cd <- c(cd, rep('>=', N), rep('<=', N))
    Amat <- rbind(Amat, diag(rep(1, N)), diag(rep(1, N)))
  }

  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
  res <- linprog::solveLP(cvec=cv, bvec=bv, Amat=Amat, maximum=FALSE, const.dir=cd, lpSolve=TRUE, solve.dual=FALSE)
  # res <- lpSolve::lp(direction='min', objective.in=cv, const.mat=Amat, const.dir=cd, const.rhs=bv)
  # res <- Rglpk::Rglpk_solve_LP(obj=cv, mat=Amat, dir=cd, rhs=bv, bounds=list(upper=rep(2, N)), types='I', max=FALSE, control=list(verbose=TRUE))
  # res <- Rsymphony::Rsymphony_solve_LP(
  #   obj=cv, mat=Amat, dir=cd, rhs=bv, bounds=list(lower=rep(0, N), upper=rep(2, N)), types='I', max=FALSE,
  #   verbosity = 0, time_limit = 20, node_limit = -1, gap_limit = -1, first_feasible = FALSE, write_lp = FALSE, write_mps = FALSE)

  Xhat <- round(pmax(pmin(as.vector(res$solution), 2), 0), 0)
  # print(table(Xhat))
  # cat('-------\n')

  return(Xhat)
}

chunk.acc <- function(s, size=100, beta=NULL, Xtrue=NULL, Y=Y, err=0.001, boundX=TRUE){
  e <- s + size - 1
  e <- min(e, ncol(beta))
  beta <- beta[, s:e]
  Xtrue <- Xtrue[, s:e]
  acc <- matrix(NA, nrow=ncol(beta), ncol=3)
  for(i in 1:ncol(beta)){
    Xhat <- solveX(beta[, i], Y=Y, err=err, boundX=boundX)
    acc[i, 1] <- sum(Xhat == Xtrue[, i], na.rm=TRUE)
    acc[i, 2] <- sum(Xhat != Xtrue[, i], na.rm=TRUE)
  }
  acc[, 3] <- acc[,1]/(acc[, 1] + acc[, 2])
  cat(s, e, '\n')
  return(acc)
}

calc.acc <- function(Y, beta, Xtrue, ncpus=10, size=100){
  # Y (sample x phe), beta (snp x phe)
  if(is.null(nrow(beta))) beta <- matrix(beta, ncol=1)
  err <- sapply(Y, function(x){2*min(abs(x))})
  Y <- t(Y) # sample x user
  N <- ncol(Y)
  M <- nrow(beta)

  idx <- c(1:M)
  # idx <- c(1:100)
  beta <- as.data.frame(t(beta[idx, ]))
  Xtrue <- Xtrue[, idx]

  # Xhat <- sapply(beta, solveX, Y=Y, err=0.001, boundX=TRUE)
  # acc <- compareX(Xhat, Xtrue)

  # acc <- chunk.acc(1, 5, beta=beta, Xtrue=Xtrue, Y=Y, err=0.001, boundX=TRUE)
  ss <- seq(1, ncol(beta), by=size)
  # acc <- do.call('rbind', lapply(ss, chunk.acc, size=size, beta=beta, Xtrue=Xtrue, Y=Y, err=0.001, boundX=TRUE))
  clust <- makeCluster(ncpus, outfile='')
  clusterExport(clust, varlist=c('solveX'))
  acc <- do.call('rbind', parLapply(clust, ss, chunk.acc, size=size, beta=beta, Xtrue=Xtrue, Y=Y, err=0.001, boundX=TRUE))
  stopCluster(clust)

  return(acc)
}

rand.phes <- function(seed, nreps, rn, np, N){
  nt <- round(N * rn)
  seeds <- get.seeds(seed, size=nreps)
  idx <- matrix(NA, nrow=nt, ncol=nreps)
  for(i in 1:nreps){
    set.seed(seeds[i])
    idx[, i] <- sort(sample(np, size=nt, replace=FALSE))
  }
  return(idx)
}

nreps <- 5
seed <- 645372
ncpus <- 40
size <- 50
overwrite <- FALSE
N <- fread('data/bed/data.fam', header=F) %>% nrow
RN <- c(1/N, 0.01, 0.05, 0.1, 0.5, 0.85, 1.0)
seeds <- get.seeds(seed, length(RN))

sim.dir <- './simulation'
for(data.file in list.files(sim.dir, pattern='*.RData')){
  data.file <- paste0(sim.dir, '/', data.file)
  out.dir <- gsub('\\.RData', '', data.file)
  load(data.file)
  for(i in 1:length(RN)){
    acc.dir <- sprintf('%s/nt-%d/', out.dir, round(RN[i]*N))
    if(!dir.exists(acc.dir)) dir.create(acc.dir, recursive=TRUE)
    idx <- rand.phes(seeds[i], nreps, RN[i], dim(y)[2], N)
    save(idx, sel.snps, file=paste0(out.dir, '/idx.RData'))
    for(irep in 1:nreps){
      acc.file <- paste0(acc.dir, '/', irep, '.tsv')
      if(file.exists(acc.file) & !overwrite) next
      phes <- idx[, irep]
      acc <- calc.acc(y[, phes], b[, phes], X0, ncpus=ncpus, size=size)
      write.table(acc, file=acc.file, sep='\t', quote=F, row.names=F, col.names=F)
    }
  }
}
