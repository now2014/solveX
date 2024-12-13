#!/usr/bin/env Rscript

wkdir <- '/share/home/lanao/projects/06-solveX/1kG'
setwd(wkdir)

read.data <- function(tsv.file, maf.cutoff=0.05){
  chrom <- gsub('\\.tsv', '', basename(tsv.file)) |> as.numeric()
  cat('Reading', tsv.file, '...\n')
  df <- fread(tsv.file) |>
    mutate(i=1:n()) |>
    # keep only biallelic SNPs
    filter(!str_detect(REF, ','), !str_detect(ALT, ',')) |>
    select(-REF, -ALT) |>
    gather(pop, maf, -POS, -i) |>
    mutate(maf = ifelse(maf > 0.5, 1-maf, maf)) |>
    filter(maf >= maf.cutoff) |>
    mutate(chrom=chrom) |>
    select(i, chrom, POS, pop, maf)
  return(df)
}

# data.dir <- 'data/1kg-extracted'
# data.files <- list.files(data.dir, pattern='tsv', full.names=T)
# df <- lapply(data.files, read.data) |> do.call(what=rbind)
# saveRDS(df, file='1kg-MAF-0.05.rds')




sample.pop <- function(cur.pop, df, nsnps=1000, nreps=30, seed=123){
  set.seed(seed)
  seeds <- sample(1:9999999, nreps)
  df <- df |> filter(pop == cur.pop)
  n <- nrow(df)
  if(n < nsnps) {
    stop('Not enough SNPs for', cur.pop)
  }
  out.df <- NULL
  for(irep in 1:nreps) {
    set.seed(seeds[irep])
    j <- sample(1:n, nsnps, replace=FALSE)
    out.df <- rbind(out.df, mutate(df[j, ], irep=irep))
  }
  return(out.df)
}

out.dir <- 'data/04-sample-identification'
nsnps <- 1000
nreps <- 30
seed <- 4754
df <- readRDS('1kg-MAF-0.05.rds') |>
  arrange(pop, chrom, i)


if(!dir.exists(out.dir)) dir.create(out.dir, recursive=TRUE)
pops <- unique(df$pop)
set.seed(seed)
seeds <- sample(1:9999999, length(pops))
snps.df <- lapply(1:length(pops), function(i) {
  sample.pop(pops[i], df, nsnps=nsnps, nreps=nreps, seed=seeds[i])
}) |> do.call(what=rbind)
row.names(snps.df) <- NULL

saveRDS(snps.df, file=paste0(out.dir, '/random-SNPs-MAF-0.05.rds'))
snps.df |>
  distinct(chrom, i, POS) |>
  arrange(chrom, i) |>
  write.table(file=paste0(out.dir, '/idx.random-SNPs-MAF-0.05.tsv'), sep='\t', quote=FALSE, row.names=FALSE)
