#!/usr/bin/env Rscript


wkdir <- '/share/home/lanao/projects/06-solveX/1kG'
setwd(wkdir)

read.geno <- function(tsv.file){
  fread(tsv.file, header=T) |>
    select(-ID, -QUAL, -FILTER, -FORMAT, -INFO, -REF, -ALT)
}

geno.files <- list.files('data/04-sample-identification/vcf-random-SNPs-MAF-0.05', pattern='tsv', full.names=T)
geno <- lapply(geno.files, read.geno) |> do.call(what=rbind)

sample2pop <- fread('./data/1kg/phase3-GRCh37/integrated_call_male_samples_v3.20130502.ALL.panel',
  header=F, skip=1, col.names=c('sample', 'pop', 'super_pop', 'gender'))
related.samples <- fread('./data/1kg/phase3-GRCh37/20140625_related_individuals.txt', header=T, select = 1) |> pull()
kept.samples <- intersect(sample2pop$sample, colnames(geno)) |> setdiff(related.samples)

geno.df <- geno |> select(all_of(c('i', 'CHROM', 'POS', kept.samples)))


saveRDS(geno.df, 'data/04-sample-identification/geno-random-SNPs-MAF-0.05.rds')
saveRDS(sample2pop |> filter(sample %in% kept.samples), 'data/04-sample-identification/sample2pop.rds')
