#!/usr/bin/env Rscript

file.list <- paste0('1kg_eur/bfile.list')
bprefix <- '1kg_eur/1000G_EUR_Phase3_plink/1000G.EUR.QC.'
paste0(bprefix, 2:22) %>% 
  write.table(., file=file.list, sep='\t', quote=F, col.names=F, row.names=F)

cmd <- paste0(
  'plink --bfile ', bprefix, '1',
  ' --merge-list ', file.list,
  ' --silent --make-bed --out 1kg_eur/data'
)
system(cmd)

cmd <- 'plink --bfile 1kg_eur/data --freq --out 1kg_eur/data-freq'
system(cmd)
