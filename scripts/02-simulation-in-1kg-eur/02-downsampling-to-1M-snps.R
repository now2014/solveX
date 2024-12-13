#!/usr/bin/env Rscript

mycols <- c('#FF6F00FF', '#C71000FF', '#008EA0FF', '#8A4198FF', '#5A9599FF', '#FF6348FF',
            '#84D7E1FF', '#FF95A8FF', '#3D3B25FF', '#ADE2D0FF', '#1A5354FF', '#3F4041FF')

vstheme <- function (w = 8, h = 6, size = 18){
    options(repr.plot.width = w, repr.plot.height = h)
    theme_bw() + theme(text=element_text(size=size), panel.grid=element_blank(),
    panel.border=element_blank(), axis.line=element_line(color='black'))
}

seed <- 5951528
size <- 1e6
out.dir <- 'data/bed'
png.dir <- 'png'

if(!dir.exists(out.dir)) dir.create(out.dir, recursive=TRUE)
if(!dir.exists(png.dir)) dir.create(png.dir, recursive=TRUE)
df <- fread('1kg_eur/data-freq.frq')

set.seed(seed)
snps <- sort(sample(1:nrow(df), size=size, replace=F))
df <- df[snps, ]
write.table(df$SNP, file=paste0(out.dir, '/kept.snps'), col.names=F, row.names=F, sep='\t', quote=F)

title <- paste0('Number of total SNPs = ', comma(nrow(df)))
title <- NULL
p <- ggplot(df) + vstheme(6, 4) +
  geom_histogram(aes(x=MAF), bins=30, fill=mycols[2], col='black', alpha=0.5) +
  labs(x='MAF', y='Number of SNPs', title=title)

png(paste0(png.dir, '/00-MAF.png'), res=100, units='in', w=6, h=4)
print(p)
tmp <- dev.off()

system(paste0(
  'plink --bfile 1kg_eur/data --extract ', paste0(out.dir, '/kept.snps'),
  ' --silent --make-bed --out ', out.dir, '/data'
))
