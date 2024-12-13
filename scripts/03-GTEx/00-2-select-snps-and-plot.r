#!/usr/bin/env Rscript

# library(BEDMatrix)
# geno <- BEDMatrix('../genotype/Whole_Blood_Genotype.bed', simple_names = T)
# maf <- apply(geno, 2, mean, na.rm=T)
# maf <- maf / 2


setwd('/share/home/lanao/projects/06-solveX/GTEx/v7/04-solve')
df.cnt <- fread('../data/snp-cnt.tsv')
df.maf <- fread('../data/snp-maf.tsv')
N <- 369
dpi <- 300
font.size <- 48

# bim.snps <- fread('genotype/Whole_Blood_Genotype.bim', select=2, header=F) |>
#   pull()

df <- df.cnt |>
  # filter(SNP %in% bim.snps) |>
  left_join(df.maf, by='SNP') |>
  mutate(RN=n.genes/N)


df.p <- df |>
  filter(MAF > 0) |>
  mutate(MAF.bin=cut(
    MAF, include.lowest=F, right=T,
    breaks=c(0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5))) |>
  # mutate(RN.bin=cut(
  #   RN, breaks=seq(0, 0.3, 0.05),
  #   labels=c('0~0.05', '0.05~0.1', '0.1~0.15', '0.15~0.2', '0.2~0.25', '0.25~0.3'))) |>
  mutate(RN.bin=cut(
    RN, include.lowest=F, right=T,
    breaks=seq(0, 0.3, 0.05))) |>
  count(MAF.bin, RN.bin, name='cnt') |>
  mutate(label=scales::comma(cnt)) |>
  mutate(label.color=ifelse(log10(cnt) > 5, '1', '2'))


p <- df.p |>
  ggplot() + vstheme(size=font.size) +
  geom_tile(aes(x=MAF.bin, y=RN.bin, fill=log10(cnt)), color='white', lwd=2) +
  geom_text(aes(x=MAF.bin, y=RN.bin, label=label, color=label.color), size=10) +
  scale_fill_distiller(palette='Blues', direction=1) +
  scale_color_manual(values=c('1'='white', '2'='black')) +
  guides(fill='none', color='none') +
  labs(x='MAF', y='Number of traits / N')

save.png(p, 'png/MAF-RN.png', w=6, h=3, dpi=dpi)


# gap.size <- 10

# df.chrom <- fread('../genotype/Whole_Blood_Genotype.bim', select=c(1, 4), col.names=c('chr', 'pos')) |>
#   slice_max(pos, n=1, by=chr) |>
#   mutate(xe=pos/1e6) |>
#   select(-pos) |>
#   mutate(xe=cumsum(xe) + gap.size) |>
#   mutate(xs=lag(xe, default=0)) |>
#   mutate(ann.pos=(xe+xs)/2) |>
#   mutate(ann=ifelse(chr %in% c(19, 21), '', as.character(chr)))



p <- df |>
  filter(MAF <= 0.05, MAF > 0) |>
  filter(RN > 0.25) |>
  separate(SNP, into=c('chr', 'pos', 'a1', 'a2', 'v'), sep='_') |>
  mutate(pos=as.numeric(pos)/1e6) |>
  mutate(chr=as.numeric(chr)) |>
  arrange(chr, pos) |>
  mutate(chr=paste0('chr', chr)) |>
  mutate(chr=factor(chr, levels=unique(chr))) |>


  ggplot() + vstheme(size=font.size) +
  geom_point(aes(x=pos, y=RN, color=MAF), size=1, alpha=0.1) +
  facet_wrap(.~chr, scales='free_x', ncol=2, strip.position = 'top') +
  scale_color_distiller(palette='Reds', direction=1) +
  labs(x='Position (Mb)', y='Number of traits / N', color='MAF') +
  scale_y_continuous(breaks=c(0.26, 0.28)) +
  theme(
    # strip.background = element_blank(),
    strip.placement = 'outside',
    legend.position = c(0.9, 0),
    legend.justification = c(1, 0.1),
    legend.background = element_blank(),
    legend.key.size = unit(12, 'pt')
  )
  
save.png(
  p, 'png/high-risk-regions.png',
  w=6, h=6, dpi=dpi
)

df |>
  filter(MAF <= 0.05) |>
  filter(RN >= 0.25) |>
  separate(SNP, into=c('chr', 'pos', 'a1', 'a2', 'v'), sep='_') |>
  mutate(pos=as.numeric(pos)) |>
  mutate(chr=as.numeric(chr)) |>
  arrange(chr, pos) |>

  select(chr, pos) |>
  group_by(chr) |>
  mutate(x=pos-min(pos)) |>
  ungroup() |>
  mutate(grp=round(x/10e6)) |>
  group_by(chr, grp) |>
  summarize(
    start=min(pos), end=max(pos), cnt=n(),
    .groups='drop'
  ) |>
  mutate(start=scales::comma(start), end=scales::comma(end)) |>
  unite(pos, start, end, sep='-') |>
  unite(pos, chr, pos, sep=':') |>
  select(pos, cnt) |>
  print()









out.dir <- './data'
if(!dir.exists(out.dir)) dir.create(out.dir, recursive=TRUE)

df |>
  filter(MAF <= 0.1) |>
  select(SNP) |>
  write.table(
    file=paste0(out.dir, '/selected-snps.tsv'), sep='\t',
    quote=FALSE, row.names=FALSE, col.names=FALSE
  )
