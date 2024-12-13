#!/usr/bin/env Rscript

wkdir <- '/share/home/lanao/projects/06-solveX/1kG'
setwd(wkdir)

mycols <- c('#FF6F00FF', '#C71000FF', '#008EA0FF', '#8A4198FF', '#5A9599FF', '#FF6348FF',
            '#84D7E1FF', '#FF95A8FF', '#3D3B25FF', '#ADE2D0FF', '#1A5354FF', '#3F4041FF')
col11 <- c('#67001F', '#B2182B', '#D6604D', '#F4A582', '#FDDBC7', '#FFFFBF',
           '#D1E5F0', '#92C5DE', '#4393C3', '#2166AC', '#053061')
col5 <- c('#FF6F00FF', '#C71000FF', '#008EA0FF', '#8A4198FF', '#053061')
mytheme <- theme_bw() + theme(text=element_text(size=18))



#     grp        b      b.se
# 1 ≥0.99 5.362494 0.3915804
# 2 <0.99 5.747249 0.4844968

maf2rn <- function(maf, b=5.362494, b.se=0.3915804, ci=F) {
  # maf <- ifelse(maf > 0.5, 1-maf, maf)
  RN <- 2 / (1 + exp(-b*maf)) - 1
  if(ci) {
    b.lwr <- b - 1.96*b.se
    b.upr <- b + 1.96*b.se
    RN.lwr <- 2 / (1 + exp(-b.lwr*maf)) - 1
    RN.upr <- 2 / (1 + exp(-b.upr*maf)) - 1
  }
  if(ci) {
    return(list(RN=RN, RN.lwr=RN.lwr, RN.upr=RN.upr))
  } else {
    return(RN)
  }
}


read.data <- function(tsv.file, maf.cutoff=0.05){
  chrom <- gsub('\\.tsv', '', basename(tsv.file)) |> as.numeric()
  df <- fread(tsv.file) |>
    select(-REF, -ALT) |>
    gather(pop, maf, -POS) |>
    mutate(maf = ifelse(maf > 0.5, 1-maf, maf)) |>
    filter(maf >= maf.cutoff) |>
    mutate(RN=maf2rn(maf)) |>
    mutate(chrom=chrom) |>
    select(chrom, POS, pop, RN)
  return(df)
}

data.dir <- 'data/1kg-extracted'
data.files <- list.files(data.dir, pattern='tsv', full.names=T)

df <- lapply(data.files, read.data) |> do.call(what=rbind)
saveRDS(df, file='1kg-RN.rds')




p <- df |>
  sample_frac(0.01, replace=FALSE) |>
  ggplot() + vstheme(size=48) +
  geom_point(aes(x=POS/1e6, y=RN, color=pop), alpha=0.1, size=0.1) +
  facet_grid(pop~chrom, scales='free_x', space='free_x', switch='x') +
  labs(x='Chromosome', y='R / N', color='Population') +
  scale_color_manual(values=col5) +
  guides(color=guide_legend(override.aes = list(alpha=1))) +
  theme(
    # legend.position = 'inside',
    # legend.position.inside = c(0, 0),
    # legend.justification = c(0, 0),
    legend.position = 'none',
    legend.background = element_rect(fill='transparent'),
    strip.background = element_blank(),
    strip.placement = 'outside',
    strip.text.x = element_text(size=15),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

save.png(p, '1kg-RN.png', 20, 10)



df.cnt <- df |>
  mutate(RN=cut(RN, breaks=seq(0, 1, by=0.1))) |>
  count(pop, chrom, RN, name='cnt')

saveRDS(df.cnt, file='1kg-RN-count.rds')

df.cnt |>
  group_by(pop, chrom) |>
  mutate(total=sum(cnt), ratio=cnt/total) |>
  ungroup() |>
  ggplot() + vstheme(size=18) +
  geom_bar(aes(x=factor(chrom), y=ratio, fill=RN), color='black', stat='identity') +
  facet_grid(pop~.) +
  labs(x='Chromosome', y='Proportion of variants', fill='R / N') +
  ggsci::scale_fill_jama() +
  # scale_fill_distiller(palette='RdBu', direction=-1) +
  theme(
    # legend.position = 'bottom',
    # legend.direction = 'horizontal',
    strip.background = element_blank()
  )



col8 <- c('#67001F', '#B2182B', '#F4A582', 'grey',
          '#92C5DE', '#4393C3', '#2166AC', '#053061')
df.cnt |>
  count(pop, RN, name='cnt', wt=cnt) |>
  arrange(pop, RN) |>
  add_count(pop, wt=cnt, name='total') |>
  mutate(ratio=cnt/total) |>
  group_by(pop) |>
  mutate(y=cumsum(ratio) - 0.5*ratio, ys=y-0.5*ratio, ye=y+0.5*ratio) |>
  ungroup() |>
  mutate(x=as.numeric(factor(pop))) |>
  mutate(label=scales::percent(ratio, accuracy=0.01)) |>
  mutate(label=paste0(label, '\n', scales::comma(cnt), '')) |>
  ggplot() + vstheme(size=18) +
  # geom_bar(aes(x=pop, y=cnt, fill=RN), color='black', stat='identity', position='fill') +
  geom_rect(aes(xmin=x-0.48, xmax=x+0.48, ymin=ys, ymax=ye, fill=RN), color='black', lwd=0.05) +
  geom_text(aes(x=x, y=y, label=label), color='white', size=5, fontface='bold') +
  labs(x='Population', y='Proportion of variants', fill='R / N') +
  scale_x_continuous(breaks=1:length(unique(df.cnt$pop)), labels=sort(unique(df.cnt$pop))) +
  scale_fill_manual(values=col8) +
  theme(
    # legend.position = 'bottom',
    # legend.direction = 'horizontal',
    strip.background = element_blank()
  )
