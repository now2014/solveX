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
  maf <- ifelse(maf > 0.5, 1-maf, maf)
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


read.data <- function(tsv.file, w=1e6){
  chrom <- gsub('\\.tsv', '', basename(tsv.file)) |> as.numeric()
  df <- fread(tsv.file) |>
    mutate(POS=floor(POS/w) + 0.5) |>
    group_by(POS) |>
    summarise(
      AFR=mean(maf2rn(AFR)),
      AMR=mean(maf2rn(AMR)),
      EAS=mean(maf2rn(EAS)),
      EUR=mean(maf2rn(EUR)),
      SAS=mean(maf2rn(SAS)),
      .groups='drop'
    ) |>
    mutate(POS=POS*(w/1e6)) |>
    gather(pop, RN, -POS) |>
    mutate(chrom=chrom)
  return(df)
}

data.dir <- 'data/1kg-extracted'
data.files <- list.files(data.dir, pattern='tsv', full.names=T)

df <- lapply(data.files, read.data) |> do.call(what=rbind)
saveRDS(df, file='1kg-RN.rds')

p <- df |>
  ggplot() + vstheme() +
  geom_point(aes(x=POS, y=1-RN, color=pop), alpha=0.3) +
  facet_grid(.~chrom, scales='free_x', space='free_x', switch='x') +
  labs(x='Chromosome', y='Mean leakage risk', color='Population') +
  scale_color_manual(values=col5) +
  guides(color=guide_legend(override.aes = list(alpha=1))) +
  theme(
    legend.position = c(0, 0),
    legend.justification = c(0, 0),
    legend.background = element_rect(fill='transparent'),
    strip.background = element_blank(),
    strip.placement = 'outside',
    strip.text = element_text(size=15),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
