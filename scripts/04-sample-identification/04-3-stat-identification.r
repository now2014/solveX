
wkdir <- '/share/home/lanao/projects/06-solveX/1kG'
setwd(wkdir)

#     grp        b      b.se
# 1 >0.99 5.362494 0.3915804
# 2 <=0.99 5.747249 0.4844968
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

calc.rate <- function(cur.pop, snps.df, geno.df, sample.df){
  kept.samples <- sample.df |> filter(super_pop==cur.pop) |> pull(sample)
  snps.df <- snps.df |>
    filter(pop==cur.pop) |>
    arrange(irep, RN, chrom, i) |>
    unique()
  geno.df <- select(geno.df, all_of(c('i', 'CHROM', 'POS', kept.samples)))
  geno.df <- left_join(snps.df, geno.df, by=c('chrom'='CHROM', 'POS', 'i')) |>
    select(-chrom, -POS, -i, -pop) |>
    relocate(irep, maf, RN)
  N <- dim(geno.df)[2] - 3
  rn.cuts <- seq(0, 0.2, by=0.002)
  df <- lapply(rn.cuts,
    function(rn.cut, geno.df){
      cnt.df <- filter(geno.df, RN <= rn.cut) |>
        select(-maf, -RN)
      if(nrow(cnt.df) == 0) return(NULL)
      cnt.df <- cnt.df |>
        group_by(irep) |>
        mutate(snp=1:n()) |>
        ungroup() |>
        gather(sample, geno, -irep, -snp) |>
        arrange(irep, sample, snp) |>
        group_by(irep, sample) |>
        summarise(geno=paste0(geno, collapse=''), .groups='drop') |>
        add_count(irep, name='total') |>
        count(irep, total, geno, name='cnt') |>
        select(-geno) |>
        filter(cnt==1) |>
        count(irep, total, name='cnt') |>
        mutate(rate=cnt/total) |>
        mutate(rn.cut=rn.cut)
      return(cnt.df)
    }, geno.df) |>
    do.call(what=rbind) |>
    select(irep, rn.cut, rate)
  
  df <- expand.grid(rn.cut=rn.cuts, irep=unique(df$irep)) |>
    as.data.frame() |>
    left_join(df, by=c('rn.cut', 'irep')) |>
    mutate(rate=ifelse(is.na(rate), 0, rate)) |>
    mutate(pop=cur.pop, N=N)
  return(df)
}


geno.df <- readRDS('data/04-sample-identification/geno-random-SNPs-MAF-0.05.rds')
snps.df <- readRDS('data/04-sample-identification/random-SNPs-MAF-0.05.rds')
sample.df <- readRDS('data/04-sample-identification/sample2pop.rds')

snps.df$RN <- maf2rn(snps.df$maf)
pops <- unique(sample.df$super_pop)


df <- lapply(pops, calc.rate, snps.df=snps.df, geno.df=geno.df, sample.df=sample.df) |>
  do.call(what=rbind)

p <- df |>
  group_by(pop, N, rn.cut) |>
  summarise(sd.rate=sd(rate), rate=mean(rate), .groups='drop') |>
  filter(rn.cut > 0.12) |>
  mutate(pop=paste0(pop, ' (N = ', N, ')')) |>
  # filter(rn.cut < 0.2) |>
  ggplot() + vstheme() +
  aes(x=rn.cut, y=rate, color=pop, group=pop) +
  geom_line() +
  # geom_errorbar(aes(ymin=rate-sd.rate, ymax=rate+sd.rate),
  #   width=0.001, alpha=0.3) +
  scale_y_continuous(labels=scales::percent_format(accuracy=1)) +
  labs(x='R / N', y='Sample identification rate\nbased on 1,000 random SNPs', color='Population') +
  ggsci::scale_color_jco() +
  theme(
    legend.position.inside = T,
    legend.background = element_blank(),
    legend.position = c(1, 0),
    legend.justification = c(1, 0)
  )


saveRDS(df, 'data/04-sample-identification/sample-identification-rate.rds')
