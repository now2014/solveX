#!/usr/bin/env Rscript

wkdir <- '~/Desktop/solve-X/Figures'
setwd(wkdir)

mycols <- c('#FF6F00FF', '#C71000FF', '#008EA0FF', '#8A4198FF', '#5A9599FF', '#FF6348FF',
            '#84D7E1FF', '#FF95A8FF', '#3D3B25FF', '#ADE2D0FF', '#1A5354FF', '#3F4041FF')
col11 <- c('#67001F', '#B2182B', '#D6604D', '#F4A582', '#FDDBC7', '#FFFFBF',
           '#D1E5F0', '#92C5DE', '#4393C3', '#2166AC', '#053061')
mytheme <- theme_bw() + theme(text=element_text(size=18))
save.png <- function(p, out.png='png/test.png', w=8, h=6, res=100, units='in'){
  out.dir <- dirname(out.png)
  out.dir <- ifelse(out.dir=='', '.', out.dir)
  if(!dir.exists(out.dir)) dir.create(out.dir, recursive=TRUE)
  png(out.png, width=w, height=h, res=res, units=units)
  print(p)
  tmp <- dev.off()
}


fit.RN <- function(df, is.min=TRUE){
  df <- mutate(df, MAF=round(MAF, 3))
  if(is.min){
    df <- df |> slice_min(RN, n=1, by='MAF', with_ties = FALSE)
  }else{
    df <- df |> slice_max(RN, n=1, by='MAF', with_ties = FALSE)
  }
  mymodel <- function(x) {
    ### erf function
    # 2 * pnorm(x * sqrt(2)) - 1
    ### tanh function
    # tanh(x)
    ### logit function
    2 / (1 + exp(-x)) - 1
  }
  model <- nls(RN ~ mymodel(b * MAF),
              data = df, 
              start = list(b=1))
  b <- coef(model)
  b.se <- summary(model)$coefficients[1, 'Std. Error']
  y <- mymodel(b * df$MAF)

  ## calculate R^2
  SSE <- sum((df$RN - y)^2)
  SST <- sum((df$RN - mean(df$RN))^2)
  R2 <- 1 - SSE / SST

  df <- data.frame(MAF=seq(0, 0.5, by=0.01))
  df$RN <- mymodel(b * df$MAF)

  b.lwr <- b - 1.96 * b.se
  b.upr <- b + 1.96 * b.se
  df$RN.low <- mymodel(b.lwr * df$MAF)
  df$RN.high <- mymodel(b.upr * df$MAF)
  df$R2 <- R2
  df$b <- b
  df$b.se <- b.se
  return(df)
}


df.sim <- fread('~/Desktop/solve-X/sim-x02.tsv') |>
  filter(N==1000) |>
  select(acc, maf, Tratio) |>
  rename('RN'='Tratio', 'MAF'='maf') |>
  filter(is.finite(acc)) |>
  mutate(grp=ifelse(acc > 0.99, '>0.99', '≤0.99'))
df.fit <- rbind(
  fit.RN(df.sim |> filter(acc > 0.99)) |> mutate(grp='>0.99'),
  fit.RN(df.sim |> filter(acc <= 0.99), is.min=F) |> mutate(grp='≤0.99')
)

maf2rn <- function(MAF, coef, acc=0.99){
  y <- log(acc / (1 - acc))
  # y = c0 * MAF + c1 * RN + c2 * MAF * RN
  RN <- (y - coef[1] * MAF) / (coef[2] + coef[3] * MAF)
  return(RN)
}



col11 <- c('#67001F', '#B2182B', '#D6604D', '#F4A582', '#FDDBC7', '#FFFFBF',
           '#D1E5F0', '#92C5DE', '#4393C3', '#2166AC', '#053061')
acc.breaks <- c(seq(0, 0.9, by=0.1), 0.95, 0.99, 1)

create.pie.polygon <- function(cur.MAF, cur.RN, df.sim, r=0.02,
  acc.breaks=c(seq(0, 0.9, by=0.1), 0.95, 0.99, 1), nsegs=100){
  ratios <- seq(0, 1, length.out=nsegs)
  # x <- cos(theta)
  # y <- sin(theta)
  df <- df.sim |> filter(RN==cur.RN, MAF==cur.MAF) |>
    mutate(acc.bin=cut(acc, breaks=acc.breaks)) |>
    filter(!is.na(acc.bin)) |>
    mutate(total=n()) |>
    count(total, acc.bin, name='cnt') |>
    mutate(ratio=cnt/total) |>
    arrange(desc(acc.bin)) |>
    mutate(ratio.s=cumsum(ratio) - ratio) |>
    mutate(ratio.e=cumsum(ratio))
  df.poly <- lapply(seq_len(nrow(df)), function(i, df, ratios, r){
    s <- df$ratio.s[i]
    e <- df$ratio.e[i]
    acc <- df$acc.bin[i]
    ratios <- ratios[which(ratios >= s & ratios <= e)]
    ratios <- c(s, ratios, e) |> unique()
    theta <- 2 * pi * ratios
    dx <- r * cos(theta)
    dy <- r * sin(theta)
    dx <- c(0, dx)
    dy <- c(0, dy)
    return(data.frame(dx=dx, dy=dy, acc=acc))
  }, df=df, ratios=ratios, r=r) |>
    do.call(what=rbind) |>
    mutate(x=cur.MAF + dx, y=cur.RN + dy) |>
    select(acc, x, y) |>
    mutate(MAF=cur.MAF, RN=cur.RN)
  return(df.poly)
}

MAF.RN.pairs <- distinct(df.sim, MAF, RN) |>
  arrange(MAF, RN)

df.poly <- lapply(seq_len(nrow(MAF.RN.pairs)), function(i, MAF.RN.pairs, df.sim){
  cur.MAF <- MAF.RN.pairs$MAF[i]
  cur.RN <- MAF.RN.pairs$RN[i]
  df.poly <- create.pie.polygon(cur.MAF, cur.RN, df.sim)
  return(df.poly)
}, MAF.RN.pairs=MAF.RN.pairs, df.sim=df.sim) |>
  do.call(what=rbind)



pA <- df.poly |>
  mutate(grp=paste0(MAF, '-', RN, '-', acc)) |>
  ggplot() + mytheme +
  geom_polygon(aes(x=x, y=y, fill=acc, group=grp), color='black', lwd=0.01) +
  geom_line(data=df.fit |> filter(grp=='>0.99'), aes(x=MAF, y=RN),
    lwd=1, color='red') +
  geom_ribbon(data=df.fit |> filter(grp=='>0.99'), aes(x=MAF, y=RN, ymin=RN.low, ymax=RN.high),
    alpha=0.2, color=NA, fill='red') +
  scale_fill_manual(values=col11) +
  labs(x='MAF', y='R / N', fill='Accuracy') +
  # coord_fixed(ratio=1) +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = 'black', linewidth = rel(1))
  )



pB <- readRDS('~/Desktop/solve-X/Figures/fit-curve-data/sample-identification-rate.rds') |>
  group_by(pop, N, rn.cut) |>
  summarise(sd.rate=sd(rate), rate=mean(rate), .groups='drop') |>
  filter(rn.cut >= 0.12, rn.cut <= 0.18) |>
  mutate(pop=paste0(pop, ' (N = ', N, ')')) |>
  mutate(ci.lwr=rate-1.96*sd.rate, ci.upr=rate+1.96*sd.rate) |>
  # filter(rn.cut < 0.2) |>
  ggplot() + mytheme +
  aes(x=rn.cut, y=rate, color=pop, group=pop) +
  geom_line(lwd=1) +
  # geom_errorbar(aes(ymin=rate-sd.rate, ymax=rate+sd.rate),
  #   width=0.001, alpha=0.3) +
  geom_ribbon(aes(ymin=ci.lwr, ymax=ci.upr, fill=pop),
    alpha=0.05, color=NA) +
  scale_y_continuous(breaks=seq(0, 1, by=0.2), labels=scales::percent_format(accuracy=1)) +
  labs(x='R / N', y='Mean sample identification rate\nbased on 1,000 randomly selected SNPs',
    color='Population', fill='Population') +
  scale_color_manual(values=cols) +
  scale_fill_manual(values=cols) +
  theme(
    axis.line = element_line(color = 'black', linewidth = rel(1)),
    axis.text = element_text(color='black'),
    panel.border = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    legend.position.inside = T,
    legend.background = element_blank(),
    legend.position = c(1, 0),
    legend.justification = c(1, 0)
  )



# p <- wrap_elements(pA + pB + plot_layout(guides='collect')) + pC + plot_layout(design='A\nB', heights=c(3, 2))
p <- wrap_elements(pA) + wrap_elements(pB) + plot_layout(design='AB')

p <- pA + pB + plot_layout(design='AB', widths=c(1, 2))
save.png(p, 'Fig4.png', w=12, h=6, units='in')


