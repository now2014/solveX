#!/usr/bin/env Rscript

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