#!/usr/bin/env Rscript

setwd('/share/home/lanao/projects/06-solveX/GTEx/v7/04-solve')

mycols <- c('#FF6F00FF', '#C71000FF', '#008EA0FF', '#8A4198FF', '#5A9599FF', '#FF6348FF',
            '#84D7E1FF', '#FF95A8FF', '#3D3B25FF', '#ADE2D0FF', '#1A5354FF', '#3F4041FF')
acc.cols <- c('#8A4198', '#008EA0', '#FF6F00', '#0846aa', '#FF95A8')
dpi <- 300
font.size <- 48

N <- 369

df <- readRDS('data/acc.rds') |>
  rename('MAF.X'='MAF') |>
  left_join(fread('../data/snp-cnt.tsv'), by='SNP') |>
  left_join(fread('../data/snp-maf.tsv'), by='SNP') |>
  mutate(RN=n.genes/N)

saveRDS(df, file='GTEx-acc-MAF-RN.rds')

# p <- df |>
#   # filter(MAF <= 0.05) |>
#   filter(MAF <= 0.1) |>
#   filter(RN > 0.2) |>
#   # filter(acc > 0.9) |>
#   count(RN, MAF, acc, name='cnt') |>
#   mutate(
#     acc.bin=cut(acc, breaks=c(0.6, 0.8, 0.9, 0.95, 0.99, 1))) |>
#   arrange(desc(acc.bin)) |>
  
#   ggplot() + vstheme(size=font.size) +
#   geom_point(aes(x=MAF, y=RN, color=acc.bin, size=cnt), alpha=0.2) +
#   geom_hline(yintercept=0.25, lty=2, lwd=0.5, color='black') +
#   geom_vline(xintercept=0.05, lty=2, lwd=0.5, color='black') +
#   guides(color=guide_legend(override.aes=list(alpha=0.4, size=2))) +
#   # facet_wrap(~acc.bin, ncol=2) +
#   labs(x='MAF', y='Number of traits / N', color='Accuracy bin', size='Count') +
#   scale_color_manual(values=rev(acc.cols))

# save.png(
#   p, 'png/GTEx-v7-acc.png', 6, 3, dpi=dpi
# )



## 2D density plot
p <- df |>
  # filter(MAF <= 0.05) |>
  filter(MAF <= 0.1) |>
  filter(RN > 0.2) |>
  # filter(acc > 0.9) |>
  mutate(
    acc.bin=cut(acc, breaks=c(0.6, 0.8, 0.9, 0.95, 0.99, 1))) |>
  arrange(desc(acc.bin)) |>

  ggplot() + vstheme(size=font.size) +
  geom_density_2d(aes(x=MAF, y=RN, color=acc.bin), alpha=0.5, lwd=0.5) +
  geom_point(aes(x=MAF, y=RN, color=acc.bin), alpha=0.01) +
  # aes(x=MAF, y=RN, group=acc.bin) +
  # stat_density_2d(
  #   geom = 'polygon',
  #   aes(alpha = after_stat(level), fill = acc.bin),
  #   bins = 6
  # ) +
  geom_hline(yintercept=0.25, lty=2, lwd=0.5, color='black') +
  geom_vline(xintercept=0.05, lty=2, lwd=0.5, color='black') +
  # guides(color=guide_legend(override.aes=list(linewidth=1))) +
  # facet_wrap(~acc.bin, ncol=2) +
  labs(x='MAF', y='Number of traits / N', color='Accuracy bin', size='Count') +
  scale_color_manual(values=rev(acc.cols)) +
  scale_fill_manual(values=rev(acc.cols))

save.png(
  p, 'png/GTEx-v7-acc.png', 6, 3, dpi=dpi
)
