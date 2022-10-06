library("dplyr")
library(ggplot2)
library(cowplot)
library("data.table")
library('showtext')
options(warn=-1)


d= fread(snakemake@input[[1]], h= T, select= c('LOG10P', 'ID'))


colorBlindBlack8= c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


showtext_opts(dpi = 300)
showtext_auto(enable = TRUE)


d= arrange(d, LOG10P)
d= d[!duplicated(d$ID), ]


df= arrange(d, LOG10P) %>% mutate(exp1= -log10(1:length(LOG10P)/length(LOG10P)))

p1= ggplot(df, aes(exp1, LOG10P)) +
  geom_point(size= 0.4, color= colorBlindBlack8[2]) +
  geom_abline(intercept = 0, slope = 1, alpha = .5) +
labs(colour="") +
theme_cowplot(font_size= 12) +
xlab('Expected (-log10(p-value))') +
ylab('Observed (-log10(p-value))')

ggsave(snakemake@output[[1]], plot= p1, width= 60, height= 60, units= 'mm', dpi= 300)
