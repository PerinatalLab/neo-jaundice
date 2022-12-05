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


d= arrange(d, desc(LOG10P))
d= d[!duplicated(d$ID), ]

df= mutate(d, exp1= -log10(1:length(LOG10P)/length(LOG10P)))

chisq= qchisq(1-10**-df$LOG10P, 1)

lambda_gc= median(chisq)/qchisq(0.5,1)

p1= ggplot(df, aes(exp1, LOG10P)) +
  geom_point(size= 0.4, color= colorBlindBlack8[2]) +
  geom_abline(intercept = 0, slope = 1, alpha = .5) +
labs(colour="") +
theme_cowplot(font_size= 10) +
xlab(expression(Expected~-log[10]~pvalue)) +
ylab(expression(Observed~-log[10]~pvalue)) +
geom_text(aes(6, 0), label= paste("lambda", "==", round(lambda_gc, 2)), size= 10/.pt, parse= T)

ggsave(snakemake@output[[1]], plot= p1, dpi = 300, width= 60, height= 60, units= 'mm')
