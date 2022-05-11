library(ggplot2)
library(dplyr)
library(data.table)
library(cowplot)
library(ggrepel)

d= fread(snakemake@input[[1]], header= T, sep= '\t', select= c('CHR', 'POS', 'ID', 'LOG10P'))

topvars= fread(snakemake@input[[2]], header= T, sep= '\t', select= c('ID', 'nearestGene'))

d= filter(d, !duplicated(ID), LOG10P> -log10(0.05))

colorBlindBlack8= c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

d= left_join(d, topvars, by= 'ID')

don <- d %>%
    arrange(CHR, POS) %>%
    group_by(CHR)      %>%
    summarise(chr_len= max(POS)) %>%
    mutate(tot= cumsum(as.numeric(chr_len))-chr_len) %>% # Calculate cumulative position of each chromosome
    select(-chr_len) %>%
    left_join(d, ., by= 'CHR') %>%
    arrange(CHR, POS) %>% # Add a cumulative position of each SNP
    mutate(BPcum=POS+tot) %>%
         ungroup()

print(table(don$nearestGene))

axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  names(axisdf)= c('CHR', 'center')
HC= -log10(5*10**-8)

print(table(don$nearestGene))

print(length(axisdf$center))
print(length(c(1:19, '', 21,'', 'X')))

p1= ggplot(data= don, aes(x= BPcum, y= LOG10P, colour= factor(CHR))) +
  geom_point(size= 0.07) +   # Show all points
  theme_cowplot(font_size= 9) +
  scale_colour_manual(values= c(rep(c(colorBlindBlack8[8], colorBlindBlack8[3]), 21), colorBlindBlack8[8]), guide= F) +
  scale_x_continuous(label = c(1:19, '', 21,'', 'X'), breaks= axisdf$center, expand= c(0.03, 0.03)) + # label = ifelse(axisdf$CHR== 23, 'X', axisdf$CHR)
#  scale_y_continuous(expand= c(0, 0), limits= c(min(don$logpval) - 2, max(don$logpval) + 2), breaks= seq(-30, 45, 10), labels= c(abs(seq(-30, 45, 10)))) + # , sec.axis = sec_axis(~ ., name = derive())) +
  ylab('-log10(pvalue)') +
  xlab('Chromosome') +
  geom_hline(yintercept= 0, size= 0.25, colour= 'black') +
  geom_hline(yintercept= HC, size= 0.2, linetype= 2, colour= '#878787') +
  coord_cartesian(clip = "off") +
  geom_text_repel(aes(label= nearestGene),
                  size= 10/ .pt)
#                  force_pull= 0, # do not pull toward data points
#                  force= 0.1,
#                  nudge_y      =  ifelse(filter(don, nearestGene!= '') %>% pull(LOG10P)>0, 1, -1), #43 - ((-log10(filter(don, GENE!= '')$pvalue))),
#                  direction    = "both",
#                  hjust        = 0,
#                  vjust=  0.5,
#                  box.padding= 0.1,
#                  angle= 0,
#                  segment.size = 0.1,
#                  segment.square= TRUE,
#                  segment.inflect= FALSE,
#                  segment.colour= colorBlindBlack8[8],
#                  segment.linetype = 4) +
 #                 ylim = c(0, 50),
 #                 xlim = c(-Inf, Inf)) +
#  theme(legend.position= 'none',
#        plot.margin = unit(c(t= 0, r=0, b= 0, l=0), 'cm'),
#        text= element_text(family="arial", size= 9),
#        axis.line= element_line(size= 0.1))

save_plot(snakemake@output[[1]], plot= p1, base_height= 90, base_width= 185, units= 'mm', dpi= 300)

