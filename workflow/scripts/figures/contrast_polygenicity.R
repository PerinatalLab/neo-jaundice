library(data.table)
library(dplyr)
library(broom)
library(ggplot2)
library(showtext)
library(cowplot)
library(ggrepel)

colorBlindBlack8= c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

showtext_opts(dpi = 300)
showtext_auto(enable = TRUE)

d= fread(snakemake@input[[1]], h= T)
d$pheno= 'child'

x= fread(snakemake@input[[2]], h= T)

x$pheno= 'height'

#d= fread('/mnt/work/pol/neo-jaundice/results/HESS/contrast-polygenicity/jaundice-fets.txt', h= T)

d= rbind(d, x)

d$pheno= factor(d$pheno, levels= c('child', 'height'))

pheno_colours= unlist(ifelse(grepl('fets', snakemake@input[[1]]), list(colorBlindBlack8[c(2, 4)]), list(colorBlindBlack8[c(6,4)])))
print(pheno_colours)

p1= ggplot(d, aes(nsnp_sorted * 100, local_h2 * 100)) +
geom_ribbon(aes(ymin= (local_h2 - se) * 100, ymax= (local_h2 + se) * 100, colour= pheno, fill= pheno), alpha= 0.2, size= 0.1) +
  geom_line(aes(colour= pheno), size= 0.4) +
  scale_colour_manual(values= pheno_colours, guide= 'none') +
  scale_fill_manual(values= pheno_colours, guide= 'none') +
theme_cowplot(font_size= 10) + 
  xlab('Proportion of genome, %') + 
  ylab('Proportion of SNP-heritability, %') +
  scale_y_continuous(limits= c(0, 100.01), expand= expansion(add= c(0, 0))) +
  scale_x_continuous(limits= c(0, 100.01), expand= expansion(add= c(0, 0))) +
  theme(axis.text.x= element_text(size= 8),
        axis.ticks.x= element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2, lineend = "square"),
        axis.ticks.y = element_line(color = "black", size = 0.2),
        plot.margin = margin(t=5,r= 5, b=5, l=5, "mm")) +
  geom_abline(slope= 1, intercept= 0, size= 0.2)

save_plot(snakemake@output[[1]], plot= p1, base_height= 65, base_width= 65, units= 'mm', dpi= 300)


