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

#d= fread('/mnt/work/pol/neo-jaundice/results/UGT-missense/delivery/jaundice.txt')

# Plot density of Polygenic scores

p1= ggplot() +
  geom_density(data= d, aes(fets_jaundice), fill= colorBlindBlack8[6], size= 0.1, alpha= 0.6, color= 'black') +
  geom_density(data= d, aes(fets_jaundice_nochr2), fill= colorBlindBlack8[2], size= 0.1, alpha= 0.6, color= 'black') +
  theme_cowplot(font_size= 10) + 
  xlab('Polygenic score of adult bilirubin levels') + 
  ylab('Denstiy') +
  scale_y_continuous(limits= c(0, 5.5), expand= expansion(add= c(0, 0.0002))) +
  theme(axis.text.x= element_text(size= 8),
        axis.ticks.x= element_blank(),
        axis.line = element_line(color = "black", size = 0.2, lineend = "square"),
        axis.ticks.y = element_line(color = "black", size = 0.2))

save_plot(snakemake@output[[1]], plot= p1, base_height= 60, base_width= 110, units= 'mm', dpi= 300)

# Plot frequency of jaundice and 95%CI by polygenic score

d$PGS_cat= ntile(d$fets_jaundice, 10)
d$PGS_cat_nochr2= ntile(d$fets_jaundice_nochr2, 10)

x= group_by(d, PGS_cat) %>% summarize(p= mean(jaundice, na.rm= T), lo95= binom.test(sum(jaundice), n())$conf.int[1], 
                                                   up95= binom.test(sum(jaundice), n())$conf.int[2])

names(x)[1]= 'PGS'

x$PGS= factor(x$PGS)
x$chr= 'Full'
x$PGS= factor(x$PGS)

x2= group_by(d, PGS_cat_nochr2) %>% summarize(p= mean(jaundice, na.rm= T), lo95= binom.test(sum(jaundice), n())$conf.int[1], 
                                             up95= binom.test(sum(jaundice), n())$conf.int[2])

names(x2)[1]= 'PGS'
x2$chr= 'No chr2'

x2$PGS= factor(x2$PGS)

x= rbind(x, x2)

p2= ggplot(data= x, aes(PGS, p*100, colour= chr, group= chr)) +
      geom_point(position= position_dodge2(width= 0.5)) +
  geom_linerange( aes(x= PGS, ymin= lo95*100, ymax= up95*100, colour= chr), alpha= 1, size= 0.7, 
                  position= position_dodge2(width= 0.5)) +
  scale_color_manual(values= colorBlindBlack8[c(6, 2)], guide= 'none') +
    theme_cowplot(font_size= 10) + 
  xlab('Adult bilirubin polygenic score') + 
  ylab('Jaundice prevalence') +
  scale_y_continuous(limits= c(0, 12), expand= expansion(add= c(0, 0.0002))) +
  theme(axis.text.x= element_text(size= 8),
        axis.ticks.x= element_blank(),
        axis.line = element_line(color = "black", size = 0.2, lineend = "square"),
        axis.ticks.y = element_line(color = "black", size = 0.2))

save_plot(snakemake@output[[2]], plot= p2, base_height= 60, base_width= 90, units= 'mm', dpi= 300)

# Plot effects of parental transmitted non-transmitted PGS

d= fread(snakemake@input[[2]], h= T)

#d= fread('/mnt/work/pol/neo-jaundice/results/UGT-missense/delivery/jaundice-transmitted.txt')

m1= glm(jaundice~ h1_jaundice + h2_jaundice + h3_jaundice + h4_jaundice + cohort + KJONN + PC1 + PC2 + PC3 + PC4 + PC5 + PC6, 
        d, family= 'binomial')

ci= data.frame(confint(m1))
ci$term= row.names(ci)
names(ci)= c('lo95', 'up95', 'term')

m1= tidy(m1) %>% filter(grepl('jaundice', term)) %>% inner_join(., ci, by= 'term')

m1$term= factor(m1$term, levels= rev(c("h1_jaundice", "h2_jaundice", "h3_jaundice",  
                                       "h4_jaundice")), labels= rev(c('Maternal\ntransmitted', 'Maternal\nnon-transmitted', 
                                                                          'Paternal\ntransmitted', 'Paternal\nnon-transmitted')))

p3= ggplot(m1, aes(x = term, y = estimate)) + 
  geom_hline(aes(yintercept = 0), size = .2, linetype = "dashed") +
  geom_hline(yintercept = log(setdiff(seq(0.6, 3.4, 0.2), 1)), size = .1, linetype = "dashed", colour= 'grey') + 
  geom_errorbar(aes(ymin = lo95, ymax = up95), size = .5, width = 0, color = colorBlindBlack8[2]) +
  theme_cowplot(font_size= 10) + 
  geom_point(size = 0.7, color = colorBlindBlack8[2]) +
  coord_trans(y = scales:::exp_trans()) +
  scale_y_continuous(breaks = log(seq(0.6000001, 3.4000002, 0.2)), labels = seq(0.6, 3.4, 0.2), limits = log(c(0.6, 3.4000002)),
                     expand= expansion(add=0)) +
  ylab('Odds Ratio') +
  theme(axis.title.x= element_blank(),
        axis.ticks.x= element_blank(),
        axis.text.x= element_text(size= 8),
        axis.line = element_line(color = "black", size = 0.2, lineend = "square"),
        axis.ticks.y = element_line(color = "black", size = 0.2))

save_plot(snakemake@output[[3]], plot= p3, base_height= 60, base_width= 90, units= 'mm', dpi= 300)