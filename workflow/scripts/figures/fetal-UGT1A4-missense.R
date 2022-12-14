library("dplyr")
library("cowplot")
library("data.table")
library('showtext')
library(ggplot2)
library(broom)
options(warn=-1)

colorBlindBlack8= c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#font_add("arial", "arial.ttf", bold= 'arial_bold.ttf')

#font_add_google("Roboto", "Roboto")

showtext_opts(dpi = 300)
showtext_auto(enable = TRUE)

d= fread(snakemake@input[[1]], h= T)

x= group_by(d, round(fets_UGT_P24T)) %>% summarize(p= mean(jaundice, na.rm= T), lo95= binom.test(sum(jaundice), n())$conf.int[1], 
                                                   up95= binom.test(sum(jaundice), n())$conf.int[2])

names(x)[1]= 'UGT_missense'

x$UGT_missense= factor(x$UGT_missense, levels= c('0', '1', '2'), labels= c('CC\n(n = 20587)', 'CA\n(n = 2779)', 'AA\n(n = 95)'))

p1= ggplot(data= x, aes(UGT_missense, p*100, alpha= UGT_missense)) +
  geom_col(fill = colorBlindBlack8[2], colour= colorBlindBlack8[1]) +
  geom_errorbar( aes(x= UGT_missense, ymin= lo95*100, ymax= up95*100), colour= colorBlindBlack8[7], alpha=1, width= 0.08, 
                 linewidth= 0.7) +
  theme_cowplot(font_size= 10) + 
  scale_alpha_manual(values= c(0.5, 0.5, 1) , guide= 'none') +
    xlab('rs6755571 genotype') + 
  ylab('Jaundice risk, %') +
  scale_y_continuous(limits= c(0, 10), expand= expansion(add= c(0, 0.0002))) +
  theme(axis.text.x= element_text(size= 8, margin(t= 0, r= 0, b= 2, l= 0, unit= 'mm')),
        axis.title.x= element_text(margin(t= 2, r= 0, b= 0, l= 0, unit= 'mm')),
        axis.ticks.x= element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.2, lineend = "square"),
        axis.ticks.y = element_line(color = "black", size = 0.2),
        plot.background = element_rect(colour = colorBlindBlack8[1], fill= NA, linewidth= 0.2))
  


save_plot(snakemake@output[[1]], plot= p1, base_height= 60, base_width= 70, units= 'mm', dpi= 300)

eaf=  group_by(d, jaundice) %>% summarize(EAF= sum(fets_UGT_P24T) / (n() * 2), 
                                     lo95= binom.test(as.integer(sum(fets_UGT_P24T)), (n() * 2))$conf.int[1],
                                     up95= binom.test(as.integer(sum(fets_UGT_P24T)), (n() * 2))$conf.int[2])

eaf$jaundice= factor(eaf$jaundice, levels= c('0', '1'), labels = c('Controls\n(n = 21828)', 'Cases\n(n = 1633)'))

p2= ggplot(data= eaf, aes(factor(jaundice), EAF * 100)) +
  geom_pointrange(aes(ymin= lo95*100, ymax= up95*100), colour = colorBlindBlack8[2], size= 0.7, fatten= 0.9) +
    theme_cowplot(font_size= 10) + 
    xlab('Jaundice') + 
  ylab('Effect allele frequency') +
  scale_y_continuous(limits= c(0, 8), expand= expansion(add= c(0, 0.0002)), position = "right") +
  theme(axis.text.x= element_text(size= 8, margin(t= 0, r= 0, b= 2, l= 0, unit= 'mm')),
        axis.title.x= element_text(margin(t= 2, r= 0, b= 0, l= 0, unit= 'mm')),
        axis.ticks.x= element_blank(),
        axis.line = element_line(color = "black", size = 0.2, lineend = "square"),
        axis.ticks.y = element_line(color = "black", size = 0.2))


save_plot(snakemake@output[[2]], plot= p2, base_height= 60, base_width= 50, units= 'mm', dpi= 300)


d= fread(snakemake@input[[2]], header = T)

m1= glm(jaundice~ chr2_234627536_C_A_MT + chr2_234627536_C_A_MnT + chr2_234627536_C_A_PT + chr2_234627536_C_A_PnT + cohort + KJONN, 
        d, family= 'binomial')

ci= data.frame(confint(m1))
ci$term= row.names(ci)
names(ci)= c('lo95', 'up95', 'term')

m1= tidy(m1) %>% filter(grepl('chr2', term)) %>% inner_join(., ci, by= 'term')

m1$term= factor(m1$term, levels= rev(c("chr2_234627536_C_A_MT", "chr2_234627536_C_A_MnT", "chr2_234627536_C_A_PT",  
                                   "chr2_234627536_C_A_PnT")), labels= rev(c('Maternal\ntransmitted', 'Maternal\nnon-transmitted', 
                                                                        'Paternal\ntransmitted', 'Paternal\nnon-transmitted')))

p3= ggplot(m1, aes(x = term, y = estimate)) + 
  geom_hline(aes(yintercept = 0), size = .2, linetype = "dashed") +
  geom_hline(yintercept = log(setdiff(seq(0.0, 1.3, 0.2), 1)), size = .1, linetype = "dashed", colour= 'grey') + 
  geom_errorbar(aes(ymin = lo95, ymax = up95), size = .5, width = 0, color = colorBlindBlack8[1]) +
  theme_cowplot(font_size= 10) + 
  geom_point(size = 1, color = colorBlindBlack8[1]) +
  coord_trans(y = scales:::exp_trans()) +
  scale_y_continuous(breaks = log(seq(0.0000001, 1.4000002, 0.2)), labels = seq(0.0, 1.4, 0.2), limits = log(c(0.070, 1.4000002)),
                     expand= expansion(add=0)) +
  ylab('Odds Ratio') +
  theme(axis.title.x= element_blank(),
        axis.ticks.x= element_blank(),
        axis.text.x= element_text(size= 8),
        axis.line = element_line(color = "black", size = 0.2, lineend = "square"),
        axis.ticks.y = element_line(color = "black", size = 0.2))

save_plot(snakemake@output[[3]], plot= p3, base_height= 60, base_width= 90, units= 'mm', dpi= 300)
