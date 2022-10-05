library("dplyr")
library("cowplot")
library("data.table")
library('showtext')
library(ggplot2)
library(broom)
options(warn=-1)

colorBlindBlack8= c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

showtext_opts(dpi = 300)
showtext_auto(enable = TRUE)

d= fread(snakemake@input[[1]], header = T)

m1= glm(jaundice~ chr2_234638006_A_G_h1 + chr2_234638006_A_G_h2 + chr2_234638006_A_G_h3 + chr2_234638006_A_G_h4 + cohort + KJONN, 
        d, family= 'binomial')

ci= data.frame(confint(m1))
ci$term= row.names(ci)
names(ci)= c('lo95', 'up95', 'term')

m1= tidy(m1) %>% filter(grepl('chr2', term)) %>% inner_join(., ci, by= 'term')

m1$term= factor(m1$term, levels= rev(c("chr2_234638006_A_G_h1", "chr2_234638006_A_G_h2", "chr2_234638006_A_G_h3",  
                                       "chr2_234638006_A_G_h4")), labels= rev(c('Maternal\ntransmitted', 'Maternal\nnon-transmitted', 
                                                                                 'Paternal\ntransmitted', 'Paternal\nnon-transmitted')))

p1= ggplot(m1, aes(x = term, y = estimate)) + 
  geom_hline(aes(yintercept = 0), size = .2, linetype = "dashed") +
  geom_hline(yintercept = log(setdiff(seq(0.0, 1.3, 0.2), 1)), size = .1, linetype = "dashed", colour= 'grey') + 
  geom_errorbar(aes(ymin = lo95, ymax = up95), size = .5, width = 0, color = colorBlindBlack8[2]) +
  theme_cowplot(font_size= 10) + 
  geom_point(size = 1, color = colorBlindBlack8[2]) +
  coord_trans(y = scales:::exp_trans()) +
  scale_y_continuous(breaks = log(seq(0.0000001, 1.4000002, 0.2)), labels = seq(0.0, 1.4, 0.2), limits = log(c(0.040, 1.4000002)),
                     expand= expansion(add=0)) +
  ylab('Odds Ratio') +
  theme(axis.title.x= element_blank(),
        axis.ticks.x= element_blank(),
        axis.text.x= element_text(size= 8),
        axis.line = element_line(color = "black", size = 0.2, lineend = "square"),
        axis.ticks.y = element_line(color = "black", size = 0.2))

save_plot(snakemake@output[[1]], plot= p1, base_height= 60, base_width= 90, units= 'mm', dpi= 300)

m1= glm(jaundice~ chr2_234522619_A_C_h1 + chr2_234522619_A_C_h2 + chr2_234522619_A_C_h3 + chr2_234522619_A_C_h4 + cohort + KJONN, 
        d, family= 'binomial')

ci= data.frame(confint(m1))
ci$term= row.names(ci)
names(ci)= c('lo95', 'up95', 'term')

m1= tidy(m1) %>% filter(grepl('chr2', term)) %>% inner_join(., ci, by= 'term')

m1$term= factor(m1$term, levels= rev(c("chr2_234522619_A_C_h1", "chr2_234522619_A_C_h2", "chr2_234522619_A_C_h3",  
                                       "chr2_234522619_A_C_h4")), labels= rev(c('Maternal\ntransmitted', 'Maternal\nnon-transmitted', 
                                                                                'Paternal\ntransmitted', 'Paternal\nnon-transmitted')))

p2= ggplot(m1, aes(x = term, y = estimate)) + 
  geom_hline(aes(yintercept = 0), size = .2, linetype = "dashed") +
  geom_hline(yintercept = log(setdiff(seq(0.0, 1.3, 0.2), 1)), size = .1, linetype = "dashed", colour= 'grey') + 
  geom_errorbar(aes(ymin = lo95, ymax = up95), size = .5, width = 0, color = colorBlindBlack8[2]) +
  theme_cowplot(font_size= 10) + 
  geom_point(size = 1, color = colorBlindBlack8[2]) +
  coord_trans(y = scales:::exp_trans()) +
  scale_y_continuous(breaks = log(seq(0.0000001, 1.4000002, 0.2)), labels = seq(0.0, 1.4, 0.2), limits = log(c(0.070, 1.4000002)),
                     expand= expansion(add=0)) +
  ylab('Odds Ratio') +
  theme(axis.title.x= element_blank(),
        axis.ticks.x= element_blank(),
        axis.text.x= element_text(size= 8),
        axis.line = element_line(color = "black", size = 0.2, lineend = "square"),
        axis.ticks.y = element_line(color = "black", size = 0.2))

save_plot(snakemake@output[[2]], plot= p2, base_height= 60, base_width= 90, units= 'mm', dpi= 300)