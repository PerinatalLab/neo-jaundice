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

d= fread(snakemake@input[[1]], header = T)

m1= glm(jaundice~ chrX_109792100_C_T_h1 + chrX_109792100_C_T_h2 + chrX_109792100_C_T_h3 + cohort + PC1 + PC2 + PC3 + PC4 + PC5 + PC6  + PC7 + PC8 + PC9 + PC10, filter(d, KJONN== 0), family= 'binomial')

ci= data.frame(confint(m1))
ci$term= row.names(ci)
names(ci)= c('lo95', 'up95', 'term')

m1= tidy(m1) %>% filter(grepl('chrX', term)) %>% inner_join(., ci, by= 'term')

m1$sex= 'Girls'

m2= glm(jaundice~ chrX_109792100_C_T_h1 + chrX_109792100_C_T_h2 + chrX_109792100_C_T_h4 + cohort + PC1 + PC2 + PC3 + PC4 + PC5 + PC6  + PC7 + PC8 + PC9 + PC10, filter(d, KJONN== 1), family= 'binomial')

ci= data.frame(confint(m2))
ci$term= row.names(ci)
names(ci)= c('lo95', 'up95', 'term')

m2= tidy(m2) %>% filter(grepl('chrX', term)) %>% inner_join(., ci, by= 'term')

m2$sex= 'Boys'

m1= rbind(m1, m2)

m1$term= factor(m1$term, levels= rev(c("chrX_109792100_C_T_h1", "chrX_109792100_C_T_h2", "chrX_109792100_C_T_h3",
                                   "chrX_109792100_C_T_h4")), labels= rev(c('Maternal\ntransmitted', 'Maternal\nnon-transmitted',
                                                                        'Paternal\ntransmitted', 'Paternal\nnon-transmitted')))

p1= ggplot(m1, aes(x = term, y = estimate, colour= sex)) +
  geom_hline(aes(yintercept = 0), size = .2, linetype = "dashed") +
  geom_hline(yintercept = log(setdiff(seq(0.4, 1.3, 0.1), 1)), size = .1, linetype = "dashed", colour= 'grey') +
  geom_errorbar(aes(ymin = lo95, ymax = up95), size = .5, width = 0, position = position_dodge(width=0.3)) +
  theme_cowplot(font_size= 10) +
  geom_point(size = 1, position = position_dodge(width=0.3)) +
scale_color_manual(values= colorBlindBlack8[c(2, 6)], name= "Sex") +
  coord_trans(y = scales:::exp_trans()) +
  scale_y_continuous(breaks = log(seq(0.4, 1.3000002, 0.2)), labels = seq(0.4, 1.3, 0.2), limits = log(c(0.4, 1.3000002)),
                     expand= expansion(add=0)) +
  ylab('Odds Ratio') +
  theme(axis.title.x= element_blank(),
        axis.ticks.x= element_blank(),
        axis.text.x= element_text(size= 8),
        axis.line = element_line(color = "black", size = 0.2, lineend = "square"),
        axis.ticks.y = element_line(color = "black", size = 0.2),
        legend.position = 'bottom',
        legend.box = "horizontal")

save_plot(snakemake@output[[1]], plot= p1, base_height= 60, base_width= 90, units= 'mm', dpi= 300)

x= fread(snakemake@input[[2]])

x= inner_join(x, d, by= 'PREG_ID_1724')

names(x)= gsub(':', '_', names(x))

x$dads_X_109792100_C_T= x$dads_X_109792100_C_T * 2
x$fets_X_109792100_C_T= ifelse(x$KJONN== 1, x$fets_X_109792100_C_T * 2, x$fets_X_109792100_C_T)

m1= glm(jaundice~  fets_X_109792100_C_T + moms_X_109792100_C_T + dads_X_109792100_C_T + cohort + KJONN +
          PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
        x, family= 'binomial')

ci= data.frame(confint(m1))
ci$term= row.names(ci)
names(ci)= c('lo95', 'up95', 'term')

m1= tidy(m1) %>% filter(grepl('_X', term)) %>% inner_join(., ci, by= 'term')

m1$term= factor(m1$term, levels= (c("fets_X_109792100_C_T", "moms_X_109792100_C_T", "dads_X_109792100_C_T")), 
                labels= (c('Neonate', 'Maternal', 'Paternal')))

p1= ggplot(m1, aes(x = term, y = estimate)) +
  geom_hline(aes(yintercept = 0), size = .2, linetype = "dashed") +
  geom_hline(yintercept = log(setdiff(seq(0.6, 1.2, 0.1), 1)), size = .1, linetype = "dashed", colour= 'grey') +
  geom_errorbar(aes(ymin = lo95, ymax = up95), size = .5, width = 0, color = colorBlindBlack8[1]) +
  theme_cowplot(font_size= 10) +
  geom_point(size = 1, color = colorBlindBlack8[1]) +
  coord_trans(y = scales:::exp_trans()) +
  scale_y_continuous(breaks = log(seq(0.6, 1.2000002, 0.2)), labels = seq(0.6, 1.2, 0.2), limits = log(c(0.6, 1.2000002)),
                     expand= expansion(add=0)) +
  ylab('Odds Ratio') +
  theme(axis.title.x= element_blank(),
        axis.ticks.x= element_blank(),
        axis.text.x= element_text(size= 8),
        axis.line = element_line(color = "black", size = 0.2, lineend = "square"),
        axis.ticks.y = element_line(color = "black", size = 0.2))

save_plot(snakemake@output[[2]], plot= p1, base_height= 60, base_width= 90, units= 'mm', dpi= 300)
