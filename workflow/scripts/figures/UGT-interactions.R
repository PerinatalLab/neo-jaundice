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
d= select(d, -names(d)[which(duplicated(names(d)))])

d= filter(d, SVLEN_UL_DG<308, SVLEN_UL_DG> 154)

d$cat_GA= with(d, ifelse(SVLEN_UL_DG< 259, 'Preterm',
                         ifelse(SVLEN_UL_DG< 273, 'Early term',
                                ifelse(SVLEN_UL_DG< 280, 'Term', 
                                       'Post-term'))))

fitted_models= d %>% nest_by(cat_GA) %>% mutate(model= list(glm(jaundice ~ fets_UGT_P24T, data = data, family= 'binomial')))

betas= fitted_models %>% summarize(tidy(model))
cis= (fitted_models %>% summarize(confint(model)))
cis= cbind(data.frame(cis[[1]]), data.frame(term = row.names(cis[[2]]), cis[[2]]))

names(cis)= c('cat_GA', 'term', 'lo95', 'up95')

models= inner_join(betas, cis, by= c('cat_GA', 'term')) %>% filter(grepl('UGT', term))

names(models) = c('cat_GA', 'term', 'estimate', 'se', 'stat', 'pvalue', 'lo95', 'up95')

models$cat_GA= factor(models$cat_GA, levels= c('Preterm', 'Early term', 'Term', 'Post-term'))
  
p1= ggplot(models, aes(x = cat_GA, y = estimate)) +
    geom_hline(yintercept = log(setdiff(seq(0.0, 1.3, 0.2), 1)), size = .1, linetype = "dashed", colour= 'grey') +
  geom_hline(aes(yintercept = 0), size = .2, linetype = "dashed") +
  geom_errorbar(aes(ymin = lo95, ymax = up95), size = .2, width = 0, color = colorBlindBlack8[2]) +
  theme_cowplot(font_size= 10) +
  geom_point(size = 1.2, color = colorBlindBlack8[2]) +
  coord_trans(y = scales:::exp_trans()) +
  scale_y_continuous(breaks = log(seq(0.0000001, 1.2000002, 0.2)), labels = seq(0.0, 1.2, 0.2), limits = log(c(0.001, 1.2000002)),
                     expand= expansion(add=0)) +
  ylab('Odds Ratio') +
  theme(axis.title.x= element_blank(),
        axis.ticks.x= element_blank(),
        axis.text.x= element_text(size= 8),
        axis.line = element_line(color = "black", size = 0.2, lineend = "square"),
        axis.ticks.y = element_line(color = "black", size = 0.2))

save_plot(snakemake@output[[1]], plot= p1, base_height= 60, base_width= 90, units= 'mm', dpi= 300)

d$fets_UGT_P24T_cat= factor(round(d$fets_UGT_P24T), levels= rev(c(0, 1, 2)), labels= rev(c('CC', 'CA', 'AA')))

dat_text= data.frame(label= c('CC', 'CA', 'AA'), fets_UGT_P24T_cat= c('CC', 'CA', 'AA'), SVLEN_UL_DG= c(156, 156, 156))

p2= ggplot(data= d, aes(SVLEN_UL_DG, group= factor(jaundice),  fill= factor(jaundice))) + 
    scale_colour_manual(values= colorBlindBlack8[c(2, 6)], guide= 'none') +
  scale_fill_manual(values= colorBlindBlack8[c(2, 6)], guide= 'none') +
  facet_grid(vars(fets_UGT_P24T_cat)) +
  geom_density(size= 0.1, alpha= 0.6, color= 'black') + 
  geom_hline(aes(yintercept = 0), size = .2,  colour= 'black') +
  theme_cowplot(font_size= 10) +
  ylab('Density') +
  xlab('Gestational duration, days') +
    scale_x_continuous(limits= c(154, 308), expand= expansion(add= 1)) +
  scale_y_continuous(limits= c(0, 0.055), expand= expansion(add= 0)) +
    theme(axis.title.x= element_blank(),
        axis.ticks.x= element_blank(),
        axis.text.x= element_text(size= 8),
        axis.text.y= element_text(size= 8),
        axis.line = element_line(color = "black", size = 0.2, lineend = "square"),
        axis.ticks.y = element_line(color = "black", size = 0.2),
        strip.background = element_blank(),
        strip.text.y = element_blank()) +
  geom_text_repel(data= dat_text, aes(x= -Inf, y= Inf, label= label),  inherit.aes = FALSE, hjust= 1)

save_plot(snakemake@output[[2]], plot= p2, base_height= 80, base_width= 110, units= 'mm', dpi= 300)

fitted_models= d %>% nest_by(ABO_incompatibility) %>% mutate(model= list(glm(jaundice ~ fets_UGT_P24T, data = data, family= 'binomial')))

betas= fitted_models %>% summarize(tidy(model)) 
cis= (fitted_models %>% summarize(confint(model)))
cis= cbind(data.frame(cis[[1]]), data.frame(term = row.names(cis[[2]]), cis[[2]]))

names(cis)= c('ABO_incompatibility', 'term', 'lo95', 'up95')

models= inner_join(betas, cis, by= c('ABO_incompatibility', 'term')) %>% filter(grepl('UGT', term))

names(models) = c('ABO', 'term', 'estimate', 'se', 'stat', 'pvalue', 'lo95', 'up95')

models$ABO= factor(models$ABO, levels = c(0, 1), labels= c('ABO\ncompatible', 'ABO\nincompatible'))

p3= ggplot(models, aes(x = ABO, y = estimate)) +
  geom_hline(yintercept = log(setdiff(seq(0.0, 1.3, 0.2), 1)), size = .1, linetype = "dashed", colour= 'grey') +
  geom_hline(aes(yintercept = 0), size = .2, linetype = "dashed") +
  geom_errorbar(aes(ymin = lo95, ymax = up95), size = .2, width = 0, color = colorBlindBlack8[2]) +
  theme_cowplot(font_size= 10) +
  geom_point(size = 1.2, color = colorBlindBlack8[2]) +
  coord_trans(y = scales:::exp_trans()) +
  scale_y_continuous(breaks = log(seq(0.0000001, 1.2000002, 0.2)), labels = seq(0.0, 1.2, 0.2), limits = log(c(0.001, 1.2000002)),
                     expand= expansion(add=0)) +
  ylab('Odds Ratio') +
  theme(axis.title.x= element_blank(),
        axis.ticks.x= element_blank(),
        axis.text.x= element_text(size= 8),
        axis.line = element_line(color = "black", size = 0.2, lineend = "square"),
        axis.ticks.y = element_line(color = "black", size = 0.2))

save_plot(snakemake@output[[3]], plot= p3, base_height= 60, base_width= 90, units= 'mm', dpi= 300)
