
library(broom)
library(dplyr)
library(data.table)
library(cowplot)
library(ggplot2)
library(ggrepel)

d= fread('results/merge_data/delivery/jaundice.txt')


d$GestAge_cat= with(d, ifelse(is.na(SVLEN_UL_DG), NA, 
ifelse(SVLEN_UL_DG< 259, 'PTD', 
ifelse(SVLEN_UL_DG>= 259 & SVLEN_UL_DG < 273, 'Early term', 
ifelse(SVLEN_UL_DG>= 273 & SVLEN_UL_DG < 287, 'Full term',
ifelse(SVLEN_UL_DG>= 287 & SVLEN_UL_DG< 308, 'Post-term', NA))))))

d$GestAge_cat= with(d, ifelse(FSTART== 1, GestAge_cat, NA))

m_non_adjusted= glm(jaundice~ h1_jaundice +  h2_jaundice + h3_jaundice + h4_jaundice + as.numeric(is.na(MISD)) + as.numeric(PARITET_5==0) + KJONN + cohort, d, family= 'binomial')
m_adjusted= glm(jaundice~ h1_jaundice +  h2_jaundice + h3_jaundice + h4_jaundice + chr2_234649665_C_T_h1  + chr2_234649665_C_T_h2 + chr2_234649665_C_T_h3 + chr2_234649665_C_T_h4 + as.numeric(is.na(MISD)) + as.numeric(PARITET_5==0) + KJONN + cohort, d, family= 'binomial')

m_non_adjusted= tidy(m_non_adjusted)
m_adjusted= tidy(m_adjusted)

m_non_adjusted$model= 'Non-adjusted'
m_adjusted$model= 'Adjusted'

m1= rbind(m_non_adjusted, m_adjusted)
m1$model= levels(m1$model, c('Non-adjusted', 'adjusted'))

m1$lo_ci= m1$estimate - 1.96 * m1$std.error
m1$up_ci= m1$estimate + 1.96 * m1$std.error
m1$OR= exp(m1$estimate)
m1$OR_lo95= exp(m1$lo_ci)
m1$OR_up95= exp(m1$up_ci)

m1= filter(m1, grepl('jaundice', term))

m1$term= with(m1, ifelse(term== 'h1_jaundice', 'MT', ifelse(term== 'h2_jaundice', 'MnT', ifelse(term== 'h3_jaundice', 'PT', ifelse(term== 'h4_jaundice', 'PnT', term)))))


ggplot(filter(m1, model== 'Non-adjusted'), aes(term, OR, colour= model)) + 
geom_pointrange(aes(ymin= OR_lo95, ymax= OR_up95 ), position=position_dodge(width=0.3)) +
geom_hline(yintercept= 1, colour= 'grey', linetype= 'dashed') +
scale_y_log10() +
theme_cowplot()


ggplot(m1, aes(term, OR, colour= model)) + 
geom_pointrange(aes(ymin= OR_lo95, ymax= OR_up95 ), position=position_dodge(width=0.3)) +
geom_hline(yintercept= 1, colour= 'grey', linetype= 'dashed') +
scale_y_log10() +
theme_cowplot() +
geom_text_repel(aes(term, OR, label= as.character(signif(p.value, 2))), hjust= 0)


