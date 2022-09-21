
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


fitted_models= group_by(d, GestAge_cat) %>% do(model= glm(as.numeric(OVERFLYTTET== 1 & jaundice== 1)~ I(h1_jaundice + h3_jaundice) +  h2_jaundice  + h4_jaundice + ABO_incompatibility  + chr2_234649665_C_T_h2 + I(chr2_234649665_C_T_h1 + chr2_234649665_C_T_h3) + chr2_234649665_C_T_h4  + as.numeric(is.na(MISD)) + as.numeric(PARITET_5==0) + KJONN + cohort, ., family= 'binomial'))

tmodel= fitted_models %>% tidy(model) %>% filter(GestAge_cat != 'NA')
tmodel$GestAge_cat= factor(tmodel$GestAge_cat, levels= c('PTD', 'Early term', 'Full term', 'Post-term'))

tmodel$lo_ci= tmodel$estimate - 1.96 * tmodel$std.error
tmodel$up_ci= tmodel$estimate + 1.96 * tmodel$std.error
tmodel$OR= exp(tmodel$estimate)
tmodel$OR_lo95= exp(tmodel$lo_ci)
tmodel$OR_up95= exp(tmodel$up_ci)

tmodel %>% filter(term== 'I(h1_jaundice + h3_jaundice)') %>%
ggplot(., aes(GestAge_cat, OR, colour= term)) + 
geom_pointrange(aes(ymin= OR_lo95, ymax= OR_up95)) +
geom_hline(yintercept= 1, colour= 'grey', linetype= 'dashed') +
scale_y_log10() +
theme_cowplot() +
geom_text_repel(aes(GestAge_cat, OR, label= as.character(signif(p.value, 2))), hjust= 0)


