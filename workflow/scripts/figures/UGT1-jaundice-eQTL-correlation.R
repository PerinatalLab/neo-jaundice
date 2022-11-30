library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)

colorBlindBlack8= c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


d= fread('results/eQTL_catalogue/delivery/SNP-jaundice-fets.txt')
pph= fread('results/eQTL_catalogue/delivery/pph-jaundice-fets.txt')
eqtls= filter(pph, PP.H4.abf > 0.9) %>% pull(eqtl_data)
eqtls= c(eqtls, 'GTEx_ge_liver')

d= filter(d, gene== 'ENSG00000241635', eqtl_data %in% eqtls)

link= fread('results/eQTL_catalogue/jaundice/temp/hg38/fets-GWAS.txt') %>% select(ID, POS)

link= separate(link, ID, into= c('CHR', 'POS2', 'REF', 'EFF'), sep= ':')

d= inner_join(d, link, by= c('snp'= 'POS'))

d$POS2= as.numeric(d$POS2)
#d= filter(d, POS2 > 234.45e6, POS2 < 234.75e6)

ld = read.table("resources/LD_variants_rs6755571.txt", h=T)
ld = separate(ld, Coord, c("CHR", "POS2"), sep=":")
ld$POS2 = as.numeric(ld$POS2)

pall = left_join(d, ld[,c("POS2", "Distance","R2")])

ggplot(pall, aes(z.df1, z.df2, size= SNP.PP.H4))+
  facet_grid(vars(eqtl_data)) +
  geom_point()

ggplot(pall, aes(z.df1, z.df2, size= SNP.PP.H4))+
  facet_grid(vars(eqtl_data)) +
  geom_point() +
geom_label(filter(pall, POS2== 234627536),  aes(z.df1, z.df2, label= 'rs6755571'))
