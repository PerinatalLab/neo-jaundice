library(dplyr)
library(data.table)


h1= fread(snakemake@input[[1]])
h2= fread(snakemake@input[[2]])
h3= fread(snakemake@input[[3]])
h4= fread(snakemake@input[[4]])

pheno= fread(snakemake@input[[5]], select= c('PREG_ID', 'jaundice', 'FID_y', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10', 'cohort', 'KJONN', 'Child', 'Father', 'Mother'))

pheno= inner_join(pheno, h1, by= 'PREG_ID') %>% inner_join(., h2, by= c('PREG_ID')) %>% inner_join(., h3, by= c('PREG_ID')) %>% inner_join(., h4, by= c('PREG_ID'))

m1= glm(jaundice~ h1_jaundice + h2_jaundice + h3_jaundice + h4_jaundice + cohort + KJONN + PC1 + PC2 + PC3 + PC4 + PC5 + PC6, pheno, family= binomial)

n= length(resid(m1))
coefs= data.frame(summary(m1)$coefficients[2:5,])
names(coefs)= c('beta', 'se', 'z', 'pvalue')

coefs$haplotype= rownames(coefs)
coefs$n= n

fwrite(coefs, snakemake@output[[1]], sep= '\t')

