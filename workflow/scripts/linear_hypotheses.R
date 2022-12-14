library(data.table)
library(dplyr)
library(tidyr)

format_haps= function(hap){
variants= paste(hap$chr, hap$pos, hap$ref, hap$eff, sep =':')
ids= names(hap)[5:ncol(hap)]
hap= as.data.frame(t(hap[, 5:ncol(hap)]))
names(hap)= variants
hap$PREG_ID_1724= as.numeric(ids)
return(hap)
}

h1= fread(snakemake@input[[1]])
h2= fread(snakemake@input[[2]])
h3= fread(snakemake@input[[3]])
h4= fread(snakemake@input[[4]])

h1= format_haps(h1)
h2= format_haps(h2)
h3= format_haps(h3)
h4= format_haps(h4)

pheno= fread(snakemake@input[[5]])

pheno$PREG_ID_1724= as.numeric(pheno$PREG_ID_1724)

print(nrow(pheno))
write( paste('snp', 'n', 'beta_h1', 'se_h1', 'pvalue_h1', 'beta_h2', 'se_h2', 'pvalue_h2', 'beta_h3', 'se_h3', 'pvalue_h3', 'beta_h4', 'se_h4', 'pvalue_h4', sep= '\t'), snakemake@output[[1]], append= T)

results_list= lapply(names(h1)[1:(length(names(h1))-1)], function(snp) {

print(snp)

h1_temp= h1[, c('PREG_ID_1724', snp)]
h2_temp= h2[, c('PREG_ID_1724', snp)]
h3_temp= h3[, c('PREG_ID_1724', snp)]
h4_temp= h4[, c('PREG_ID_1724', snp)]

names(h1_temp)= c('PREG_ID_1724', 'h1')
names(h2_temp)= c('PREG_ID_1724', 'h2')
names(h3_temp)= c('PREG_ID_1724', 'h3')
names(h4_temp)= c('PREG_ID_1724', 'h4')

d= inner_join(pheno, h1_temp, by= 'PREG_ID_1724') %>% inner_join(., h2_temp, by= 'PREG_ID_1724') %>% inner_join(., h3_temp, by= 'PREG_ID_1724') %>% inner_join(., h4_temp, by= 'PREG_ID_1724')

if (grepl('X', snp)) {

df= filter(d, KJONN== 0)

m1= glm(jaundice~ h1 + h2 + h3 + cohort + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, df, family= binomial)

n= length(resid(m1))
coefs= summary(m1)$coefficients[2:5,]
beta_h1= coefs[1,1]
se_h1= coefs[1,2]
pvalue_h1= coefs[1,4]
beta_h2= coefs[2,1]
se_h2= coefs[2,2]
pvalue_h2= coefs[2,4]
beta_h3= coefs[3,1]
se_h3= coefs[3,2]
pvalue_h3= coefs[3,4]
beta_h4= ''
se_h4= ''
pvalue_h4= ''

results= paste(snp, n, beta_h1, se_h1, pvalue_h1, beta_h2, se_h2, pvalue_h2, beta_h3, se_h3, pvalue_h3, beta_h4, se_h4, pvalue_h4, sep= '\t')

write(results, file= snakemake@output[[1]], append=TRUE)

df= filter(d, KJONN== 1)

m1= glm(jaundice~ h1 + h2 + h4 + cohort + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, df, family= binomial)

n= length(resid(m1))
coefs= summary(m1)$coefficients[2:5,]
beta_h1= coefs[1,1]
se_h1= coefs[1,2]
pvalue_h1= coefs[1,4]
beta_h2= coefs[2,1]
se_h2= coefs[2,2]
pvalue_h2= coefs[2,4]
beta_h3= ''
se_h3= ''
pvalue_h3= '' 
beta_h4= coefs[3, 1]
se_h4= coefs[3, 2]
pvalue_h4= coefs[3,4]

results= paste(snp, n, beta_h1, se_h1, pvalue_h1, beta_h2, se_h2, pvalue_h2, beta_h3, se_h3, pvalue_h3, beta_h4, se_h4, pvalue_h4, sep= '\t')

write(results, file= snakemake@output[[1]], append=TRUE)


} else {

m1= glm(jaundice~ h1 + h2 + h3 + h4 + KJONN + cohort + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, d, family= binomial)

n= length(resid(m1))
coefs= summary(m1)$coefficients[2:5,]
beta_h1= coefs[1,1]
se_h1= coefs[1,2]
pvalue_h1= coefs[1,4]
beta_h2= coefs[2,1]
se_h2= coefs[2,2]
pvalue_h2= coefs[2,4]
beta_h3= coefs[3,1]
se_h3= coefs[3,2]
pvalue_h3= coefs[3,4]
beta_h4= coefs[4,1]
se_h4= coefs[4,2]
pvalue_h4= coefs[4,4]

results= paste(snp, n, beta_h1, se_h1, pvalue_h1, beta_h2, se_h2, pvalue_h2, beta_h3, se_h3, pvalue_h3, beta_h4, se_h4, pvalue_h4, sep= '\t')
write(results, file= snakemake@output[[1]], append=TRUE)


}


}

)



