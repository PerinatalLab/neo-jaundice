library(data.table)
library(dplyr)
library(tidyr)

format_haps= function(hap){
variants= paste(hap$chr, hap$pos, hap$ref, hap$eff, sep =':')
ids= names(hap)[5:ncol(hap)]
hap= as.data.frame(t(hap[, 5:ncol(hap)]))
names(hap)=  variants
hap$IID= ids

return(hap)
}

fets= fread(snakemake@input[[1]])
moms= fread(snakemake@input[[2]])
dads= fread(snakemake@input[[3]])

fets= format_haps(fets)
moms= format_haps(moms)
dads= format_haps(dads)

pheno= fread(snakemake@input[[4]])
covar= fread(snakemake@input[[5]])
trios= fread(snakemake@input[[6]])
trio_ids= readLines(snakemake@input[[7]])
trio_ids= as.numeric(trio_ids)

pheno= inner_join(pheno, covar, by= 'IID') %>% inner_join(., trios, by= c('IID'= 'Child')) 

moms= inner_join(moms, trios, by= c('IID'= 'Mother')) 
dads= inner_join(dads, trios, by= c('IID'= 'Father')) 
fets= inner_join(fets, trios, by= c('IID'= 'Child'))

pheno= filter(pheno, PREG_ID_1724 %in% trio_ids)
print(nrow(pheno))
write( paste('snp', 'n', 'beta_fets', 'se_fets', 'pvalue_fets', 'beta_moms', 'se_moms', 'pvalue_moms', 'beta_dads', 'se_dads', 'pvalue_dads', sep= '\t'), snakemake@output[[1]], append= T)

results_list= lapply(names(fets)[grepl(':', names(fets))], function(snp) {

print(snp)
moms_temp= moms[, c('PREG_ID_1724', snp)]
fets_temp= fets[, c('PREG_ID_1724', snp)]
dads_temp= dads[, c('PREG_ID_1724', snp)]

names(moms_temp)= c('PREG_ID_1724', 'moms')
names(fets_temp)= c('PREG_ID_1724', 'fets')
names(dads_temp)= c('PREG_ID_1724', 'dads')

d= inner_join(pheno, moms_temp, by= 'PREG_ID_1724') %>% inner_join(., fets_temp, by= 'PREG_ID_1724') %>% inner_join(., dads_temp, by= 'PREG_ID_1724')

if (grepl('X', snp)){

d$dads= d$dads * 2
d$fets= ifelse(d$KJONN== 1, d$fets * 2, d$fets)

}
m1= glm(jaundice~ fets + moms + dads + KJONN + cohort + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, d, family= binomial)

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



results= paste(snp, n, beta_h1, se_h1, pvalue_h1, beta_h2, se_h2, pvalue_h2, beta_h3, se_h3, pvalue_h3, sep= '\t')
write(results, file= snakemake@output[[1]], append=TRUE)

}

)

print('Analyses performed, saving data.')

names(fets)[grepl(':', names(fets))]= paste0('fets_'  ,names(fets)[grepl(':', names(fets))])
names(moms)[grepl(':', names(moms))]= paste0('moms_'  ,names(moms)[grepl(':', names(moms))])
names(dads)[grepl(':', names(dads))]= paste0('dads_'  ,names(dads)[grepl(':', names(dads))])

fets= select(fets, -c(IID, Mother, Father))
moms= select(moms, -c(IID, Child, Father))
dads= select(dads, -c(IID, Mother, Child))

x= inner_join(fets, moms, by= 'PREG_ID_1724') %>% inner_join(., dads, by= 'PREG_ID_1724')

fwrite(x, snakemake@output[[2]], sep= '\t')
