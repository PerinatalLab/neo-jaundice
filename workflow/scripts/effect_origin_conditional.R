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

pheno= inner_join(pheno, covar, by= 'IID') %>% inner_join(., trios, by= c('IID'= 'Child')) 

moms= inner_join(moms, trios, by= c('IID'= 'Mother')) 
dads= inner_join(dads, trios, by= c('IID'= 'Father')) 
fets= inner_join(fets, trios, by= c('IID'= 'Child'))

pheno$PREG_ID= as.character(pheno$PREG_ID)
moms$PREG_ID= as.character(moms$PREG_ID)
dads$PREG_ID= as.character(dads$PREG_ID)
fets$PREG_ID= as.character(fets$PREG_ID)

print(nrow(pheno))
write( paste('snp', 'n', 'beta_fets', 'se_fets', 'pvalue_fets', 'beta_moms', 'se_moms', 'pvalue_moms', 'beta_dads', 'se_dads', 'pvalue_dads', sep= '\t'), snakemake@output[[1]], append= T)

results_list= lapply(names(fets)[1:(length(names(fets))-4)], function(snp) {

print(snp)
moms_temp= moms[, c('PREG_ID', snp)]
fets_temp= fets[, c('PREG_ID', snp)]
dads_temp= dads[, c('PREG_ID', snp)]

names(moms_temp)= c('PREG_ID', 'moms')
names(fets_temp)= c('PREG_ID', 'fets')
names(dads_temp)= c('PREG_ID', 'dads')

d= inner_join(pheno, moms_temp, by= 'PREG_ID') %>% inner_join(., fets_temp, by= 'PREG_ID') %>% inner_join(., dads_temp, by= 'PREG_ID')

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

names(fets)[1:(length(names(fets))-4)]= paste0('fets_'  ,names(fets)[1:(length(names(fets))-4)])
names(moms)[1:(length(names(moms))-4)]= paste0('moms_'  ,names(moms)[1:(length(names(moms))-4)])
names(dads)[1:(length(names(dads))-4)]= paste0('dads_'  ,names(dads)[1:(length(names(dads))-4)])

fets= select(fets, -c(IID, Mother, Father))
moms= select(moms, -c(IID, Child, Father))
dads= select(dads, -c(IID, Mother, Child))

x= inner_join(fets, moms, by= 'PREG_ID') %>% inner_join(., dads, by= 'PREG_ID')

fwrite(x, snakemake@output[[2]], sep= '\t')
