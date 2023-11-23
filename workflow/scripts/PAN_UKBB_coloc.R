library(data.table)
library(parallel)
library(dplyr)
library(coloc)


pph_outfile= snakemake@output[[1]]
results_outfile= snakemake@output[[2]]

cat('nsnps\tPP.H0.abf\tPP.H1.abf\tPP.H2.abf\tPP.H3.abf\tPP.H4.abf\tpreg_trait\tpheno_PAN_UKBB\n', file = pph_outfile)
cat('snp\tV.df\tz.df1\tr.df1\tlABF.df1\tV.df2\tz.df2\tr.df2\tlABF.df2\tinternal.sum.lABF\tSNP.PP.H4\tpreg_trait\tpheno_PAN_UKBB\n', file= results_outfile)

prior1= 1 * 10**-4
prior2= 1 * 10**-4
prior12= 5 * 10**-6


format_PANUKBB= function(df){

pheno= as.character(gsub('.txt.gz', '', unlist(strsplit(snakemake@input[[3]], '/'))[6]))
trait_type= mani[mani$trait== pheno, ]$trait_type
print(paste('This is the trait type:', trait_type))
if (trait_type %in% c('categorical', 'prescriptions', 'phecode', 'icd10')) {
cols= c('chr', 'pos', 'ref', 'alt', 'af_controls_EUR', 'beta_EUR', 'se_EUR', 'neglog10_pval_EUR', 'low_confidence_EUR')

} else {
cols= c('chr', 'pos', 'ref', 'alt', 'af_EUR', 'beta_EUR', 'se_EUR', 'neglog10_pval_EUR', 'low_confidence_EUR')
}
df= data.frame(df)
df= select(df, all_of(cols))

print(pheno)
names(df)= c('chr', 'pos', 'ref', 'alt', 'eaf', 'beta', 'se', 'pvalue', 'low_confidence')

df= filter(df, low_confidence== FALSE)

df$beta= ifelse(df$ref> df$alt, -1 * df$beta, df$beta)

df$ID= ifelse(df$ref> df$alt, paste(df$chr, df$pos, df$alt, df$ref, sep= ':'), paste(df$chr, df$pos, df$ref, df$alt, sep= ':'))

df$maf= ifelse(df$eaf> 0.5, 1 - df$eaf, df$eaf)

df$N= mani[mani$trait== pheno, ]$N

s_pheno= mani[mani$trait== pheno, ]$s_pheno

z= coloc_FINNGEN(d, df, 'Neonatal jaundice', pheno, s_pheno, trait_type)
z= coloc_FINNGEN(bil, df, 'Adult bilirubin', pheno, s_pheno, trait_type)
}

coloc_FINNGEN= function(trait1df, trait2df, trait, phenotype, prop_cases, trait_type){

temp_df= inner_join(trait1df, trait2df, by= 'ID')

temp_df= filter(temp_df, SE>0, se>0)

if (!('MAF' %in% names(temp_df))){
temp_df$MAF= temp_df$maf
}
data1= list(beta= temp_df$BETA, varbeta= temp_df$SE**2, N=temp_df$TOTALSAMPLESIZE, type= 'quant', snp= temp_df$ID, MAF= temp_df$MAF)

if (trait_type %in% c('categorical', 'prescriptions', 'phecode', 'icd10')) {
data2= list(beta= temp_df$beta, varbeta= temp_df$se**2, N= temp_df$N, type= 'cc', snp= temp_df$ID, s= prop_cases, MAF= temp_df$maf)

} else {

data2= list(beta= temp_df$beta, varbeta= temp_df$se**2, N=temp_df$N, type= 'quant', snp= temp_df$ID, MAF= temp_df$maf)

}
print(phenotype)

myres= tryCatch({(coloc.abf(data1, data2, p1= prior1, p2= prior2, p12= prior12))}, error= function(e) { return(0)}

)
        if (length(myres)==1 ) {
        PPH= data.frame(nsnps= 0, PP.H0.abf= 0, PP.H1.abf= 0, PP.H2.abf= 0, PP.H3.abf= 0, PP.H4.abf= 0, preg_trait= trait, pheno_PAN_UKBB= phenotype)
        fwrite(PPH, pph_outfile, sep= '\t', row.names=F, col.names= F, quote=F, append= T)
        res= data.frame(snp= 'none', V.df1= 0, z.df1= 0, r.df1= 0, lABF.df1= 0, V.df2= 0, z.df2= 0, r.df2= 0, lABF.df2= 0, internal.sum.lABF= 0, SNP.PP.H4= 0,  preg_trait= trait, pheno_PAN_UKBB= phenotype)
        fwrite(res, results_outfile, sep= '\t', row.names=F, col.names= F, quote=F, append= T)
        print('next')
        } else {
        PPH= data.frame(t(myres[[1]]))
        PPH$preg_trait= trait
	PPH$pheno_PAN_UKBB= phenotype
        if ((PPH$PP.H3.abf + PPH$PP.H4.abf) >= 0.01) {
        fwrite(PPH, pph_outfile, sep= '\t', row.names=F, col.names= F, quote=F, append= T)
        res= myres[[2]]
        res$preg_trait= trait
        res$pheno_PAN_UKBB= phenotype
	if ((PPH$PP.H3.abf + PPH$PP.H4.abf) >= 0.5) {
        fwrite(res, results_outfile, sep= '\t', row.names=F, col.names= F, quote=F, append= T)
	}
        } else {
        print('Not enough power')
        }
        }


}

d= fread(snakemake@input[[1]])
bil= fread(snakemake@input[[5]])
mani= fread(snakemake@input[[2]])

mani= data.frame(mani)


mani$trait= paste(mani$trait_type, mani$phenocode, mani$pheno_sex, mani$coding, sep= '-')

mani$trait= ifelse(substr(mani[,'trait'], nchar(mani[,'trait']), nchar(mani[,'trait']))== '-',substr(mani$trait, 1, nchar(mani$trait)-1), mani$trait)
mani$trait= gsub(' ', '_', mani$trait)
mani$trait= gsub('/', '_', mani$trait)

mani$s_pheno= with(mani, ifelse(mani$trait_type %in% c('categorical', 'prescriptions', 'phecode', 'icd10'), n_cases_EUR / (n_controls_EUR + n_cases_EUR), 0))
mani$N= with(mani, ifelse(mani$trait_type %in% c('categorical', 'prescriptions', 'phecode', 'icd10'), n_controls_EUR + n_cases_EUR, n_cases_EUR))

df= fread(snakemake@input[[3]], h= F)
header_line= scan(snakemake@input[[4]], what= 'character', sep= '\t')

names(df)= header_line


z= format_PANUKBB(df)

