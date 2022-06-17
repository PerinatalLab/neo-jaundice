library(data.table)
library(dplyr)
library(coloc)
library(parallel)

prior1= 1 * 10**-4
prior2= 1 * 10**-4
prior12= 5 * 10**-6

cat('nsnps\tPP.H0.abf\tPP.H1.abf\tPP.H2.abf\tPP.H3.abf\tPP.H4.abf\tlocus\n', file = snakemake@output[[1]])

cat('snp\tV.df\tz.df1\tr.df1\tlABF.df1\tV.df2\tz.df2\tr.df2\tlABF.df2\tinternal.sum.lABF\tSNP.PP.H4\tlocus\n', file= snakemake@output[[2]])


d= fread(snakemake@input[[1]], select= c('ID', 'CHR', 'POS', 'N', 'BETA', 'SE', 'LOG10P', 'EAF'))

d= filter(d, !duplicated(ID))

d$MAF= ifelse(d$EAF>0.5, 1 - d$EAF, d$EAF)

x= fread(snakemake@input[[2]])

x$BETA= ifelse(x$REF > x$EFF, -1 * x$BETA, x$BETA)

x= select(x, CHROM, POS, EFF, REF, BETA, SE, pvalue)

x$TOTALSAMPLESIZE= 363228

x$EFF= toupper(x$EFF)
x$REF= toupper(x$REF)

x$BETA= ifelse(x$REF> x$EFF, -1 * x$BETA, x$BETA)
x$ID= with(x, ifelse(REF> EFF, paste(CHROM, POS, EFF, REF, sep= ':'), paste(CHROM, POS, REF, EFF, sep= ':')))

x= select(x, ID, TOTALSAMPLESIZE, BETA, SE, pvalue)
names(x)= c('ID', 'TOTALSAMPLESIZE', 'beta', 'se', 'p')

d= inner_join(d, x, by= 'ID')

ld_indep= fread(snakemake@input[[3]])

ld_indep$chr= gsub('chr', '', ld_indep$chr)
ld_indep$chr= as.numeric(ld_indep$chr)

ld_indep$locus= paste(ld_indep$chr, ld_indep$start, ld_indep$stop, sep= '-')


funk= function(i) {
        row= ld_indep[i, ]
        locus1= row[, 'locus']
        temp_df= filter(d, CHR== as.integer(row[, 'chr']), POS >= as.integer(row[, 'start']), POS<= as.integer(row[, 'stop']))
	print(nrow(temp_df))

        if (nrow(temp_df)== 0) {
        PPH= data.frame(nsnps= 0, PP.H0.abf= 0,PP.H1.abf= 0,  PP.H2.abf= 0, PP.H3.abf= 0, PP.H4.abf= 0, locus= locus1)
        res= data.frame(snp= 'none', V.df1= 0, z.df1= 0, r.df1= 0, lABF.df1= 0, V.df2= 0, z.df2= 0, r.df2= 0, lABF.df2= 0, internal.sum.lABF= 0, SNP.PP.H4= 0, locus= locus1)
        fwrite(PPH, snakemake@output[[1]], sep= '\t', row.names=F, col.names= F, quote=F, append= T)
        res= data.frame(snp= 'none', V.df1= 0, z.df1= 0, r.df1= 0, lABF.df1= 0, V.df2= 0, z.df2= 0, r.df2= 0, lABF.df2= 0, internal.sum.lABF= 0, SNP.PP.H4= 0, locus= locus1)
        fwrite(res, snakemake@output[[2]], sep= '\t', row.names=F, col.names= F, quote=F, append= T)
        print('next')

        } else {
        temp_df= filter(temp_df, SE>0, se>0)

        data1= list(beta= temp_df$BETA, varbeta= temp_df$SE**2, N=temp_df$N, type= 'cc', snp= temp_df$ID, MAF= temp_df$MAF, s= 0.069)

        data2= list(beta= temp_df$beta, varbeta= temp_df$se**2, N=temp_df$TOTALSAMPLESIZE, type= 'quant', snp= temp_df$ID, MAF= temp_df$MAF)
        myres= tryCatch({suppressWarnings(coloc.abf(data1, data2, p1= prior1, p2= prior2, p12= prior12))}, error= function(e) { return(0)}
)
        if (length(myres)==1 ) {
        PPH= data.frame(nsnps= 0, PP.H0.abf= 0, PP.H1.abf= 0, PP.H2.abf= 0, PP.H3.abf= 0, PP.H4.abf= 0, locus= locus1)
        res= data.frame(snp= 'none', V.df1= 0, z.df1= 0, r.df1= 0, lABF.df1= 0, V.df2= 0, z.df2= 0, r.df2= 0, lABF.df2= 0, internal.sum.lABF= 0, SNP.PP.H4= 0, locus= locus1)
        fwrite(res, snakemake@output[[2]], sep= '\t', row.names=F, col.names= F, quote=F, append= T)
        print('next')
        next
        } else {
	print(str(myres))
        PPH= as.data.frame(t(myres[[1]]))
        PPH$locus= locus1
        fwrite(PPH, snakemake@output[[1]], sep= '\t', row.names=F, col.names= F, quote=F, append= T)
        res= myres[[2]]
        res$locus= rep(locus1, nrow(res))
        fwrite(res, snakemake@output[[2]], sep= '\t', row.names=F, col.names= F, quote=F, append= T)

}
}
}




lapply(1:nrow(ld_indep), funk)

