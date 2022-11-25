library(data.table)
library(dplyr)
library(coloc)
library(parallel)

pph_outfile= snakemake@output[[1]]
results_outfile= snakemake@output[[2]]

cat('nsnps\tPP.H0.abf\tPP.H1.abf\tPP.H2.abf\tPP.H3.abf\tPP.H4.abf\tgene\teqtl_data\n', file = snakemake@output[[1]])

cat('snp\tV.df\tz.df1\tr.df1\tlABF.df1\tV.df2\tz.df2\tr.df2\tlABF.df2\tinternal.sum.lABF\tSNP.PP.H4\tgene\teqtl_data\n', file= snakemake@output[[2]])

prior1= 1 * 10**-4
prior2= 1 * 10**-4
prior12= 5 * 10**-6

s_cc= ifelse(grepl('fets', snakemake@input[[1]]), 0.069, ifelse(grepl('moms', snakemake@input[[1]]),  0.085, 0.086))

eqtl_coloc= function(temp_df, gene, eqtl_data){
if (nrow(temp_df)== 0) {
        
        PPH= data.frame(nsnps= 0, PP.H0.abf= 0,PP.H1.abf= 0, PP.H2.abf= 0, PP.H3.abf= 0, PP.H4.abf= 0, geneid= gene, eqtl_data_n= eqtl_data)
        fwrite(PPH, pph_outfile, sep= '\t', row.names=F, col.names= F, quote=F, append= T)
	res= data.frame(snp= 'none', V.df1= 0, z.df1= 0, r.df1= 0, lABF.df1= 0, V.df2= 0, z.df2= 0, r.df2= 0, lABF.df2= 0, internal.sum.lABF= 0, SNP.PP.H4= 0,  geneid= gene, eqtl_data_n= eqtl_data)
        fwrite(res, results_outfile, sep= '\t', row.names=F, col.names= F, quote=F, append= T)    
        print('No data available.')

        } else {

        data1= list(beta= temp_df$BETA, varbeta= temp_df$SE**2, N= temp_df$TOTALSAMPLESIZE, type= 'cc', snp= temp_df$POS, s= s_cc, MAF= temp_df$MAF)

        data2= list(beta= temp_df$beta, varbeta= temp_df$se**2, N= temp_df$N, type= 'quant', snp= temp_df$POS, MAF= temp_df$maf)

        myres= tryCatch({suppressWarnings(coloc.abf(data1, data2, p1= prior1, p2= prior2, p12= prior12))}, error= function(e) { return(0)}
)       
        if (length(myres)==1 ) {  
        PPH= data.frame(nsnps= 0, PP.H0.abf= 0, PP.H1.abf= 0, PP.H2.abf= 0, PP.H3.abf= 0, PP.H4.abf= 0, geneid= gene, eqtl_data_n= eqtl_data)
        
fwrite(PPH, pph_outfile, sep= '\t', row.names=F, col.names= F, quote=F, append= T)
        
res= data.frame(snp= 'none', V.df1= 0, z.df1= 0, r.df1= 0, lABF.df1= 0, V.df2= 0, z.df2= 0, r.df2= 0, lABF.df2= 0, internal.sum.lABF= 0, SNP.PP.H4= 0, geneid= gene, eqtl_data_n= eqtl_data)
        
	fwrite(res, results_outfile, sep= '\t', row.names=F, col.names= F, quote=F, append= T)

        print(paste0('Error in coloc.abf. Nrow temp_df:', nrow(temp_df)))

        } else {

        PPH= data.frame(t(myres[[1]]))
        PPH$geneid= gene
	PPH$eqtl_data_n= eqtl_data
        if ((PPH$PP.H3.abf + PPH$PP.H4.abf) >= 0.01) {
        fwrite(PPH, pph_outfile, sep= '\t', row.names=F, col.names= F, quote=F, append= T)
        res= myres[[2]]
        res$geneid= gene
	res$eqtl_data_n= eqtl_data
        fwrite(res, results_outfile, sep= '\t', row.names=F, col.names= F, quote=F, append= T)
        } else {
        print('Not enough power')
        }

}
}
}

format_eqtl= function(temp_df){
	gene= unique(temp_df$gene_id)
	eqtl_data= unique(temp_df$eqtl_data)
	temp_df = filter(temp_df, SE>0, se> 0)
	print(nrow(temp_df))
	eqtl_coloc(temp_df, gene, eqtl_data)
	
}


d= fread(snakemake@input[[1]])
d= filter(d, !duplicated(POS))

df= fread(snakemake@input[[2]], h= T)
df= group_by(df, gene_id) %>% filter(!duplicated(position)) %>% ungroup()

if (nrow(df)== 0) {

print('No data for the genes selected.')
} else {

eqtl_data_name= snakemake@wildcards[[1]]

df$eqtl_data= eqtl_data_name

d= inner_join(d, df, by= c('POS'= 'position'))

d= filter(d, EFF == alt | EFF == ref)

(print(paste('# rows:', nrow(d))))

z= mclapply(split(d, d$gene_id), format_eqtl, mc.cores= 3)
 
} 
