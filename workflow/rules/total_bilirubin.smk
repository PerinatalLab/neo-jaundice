rule merge_array_imputed_total_bilirubin:
	'Merge Array- and imputed-based GWAS of bilirubin in UK Biobank. doi: 10.1038/s41588-020-00757-z'
	input:
		'resources/Total_bilirubin.array.gz',
		'resources/Total_bilirubin.imp.gz'
	output:
		'resources/bilirubin/delivery/total-bilirubin-GWAS.txt.gz'
	run:
		array_res= pd.read_csv(input[0], sep= '\t', header= 0, usecols= ['#CHROM', 'POS', 'MarkerName', 'ALT', 'REF', 'Effect', 'StdErr', 'P-value'])
		imp_res= pd.read_csv(input[1], sep= '\t', header= 0, usecols= ['#CHROM', 'POS', 'MarkerName', 'ALT', 'REF', 'Effect', 'StdErr', 'P-value', 'MACH_R2', 'MAF'])
		imp_res= imp_res.loc[imp_res.MAF>= 0.01, :]
		imp_res= imp_res.loc[imp_res.MACH_R2>= 0.7, :]
		d= pd.concat([array_res, imp_res])
		print(d.columns)
		d['ALT']= d.ALT.str.upper()
		d['REF']= d.REF.str.upper()
		d= d[['#CHROM', 'POS', 'MarkerName', 'ALT', 'REF', 'Effect', 'StdErr', 'P-value']]
		d.columns= ['CHROM', 'POS', 'RSID', 'EFF', 'REF', 'BETA', 'SE', 'pvalue']
		d.to_csv(output[0], sep= '\t', header= True, index= False, compression= 'gzip')

rule independent_bilirubin_regions:
        'Obtain a file with independent regions for top loci with a radius of 1.5 Mb.'
        input:
                'resources/bilirubin/delivery/total-bilirubin-GWAS.txt.gz',
		'results/GWAS/delivery/MoBa-GWAS-jaundice-fets.txt.gz'
        output:
                'resources/bilirubin/delivery/top-variants/total-bilirubin-top.txt',
		'resources/bilirubin/aux/regions_to_extract.txt'
        run:
                df= pd.read_csv(input[0], sep= '\t')
		df['pvalue']= -np.log10(st.norm.sf(abs(df.BETA / df.SE))*2)
		df= df.loc[df.pvalue> -np.log10(5*10**-8), :]
		df['ID']= np.where(df.REF> df.EFF, df.CHROM.apply(str) + ':' + df.POS.apply(str) + ':' + df.EFF + ':' + df.REF,  df.CHROM.apply(str) + ':' + df.POS.apply(str) + ':' + df.REF + ':' + df.EFF)
		x= pd.read_csv(input[1], sep= '\t', header= 0, usecols= ['ID'])
		df= pd.merge(df, x, on= 'ID', how= 'inner')
		df.sort_values(by= 'pvalue', ascending= False, inplace= True)
                df.drop_duplicates(subset= ['CHROM', 'POS'], keep= 'first', inplace= True)
                df_list= list()
                for chrom in set(df.CHROM):
                        d_temp= df.loc[df.CHROM== chrom, :]
                        positions= d_temp.POS.values
                        for pos in positions:
                                if pos in d_temp.POS.values:
                                        df_list.append(d_temp.loc[d_temp.POS== pos, :])
                                        d_temp= d_temp.loc[(d_temp.POS < pos - (1.5*10**6)) | (d_temp.POS> pos + (1.5 * 10**6)), :]
                                else:
                                        continue
                x= pd.concat(df_list)
                x['CHROM']= x.CHROM.astype(str)
                x['CHROM']= np.where(x.CHROM== '23', 'X', x.CHROM)
		x.to_csv(output[0], sep='\t', header= True, index= False)
		x['POS']= x.POS.apply(int).apply(str)
		x['pos2']= x.POS
		x['CHROM']= x.CHROM.apply(str)
		x.sort_values(['CHROM', 'POS'], inplace= True, ascending= True)
		x.to_csv(output[1], header= False, index= False, sep= '\t', columns= ['CHROM', 'POS', 'pos2'])


rule get_GT_effect_origin_total_bilirubin:
        'Extract GT from VCF file for a subset of genetic variants.'
        input:
                'resources/bilirubin/aux/regions_to_extract.txt',
                'results/effect_origin/aux/ids/{sample}_toextract.txt',
                '/mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/imputed/all/vcf/{CHR}.vcf.gz'
        output:
                temp('results/bilirubin/aux/GT/temp/{sample}_gt{CHR}')
        run:
                shell("bcftools query -S {input[1]} -R {input[0]} -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' {input[2]} -o {output[0]}")

rule add_header_GT_effect_origin_total_bilirubin:
        'Add header to genotype files.'
        input:
                'results/effect_origin/aux/ids/{sample}_toextract.txt',
                'results/bilirubin/aux/GT/temp/{sample}_gt{CHR}'
        output:
                temp('results/bilirubin/GT/{sample}_GT{CHR}')
        run:
                cols= ['chr','pos','ref','eff'] + [line.strip() for line in open(input[0], 'r')]
                d= pd.read_csv(input[1], header= None, names= cols, sep= '\t')
                d.drop_duplicates(['chr', 'pos'], keep=False, inplace= True)
                d.to_csv(output[0], sep= '\t', header= True, index= False)

rule concat_GT_effect_origin_total_bilirubin:
        'Collect GT from all CHR.'
        input:
                expand('results/bilirubin/GT/{{sample}}_GT{CHR}', CHR= CHROM)
        output:
                'results/bilirubin/GT/allchr/{sample}_GT.txt'
        shell:
                '''
                set +o pipefail;
                head -1 {input[0]} > {output[0]}
                cat {input} | grep -v 'chr' >> {output[0]}
                '''

rule get_allele_transmission_effect_origin_total_bilirubin:
        'Retrieve allele transmission from family trios (after phasing).'
        input:
                'results/bilirubin/GT/allchr/fets_GT.txt',
                'results/bilirubin/GT/allchr/moms_GT.txt',
                'results/effect_origin/aux/ids/parent_offspring_trios.txt',
                'results/bilirubin/GT/allchr/dads_GT.txt'
        output:
                'results/bilirubin/haplotypes/bilirubin-h1_PREG_ID',
                'results/bilirubin/haplotypes/bilirubin-h2_PREG_ID',
                'results/bilirubin/haplotypes/bilirubin-h3_PREG_ID',
                'results/bilirubin/haplotypes/bilirubin-h4_PREG_ID'
        script:
                '../scripts/phase_by_transmission.py'

rule merge_haplotype_pheno_total_bilirubin:
        'Merge each haplotype and the pheno file.'
        input:
                'results/pheno/fets_pheno_bin.txt',
                'results/pheno/fets_covars.txt',
                'results/effect_origin/aux/ids/parent_offspring_trios.txt',
                'results/bilirubin/haplotypes/bilirubin-h1_PREG_ID',
                'results/bilirubin/haplotypes/bilirubin-h2_PREG_ID',
                'results/bilirubin/haplotypes/bilirubin-h3_PREG_ID',
                'results/bilirubin/haplotypes/bilirubin-h4_PREG_ID'
        output:
                temp('results/bilirubin/pheno/temp/all_subjects.txt')
        run:
                d= pd.read_csv(input[0], sep= '\t', header= 0)
                covar= pd.read_csv(input[1], sep= '\t', header= 0)
                d= pd.merge(d, covar, on= 'IID')
                ids= pd.read_csv(input[2], sep= '\t', header= 0)
                d= pd.merge(d, ids, left_on= 'IID', right_on= 'Child')
                df_list= list()
                for i in range(3, len(input)):
                        x= pd.read_csv(input[i], sep= '\t', header= 0)
                        varnames= ('chr' + x.chr.apply(str) + '_' + x.pos.apply(str) + '_' + x.ref + '_' + x.eff).values.tolist()
                        x= pd.DataFrame(x.iloc[:, 4:].T)
                        haplo= input[i].split('-')[1].replace('_PREG_ID', '')
                        x.columns= [i + '_' + haplo for i in varnames]
                        x['PREG_ID']= x.index
                        df_list.append(x)
                x= reduce(lambda x, y: pd.merge(x, y, on = 'PREG_ID', how = 'inner'), df_list)
                print(x)
                print('Now d')
                print(d)
                x['PREG_ID']= x.PREG_ID.apply(str)
                d['PREG_ID']= d.PREG_ID.apply(str)
                x= pd.merge(x, d, on= 'PREG_ID')
                print(x.columns)
                x.to_csv(output[0], sep= '\t', header= True, index= False)

rule remove_related_effect_origin_total_bilirubin:
        'Remove related individuals'
        input:
                'results/bilirubin/pheno/temp/all_subjects.txt',
                '/mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/aux/pedigree/mobagen-ethnic-core-samples.kin0'
        output:
                'results/bilirubin/delivery/total-bilirubin.txt'
        run:
                d= pd.read_csv(input[0], sep= '\t', header= 0)
                remove= selectUnrelated(input[1], d, d.Child)
                d= d.loc[~d.Child.isin(remove), :]
                remove= selectUnrelated(input[1], d, d.Mother)
                d= d.loc[~d.Mother.isin(remove), :]
                remove= selectUnrelated(input[1], d, d.Father)
                d= d.loc[~d.Father.isin(remove), :]
                d.to_csv(output[0], sep= '\t', header= True, index= False)

rule linear_hypotheses_total_bilirubin:
        ''
        input:
                'results/bilirubin/haplotypes/bilirubin-h1_PREG_ID',
                'results/bilirubin/haplotypes/bilirubin-h2_PREG_ID',
                'results/bilirubin/haplotypes/bilirubin-h3_PREG_ID',
                'results/bilirubin/haplotypes/bilirubin-h4_PREG_ID',
                'results/bilirubin/delivery/total-bilirubin.txt'
        output:
                'results/bilirubin/delivery/lh/total-bilirubin.txt'
	conda:
		'../envs/plots.yml'
	script:
		'../scripts/linear_hypotheses.R'

