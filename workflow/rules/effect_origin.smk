

rule list_vcf_ids:
	'Obtain list of IID present in each chromosome.'
	input:
		'/mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/imputed/all/vcf/{CHR}.vcf.gz'
	output:
		temp('results/effect_origin/aux/vcf_ids/temp/{CHR}-ids.txt')
	shell:
		'bcftools query -l {input[0]} > {output[0]}'

rule merge_vcf_ids:
	'Keep only IIDs present in all chromosomes.'
	input:
		expand('results/GWAS/effect_origin/aux/vcf_ids/temp/{CHR}-ids.txt', CHR= CHROM)
	output:
		'results/effect_origin/aux/vcf_ids/allchr-ids.txt'
	run:
		df_list= list()
		for infile in input:
			d= pd.read_csv(infile, header= None, names= ['IID'])
			df_list.append(d)
		d= reduce(lambda x, y: pd.merge(x, y, on = 'IID', how = 'inner'), df_list)
		d.to_csv(output[0], sep= '\t', header= True, index= False)
	

rule list_trio_ids:
        'Make a list of trio IDs with genotype data.'
        input:
                '/mnt/work/pol/MoBaGenetics-1.0/delivery/sample-ids.txt',
                '/mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/aux/flaglist-merged/mobagen-flaglist-n99259.txt',
                '/mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/aux/pca/ethnic-core-samples',
                'results/effect_origin/aux/vcf_ids/allchr-ids.txt',
		'/mnt/work/pol/MoBaGenetics-1.0/delivery/trio-ids.txt'
        output:
                'results/effect_origin/aux/ids/fets_toextract.txt',
                'results/effect_origin/aux/ids/moms_toextract.txt',
                'results/effect_origin/aux/ids/dads_toextract.txt',
                'results/effect_origin/aux/ids/parent_offspring_trios.txt'
        run:
                d= pd.read_csv(input[0], sep= '\t', header= 0)
                flag= pd.read_csv(input[1], sep= '\t', header= 0)
                flag= flag.loc[flag.genotypesOK== True, :]
                flag= flag.loc[flag.phenoOK== True, :]
                pcs= [line.strip() for line in open(input[2], 'r')]
                x= pd.read_csv(input[3], sep= '\t', header= 0)
		fets= format_trios(d, 'Child', flag, pcs, x)
		print(fets.shape)
		moms= format_trios(d, 'Mother', flag, pcs, x)
		dads= format_trios(d, 'Father', flag, pcs, x)
                fets.to_csv(output[0], columns= ['IID'], sep= '\t', header= False, index= False)
                moms.to_csv(output[1], columns= ['IID'], sep= '\t', header= False, index= False)
                dads.to_csv(output[2], columns= ['IID'], sep= '\t', header= False, index= False)
		df= pd.read_csv(input[4], sep= '\t', header= 0)
		df= df.loc[df.Child.isin(pcs) & (df.Mother.isin(pcs)) & (df.Father.isin(pcs)), :]
		df= df.loc[(df.Child.isin(flag.IID.values)) & (df.Mother.isin(flag.IID.values)) & (df.Father.isin(flag.IID.values)), :]
		df= df.loc[df.BATCH != 'TED', ]
		df.drop_duplicates('PREG_ID_1724', inplace= True, keep= 'first')
		df.drop_duplicates('Mother', inplace= True, keep= 'first')
		df.drop_duplicates('Father', inplace= True, keep= 'first')
                df.to_csv(output[3], sep= '\t', header= True, index= False)

rule format_sumstats:
	'Remove non-necessary rows from summary statistics.'
	input:
		expand('results/topregions/delivery/loci-{{pheno}}-{sample}.txt', sample= fam_ids['fam_id'])
	output:
		temp('results/effect_origin/aux/top_signals/{pheno}-regions_to_extract.txt'),
	run:
		df_list= list()
		for infile in input:
			d= pd.read_csv(infile, sep = '\t', header= 0)
			d['POS']= d.POS.apply(int).apply(str)
			d['pos2']= d.POS
			d['CHR']= d.CHR.apply(str)
			d['CHR']= np.where(d.CHR== '23', 'X', d.CHR)
			df_list.append(d)
		d= pd.concat(df_list)
		d.sort_values(['CHR', 'POS'], inplace= True, ascending= True)
		d.to_csv(output[0], header= False, index= False, sep= '\t', columns= ['CHR', 'POS', 'pos2'])

rule get_GT_effect_origin:
        'Extract GT from VCF file for a subset of genetic variants.'
        input:
                'results/effect_origin/aux/top_signals/{pheno}-regions_to_extract.txt',
                'results/effect_origin/aux/ids/{sample}_toextract.txt',
                '/mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/imputed/all/vcf/{CHR}.vcf.gz'
        output:
                temp('results/effect_origin/aux/GT/temp/{pheno}/{sample}_gt{CHR}')
        run:
                shell("bcftools query -S {input[1]} -R {input[0]} -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' {input[2]} -o {output[0]}")

rule add_header_GT_effect_origin:
        'Add header to genotype files.'
        input:
                'results/effect_origin/aux/ids/{sample}_toextract.txt',
                'results/effect_origin/aux/GT/temp/{pheno}/{sample}_gt{CHR}'
        output:
                temp('results/effect_origin/GT/{pheno}/{sample}_GT{CHR}')
        run:
                cols= ['chr','pos','ref','eff'] + [line.strip() for line in open(input[0], 'r')]
                d= pd.read_csv(input[1], header= None, names= cols, sep= '\t')
                d.drop_duplicates(['chr', 'pos'], keep=False, inplace= True)
                d.to_csv(output[0], sep= '\t', header= True, index= False)

rule concat_GT_effect_origin:
        'Collect GT from all CHR.'
        input:
                expand('results/effect_origin/GT/{{pheno}}/{{sample}}_GT{CHR}', CHR= CHROM)
        output:
                'results/effect_origin/GT/allchr/{pheno}-{sample}_GT.txt'
        shell:
                '''
                set +o pipefail;
                head -1 {input[0]} > {output[0]}
                cat {input} | grep -v 'chr' >> {output[0]}
                '''

rule get_allele_transmission_effect_origin:
        'Retrieve allele transmission from family trios (after phasing).'
        input:
                'results/effect_origin/GT/allchr/{pheno}-fets_GT.txt',
                'results/effect_origin/GT/allchr/{pheno}-moms_GT.txt',
                'results/effect_origin/aux/ids/parent_offspring_trios.txt',
                'results/effect_origin/GT/allchr/{pheno}-dads_GT.txt'
        output:
                'results/effect_origin/haplotypes/{pheno}-h1_PREG_ID',
                'results/effect_origin/haplotypes/{pheno}-h2_PREG_ID',
                'results/effect_origin/haplotypes/{pheno}-h3_PREG_ID',
                'results/effect_origin/haplotypes/{pheno}-h4_PREG_ID'
	conda:
		'../envs/plots.yml'
	script:
		'../scripts/phase_by_transmission.py'

rule merge_haplotype_pheno:
	'Merge each haplotype and the pheno file.'
	input:
		'results/pheno/fets_pheno_bin.txt',
                'results/pheno/fets_covars.txt',
		'results/effect_origin/aux/ids/parent_offspring_trios.txt',
		'results/effect_origin/haplotypes/{pheno}-h1_PREG_ID',
                'results/effect_origin/haplotypes/{pheno}-h2_PREG_ID',
                'results/effect_origin/haplotypes/{pheno}-h3_PREG_ID',
                'results/effect_origin/haplotypes/{pheno}-h4_PREG_ID'
	output:
		temp('results/effect_origin/pheno/temp/{pheno}-all_subjects.txt')
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
			x['PREG_ID_1724']= x.index
			df_list.append(x)
		x= reduce(lambda x, y: pd.merge(x, y, on = 'PREG_ID_1724', how = 'inner'), df_list)
		print(x)
		print('Now d')
		print(d)
		x['PREG_ID_1724']= x.PREG_ID_1724.apply(str)
		d['PREG_ID_1724']= d.PREG_ID_1724.apply(str)
		x= pd.merge(x, d, on= 'PREG_ID_1724')
		print(x.columns)
		x.to_csv(output[0], sep= '\t', header= True, index= False)


rule remove_related_effect_origin:
	'Remove related individuals'
	input:
		'results/effect_origin/pheno/temp/{pheno}-all_subjects.txt',
		'/mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/aux/pedigree/mobagen-ethnic-core-samples.kin0'
	output:
		'results/effect_origin/delivery/{pheno}.txt',
		'results/effect_origin/ids/PREG_ID_{pheno}.txt'
	run:
		d= pd.read_csv(input[0], sep= '\t', header= 0)
		remove= selectUnrelated(input[1], d, d.Child)
		d= d.loc[~d.Child.isin(remove), :]
		remove= selectUnrelated(input[1], d, d.Mother)
                d= d.loc[~d.Mother.isin(remove), :]
		remove= selectUnrelated(input[1], d, d.Father)
                d= d.loc[~d.Father.isin(remove), :]
		d= d.sample(frac=1).reset_index(drop=True)
		d.drop_duplicates(subset= ['Mother'], keep= 'first', inplace= True)
		d.drop_duplicates(subset= ['Father'], keep= 'first', inplace= True)
		d.to_csv(output[0], sep= '\t', header= True, index= False)
		d.to_csv(output[1], sep= '\t', header= False, index= False, columns= ['PREG_ID_1724'])

rule linear_hypotheses:
	''
	input:
		'results/effect_origin/haplotypes/{pheno}-h1_PREG_ID',
                'results/effect_origin/haplotypes/{pheno}-h2_PREG_ID',
                'results/effect_origin/haplotypes/{pheno}-h3_PREG_ID',
                'results/effect_origin/haplotypes/{pheno}-h4_PREG_ID',
		'results/effect_origin/delivery/{pheno}.txt'
	output:
		'results/effect_origin/delivery/lh/{pheno}-results.txt'
	conda:
		'../envs/plots.yml'
	script:
		'../scripts/linear_hypotheses.R'

rule get_DS_effect_origin:
        'Extract DS from VCF file for a subset of genetic variants.'
        input:
                'results/effect_origin/aux/top_signals/{pheno}-regions_to_extract.txt',
                'results/effect_origin/aux/ids/{sample}_toextract.txt',
                '/mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/imputed/all/vcf/{CHR}.vcf.gz'
        output:
                temp('results/effect_origin/aux/DS/temp/{pheno}/{sample}_ds{CHR}')
        run:
                shell("bcftools query -S {input[1]} -R {input[0]} -f '%CHROM\t%POS\t%REF\t%ALT[\t%DS]\n' {input[2]} -o {output[0]}")

rule add_header_DS_effect_origin:
        'Add header to genotype files.'
        input:
                'results/effect_origin/aux/ids/{sample}_toextract.txt',
                'results/effect_origin/aux/DS/temp/{pheno}/{sample}_ds{CHR}'
        output:
                temp('results/effect_origin/DS/{pheno}/{sample}_DS{CHR}')
        run:
                cols= ['chr','pos','ref','eff'] + [line.strip() for line in open(input[0], 'r')]
                d= pd.read_csv(input[1], header= None, names= cols, sep= '\t')
                d.drop_duplicates(['chr', 'pos'], keep=False, inplace= True)
                d.to_csv(output[0], sep= '\t', header= True, index= False)

rule concat_DS_effect_origin:
        'Collect DS from all CHR.'
        input:
                expand('results/effect_origin/DS/{{pheno}}/{{sample}}_DS{CHR}', CHR= CHROM)
        output:
                'results/effect_origin/DS/allchr/{pheno}-{sample}_DS.txt'
        shell:
                '''
                set +o pipefail;
                head -1 {input[0]} > {output[0]}
                cat {input} | grep -v 'chr' >> {output[0]}
                '''

rule conditional_analysis_effect_origin:
	''
	input:
		'results/effect_origin/DS/allchr/{pheno}-fets_DS.txt',
		'results/effect_origin/DS/allchr/{pheno}-moms_DS.txt',
		'results/effect_origin/DS/allchr/{pheno}-dads_DS.txt',
		'results/pheno/fets_pheno_bin.txt',
                'results/pheno/fets_covars.txt',
                'results/effect_origin/aux/ids/parent_offspring_trios.txt',
		'results/effect_origin/ids/PREG_ID_{pheno}.txt'
	output:
		'results/effect_origin/delivery/conditional/{pheno}.txt',
		'results/effect_origin/delivery/conditional/dosage-{pheno}.txt'
	conda:
                '../envs/plots.yml'
	script:
		'../scripts/effect_origin_conditional.R'

rule check_pheno_effect_origin:
	''
	input:
		'results/effect_origin/delivery/lh/jaundice-results.txt',
		'results/effect_origin/delivery/conditional/jaundice.txt'
	output:
		'results/effect_origin/delivery/checks/effect_origin_performed.txt'
	shell:
		'touch {output[0]}'

