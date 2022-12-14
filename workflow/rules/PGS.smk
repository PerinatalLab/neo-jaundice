
rule format_sumstats_PGS_catalog:
	'Format summary statistics from GWAS Catalog.'
	input:
		'resources/PGS002160.txt.gz',
		'results/GWAS/delivery/MoBa-GWAS-{pheno}-{sample}.txt.gz'
	output:
		'results/PGS/aux/weights/betas-{pheno}-{sample}.txt',
		'results/PGS/aux/variant_id/{pheno}-{sample}.txt'
	run:
		d= pd.read_csv(input[0], sep= '\t', header= 0, comment= '#')
		d['ID']= np.where(d.effect_allele < d.other_allele, d.chr_name.apply(str) + ':' + d.chr_position.apply(str) + ':' + d.effect_allele + ':' + d.other_allele, d.chr_name.apply(str) + ':' + d.chr_position.apply(str) + ':' + d.other_allele + ':' + d.effect_allele)
		x= pd.read_csv(input[1], sep= '\t', header= 0, usecols= ['POS', 'ID'])
		x.drop_duplicates('ID', inplace= True)
		d= d.loc[d.ID.isin(x.ID.values), :]
		d= d[['chr_name', 'chr_position', 'other_allele', 'effect_allele', 'effect_weight']]
		d.columns= ['chr', 'pos', 'REF', 'EFF', 'beta']
		d.to_csv(output[0], sep= '\t', header= True, index= False)
		d.to_csv(output[1], sep= '\t', columns= ['chr', 'pos', 'pos'], header= False, index= False)

rule extract_DS_PGS:
	''
	input:
		'results/PGS/aux/variant_id/{pheno}-{sample}.txt',
		'results/effect_origin/aux/ids/{sample}_toextract.txt',
		'/mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/imputed/all/vcf/{autoCHR}.vcf.gz'
	output:
		temp('results/PGS/aux/DS/temp/{pheno}/{sample}_ds{autoCHR}')
	shell:
		"bcftools query -S {input[1]} -R {input[0]} -f '%CHROM\t%POS\t%REF\t%ALT[\t%DS]\n' {input[2]} -o {output[0]}"

rule calculate_PGS:
	'Calculate PGS using dosage.'
	input:
		'results/effect_origin/aux/ids/{sample}_toextract.txt',
		'results/PGS/aux/DS/temp/{pheno}/{sample}_ds{autoCHR}',
		'results/PGS/aux/weights/betas-{pheno}-{sample}.txt'
	output:
		'results/PGS/aux/DS/temp/{autoCHR}-{pheno}-{sample}.txt'
	script:
		'../scripts/calculate_PGS.py'

rule concat_CHR_PGS:
        'Concat PGS from all CHR for each sample.'
        input:
                expand('results/PGS/aux/DS/temp/{autoCHR}-{{pheno}}-{{sample}}.txt', autoCHR= autosomal_CHR)
        output:
                temp('results/PGS/temp/PGS-{sample}-{pheno}-tosum.txt')
        shell:
                '''
                head -1 {input[0]} > {output[0]}
                tail -n +2 -q {input} >> {output[0]}
                '''

rule sum_PGS:
        'Sum chromosome-based GRS for each sample.'
        input:
                'results/PGS/temp/PGS-{sample}-{pheno}-tosum.txt'
        output:
                'results/PGS/delivery/{sample}-{pheno}-PGS.txt'
        run:
                df= pd.read_csv(input[0], sep= '\t', header= 0)
		df.columns= ['IID', df.columns[0]]
                df= df.groupby('IID').sum().reset_index()
                df.to_csv(output[0], sep= '\t', header= True, index= False)


rule overlapping_variants:
	''
	input:
		expand('results/PGS/aux/weights/betas-{{pheno}}-{sample}.txt', sample= fam_ids['fam_id']),
                expand('results/PGS/aux/variant_id/{{pheno}}-{sample}.txt', sample= fam_ids['fam_id'])
	output:
		'results/PGS/aux/weights/betas-{pheno}-all-samples.txt',
		'results/PGS/aux/variant_id/{pheno}-all-samples.txt'
	run:
		df_list= list()
		beta_files= [i for i in input if 'betas' in i]
		for infile in beta_files:
			d= pd.read_csv(infile, sep= '\t', header= 0)
			df_list.append(d)
		d= reduce(lambda df1, df2: pd.merge(df1, df2, on= ['chr', 'pos', 'REF', 'EFF', 'beta']), df_list)
		d.to_csv(output[0], sep= '\t', header= True, index= False)
		df_list= list()
		variant_ids= [i for i in input if i not in beta_files]
		for infile in variant_ids:
			d= pd.read_csv(infile, sep= '\t', header= None, names= ['chr', 'pos1', 'pos2'])
                        df_list.append(d)
		d= reduce(lambda df1, df2: pd.merge(df1, df2, on= ['chr', 'pos1', 'pos2']), df_list)
                d.to_csv(output[1], sep= '\t', header= False, index= False)


rule calculate_PGS_nochr2:
        'Calculate PGS using dosage.'
        input:
                'results/effect_origin/aux/ids/{sample}_toextract.txt',
                'results/PGS/aux/DS/temp/{pheno}/{sample}_ds{autoCHR}',
                'results/PGS/aux/weights/betas-{pheno}-{sample}.txt'
        output:
                'results/PGS/aux/DS/temp/nochr2/{autoCHR}-{pheno}-{sample}.txt'
        script:
                '../scripts/calculate_PGS.py'

rule concat_CHR_PGS_no_chr2:
        'Concat PGS from all CHR for each sample. Exclude CHROM 2'
        input:
                expand('results/PGS/aux/DS/temp/nochr2/{autoCHR}-{{pheno}}-{{sample}}.txt', autoCHR= [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22])
        output:
                temp('results/PGS/temp/PGS-{sample}-{pheno}-tosum-nochr2.txt')
        shell:
                '''
                head -1 {input[0]} > {output[0]}
                tail -n +2 -q {input} >> {output[0]}
                '''

rule sum_PGS_no_chr2:
        'Sum chromosome-based GRS for each sample. Exclude CHROM 2'
        input:
                'results/PGS/temp/PGS-{sample}-{pheno}-tosum-nochr2.txt'
        output:
                'results/PGS/delivery/{sample}-{pheno}-PGS-nochr2.txt'
        run:
                df= pd.read_csv(input[0], sep= '\t', header= 0)
		df.columns= ['IID', df.columns[0]]
                df= df.groupby('IID').sum().reset_index()
                df.to_csv(output[0], sep= '\t', header= True, index= False)

rule extract_GT_PGS:
	''
	input:
		'results/PGS/aux/variant_id/{pheno}-all-samples.txt',
		'results/effect_origin/aux/ids/{sample}_toextract.txt',
		'/mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/imputed/all/vcf/{autoCHR}.vcf.gz'
	output:
		temp('results/PGS/aux/GT/temp/{pheno}/vcf/{sample}_gt{autoCHR}.vcf'),
		temp('results/PGS/aux/GT/temp/{pheno}/{sample}_gt{autoCHR}')
	shell:
		"""
		tabix -h -R {input[0]} {input[2]} > {output[0]}
		bcftools query -S {input[1]} -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' {output[0]} -o {output[1]}
		"""

rule add_header_GT_PGS:
        'Add header to genotype files.'
        input:
                'results/effect_origin/aux/ids/{sample}_toextract.txt',
                'results/PGS/aux/GT/temp/{pheno}/{sample}_gt{autoCHR}'
        output:
                temp('results/PGS/aux/temp/with_header/{pheno}-{sample}-{autoCHR}')
        run:
                cols= ['chr','pos','ref','eff'] + [line.strip() for line in open(input[0], 'r')]
		d= pd.DataFrame(columns= cols)
#                d= pd.read_csv(input[1], header= None, names= cols, sep= '\t')
#                d.drop_duplicates(['chr', 'pos'], keep=False, inplace= True)
                d.to_csv(output[0], sep= '\t', header= True, index= False)
		shell('cat {input[1]} >> {output[0]} ')

rule get_allele_transmission_GT_PGS:
        'Retrieve allele transmission from family trios (after phasing).'
        input:
                'results/PGS/aux/temp/with_header/{pheno}-fets-{autoCHR}',
                'results/PGS/aux/temp/with_header/{pheno}-moms-{autoCHR}',
                'results/effect_origin/aux/ids/parent_offspring_trios.txt',
                'results/PGS/aux/temp/with_header/{pheno}-dads-{autoCHR}'
        output:
                temp('results/PGS/aux/temp/haplotypes/CHR/h1_PREG_ID_{autoCHR}_{pheno}'),
                temp('results/PGS/aux/temp/haplotypes/CHR/h2_PREG_ID_{autoCHR}_{pheno}'),
                temp('results/PGS/aux/temp/haplotypes/CHR/h3_PREG_ID_{autoCHR}_{pheno}'),
                temp('results/PGS/aux/temp/haplotypes/CHR/h4_PREG_ID_{autoCHR}_{pheno}')
	threads: 6
        script:
                '../scripts/allele_transmission.py'

rule calculate_haplotype_PGS:
        'Calculate PGS using dosage.'
        input:
                'results/PGS/aux/temp/haplotypes/CHR/{haplo}_PREG_ID_{autoCHR}_{pheno}',
                'results/PGS/aux/weights/betas-{pheno}-all-samples.txt'
        output:
                'results/PGS/aux/GT/haplotypes/temp/{autoCHR}-{pheno}-{haplo}.txt'
        script:
                '../scripts/calculate_PGS.py'

rule concat_CHR_haplotype_PGS:
        'Concat PGS from all CHR for each sample.'
        input:
                expand('results/PGS/aux/GT/haplotypes/temp/{autoCHR}-{{pheno}}-{{haplo}}.txt', autoCHR= autosomal_CHR)
        output:
                temp('results/PGS/temp/haplotypes/{haplo}-{pheno}-tosum.txt')
        shell:
                '''
                head -1 {input[0]} > {output[0]}
                tail -n +2 -q {input} >> {output[0]}
                '''

rule sum_haplotype_PGS:
        'Sum chromosome-based GRS for each sample.'
        input:
                'results/PGS/temp/haplotypes/{haplo}-{pheno}-tosum.txt'
        output:
                'results/PGS/delivery/{haplo}-{pheno}-PGS.txt'
        run:
                df= pd.read_csv(input[0], sep= '\t', header= 0)
                df= df.groupby('PREG_ID').sum().reset_index()
                df.to_csv(output[0], sep= '\t', header= True, index= False)

rule check_PGS:
	'Check all PGS files are created.'
	input:
		expand('results/PGS/delivery/{sample}-{pheno}-PGS.txt', sample= 'fets', pheno= pheno_file['phenotypes']),
#		expand('results/PGS/delivery/{haplo}-{pheno}-PGS.txt', haplo= haplotypes, pheno= pheno_file['phenotypes'])
	output:
		'results/PGS/delivery/checks/PGS_performed.txt'
	shell:
		'touch {output[0]}'

rule PGS_haplotype_analysis:
	''
	input:
		'results/PGS/delivery/h1-jaundice-PGS.txt',
		'results/PGS/delivery/h2-jaundice-PGS.txt',
		'results/PGS/delivery/h3-jaundice-PGS.txt',
		'results/PGS/delivery/h4-jaundice-PGS.txt',
		'results/effect_origin/delivery/jaundice.txt'
	output:
		'results/PGS/delivery/glm/PGS_haplotype.txt'
	conda:
		'../envs/plots.yml'
	script:
		'../scripts/PGS_haplotype.R'

