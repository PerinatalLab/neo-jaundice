
rule extract_top_fetal_SNP:
	''
	input:
		'results/effect_origin/aux/ids/fets_toextract.txt',
		'/mnt/archive/MOBAGENETICS/genotypes-base/imputed/all/vcf/2.vcf.gz'
	output:
		'results/conditional_GWAS/aux/temp/fetal_SNP.txt'
	shell:
		"bcftools query -S {input[0]} -r 2:234649665 -f '%CHROM\t%POS\t%REF\t%ALT[\t%DS]\n' {input[1]} -o {output[0]}"

rule covariate_file_with_fetal_SNP:
	''
	input:
		'results/conditional_GWAS/aux/temp/fetal_SNP.txt',
		'results/effect_origin/aux/ids/fets_toextract.txt',
		'results/effect_origin/aux/ids/parent_offspring_trios.txt',
		'results/pheno/{sample}_pheno_bin.txt',
		'results/pheno/{sample}_covars.txt'
	output:
		'results/pheno/conditional_analysis/{sample}_pheno_bin.txt',
                'results/pheno/conditional_analysis/{sample}_covars.txt'
	run:
		with open(input[0]) as file:
                        ds= [float(i) for i in [line.strip().split('\t') for line in file][0][4:]]
		var= np.where(wildcards.sample == 'moms', 'Mother', np.where(wildcards.sample== 'fets', 'Child', 'Father'))
		with open(input[1]) as file:
			fet_ids = [line.strip() for line in file]
		d= pd.DataFrame({'top_fetal_SNP': ds, 'Child': fet_ids})
		ids= pd.read_csv(input[2], sep= '\t', header= 0)
		ids= ids[['Child', str(var)]]
		ids.columns= ['Child', 'IID']
		d= pd.merge(d, ids, on= 'Child')
		pheno= pd.read_csv(input[3], sep= '\t', header= 0)
		covars= pd.read_csv(input[4], sep= '\t', header= 0)
		covars= pd.merge(covars, d, on= 'IID')
		pheno= pheno.loc[pheno.IID.isin(covars.IID.values), :]
		pheno.drop_duplicates('IID', inplace= True)
		covars.dropna(inplace= True)
		covars.drop_duplicates('IID', inplace= True)
		pheno.to_csv(output[0], sep= '\t', header= True, index= False)
		covars.to_csv(output[1], sep= '\t', header= True, index= False, columns= ['FID', 'IID', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10', 'cohort', 'KJONN', 'top_fetal_SNP'])
		

rule REGENIE_step1_conditional:
        'Whole genome regression model is fit to the traits.'
        input:
                '/mnt/archive/MOBAGENETICS/genotypes-base/imputed/subset/grm-high-quality-pruned/grm-high-quality-pruned.bed',
                'results/pheno/conditional_analysis/{sample}_pheno_bin.txt',
                'results/pheno/conditional_analysis/{sample}_covars.txt',
                'results/aux/ids/samples/{sample}_ids.txt',
                'results/GWAS/regenie/step1/snp_to_filter/{sample}.snplist'
        output:
                temp('results/conditional_GWAS/regenie/step1/results/{sample}_1.loco.gz'),
                temp('results/conditional_GWAS/regenie/step1/results/{sample}_pred.list')
        params:
                '/mnt/archive/MOBAGENETICS/genotypes-base/imputed/subset/grm-high-quality-pruned/grm-high-quality-pruned',
                'results/conditional_GWAS/regenie/step1/results/{sample}',
                'results/conditional_GWAS/regenie/step1/results/{sample}_temp'
        threads: 30
        shell:
                '''
                /home/pol.sole.navais/soft/regenie_v3.1.gz_x86_64_Linux \
                --step 1 \
                --threads {threads} \
                --gz \
                --bed {params[0]} \
                --covarFile {input[2]} \
                --phenoFile {input[1]} \
                --keep {input[3]} \
                --extract {input[4]} \
                --bsize 1000 \
                --bt --lowmem \
                --lowmem-prefix {params[2]} \
                --catCovarList cohort \
                --out {params[1]}
                '''

rule REGENIE_step2_conditional:
        'Whole genome regression model is fit to the traits.'
        input:
                '/mnt/archive/MOBAGENETICS/genotypes-base/imputed/all/bgen/{CHR}.bgen',
                'results/pheno/conditional_analysis/{sample}_pheno_bin.txt',
                'results/pheno/conditional_analysis/{sample}_covars.txt',
                'results/aux/ids/samples/{sample}_ids.txt',
                'results/conditional_GWAS/regenie/step1/results/{sample}_pred.list',
                'results/aux/ids/samples/bgen/{CHR}_samples.txt',
                'results/conditional_GWAS/regenie/step1/results/{sample}_1.loco.gz'
        output:
                temp(expand('results/conditional_GWAS/regenie/step2/temp/{{sample}}/{{CHR}}_{pheno}.regenie', pheno= pheno_file['phenotypes']))
        params:
                'results/conditional_GWAS/regenie/step2/temp/{sample}/{CHR}',
                'results/GWAS/pgen/X'
        threads: 5
        run:
                shell('''
                /home/pol.sole.navais/soft/regenie_v3.1.gz_x86_64_Linux \
                --step 2 \
                --bgen {input[0]} \
                --covarFile {input[2]} \
                --keep {input[3]} \
                --phenoFile {input[1]} \
                --bsize 400 \
                --bt \
                --firth --approx \
                --minINFO 0.6 \
                --threads {threads} \
                --sample {input[5]} \
                --pred {input[4]} \
                --out {params[0]} \
                --af-cc \
                --catCovarList cohort \
                --verbose
                ''')

rule concat_GWAS_results_conditional:
        'Concatenate results from regenie'
        input:
                expand('results/conditional_GWAS/regenie/step2/temp/{{sample}}/{CHR}_{{pheno}}.regenie', CHR= CHROM)
        output:
                temp('results/conditional_GWAS/sumstats/GWAS-{pheno}/temp/{sample}.allchr.txt')
        shell:
                '''
                head -1 {input[0]} > {output[0]}
                tail -n +2 -q {input} >> {output[0]}
                '''

rule gzip_results_conditional_GWAS:
        'Gzip results.'
        input:
                'results/conditional_GWAS/sumstats/GWAS-{pheno}/temp/{sample}.allchr.txt'
        output:
                'results/conditional_GWAS/sumstats/GWAS-{pheno}/{sample}.allchr.txt.gz'
        shell:
                'gzip -c {input[0]} > {output[0]}'

rule check_results_conditional_GWAS:
        ''
        input:
                expand('results/conditional_GWAS/sumstats/GWAS-{pheno}/{sample}.allchr.txt.gz', pheno= pheno_file['phenotypes'], sample= fam_ids['fam_id'])
        output:
                'results/conditional_GWAS/checks/GWAS_performed.txt'
        shell:
                'touch {output[0]}'

