
rule top_SNP_list:
	'List of top rsid to condition on.'
	output:
		'results/conditional_GWAS/regenie/SNP_list/rs6755571.txt'
	shell:
		'echo rs6755571 > {output[0]}'


rule extract_variant:
	'Extract rs6755571 from child bgen.'
	input:
		'/mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/imputed/all/bgen/2.bgen',
		'results/conditional_GWAS/regenie/SNP_list/rs6755571.txt'
	output:
		temp('results/conditional_GWAS/regenie/BGEN/fets_rs6755571.bgen')
	shell:
		'/home/pol.sole.navais/soft/qctool_v2.2.0/qctool -g {input[0]} -incl-rsids {input[1]} -og {output[0]}'

rule REGENIE_step1_conditional:
        'Whole genome regression model is fit to the traits.'
        input:
                '/mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/imputed/subset/grm-high-quality-pruned/grm-high-quality-pruned.bed',
                'results/pheno/{sample}_pheno_bin.txt',
                'results/pheno/{sample}_covars.txt',
                'results/aux/ids/samples/{sample}_ids.txt',
                'results/GWAS/regenie/step1/snp_to_filter/{sample}.snplist',
		'results/conditional_GWAS/regenie/BGEN/fets_rs6755571.bgen',
		'results/conditional_GWAS/regenie/SNP_list/rs6755571.txt',
		'results/aux/ids/samples/bgen/2_samples.txt'
        output:
                temp('results/conditional_GWAS/regenie/step1/results/{sample}_1.loco.gz'),
                temp('results/conditional_GWAS/regenie/step1/results/{sample}_pred.list')
	params:
		'/mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/imputed/subset/grm-high-quality-pruned/grm-high-quality-pruned',
		'results/conditional_GWAS/regenie/step1/results/{sample}',
		'results/conditional_GWAS/regenie/step1/results/{sample}_temp'
	threads: 30
	shell:
                '''
                /home/pol.sole.navais/soft/regenie_v3.2.1.gz_x86_64_Linux \
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
		--condition-list {input[6]} \
		--condition-file bgen,{input[5]} \
		--condition-file-sample {input[7]} \
                --out {params[1]}
                '''

rule REGENIE_step2_conditional:
        'Whole genome regression model is fit to the traits.'
        input:
                '/mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/imputed/all/bgen/{CHR}.bgen',
                'results/pheno/{sample}_pheno_bin.txt',
                'results/pheno/{sample}_covars.txt',
                'results/aux/ids/samples/{sample}_ids.txt',
                'results/conditional_GWAS/regenie/step1/results/{sample}_pred.list',
                'results/aux/ids/samples/bgen/{CHR}_samples.txt',
                'results/conditional_GWAS/regenie/step1/results/{sample}_1.loco.gz',
		'results/conditional_GWAS/regenie/BGEN/fets_rs6755571.bgen',
                'results/conditional_GWAS/regenie/SNP_list/rs6755571.txt',
		'results/aux/ids/samples/bgen/2_samples.txt'
        output:
                temp(expand('results/conditional_GWAS/regenie/step2/temp/{{sample}}/{{CHR}}_{pheno}.regenie', pheno= pheno_file['phenotypes']))
        params:
                'results/conditional_GWAS/regenie/step2/temp/{sample}/{CHR}',
                'results/GWAS/pgen/X'
        threads: 5
        shell:
                """
		/home/pol.sole.navais/soft/regenie_v3.2.1.gz_x86_64_Linux \
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
		--condition-list {input[8]} \
                --condition-file bgen,{input[7]} \
		--condition-file-sample {input[9]} \
                --verbose
                """

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
                expand('results/conditional_GWAS/sumstats/GWAS-{pheno}/{sample}.allchr.txt.gz', pheno= pheno_file['phenotypes'], sample= ['fets'])
        output:
                'results/conditional_GWAS/checks/GWAS_performed.txt'
        shell:
                'touch {output[0]}'

