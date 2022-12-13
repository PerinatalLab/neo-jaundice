

rule filter_mac:
        'Filter mac for Step 1 of REGENIE'
        input:
                '/mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/imputed/subset/grm-high-quality-pruned/grm-high-quality-pruned.bed',
                'results/aux/ids/samples/{sample}_ids.txt'
        output:
                'results/GWAS/regenie/step1/snp_to_filter/{sample}.snplist'
        params:
                '/mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/imputed/subset/grm-high-quality-pruned/grm-high-quality-pruned',
                'results/GWAS/regenie/step1/snp_to_filter/{sample}'
        shell:
                '''
                /home/pol.sole.navais/soft/plink2 \
                --bfile {params[0]} \
                --mac 100 \
                --write-snplist \
                --keep {input[1]} \
                --out {params[1]}
                '''

rule list_bgen_samples:
        'List of bGEN samples.'
        input:
                '/mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/imputed/all/bgen/{CHR}.bgen'
        output:
                'results/aux/ids/samples/bgen/{CHR}_samples.txt'
        shell:
                '/home/pol.sole.navais/soft/qctool_v2.2.0/qctool -g {input[0]} -os {output[0]}'

rule chrX_to_diploid:
	'Use PLINK2 to convert BGEN haploid males to diplod (pgen file format).'
	input:
		'/mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/imputed/all/vcf/X.vcf.gz'
	output:
		temp(multiext('results/GWAS/pgen/X', '.pgen', '.pvar', '.psam'))
	params:
		'results/GWAS/pgen/X'
	shell:
		"/home/pol.sole.navais/soft/plink2 --vcf {input[0]} dosage='DS' --double-id --make-pgen psam-cols=+fid --out {params[0]}"


rule REGENIE_step1:
        'Whole genome regression model is fit to the traits.'
        input:
                '/mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/imputed/subset/grm-high-quality-pruned/grm-high-quality-pruned.bed',
                'results/pheno/{sample}_pheno_bin.txt',
                'results/pheno/{sample}_covars.txt',
                'results/aux/ids/samples/{sample}_ids.txt',
                'results/GWAS/regenie/step1/snp_to_filter/{sample}.snplist'
        output:
                temp('results/GWAS/regenie/step1/results/{sample}_1.loco.gz'),
                temp('results/GWAS/regenie/step1/results/{sample}_pred.list')
        params:
                '/mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/imputed/subset/grm-high-quality-pruned/grm-high-quality-pruned',
                'results/GWAS/regenie/step1/results/{sample}',
                'results/GWAS/regenie/step1/results/{sample}_temp'
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
                --out {params[1]}
                '''

rule REGENIE_step2:
	'Whole genome regression model is fit to the traits.'
	input:
		'/mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/imputed/all/bgen/{CHR}.bgen',
		'results/pheno/{sample}_pheno_bin.txt',
		'results/pheno/{sample}_covars.txt',
		'results/aux/ids/samples/{sample}_ids.txt',
		'results/GWAS/regenie/step1/results/{sample}_pred.list',
		'results/aux/ids/samples/bgen/{CHR}_samples.txt',
		'results/GWAS/regenie/step1/results/{sample}_1.loco.gz',
		multiext('results/GWAS/pgen/X', '.pgen', '.pvar', '.psam')
	output:
		temp(expand('results/GWAS/regenie/step2/temp/{{sample}}/{{CHR}}_{pheno}.regenie', pheno= pheno_file['phenotypes']))
	params:
		'results/GWAS/regenie/step2/temp/{sample}/{CHR}',
		'results/GWAS/pgen/X'
	threads: 5
	run:
		if wildcards.CHR != 'X':
			shell('''
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
                --verbose
                ''')
		else:
			shell('''
                /home/pol.sole.navais/soft/regenie_v3.2.1.gz_x86_64_Linux \
                --step 2 \
                --pgen {params[1]} \
                --covarFile {input[2]} \
                --keep {input[3]} \
                --phenoFile {input[1]} \
                --bsize 400 \
                --bt \
                --firth --approx \
                --minINFO 0.6 \
                --threads {threads} \
                --pred {input[4]} \
                --out {params[0]} \
                --af-cc \
                --catCovarList cohort \
                --verbose
                ''')

rule concat_GWAS_results:
	'Concatenate results from regenie'
	input:
		expand('results/GWAS/regenie/step2/temp/{{sample}}/{CHR}_{{pheno}}.regenie', CHR= CHROM)
	output:
		temp('results/GWAS/sumstats/GWAS-{pheno}/temp/allchr/{sample}.txt')
	shell:
		'''
		head -1 {input[0]} > {output[0]}
		tail -n +2 -q {input} >> {output[0]}
		'''

rule gzip_results:
        'Gzip results.'
        input:
                'results/GWAS/sumstats/GWAS-{pheno}/temp/allchr/{sample}.txt'
        output:
                'results/GWAS/sumstats/GWAS-{pheno}/allchr-{sample}.txt.gz'
        shell:
                'gzip -c {input[0]} > {output[0]}'

rule check_results:
	''
	input:
		expand('results/GWAS/sumstats/GWAS-{pheno}/allchr-{sample}.txt.gz', pheno= pheno_file['phenotypes'], sample= fam_ids['fam_id'])
	output:
		'results/GWAS/checks/GWAS_performed.txt'
	shell:
		'touch {output[0]}'
