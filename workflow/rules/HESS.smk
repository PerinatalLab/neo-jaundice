
rule format_sumstats_HESS:
	'Format summary statistics for HESS heritability estimation.'
	input:
		'results/delivery/MoBa-GWAS-{pheno}-{sample}.txt.gz'
	output:
		'results/HESS/aux/sumstats/{pheno}-{sample}.txt'
	run:
		d= pd.read_csv(input[0], sep= '\t', header= 0, usecols= ['CHR', 'POS', 'rsid', 'REF', 'EFF', 'N', 'BETA', 'SE'])[['CHR', 'POS', 'rsid', 'REF', 'EFF', 'N', 'BETA', 'SE']]
		d= d.loc[d.rsid != '.', :]
		d= d.loc[d.rsid != '', :]
		d.drop_duplicates(['CHR', 'POS'], inplace= True)
		d['Z']= d.BETA / d.SE
		d.columns= ['CHR', 'BP', 'SNP', 'A2', 'A1', 'N', 'BETA', 'SE', 'Z']
		d.to_csv(output[0], sep= '\t', header= True, index= False, columns= ['SNP', 'CHR', 'BP', 'A1', 'A2', 'Z', 'N'])


rule HESS_step1:
	'Run HESS step1.'
	input:
		'results/HESS/aux/sumstats/{pheno}-{sample}.txt',
		'resources/ld_indep_regions.txt',
		multiext('/mnt/archive/MOBAGENETICS/genotypes-base/imputed/all/plink/{autoCHR}', '.bim', '.bed', '.fam') 
	output:
		temp(multiext('results/HESS/step1/{pheno}-{sample}_chr{autoCHR}', '.info.gz', '.eig.gz', '.prjsq.gz', '.log'))
	params:
		'/mnt/archive/MOBAGENETICS/genotypes-base/imputed/all/plink/{autoCHR}',
		'results/HESS/step1/{pheno}-{sample}'
	threads: 10
	conda:
		'../envs/HESS.yaml'
	shell:
		'''
		python2 ~/soft/hess-0.5.3-beta/hess.py \
		--local-hsqg {input[0]} \
		--chrom {wildcards.autoCHR} \
		--bfile {params[0]} \
		--partition {input[1]} \
		--out {params[1]}
		'''

rule HESS_step2:
	'Run HESS step2.'
	input:
		expand('results/HESS/step1/{{pheno}}-{{sample}}_chr{autoCHR}.info.gz', autoCHR= autosomal_CHR),
                expand('results/HESS/step1/{{pheno}}-{{sample}}_chr{autoCHR}.eig.gz', autoCHR= autosomal_CHR),
                expand('results/HESS/step1/{{pheno}}-{{sample}}_chr{autoCHR}.prjsq.gz', autoCHR= autosomal_CHR),
                expand('results/HESS/step1/{{pheno}}-{{sample}}_chr{autoCHR}.log', autoCHR= autosomal_CHR)
	output:
		multiext('results/HESS/step2/{pheno}-{sample}', '.log', '.txt')
	params:
		'results/HESS/step1/{pheno}-{sample}',
		'results/HESS/step2/{pheno}-{sample}'
	threads: 10
	conda:
		'../envs/HESS.yaml'
	shell:
		'''
		python2 ~/soft/hess-0.5.3-beta/hess.py --prefix {params[0]} --out {params[1]}
		'''

rule check_HESS:
	'Check HESS is created.'
	input:
		expand('results/HESS/step2/{pheno}-{sample}.{ext}', ext= ['log', 'txt'], pheno= pheno_file['phenotypes'], sample= fam_ids['fam_id'])
	output:
		'results/HESS/checks/h2_estimated.txt'
	shell:
		'touch {output[0]}'
