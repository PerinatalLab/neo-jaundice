
rule manhattan_plot:
	'Manhattan plot of GWAS results.'
	input:
		'results/delivery/MoBa-GWAS-{pheno}-{sample}.txt.gz',
		'results/delivery/topregions/loci-{pheno}-{sample}.txt'
	output:
		'results/plots/manhattan/{pheno}-{sample}.png'
	script:
		'../scripts/figures/manhattan.R'

rule check_manhattan:
	''
	input:
		expand('results/plots/manhattan/{pheno}-{sample}.png', pheno= pheno_file['phenotypes'], sample= fam_ids['fam_id'])
	output:
		'results/plots/manhattan/checks/manhattan_performed.txt'
	shell:
		'touch {output[0]}'
