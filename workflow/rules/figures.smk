
rule manhattan_plot_mother_child:
	'Manhattan plot of GWAS results.'
	input:
		'results/GWAS/delivery/MoBa-GWAS-jaundice-fets.txt.gz',
		'results/GWAS/delivery/MoBa-GWAS-jaundice-moms.txt.gz'
	output:
		'results/plots/manhattan/jaundice-manhattan-mother-child.png'
	conda:
		'../envs/plots.yml'
	script:
		'../scripts/figures/manhattan-mother-child.R'

rule manhattan_paternal:
	''
	input:
		expand('results/plots/manhattan/{pheno}-{sample}.png', pheno= pheno_file['phenotypes'], sample= fam_ids['fam_id'])
	output:
		'results/plots/manhattan/checks/manhattan_performed.txt'
	shell:
		'touch {output[0]}'
