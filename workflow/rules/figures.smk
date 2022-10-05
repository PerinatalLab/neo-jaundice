
rule manhattan_plot_mother_child:
	'Manhattan plot of GWAS results.'
	input:
		'results/GWAS/delivery/MoBa-GWAS-jaundice-fets.txt.gz',
		'results/GWAS/delivery/MoBa-GWAS-jaundice-moms.txt.gz'
	output:
		'results/plots/jaundice-manhattan-mother-child.png'
	conda:
		'../envs/plots.yml'
	script:
		'../scripts/figures/manhattan-mother-child.R'

rule UGT_P24T_jaundice:
	'Plot of fetal UGT1A4 missense variant and proportion of neonatal jaundice.'
	input:
		'results/UGT-missense/delivery/jaundice.txt',
		'results/UGT-missense/delivery/jaundice-transmitted.txt'
	output:
		'results/plots/fetal-UGT1A4-missense-jaundice-prevalence.pdf',
		'results/plots/fetal-UGT1A4-missense-EAF.pdf',
		'results/plots/fetal-UGT1A4-missense-effect-origin.pdf'
	conda:
                '../envs/plots.yml'
	script:
		'../scripts/figures/fetal-UGT1A4-missense.R'


rule ABO_effect:
	'Effect of rs687621 on neonatal jaundice with and without adjusting for ABO blood group incompatibility.'
	input:
		'/mnt/work/pol/neo-jaundice/results/UGT-missense/delivery/jaundice-transmitted.txt'
	output:
		'results/plots/ABO-alleles.pdf'
	conda:
		'../envs/plots.yml'
	script:
		'../scripts/figures/ABO-maternal.R'

rule UGT_P24T_interaction:
	'Plot showing the interaction between the effect of the missense variant at UGT1A4 and gestational duration and maternal-fetal ABO incompatibility.'
	input:
		'/mnt/work/pol/neo-jaundice/results/UGT-missense/delivery/jaundice.txt'
	output:
		'results/plots/fetal-UGT1A4-GA-interaction.pdf',
		'results/plots/fetal-UGT1A4-GA-interaction-density.pdf',
		'results/plots/fetal-UGT1A4-ABO-interaction.pdf'
	conda:
		'../envs/plots.yml'
	script:
		'../scripts/figures/UGT-interactions.R'

rule parental_UGT:
	'Plot of parental transmitted and non-transmitted effects on neonatal jaundice for variants identified in the maternal and paternal genome at UGT1A* gene region.'
	input:
		'/mnt/work/pol/neo-jaundice/results/UGT-missense/delivery/jaundice-transmitted.txt'
	output:
		'results/plots/maternal-UGT-alleles.pdf',
		'results/plots/paternal-UGT-alleles.pdf'
	conda:
		'../envs/plots.yml'
	script:
		'../scripts/figures/parental-UGT-variants.R'

rule PGS:
	'Plots for the adult bilirubin PGS.'
	input:
		'/mnt/work/pol/neo-jaundice/results/UGT-missense/delivery/jaundice.txt',
		'/mnt/work/pol/neo-jaundice/results/UGT-missense/delivery/jaundice-transmitted.txt'
	output:
		'results/plots/adult-bilirubin-PGS-distribution.pdf',
		'results/plots/adullt-bilirubin-PGS-jaundice.pdf',
		'results/plots/adult-bilirubin-PGS-phased-alleles.pdf'
	conda:
		'../envs/plots.yml'
	script:
		'../scripts/figures/PGS.R'
