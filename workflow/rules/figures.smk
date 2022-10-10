
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

rule QQ_plots:
	'QQ-plots for GWAS.'
	input:
		'results/GWAS/delivery/MoBa-GWAS-jaundice-{sample}.txt.gz'
	output:
		'results/plots/QQ-plot-jaundice-{sample}.png'
	conda:
		'../envs/plots.yml'
	script:
                '../scripts/figures/QQ-plot.R'

rule check_QQ_plot:
	'Rule to check that all QQ plots are done.'
	input:
		expand('results/plots/QQ-plot-jaundice-{sample}.png', sample= fam_ids['fam_id'])
	output:
		'results/plots/checks/QQ-plot.txt'
	shell:
		'touch {output[0]}'

rule manhattan:
	'Manhattan plot of Paternal GWAS.'
	input:
		'results/GWAS/delivery/MoBa-GWAS-jaundice-{sample}.txt.gz'
	output:
		'results/plots/manhattan-{sample}-jaundice.png'
	conda:
		'../envs/plots.yml'
	script:
		'../scripts/figures/manhattan.R'

rule check_manhattan:
        'Manhattan plot of GWAS.'
        input:
                expand('results/plots/manhattan-{sample}-jaundice.png', sample= fam_ids['fam_id'])
        output:
                'results/plots/checks/manhattan-jaundice.txt'
        shell:
                'touch {output[0]}'


rule contrast_polygenicity_plot:
	''
	input:
		'results/HESS/contrast-polygenicity/jaundice-{sample}.txt',
		'results/HESS/contrast-polygenicity/height.txt'
	output:
		'results/plots/polygenicity-contrast-jaundice-{sample}.pdf'
	conda:
		'../envs/plots.yml'
	script:
		'../scripts/figures/contrast_polygenicity.R'

rule check_contrast_polygenicity:
        'Rule to check that all QQ plots are done.'
        input:
                expand('results/plots/polygenicity-contrast-jaundice-{sample}.pdf', sample= fam_ids['fam_id'])
        output:
                'results/plots/checks/HESS-plot.txt'
        shell:
                'touch {output[0]}'


