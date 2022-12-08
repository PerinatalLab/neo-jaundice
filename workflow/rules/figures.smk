
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
		'results/UGT-missense/delivery/jaundice-transmitted.txt'
	output:
		'results/plots/ABO-alleles.pdf',
		'results/plots/ABO-alleles-not-adjusted.pdf'
	conda:
		'../envs/plots.yml'
	script:
		'../scripts/figures/ABO-maternal.R'

rule chrX_locus:
	'Haplotype-based analysis for the fetal chromosome X locus.'
	input:
		'results/UGT-missense/delivery/jaundice-transmitted.txt',
		'results/effect_origin/delivery/conditional/dosage-jaundice.txt'
	output:
		'results/plots/fetal-chrX-haplotype.pdf',
                'results/plots/fetal-chrX-conditional.pdf'
	conda:
                '../envs/plots.yml'
	script:
		'../scripts/figures/chrX_locus_haplotype.R'

rule UGT_P24T_interaction:
	'Plot showing the interaction between the effect of the missense variant at UGT1A4 and gestational duration and maternal-fetal ABO incompatibility.'
	input:
		'results/UGT-missense/delivery/jaundice.txt'
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
	'Plots for the adult bilirubin PGS and locus zoom plot.'
	input:
		'results/GWAS/delivery/MoBa-GWAS-jaundice-fets.txt.gz',
		'resources/bilirubin/delivery/total-bilirubin-GWAS.txt.gz',
		'results/LD/delivery/UGT1A.ld',
		'results/UGT-missense/delivery/jaundice.txt'
	output:
		'results/plots/locus-zoom-adult-bilirubin-PGS-jaundice.png',
		'results/plots/adult-bilirubin-PGS-distribution.pdf'
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

rule eqtl_colocalization:
	'Circular Column plot for the tissue-agnostic colocalization with all UGT1A* genes.'
	input:
		'results/eQTL_catalogue/delivery/pph-jaundice-fets.txt'
	output:
		'results/plots/jaundice-fets-eQTL-coloc-{UGT_genes}.pdf'
	conda:
		'../envs/plots.yml'
	script:
		'../scripts/figures/circular_eQTL_coloc.R'

rule check_eQTL_UGT_genes:
	'Rule to check that all eqtl colocalization plots are done.'
	input:
		expand('results/plots/jaundice-fets-eQTL-coloc-{UGT_genes}.pdf', UGT_genes= ['UGT1A6', 'UGT1A9', 'UGT1A1', 'UGT1A8', 'UGT1A10', 'UGT1A7', 'UGT1A4'])
	output:
		'results/plots/checks/jaundice-fets-eQTL-coloc-plot.txt'
	shell:
		'touch {output[0]}'

rule UGT1_locuszoom:
	'Locus zoom plot of UGT1A* genes region for liver and colon eQTL and neonatla jaundice GWAS.'
	input:
		'resources/Homo_sapiens.GRCh37.87.chromosome.2.gff3.gz',
		'results/GWAS/delivery/MoBa-GWAS-jaundice-fets.txt.gz',
		'results/eQTL_catalogue/data/hg38/formatted/GTEx_ge_colon_transverse.txt.gz',
		'results/eQTL_catalogue/data/hg38/formatted/GTEx_ge_liver.txt.gz',
		'results/eQTL_catalogue/jaundice/temp/hg38/fets-GWAS.txt',
		'results/LD/delivery/UGT1A.ld'
	output:
		'results/plots/UGT1A-eQTL-locus-zoom.png'
	conda:
		'../envs/plots.yml'
	script:
		'../scripts/figures/UGT-locuszoom.R'


rule PIP_UGT1A:
	''
	input:
		'results/eQTL_catalogue/delivery/SNP-jaundice-fets.txt',
		'results/eQTL_catalogue/delivery/pph-jaundice-fets.txt',
		'results/eQTL_catalogue/jaundice/temp/hg38/fets-GWAS.txt',
		'results/LD/delivery/UGT1A.ld'
	output:
		'results/plots/UGT1A-jaundice-eQTL-correlation.pdf'
	conda:
		'../envs/plots.yml'
	script:
		'../scripts/figures/UGT1-jaundice-eQTL-correlation.R'

rule replication_locus_zoom:
	'Locus zoom plot of UGT1A region for the Danish replication study.'
	input:
		'resources/Homo_sapiens.GRCh37.87.chromosome.2.gff3.gz',
                'results/replication/meta.fetal.NJ.gwas.tbl.gz',
                'results/LD/delivery/UGT1A.ld'
	output:
		'results/plots/replication-UGT1A-locus-zoom.png'
	conda:
		'../envs/plots.yml'
	script:
		'../scripts/figures/replication-UGT-locuszoom.R'
