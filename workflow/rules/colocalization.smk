
rule colocalization:
        'Run coloc using reproductive traits.'
        input:
                'results/GWAS/delivery/MoBa-GWAS-jaundice-fets.txt.gz',
                'resources/bilirubin/delivery/total-bilirubin-GWAS.txt.gz',
                'resources/ld_indep_regions.txt'
        output:
                'results/colocalization/delivery/pph_jaundice_fets.txt',
                'results/colocalization/delivery/results_jaundice_fets.txt'
        conda:
                '../envs/coloc.yml'
        script:
                '../scripts/coloc.R'


