import pandas as pd
import numpy as np
import csv
import gzip

##### load config and sample sheets #####

configfile: "config/config.yaml"

fam_ids = pd.read_csv(config["fam_member"], sep= "\t").set_index("fam_id", drop=False)
pheno_file= pd.read_csv(config["pheno"], sep= '\t')
CHROM= [i.strip() for i in open(config["CHROM"], 'r')]

rule phenofile:
        'Pheno file for Jaundice.'
        input:
                '/mnt/work/pol/MOBAGENETICS/rotterdam2_mfr.csv',
                '/mnt/work/pol/MOBAGENETICS/PREG_ID_to_IID.txt',
                '/mnt/archive/MOBAGENETICS/genotypes-base/aux/flaglist-merged/mobagen-flaglist-n99259.txt',
                '/mnt/work/pol/MOBAGENETICS/pca_out.txt'
        output:
                temp('results/aux/pheno/temp/temp_pheno_moms.txt'),
                temp('results/aux/pheno/temp/temp_pheno_fets.txt'),
                temp('results/aux/pheno/temp/temp_pheno_dads.txt')
        script:
                '../scripts/pheno_file.R'


rule concat_phenos_PCA:
        'Concat pheno files, and add PCA.'
        input:
                '/mnt/archive/MOBAGENETICS/genotypes-base/aux/pca/mobagen-total/mobagen-total-proj-pc',
                'results/aux/pheno/temp/temp_pheno_{sample}.txt',
                '/mnt/archive/MOBAGENETICS/genotypes-base/aux/pedigree/mobagen-ethnic-core-samples.kin0'
        output:
                'results/pheno/{sample}_pheno_bin.txt',
                'results/pheno/{sample}_covars.txt',
                'results/aux/ids/samples/{sample}_ids.txt',
        run:
                d= pd.read_csv(input[1], header= 0, sep= '\t')
                pca= pd.read_csv(input[0], header= 0, sep= '\t')
                d= pd.merge(d, pca, how= 'inner', on= 'IID')
                d.fillna('NA', inplace= True)
                for elem in set(d['cohort']):
                        d[str(elem)]= d['cohort'] == elem
                d[['HARVEST', 'NORMENT', 'ROTTERDAM1', 'ROTTERDAM2']] *= 1
                d['FID']= d.IID
                d.to_csv(output[0], sep= '\t', header= True, index= False, columns= ['FID', 'IID', 'jaundice'])
                d.to_csv(output[1], sep= '\t', header= True, index= False, columns= ['FID', 'IID', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10', 'cohort', 'KJONN'])
                d.to_csv(output[2], sep= '\t', header= False, index= False, columns= ['FID', 'IID'])
