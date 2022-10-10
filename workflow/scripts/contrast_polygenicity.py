import pandas as pd
import numpy as np

def contrast(input):
    all_local_hsq = pd.read_table(input, sep='\t')
    all_local_hsq.loc[all_local_hsq['local_h2g']<0, 'local_h2g'] = 0.0
    # sort local SNP-heritability by descending order and get total hsq and
    # total number of snps
    # sort local SNP-heritability by descending order
    idx = np.argsort(-all_local_hsq['local_h2g'])
    all_local_hsq_sorted= (all_local_hsq['local_h2g'][idx].reset_index(drop=True))
    all_nsnp_sorted= (all_local_hsq['num_snp'][idx].reset_index(drop=True))
    # get total SNP-heritability and total number of SNPs
    all_total_hsq= (np.sum(all_local_hsq['local_h2g']))
    all_total_nsnp= (np.float(np.sum(all_local_hsq['num_snp'])))
    # get the proportions
    all_xval = np.array([])
    all_yval = np.array([])
    all_yerr = np.array([])
    nloci = all_local_hsq_sorted.shape[0]
    # iterate through the loci
    for j in range(nloci):
        hsq_sum = np.sum(all_local_hsq_sorted[0:j+1])
        hsq_frac = hsq_sum / all_total_hsq
        snp_frac = np.sum(all_nsnp_sorted[0:j+1]) / all_total_nsnp
        all_xval = np.append(all_xval, snp_frac)
        all_yval = np.append(all_yval, hsq_frac)
        # get jack knife standard error
        jk_est = []
        for k in range(nloci):
            total_hsq_jk = all_total_hsq-all_local_hsq_sorted[k]
            hsq_sum_jk = hsq_sum
            if k <= j:
                hsq_sum_jk -= all_local_hsq_sorted[k]
            jk_est.append(hsq_sum_jk/total_hsq_jk)
        jk_est = np.array(jk_est)
        se = np.sqrt((nloci-1)*np.mean(np.square(jk_est-hsq_frac)))
        all_yerr = np.append(all_yerr, se)
    d= pd.DataFrame({'local_h2': all_yval, 'nsnp_sorted': all_xval, 'se': all_yerr})
    return(d)


d= contrast(snakemake.input[0])

d.to_csv(snakemake.output[0], header= True, index= False, sep= '\t')
