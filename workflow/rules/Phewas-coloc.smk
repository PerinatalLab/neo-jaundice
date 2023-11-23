import itertools
from urllib import request
from concurrent.futures import ThreadPoolExecutor
import glob
import os
import socket

socket.setdefaulttimeout(3) # 60 seconds

#FINNGEN_pheno_nms= pd.read_csv('/mnt/hdd/common/pol/references/FINNGEN/manifest_R5.txt', header= 0, sep= '\t').phenocode.to_list()



#### Pan UK Biobank

rule format_manifest_Pan_UKBB:
	'Format manifest files.'
	input:
		'resources/PanUKBB-h2-manifest.txt',
		'resources/PanUKBB-pheno-manifest.txt'
	output:
		'resources/Pan-UK-BB-manifest.txt'
	run:
		h2= pd.read_csv(input[0], sep= '\t', header= 0, low_memory= False)
		pheno= pd.read_csv(input[1], sep= '\t', header= 0, low_memory=False)
		h2= h2.loc[h2['pop']== 'EUR', :]
		h2['coding']= np.where(h2.coding.isnull(), h2.modifier, h2.coding)
		pheno['coding']= np.where(pheno.coding.isnull(), pheno.modifier, pheno.coding)
		h2= pd.merge(h2, pheno, on= ['trait_type', 'phenocode', 'pheno_sex', 'coding'])
		h2= h2.loc[h2['estimates.ldsc.h2_liability'] > 0.01, :]
		h2= h2.loc[h2['qcflags.GWAS_run'], :]
		h2= h2.loc[h2['qcflags.defined_h2'], :]
		h2= h2.loc[h2['qcflags.significant_z'], :]
		h2= h2.loc[h2['qcflags.in_bounds_h2'], :]
		h2= h2.loc[h2['lambda_gc_EUR']> 0.9, :]
		h2.to_csv(output[0], sep= '\t', header= True, index= False)

rule format_trait_PAN_UKBB:
        'Format neonatal jaundice summary statistics for speed.'
        input:
                'results/GWAS/delivery/MoBa-GWAS-jaundice-fets.txt.gz'
        output:
                temp('results/PheWas/PAN_UKBB/temp/traits/fets-jaundice.txt'),
        run:
                d= pd.read_csv(input[0], sep= '\t', header= 0, usecols= ['ID', 'CHR', 'POS', 'EAF', 'BETA', 'SE', 'N'])[['ID', 'CHR', 'POS', 'EAF', 'BETA', 'SE', 'N']]
                d= d.loc[((d.CHR== 2) & (d.POS>= 233127536) & (d.POS<= 236127536)), :]
                d.drop_duplicates('ID', inplace= True, keep= 'first')
		d['MAF']= np.where(d.EAF> 0.5, 1 - d.EAF, d.EAF)
		d= d[['ID', 'CHR', 'POS', 'MAF', 'BETA', 'SE', 'N']]
                d.columns= ['ID', 'CHR', 'POS', 'MAF', 'BETA', 'SE', 'TOTALSAMPLESIZE']
                d.to_csv(output[0], sep= '\t', header= True, index= False)

rule format_adult_bilirubin_PAN_UKBB:
        'Format bilirubin summary statistics for speed.'
        input:
                'resources/bilirubin/delivery/total-bilirubin-GWAS.txt.gz'
        output:
                temp('results/PheWas/PAN_UKBB/temp/traits/adult-bilirubin.txt'),
        run:
                d= pd.read_csv(input[0], sep= '\t', header= 0, usecols= ['CHROM', 'POS', 'EFF', 'REF', 'BETA', 'SE'])[['CHROM', 'POS', 'EFF', 'REF', 'BETA', 'SE']]
                d= d.loc[((d.CHROM== 2) & (d.POS>= 233127536) & (d.POS<= 236127536)), :]
		d['TOTALSAMPLESIZE']= 363228
		d['BETA']= np.where(d.REF > d.EFF, -1 * d.BETA, d.BETA)
		d['REF'], d['EFF']= np.where(d['REF'] > d['EFF'], (d['EFF'], d['REF']), (d['REF'], d['EFF']))
		d['ID']= d.CHROM.apply(str) + ':' + d.POS.apply(str) + ':' + d.REF + ':' + d.EFF
                d.drop_duplicates('ID', inplace= True, keep= 'first')
                d= d[['ID', 'CHROM', 'POS', 'BETA', 'SE', 'TOTALSAMPLESIZE']]
                d.columns= ['ID', 'CHR', 'POS', 'BETA', 'SE', 'TOTALSAMPLESIZE']
                d.to_csv(output[0], sep= '\t', header= True, index= False)


def download(url):
	pheno_name= url.split('/')[-1].replace('.tsv.bgz', '')
	temp_outfile= '/mnt/work/pol/neo-jaundice/results/PheWas/PAN_UKBB/temp/data/' + pheno_name + '-header.tsv.bgz'
	tbi_url= url.split('files')[0] + 'files_tabix' + url.split('files')[1] + '.tbi'
	tbi_file= pheno_name + '.tsv.bgz.tbi'
	print('Downloading' + pheno_name)
	#socket.setdefaulttimeout(3)
	#request.urlretrieve(url, temp_outfile)
	#request.urlretrieve(tbi_url, tbi_file)
	outfile_header= '/mnt/work/pol/neo-jaundice/results/PheWas/PAN_UKBB/temp/data/' + pheno_name + '-header.txt'
	outfile= '/mnt/work/pol/neo-jaundice/results/PheWas/PAN_UKBB/temp/data/' + pheno_name + '.txt.gz'
	shell("wget {tbi_url}")
	shell("tabix {url} -h 2:233127536-236127536 | gzip > {outfile}")
	print('This is a test')
	shell('rm {tbi_file}')
	shell("set +e; timeout --foreground 5 wget -O {temp_outfile} {url}; exitcode=$?; if [ $exitcode -eq 1 ]; then exit; else echo 'Done!'; fi")

checkpoint dl_PAN_UKBB:
	'Download summary statistics from PAN_UK Biobank.'
	input:
		'resources/Pan-UK-BB-manifest.txt'
	output:
		directory('results/PheWas/PAN_UKBB/temp/data/')
	params:
		'results/PheWas/PAN_UKBB/temp/data/'
	run:
		if not os.path.exists(params[0]):
			os.makedirs(params[0])
		d= pd.read_csv(input[0], sep= '\t', header= 0)
		urls= d.aws_link.tolist()
		with ThreadPoolExecutor(max_workers=8) as executor:
			executor.map(download, urls)

rule UGT1A4_coloc_PAN_UKBB:
        'Colocalization neonatal jaundice with other phenotypes from PAN UK Biobank.'
        input:
                'results/PheWas/PAN_UKBB/temp/traits/fets-jaundice.txt',
                'resources/Pan-UK-BB-manifest.txt',
		'results/PheWas/PAN_UKBB/temp/data/{PAN_UKBB_trait}.txt.gz',
		'results/PheWas/PAN_UKBB/temp/data/{PAN_UKBB_trait}-header.txt',
		'results/PheWas/PAN_UKBB/temp/traits/adult-bilirubin.txt'
        output:
                temp('results/PheWas/PAN_UKBB/coloc_results/{PAN_UKBB_trait}/pph.txt'),
                temp('results/PheWas/PAN_UKBB/coloc_results/{PAN_UKBB_trait}/results.txt')
        conda:
                '../envs/coloc.yml'
        script:
                '../scripts/PAN_UKBB_coloc.R'


def aggregate_output_coloc_PAN_UKBB(wildcards):
        'Aggregate the files from PAN_UKBB_trait wildcard.'
        checkpoint_output= checkpoints.dl_PAN_UKBB.get(**wildcards).output[0]
	return expand('results/PheWas/PAN_UKBB/coloc_results/{PAN_UKBB_trait}/{coloc_out}.txt', coloc_out= wildcards.coloc_out, PAN_UKBB_trait= glob_wildcards(os.path.join(checkpoint_output, '{PAN_UKBB_trait}.txt.gz')).PAN_UKBB_trait)

rule cat_coloc_PAN_UKBB:
	'Concatenate posterior probabilities.'
	input:
		aggregate_output_coloc_PAN_UKBB
	output:
		'results/PheWas/PAN_UKBB/final/{coloc_out}.txt'
	shell:
		'''
		head -1 {input[0]} > {output[0]}
		tail -n +2 -q {input} >> {output[0]}
		'''


rule clear_UGT1A4_coloc_direct_UKBB:
	'Clean temp directory from checkpoint'
	input:
		expand('results/PheWas/PAN_UKBB/final/{coloc_out}.txt', coloc_out= ['pph', 'results'])
	output:
		'results/PheWas/PAN_UKBB/flag/cleaned_file.txt'
	params:
		tmp_dir= 'results/PheWas/PAN_UKBB/temp/data/'
	shell:
                '''
		rm -r {params.tmp_dir}
		touch {output[0]}
		'''
