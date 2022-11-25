rule tabix_eQTL_Catalog:
	'Download summary stats from eQTL Catalog using tabix, for all variants in chromosome 2.'
	priority: 1
	input:
		'resources/eQTL_Catalogue_manifest.txt'
	output:
		temp('results/eQTL_catalogue/data/hg38/temp/{eqtl_catalogue}.txt')
	threads: workflow.cores * 0.25
	run:
		d= pd.read_csv(input[0], sep= '\t', header= 0)
		d= d.loc[(d.quant_method== 'ge') | (d.quant_method== 'microarray'), :]
		d['study']= d.study + '_' + d.quant_method + '_' + d.qtl_group
		d= d.loc[d.study == wildcards.eqtl_catalogue, :]
		url= 'http://' + d.ftp_path.values[0].replace('ftp://', '')
		time.sleep(2)
		print('Downloading from the following url: ' + url)
		shell("tabix {url} 2: > {output[0]}")
		tbi= url.split('/')[-1] + '.tbi'
		shell("rm {tbi}")

rule format_eQTL_catalog:
        'Format eQTL Catalogue data.'
	priority: 10
	input:
                'results/eQTL_catalogue/data/hg38/temp/{eqtl_catalogue}.txt'
        output:
                'results/eQTL_catalogue/data/hg38/formatted/{eqtl_catalogue}.txt.gz'
        run:
                cols= ['molecular_trait_id', 'chromosome', 'position', 'ref', 'alt', 'variant', 'ma_samples', 'maf', 'pvalue', 'beta', 'se', 'type', 'ac', 'an', 'r2', 'molecular_trait_object_id', 'gene_id', 'median_tpm', 'rsid']
                d= pd.read_csv(input[0], sep= '\t', header= None, names= cols)
		genes= ['ENSG00000241635', 'ENSG00000228445', 'ENSG00000288702', 'ENSG00000244474', 'ENSG00000288705', 'ENSG00000167165', 'ENSG00000244122', 'ENSG00000242366', 'ENSG00000241119', 'ENSG00000242515', 'ENSG00000227846', 'ENSG00000228949', 'ENSG00000234143']
		d= d.loc[d.gene_id.isin(genes), :]
                d['ref']= np.where(d.ref.str.len() > d.alt.str.len(), 'I', d.ref)
                d['alt']= np.where(d.ref.str.len() < d.alt.str.len(), 'I', d.alt)
                d['ref']= np.where(d.alt== 'I', 'D', d.ref)
                d['alt']= np.where(d.ref== 'I', 'D', d.alt)
		d['beta']= np.where(d.ref> d.alt, -1 * d.beta, d.beta)
                d.drop_duplicates(subset= ['rsid', 'gene_id'], inplace= True, keep= 'first')
                d['N']= d.an / 2
                d= d[['gene_id', 'rsid', 'position', 'maf', 'N', 'ref', 'alt', 'pvalue', 'beta', 'se']]
                d.to_csv(output[0], sep= '\t', header= True, index= False, compression= 'gzip')


rule format_jaundice_liftover:
	'Format neonatal jaundice summary statistics for speed.'
	input:
		'results/GWAS/delivery/MoBa-GWAS-jaundice-fets.txt.gz'
	output:
		temp('results/eQTL_catalogue/jaundice/temp/hg19/fets-GWAS.txt'),
		temp('results/eQTL_catalogue/jaundice/temp/bed/fets-GWAS.bed')
	run:
		d= pd.read_csv(input[0], sep= '\t', header= 0, usecols= ['rsid', 'CHR', 'POS', 'EAF', 'BETA', 'SE', 'N', 'REF', 'EFF', 'ID'])
		d= d.loc[((d.CHR== 2)), :]
		d.drop_duplicates('ID', keep= 'first', inplace= True)
                d['MAF']= np.where(d.EAF> 0.5, 1 - d.EAF, d.EAF)
		d= d[['rsid', 'CHR', 'POS', 'MAF', 'BETA', 'SE', 'N', 'REF', 'EFF', 'ID']]
		d.columns= ['rsid', 'CHR', 'POS', 'MAF', 'BETA', 'SE', 'TOTALSAMPLESIZE', 'REF', 'EFF', 'ID']
                d.to_csv(output[0], sep= '\t', header= True, index= False)
		d['CHR']= 'chr' + d['CHR'].apply(str)
		d['start']= d.POS.apply(int) - 1
		df= d[['CHR', 'start', 'POS', 'ID']]
		df.to_csv(output[1], sep= '\t', header= False, index= False)

rule liftOver_jaundice:
	'LiftOver from jaundice data from hg19 to hg38 using UCSC liftOver.'
	input:
		'results/eQTL_catalogue/jaundice/temp/bed/fets-GWAS.bed',
		'/home/pol.sole.navais/soft/hg19ToHg38.over.chain.gz'
	output:
		temp('results/eQTL_catalogue/jaundice/temp/fets-GWAS-hg38-keys.txt'),
                temp('results/eQTL_catalogue/jaundice/temp/fets-GWAS-hg38-keys-non-lifted.txt')
	shell:
		'/home/pol.sole.navais/soft/liftOver {input[0]} {input[1]} {output[0]} {output[1]}'

rule jaundice_hg38:
        ''
        input:
                'results/eQTL_catalogue/jaundice/temp/hg19/fets-GWAS.txt',
                'results/eQTL_catalogue/jaundice/temp/fets-GWAS-hg38-keys.txt'
        output:
                'results/eQTL_catalogue/jaundice/temp/hg38/fets-GWAS.txt'
        run:
                x= pd.read_csv(input[1], sep= '\t', header= None, names= ['CHR', 'start', 'POS', 'ID'])
                d= pd.read_csv(input[0], header= 0, sep= '\t', usecols = ['rsid', 'MAF', 'BETA', 'SE', 'TOTALSAMPLESIZE', 'REF', 'EFF', 'ID'])
                d= pd.merge(d, x[['CHR', 'POS', 'ID']], on= 'ID')
                d['CHR']= np.where(d.CHR== 'X', '23', d.CHR)
                d.to_csv(output[0], sep= '\t', header= True, index= False)

rule coloc_eQTL_catalog:
        'Colocalization jaundice with eQTL data form eQTL Catalogue.'
        input:
                'results/eQTL_catalogue/jaundice/temp/hg38/fets-GWAS.txt',
		'results/eQTL_catalogue/data/hg38/formatted/{eqtl_catalogue}.txt.gz'
        output:
                temp('results/eQTL_catalogue/delivery/temp/pph-{eqtl_catalogue}-jaundice-fets.txt'),
                temp('results/eQTL_catalogue/delivery/temp/SNP-{eqtl_catalogue}-jaundice-fets.txt'),
        conda:
                '../envs/coloc.yml'
	threads: 4
        script:
                '../scripts/coloc_eqtl_catalogue.R'

rule concat_eQTL_catalog_coloc_pph:
	'concat results from Colocalization analysis with eQTL catalogue.'
	input:
		expand('results/eQTL_catalogue/delivery/temp/pph-{eqtl_catalogue}-jaundice-fets.txt', eqtl_catalogue= eQTL_catalogue_conditions)
	output:
		'results/eQTL_catalogue/delivery/pph-jaundice-fets.txt'
	shell:
		'''
		head -1 {input[0]} > {output[0]}
		tail -n +2 -q {input} >> {output[0]}
		'''

rule concat_eQTL_catalog_coloc_results:
        'concat results from Colocalization analysis with eQTL catalogue.'
        input:
                expand('results/eQTL_catalogue/delivery/temp/SNP-{eqtl_catalogue}-jaundice-fets.txt', eqtl_catalogue= eQTL_catalogue_conditions)
        output:
                'results/eQTL_catalogue/delivery/SNP-jaundice-fets.txt'
        shell:
                '''
                head -1 {input[0]} > {output[0]}
                tail -n +2 -q {input} >> {output[0]}
                '''

