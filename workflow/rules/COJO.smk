
rule list_variants_COJO:
	'Obtain a list of variants included in the GWAS.'
	input:
		'results/GWAS/delivery/MoBa-GWAS-{pheno}-{sample}.txt.gz',
		'results/topregions/delivery/loci-{pheno}-{sample}.txt'
	output:
		'results/aux/COJO/variants/{pheno}-{sample}.txt'
	run:
		d= pd.read_csv(input[0], sep='\t', header=0, usecols= ['CHR', 'POS', 'ID'])
		x= pd.read_csv(input[1], sep= '\t', header= 0, usecols= ['CHR', 'POS'])
		x.columns= ['CHR', 'pos']
		d['CHR']= d.CHR.apply(str)
		d['CHR']= np.where(d.CHR== '23', 'X', d.CHR)
		x['CHR']= x.CHR.apply(str)
                x['CHR']= np.where(x.CHR== '23', 'X', x.CHR)
                d= pd.merge(d, x, on= 'CHR')
		d= d.loc[(d.POS> d.pos - 1.5e6) & (d.POS< d.pos + 1.5e6), :]
                d.drop_duplicates('ID', inplace= True, keep ='first')
		d['CHR']= np.where(d.CHR== 'X', '23', d.CHR)
		d.sort_values(['CHR', 'POS'], ascending= True, inplace= True)
                d.to_csv(output[0], sep= '\t', header= False, index= False, columns= ['CHR', 'POS', 'POS', 'ID'])


rule list_non_related_samples:
	'Obtain a list of non-related samples for LD calculation.'
	input:
		'results/aux/ids/samples/{sample}_ids.txt',
		'/mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/aux/pedigree/mobagen-ethnic-core-samples.kin0'
	output:
		'results/aux/COJO/non-related/{sample}-ids-plink.txt'
	run:
		d= pd.read_csv(input[0], sep= '\t', header= None)
		d.columns= ['FID', 'IID']
		remove= selectUnrelated(input[1], d, d.IID)
		d= d.loc[~d.IID.isin(remove), :]
		d.to_csv(output[0], sep= '\t', header= True, index= False)

rule filter_bed:
	'Extract from MoBaGenetics all genetic variants matching women.'
	input:
		'results/aux/COJO/non-related/{sample}-ids-plink.txt',
		'results/aux/COJO/variants/{pheno}-{sample}.txt',
		'/mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/imputed/all/plink/{CHR}.bed'
	output:
		temp('results/COJO/data/plink/temp/{pheno}/{sample}-{CHR}.bed'),
		temp('results/COJO/data/plink/temp/{pheno}/{sample}-{CHR}.bim'),
		temp('results/COJO/data/plink/temp/{pheno}/{sample}-{CHR}.fam')
	params:
		'/mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/imputed/all/plink/{CHR}',
		'results/COJO/data/plink/temp/{pheno}/{sample}-{CHR}'
	threads: 7
	run:
		d= pd.read_csv(input[1], sep= '\t', header= None, names= ['CHR', 'POS', 'POS2', 'ID'])
		d['CHR']= d['CHR'].apply(str)
		wildCHROM= np.where(str(wildcards.CHR)== 'X', '23', wildcards.CHR)
		if str(wildCHROM) not in set(d.CHR.values):
			print('CHR: ' + wildcards.CHR + ' not in ' + wildcards.sample + ' GWAS.')
			shell('touch {output[0]}')
			shell('touch {output[1]}')
			shell('touch {output[2]}')
		else:
			shell('~/soft/plink2 --bfile {params[0]} --keep {input[0]} --extract bed1 {input[1]} --memory 10000 --threads {threads} --make-bed --out {params[1]}')


rule modify_bim:
	'Modify bim so that variant id== CHR:POS:REF:EFF, where REF is alphabetically lower than EFF.'
	input:
		'results/COJO/data/plink/temp/{pheno}/{sample}-{CHR}.bim',
                'results/COJO/data/plink/temp/{pheno}/{sample}-{CHR}.bed',
                'results/COJO/data/plink/temp/{pheno}/{sample}-{CHR}.fam'
	output:
		temp('results/COJO/data/plink/temp/norsid/{sample}/{pheno}-{CHR}_norsid.bim'),
		temp('results/COJO/data/plink/temp/norsid/{sample}/{pheno}-{CHR}_norsid.bed'),
		temp('results/COJO/data/plink/temp/norsid/{sample}/{pheno}-{CHR}_norsid.fam'),
		temp('results/COJO/data/plink/temp/norsid/{sample}/{pheno}-{CHR}_norsid_duplicates.txt')
	run:
		if os.stat(input[0]).st_size == 0:
			shell('mv {input[0]} {output[0]}')
			shell('mv {input[1]} {output[1]}')
	                shell('mv {input[2]} {output[2]}')
			shell('touch {output[3]}')
		else:
			d= pd.read_csv(input[0], sep= '\t', header= None, names= ['CHR', 'SNP', 'x1', 'POS', 'A1', 'A2'])
		        d['A1']= np.where(len(d.A1) > len(d.A2), 'I', d.A1)
			d['A2']= np.where(len(d.A1) < len(d.A2), 'I', d.A2)
	                d['A1']= np.where(d.A2== 'I', 'D', d.A1)
		        d['A2']= np.where(d.A1== 'I', 'D', d.A2)
			d['CHR']= d.CHR.apply(str)
			d['CHR']= np.where(d.CHR== 'X', '23', d.CHR)
		        d['SNP']= np.where(d.A1>d.A2, d.CHR.apply(str) + ':' + d.POS.apply(str) + ':' + d.A2 + ':' + d.A1, d.CHR.apply(str) + ':' + d.POS.apply(str) + ':' + d.A1 + ':' + d.A2)
			d.to_csv(output[0], sep= '\t', header= False, index= False)
	                d= d[d.duplicated(['SNP'], keep= False)]
		        d.drop_duplicates('SNP', inplace= True, keep= 'first')
			d.to_csv(output[3], sep='\t', columns= ['SNP'])
			shell('mv {input[1]} {output[1]}')
			shell('mv {input[2]} {output[2]}')

rule bed_list:
	'File with list of bed files for COJO.'
	input:
		expand('results/COJO/data/plink/temp/norsid/{{sample}}/{{pheno}}-{CHR}_norsid.bim', CHR= CHROM)
	output:
		'results/COJO/data/plink/temp/norsid/lf/list_files_{pheno}_{sample}.txt'
	run:
		bed_list= [infile.replace('.bim', '') for infile in input]
		with open(output[0], 'w') as f:
			for item in bed_list:
				f.write("%s\n" % item)
		


rule format_sumstats_cojo:
	'Format sumstats according to CGTA cojo.'
	input:
		'results/GWAS/delivery/MoBa-GWAS-{pheno}-{sample}.txt.gz'
	output:
		'results/COJO/sumstats/sumstats-{pheno}-{sample}.txt'
	run:
		d= pd.read_csv(input[0], sep= '\t', header= 0, compression= 'gzip', usecols= ['CHR', 'POS', 'EFF', 'REF', 'N', 'EAF', 'BETA', 'SE', 'LOG10P', 'ID'])[['CHR', 'POS', 'ID', 'EFF', 'REF', 'N', 'EAF', 'BETA', 'SE', 'LOG10P']]
		d.columns= ['CHR', 'POS', 'SNP', 'A1', 'A2', 'N', 'freq', 'b', 'se', 'p']
		d['p']= 10**-d['p']
		d['CHR']= d.CHR.apply(str)
		d['CHR']= np.where(d.CHR== 'X', '23', d.CHR)
		d['SNP']= d.CHR.apply(str) + ':' + d.POS.apply(str) + ':' + d.A2 + ':' + d.A1
		d.drop_duplicates('SNP', keep= 'first', inplace= True)
		d.to_csv(output[0], sep= '\t', header= True, index= False, columns= ['SNP', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'N'])

rule COJO_slct:
	'Perform conditional analysis using COJO from GCTA'
	input:
		'results/COJO/sumstats/sumstats-{pheno}-{sample}.txt',
		'results/topregions/delivery/loci-{pheno}-{sample}.txt',
		'results/COJO/data/plink/temp/norsid/lf/list_files_{pheno}_{sample}.txt',
		expand('results/COJO/data/plink/temp/norsid/{{sample}}/{{pheno}}-{{CHR}}_norsid.{ext}', ext= ['bim', 'fam', 'bed']),
	output:
		'results/COJO/aux/{sample}/{pheno}-{CHR}.cma.cojo',
		'results/COJO/aux/{sample}/{pheno}-{CHR}.jma.cojo'
	threads: 7
	params:
		'results/COJO/data/plink/temp/norsid/{sample}/{pheno}-{CHR}_norsid',
		'results/COJO/aux/{sample}/{pheno}-{CHR}',
	run:
		cma_list= list()
		jma_list= list()
		cma= ['Chr', 'SNP', 'bp', 'refA', 'freq', 'b', 'se', 'p', 'n', 'freq_geno', 'bC', 'bC_se', 'pC', 'locus']
		jma= ['Chr', 'SNP', 'bp', 'refA', 'freq', 'b', 'se', 'p', 'n', 'freq_geno', 'bJ', 'bJ_se', 'pJ', 'LD_r', 'locus']
		d= pd.read_csv(input[1], sep= '\t', header= 0)
		d['CHR']= d['CHR'].apply(str)
		d['CHR']= np.where(d.CHR== 'X', '23', d.CHR)
		wildCHROM= np.where(str(wildcards.CHR)== 'X', '23', wildcards.CHR)
                if str(wildCHROM) not in d.CHR.values:
			cma_df= pd.DataFrame(columns= cma)
			cma_df.to_csv(output[0], sep= '\t', header= True, index= False)
			jma_df= pd.DataFrame(columns= jma)
                        jma_df.to_csv(output[1], sep= '\t', header= True, index= False)
		else:
			d= d.loc[d.CHR == str(wildCHROM), :]
			for index, row in d.iterrows():
				pos= int(row['POS'])
				region= str(row['CHR']) + ' ' + str(pos) + ' ' + '1500'
				outfile= params[0]
				shell('~/soft/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --bfile {params[0]} --maf 0.001 --extract-region-bp {region} --cojo-file {input[0]} --cojo-slct --thread-num {threads} --out {params[1]}')
				x= pd.read_csv(output[0], sep= '\t', header= 0)
				print(x.columns)
				x['locus']= row['nearestGene']
				print(x.columns)
				cma_list.append(x)
                                x= pd.read_csv(output[1], sep= '\t', header= 0)
                                x['locus']= row['nearestGene']
				jma_list.append(x)
			cma= pd.concat(cma_list)
			cma.to_csv(output[0], sep= '\t', header= True, index= False)
			print(cma)
			print(cma.columns)
			jma= pd.concat(jma_list)
			jma.to_csv(output[1], sep= '\t', header= True, index= False)

rule concat_COJO:
	'Concat results from COJO analysis.'
	input:
		expand('results/COJO/aux/{{sample}}/{{pheno}}-{CHR}.{{cojo}}.cojo', CHR= CHROM)
	output:
		'results/COJO/delivery/{sample}-{pheno}.{cojo}'
	shell:
		'''
		head -1 {input[0]} > {output[0]}
                tail -n +2 -q {input} >> {output[0]}
		'''

rule check_COJO:
	'File to check COJO.'
	input:
		expand('results/COJO/delivery/{sample}-{pheno}.{cojo}', sample= fam_ids['fam_id'], pheno= pheno_file['phenotypes'], cojo= ['jma', 'cma']),
	output:
		'results/COJO/checks/COJO_performed'
	shell:
		'touch {output[0]}'
