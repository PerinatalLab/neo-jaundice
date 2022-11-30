
rule format_GWAS:
	'Format GWAS summary statistics.'
	input:
		'results/GWAS/sumstats/GWAS-{pheno}/allchr-{sample}.txt.gz'
	output:
		'results/GWAS/delivery/tmp/GWAS-{pheno}-{sample}.allchr.txt',
		'results/GWAS/delivery/tmp/bedfile-{pheno}-{sample}.allchr.txt'
	run:
		df_list= list()
		for d in pd.read_csv(input[0], sep= ' ', header= 0, chunksize= 200000):
			d= d[['CHROM', 'GENPOS', 'ALLELE0', 'ALLELE1', 'A1FREQ', 'INFO', 'N' ,'BETA', 'SE', 'LOG10P']]
			d= d.loc[d.INFO> 0.7, :]
			d= d.loc[(d.A1FREQ> 0.001) & (d.A1FREQ< 0.991), :]
			d.columns= ['CHR', 'POS', 'REF', 'EFF', 'EAF', 'INFO', 'N', 'BETA', 'SE', 'LOG10P']
			d['CHR']= np.where(d.CHR== 'X', '23', d.CHR)
			d['BETA']= np.where(d.REF > d.EFF, -1 * d.BETA, d.BETA)
			d['EAF']= np.where(d.REF > d.EFF, 1 - d.EAF, d.EAF)
			d['REF'], d['EFF']= np.where(d['REF'] > d['EFF'], (d['EFF'], d['REF']), (d['REF'], d['EFF']))
			d['ID']= d.CHR.apply(str) + ':' + d.POS.apply(str) + ':' + d.REF + ':' + d.EFF
			df_list.append(d)
		d= pd.concat(df_list)
		d.to_csv(output[0], sep= '\t', header= True, index= False)
		d['start']= d.POS - 1
		d.to_csv(output[1], sep= '\t', header= False, index= False, columns= ['CHR', 'start', 'POS', 'ID'])
		
rule cut_HRC:
	'Keep only columns wanted.'
	input:
		'/mnt/work/pol/resources/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz'
	output:
		temp('resources/HRC/tmp/hg19_rsids.txt')
	shell:
		'gzip -cd {input[0]} | cut -f1-5 > {output[0]}'

rule format_HRC:
	'Add ID column to HRC file.'
	input:
		'resources/HRC/tmp/hg19_rsids.txt'
	output:
		'resources/HRC/hg19_rsids.txt'
	run:
		with open(input[0], 'rt') as f:
			dialect = csv.Sniffer().sniff(f.readline(), delimiters= ' \t')
			f.seek(0)
			input_file= csv.DictReader(f, delimiter= dialect.delimiter)
			df_list= list()
			with open(output[0], 'w') as csvfile:
				writer = csv.writer(csvfile, delimiter= '\t')
				writer.writerow([g for g in ['ID', 'rsid']])
			for row in input_file:
				if row['ID']== '.':
					continue
				if row['REF'] > row['ALT']:
					ID= row['#CHROM'] + ':' + row['POS'] + ':' + row['ALT'] + ':' + row['REF']
				else:
					ID= row['#CHROM'] + ':' + row['POS'] + ':' + row['REF'] + ':' + row['ALT']
				df_list.append([ID, row['ID']])
				if len(df_list)== 10**6:
					with open(output[0], 'a', newline= '') as file_handler:
						writer1= csv.writer(file_handler, delimiter= '\t')
						for item in df_list:
							writer1.writerow(item)
					df_list= list()
		with open(output[0], 'a', newline= '') as file_handler:
			writer1= csv.writer(file_handler, delimiter= '\t')
			for item in df_list:
				writer1.writerow(item)

rule format_gene_map:
	'Format gene coordinates for bedtools intersect.'
	input:
		'/mnt/work/pol/refdata/UCSC_gene_coordinates_hg19_20210131.txt',
	output:
		'resources/UCSC_coding_gene_coordinates_tx_hg19.txt',
		'resources/bedfiles/UCSC_coding_gene_coordinates_tx_hg19.txt'
	run:
		d= pd.read_csv(input[0], sep= '\t', header= 0)
		d= d.loc[d.cdsStart!= d.cdsEnd, :]
		d['chrom']= d.chrom.str.replace(' ', '')
		d['chrom']= d.chrom.str.replace('chr', '')
		d['chrom']= d.chrom.str.replace('X', '23')
		d['chrom']= pd.to_numeric(d.chrom, errors= 'coerce')
		d.dropna(subset= ['chrom'], inplace= True)
		d.loc[d.txStart > d.txEnd, ['txStart', 'txEnd']]= d.loc[d.txStart > d.txEnd, ['txEnd', 'txStart']].values
		x= d.groupby(['chrom', 'geneSymbol', 'ENSEMBLE_ID'])['txStart'].min().reset_index()
		x1= d.groupby(['chrom', 'geneSymbol', 'ENSEMBLE_ID'])['txEnd'].max().reset_index()
		df= pd.merge(x, x1, on= ['chrom', 'geneSymbol', 'ENSEMBLE_ID'])
		df.columns= ['CHR', 'geneSymbol', 'ENSEMBLE_ID', 'start', 'end']
		df= df.loc[df.start != df.end, :]
		df['start']= df.start - 1
		df.sort_values(by= ['CHR', 'start'], inplace= True)
		df= df[['CHR', 'start', 'end', 'geneSymbol', 'ENSEMBLE_ID']]
		df.to_csv(output[0], sep= '\t', header= True, index= False)
		df[['CHR', 'start', 'end']]= df[['CHR', 'start', 'end']].apply(np.int64)
		df.to_csv(output[1], sep= '\t', header= False, index= False)

rule bedtools_nearest_gene:
        'Use bedtools to add nearest protein coding gene.'
        input:
                'results/GWAS/delivery/tmp/bedfile-{pheno}-{sample}.allchr.txt',
                'resources/bedfiles/UCSC_coding_gene_coordinates_tx_hg19.txt'
        output:
                temp('resources/tmp/nearest_gene/{pheno}-{sample}.txt')
        shell:
                'bedtools closest -t all -a {input[0]} -b {input[1]} > {output[0]}'

rule add_rsid_nearestGene:
	'Add rsid and nearest protein coding gene to summary statistics.'
	input:
		'results/GWAS/delivery/tmp/GWAS-{pheno}-{sample}.allchr.txt',
		'resources/HRC/hg19_rsids.txt',
		'resources/tmp/nearest_gene/{pheno}-{sample}.txt'
	output:
		temp('results/GWAS/delivery/tmp2/MoBa-GWAS-{pheno}-{sample}.txt')
	run:
		rsid= pd.read_csv(input[1], sep= '\t', header=0)
		nearest_gene= pd.read_csv(input[2], sep= '\t', header= None, names= ['CHR', 'X', 'POS', 'ID', 'c1', 'p1', 'p2', 'nearestGene', 'Ensembl_gene'], usecols= ['ID', 'nearestGene'])
		for d in pd.read_csv(input[0], sep= '\t', header= 0, chunksize= 200000):
			d= pd.merge(d, rsid, on= 'ID', how= 'left')
			d= pd.merge(d, nearest_gene, on= 'ID', how= 'left')
			d.to_csv(output[0], sep= '\t', header= not os.path.isfile(output[0]), index= False, mode= 'a')


rule gzip_GWAS_sumstats:
	'Gzip sumstats.'
	input:
		'results/GWAS/delivery/tmp2/MoBa-GWAS-{pheno}-{sample}.txt'
	output:
		'results/GWAS/delivery/MoBa-GWAS-{pheno}-{sample}.txt.gz'
	shell:
		'gzip -c {input[0]} > {output[0]}'


rule independent_GWAS_regions:
        'Obtain a file with independent regions for top loci with a radius of 1.5 Mb.'
        input:
                'results/GWAS/delivery/MoBa-GWAS-{pheno}-{sample}.txt.gz'
        output:
                'results/topregions/delivery/loci-{pheno}-{sample}.txt'
        run:
                d= pd.read_csv(input[0], sep= '\t', compression= 'gzip')
                df= d.loc[d.LOG10P> -np.log10(5*10**-8), :]
                df.sort_values(by= 'LOG10P', ascending= False, inplace= True)
                df.drop_duplicates(subset= ['CHR', 'POS'], keep= 'first', inplace= True)
                df_list= list()
                for chrom in set(df.CHR):
                        d_temp= df.loc[df.CHR== chrom, :]
                        positions= d_temp.POS.values
                        for pos in positions:
                                if pos in d_temp.POS.values:
                                        df_list.append(d_temp.loc[d_temp.POS== pos, :])
                                        d_temp= d_temp.loc[(d_temp.POS < pos - (1.5*10**6)) | (d_temp.POS> pos + (1.5 * 10**6)), :]
                                else:
                                        continue
                x= pd.concat(df_list)
                x['CHR']= x.CHR.astype(str)
                x['CHR']= np.where(x.CHR== '23', 'X', x.CHR)
                x.to_csv(output[0], sep='\t', header= True, index= False)

rule check_files_QC_GWAS:
	'Check that the QC files are created.'
	input:
		expand('results/topregions/delivery/loci-{pheno}-{sample}.txt', pheno= pheno_file['phenotypes'], sample= fam_ids['fam_id'])
	output:
		'results/topregions/delivery/checks/QC_performed.txt'
	shell:
		'touch {output[0]}'

