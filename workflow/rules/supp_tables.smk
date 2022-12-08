rule main_results:
	'Supplementary table with main results.'
	input:
		'results/topregions/delivery/loci-jaundice-fets.txt',
		'results/topregions/delivery/loci-jaundice-moms.txt',
		'results/topregions/delivery/loci-jaundice-dads.txt',
		'results/GWAS/delivery/MoBa-GWAS-jaundice-fets.txt.gz'
	output:
		'results/supp_tables/main_results.txt'
	run:
		fets= pd.read_csv(input[0], sep= '\t', header= 0)
		moms= pd.read_csv(input[1], sep= '\t', header= 0)
		dads= pd.read_csv(input[2], sep= '\t', header= 0)
		fets['OR']= np.exp(fets.BETA)
		fets['origin']= 'fets'
		moms['OR']= np.exp(moms.BETA)
		moms['origin']= 'moms'
		dads['OR']= np.exp(dads.BETA)
		dads['origin']= 'dads'
		x= pd.read_csv(input[3], sep= '\t', header= 0)
		x= x.loc[x.rsid== 'rs6755571', :]
		x['origin']= 'fets'
		x.drop_duplicates(subset= 'ID', inplace= True, keep= 'first')
		d= pd.concat([fets, moms, dads, x])
		d.to_csv(output[0], sep= '\t', header= True, index= False)
