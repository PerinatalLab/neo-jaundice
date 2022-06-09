rule merge_PGS:
	''
	input:
		'results/PGS/delivery/moms-jaundice-PGS.txt',
		'results/PGS/delivery/dads-jaundice-PGS.txt',
		'results/PGS/delivery/fets-jaundice-PGS.txt',
		'results/effect_origin/aux/ids/parent_offspring_trios.txt'
	output:
		'results/PGS/delivery/all-jaundice-PGS.txt'
	run:
		moms= pd.read_csv(input[0], sep= '\t', header= 0)
		dads= pd.read_csv(input[1], sep= '\t', header= 0)
		fets= pd.read_csv(input[2], sep= '\t', header= 0)
		trios= pd.read_csv(input[3], sep= '\t', header= 0)
		d= pd.merge(trios, moms, left_on= 'Mother', right_on= 'IID')
		d= pd.merge(d, dads, left_on= 'Father', right_on= 'IID')
		d= pd.merge(d, fets, left_on= 'Child', right_on= 'IID')
		d.to_csv(output[0], sep= '\t', header= True, index= False)

rule merge_PGS_nochr2:
        ''
        input:
                'results/PGS/delivery/moms-jaundice-PGS-nochr2.txt',
                'results/PGS/delivery/dads-jaundice-PGS-nochr2.txt',
                'results/PGS/delivery/fets-jaundice-PGS-nochr2.txt',
                'results/effect_origin/aux/ids/parent_offspring_trios.txt'
        output:
                'results/PGS/delivery/all-jaundice-PGS-nochr2.txt'
        run:
                moms= pd.read_csv(input[0], sep= '\t', header= 0)
		moms.columns= ['IID', 'moms_jaundice_nochr2']
                dads= pd.read_csv(input[1], sep= '\t', header= 0)
		dads.columns= ['IID', 'dads_jaundice_nochr2']
                fets= pd.read_csv(input[2], sep= '\t', header= 0)
		fets.columns= ['IID', 'fets_jaundice_nochr2']
                trios= pd.read_csv(input[3], sep= '\t', header= 0)
                d= pd.merge(trios, moms, left_on= 'Mother', right_on= 'IID')
                d= pd.merge(d, dads, left_on= 'Father', right_on= 'IID')
                d= pd.merge(d, fets, left_on= 'Child', right_on= 'IID')
                d.to_csv(output[0], sep= '\t', header= True, index= False)


rule merge_all_data_phased_transmission:
	'Merge all data, including MFR.'
	input:
		'results/effect_origin/delivery/jaundice.txt',
		'results/ABO/delivery/ABO-blood-groups.txt',
		'/mnt/archive2/p1724/v12/PDB1724_MFR_541_v12.csv',
		'results/PGS/delivery/h1-jaundice-PGS.txt',
		'results/PGS/delivery/h2-jaundice-PGS.txt',
		'results/PGS/delivery/h3-jaundice-PGS.txt',
		'results/PGS/delivery/h4-jaundice-PGS.txt',
	output:
		'results/merge_data/delivery/jaundice-transmission.txt'
	run:
		snps= pd.read_csv(input[0], sep= '\t', header= 0)
		abo= pd.read_csv(input[1], sep= '\t', header= 0)
		abo= abo[['ABO_incompatibility', 'PREG_ID']]
		mfr= pd.read_csv(input[2], sep= ';', header= 0)
		mfr.drop('KJONN', axis= 1, inplace= True)
		df_list= list()
		for i in range(3, len(input)):
			h= pd.read_csv(input[i], sep= '\t', header= 0)
			df_list.append(h)
		d= reduce(lambda df1, df2: pd.merge(df1, df2, on= 'PREG_ID'), df_list)
		d= pd.merge(d, snps, on= 'PREG_ID')
		d= pd.merge(d, abo, on= 'PREG_ID')
		d= pd.merge(d, mfr, left_on= 'PREG_ID', right_on= 'PREG_ID_1724')
		d.to_csv(output[0], sep= '\t', header= True, index= False)

rule merge_all_data:
        'Merge all data, including MFR.'
        input:
                'results/ABO/delivery/ABO-blood-groups.txt',
                '/mnt/archive2/p1724/v12/PDB1724_MFR_541_v12.csv',
                'results/effect_origin/delivery/conditional/dosage-jaundice.txt',
                'results/PGS/delivery/all-jaundice-PGS.txt',
                'results/PGS/delivery/all-jaundice-PGS-nochr2.txt',
		'results/effect_origin/delivery/jaundice.txt'
        output:
                'results/merge_data/delivery/jaundice.txt'
        run:
                abo= pd.read_csv(input[0], sep= '\t', header= 0)
                abo= abo[['ABO_incompatibility', 'PREG_ID']]
                mfr= pd.read_csv(input[1], sep= ';', header= 0)
                mfr.drop('KJONN', axis= 1, inplace= True)
                ds= pd.read_csv(input[2], header= 0, sep= '\t')
                pgs= pd.read_csv(input[3], header= 0, sep= '\t')
                pgsno2= pd.read_csv(input[4], header= 0, sep= '\t')
		pheno= pd.read_csv(input[5], sep= '\t', header= 0)
		pheno= pheno[['PREG_ID', 'jaundice', 'KJONN', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6',  'PC7', 'PC8', 'PC9', 'PC10', 'cohort']]
                d= pd.merge(ds, abo, on= 'PREG_ID')
                d= pd.merge(d, mfr, left_on= 'PREG_ID', right_on= 'PREG_ID_1724')
                d= pd.merge(d, pgs, on= 'PREG_ID')
                d= pd.merge(d, pgsno2, on= 'PREG_ID')
		d= pd.merge(d, pheno, on= 'PREG_ID')
		d.columns= d.columns.str.replace(':', '_')
                d.to_csv(output[0], sep= '\t', header= True, index= False)

