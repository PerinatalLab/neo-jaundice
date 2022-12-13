
rule extract_ABO_geno:
	'''
	Rule to extract two missense variants that (largely) determine ABO group.
	# rs8176746 --> chr9:136131322
	#rs8176719 --> chr9:136132908
	# OBS! rs8176719 is a deletion, do therefore not exist in our data. Using chr9:136139265 (rs657152) instead - high LD and close to rs8176719. "A" is instead of the deletion and C instead of G
	'''
        input: 
		'/mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/imputed/all/vcf/9.vcf.gz',
		'results/effect_origin/aux/ids/{sample}_toextract.txt',
	output:
		temp('results/ABO/GT/temp/ABO_{sample}.txt')
	shell:
		"bcftools query -S {input[1]} -r 9:136131322,9:136139265 -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' {input[0]} > {output[0]}"

rule add_header_GT_blood_group:
        'Add header to genotype files.'
        input:
                'results/effect_origin/aux/ids/{sample}_toextract.txt',
                'results/ABO/GT/temp/ABO_{sample}.txt'
        output:
                'results/ABO/GT/{sample}_GT'
        run:
                cols= ['chr','pos','ref','eff'] + [line.strip() for line in open(input[0], 'r')]
                d= pd.read_csv(input[1], header= None, names= cols, sep= '\t')
                d.drop_duplicates(['chr', 'pos'], keep=False, inplace= True)
                d.to_csv(output[0], sep= '\t', header= True, index= False)

rule estimate_blood_group:
	'Estimate blood group from two missense variants.'
        input:
		'results/ABO/GT/{sample}_GT'
	output:
		'results/ABO/blood_group/ABO_{sample}.txt'
	conda:
		'../envs/plots.yml'
	script:
		"../scripts/ABO.R"

rule check_ABO_blood_group:
        ''
        input:
                expand('results/ABO/blood_group/ABO_{sample}.txt', sample= fam_ids['fam_id'])
        output:
                'results/ABO/blood_group/checks/ABO_group_determined.txt'
        shell:
                'touch {output[0]}'

rule merge_haplotype_pheno_ABO:
        'Merge each haplotype and the pheno file.'
        input:
                'results/pheno/fets_pheno_bin.txt',
                'results/pheno/fets_covars.txt',
                'results/effect_origin/aux/ids/parent_offspring_trios.txt',
                'results/ABO/blood_group/ABO_moms.txt',
                'results/ABO/blood_group/ABO_fets.txt',
                'results/ABO/blood_group/ABO_dads.txt',
		'results/effect_origin/delivery/jaundice.txt'
        output:
                temp('results/ABO/pheno/temp/all_subjects.txt')
        run:
                d= pd.read_csv(input[0], sep= '\t', header= 0)
                covar= pd.read_csv(input[1], sep= '\t', header= 0)
                d= pd.merge(d, covar, on= 'IID')
                ids= pd.read_csv(input[2], sep= '\t', header= 0)
                d= pd.merge(d, ids, left_on= 'IID', right_on= 'Child')
		moms= pd.read_csv(input[3], sep= '\t', header= 0)
		fets= pd.read_csv(input[4], sep= '\t', header= 0)
		dads= pd.read_csv(input[5], sep= '\t', header= 0)
		moms.columns= ['Mother', 'ABO_mom']
		fets.columns= ['Child', 'ABO_fet']
		dads.columns= ['Father', 'ABO_dad']
		d= pd.merge(d, moms, on= 'Mother')
		d= pd.merge(d, dads, on= 'Father')
		d= pd.merge(d, fets, on= 'Child')
                d['PREG_ID_1724']= d.PREG_ID_1724.apply(str)
		haplotypes= pd.read_csv(input[6], sep= '\t', header= 0, usecols= ['PREG_ID_1724', 'chr9_136137065_A_G_h2', 'chr9_136137065_A_G_h3'])
		haplotypes['PREG_ID_1724']= haplotypes.PREG_ID_1724.apply(str)
		d= pd.merge(d, haplotypes, on= 'PREG_ID_1724')
		d['ABO_incompatibility']= np.where((d.ABO_mom== 'O') & ((d.ABO_fet== 'A') | (d.ABO_fet== 'B') | (d.ABO_fet== 'AB')), 1, 0)
                d.to_csv(output[0], sep= '\t', header= True, index= False)

rule remove_related_effect_origin_ABO:
        'Remove related individuals'
        input:
                'results/ABO/pheno/temp/all_subjects.txt',
                '/mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/aux/pedigree/mobagen-ethnic-core-samples.kin0'
        output:
                'results/ABO/delivery/ABO-blood-groups.txt'
        run:
                d= pd.read_csv(input[0], sep= '\t', header= 0)
                remove= selectUnrelated(input[1], d, d.Child)
                d= d.loc[~d.Child.isin(remove), :]
                remove= selectUnrelated(input[1], d, d.Mother)
                d= d.loc[~d.Mother.isin(remove), :]
                remove= selectUnrelated(input[1], d, d.Father)
                d= d.loc[~d.Father.isin(remove), :]
                d.to_csv(output[0], sep= '\t', header= True, index= False)

