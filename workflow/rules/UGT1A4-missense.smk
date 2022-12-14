
rule get_dosage_missense:
	'Extract dosage for missense variant at UGT1A4.'
	input:
		'results/effect_origin/aux/ids/{sample}_toextract.txt',
		'/mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/imputed/all/vcf/2.vcf.gz'
	output:
		'results/UGT-missense/aux/DS/temp/{sample}_ds2'
	run:
		shell("bcftools query -S {input[0]} -r 2:234627536 -f '%CHROM\t%POS\t%REF\t%ALT[\t%DS]\n' {input[1]} -o {output[0]}")

rule add_header_missense_dosage:
	''
	input:
		'results/effect_origin/aux/ids/{sample}_toextract.txt',
		'results/UGT-missense/aux/DS/temp/{sample}_ds2'
	output:
		'results/UGT-missense/DS/{sample}_DS2'
	run:
		cols= ['chr','pos','ref','eff'] + [line.strip() for line in open(input[0], 'r')]
                d= pd.read_csv(input[1], header= None, names= cols, sep= '\t')
                d.drop_duplicates(['chr', 'pos'], keep=False, inplace= True)
                d.to_csv(output[0], sep= '\t', header= True, index= False)


rule merge_missense_data:
	''
	input:
		'results/merge_data/delivery/jaundice.txt',
		'results/effect_origin/aux/ids/parent_offspring_trios.txt',
		'results/UGT-missense/DS/moms_DS2',
		'results/UGT-missense/DS/fets_DS2',
		'results/UGT-missense/DS/dads_DS2',
		'results/effect_origin/ids/PREG_ID_jaundice.txt'
	output:
		'results/UGT-missense/delivery/jaundice.txt'
	run:
		alld= pd.read_csv(input[0], sep= '\t', header= 0)
		trios= pd.read_csv(input[1], sep= '\t', header= 0)
		moms= pd.read_csv(input[2], sep= '\t', header= 0)
		moms= moms.iloc[:, 4:].T.reset_index()
                moms.columns= ['Mother', 'moms_UGT_P24T']
		fets= pd.read_csv(input[3], sep= '\t', header= 0)
		fets= fets.iloc[:, 4:].T.reset_index()
		fets.columns= ['Child', 'fets_UGT_P24T']
		dads= pd.read_csv(input[4], sep= '\t', header= 0)
		dads= dads.iloc[:, 4:].T.reset_index()
                dads.columns= ['Father', 'dads_UGT_P24T']
		d= pd.merge(trios, moms, on= 'Mother')
		d= pd.merge(d, dads, on= 'Father')
		d= pd.merge(d, fets, on= 'Child')
		d['PREG_ID_1724']= d.PREG_ID_1724.apply(int).apply(str)
                alld['PREG_ID_1724']= alld.PREG_ID_1724.apply(int).apply(str)
		d= pd.merge(alld, d, on= 'PREG_ID_1724')
		with open(input[5]) as f:
			trio_ids = [line.rstrip('\n').replace('.0', '') for line in f]
		d= d.loc[d.PREG_ID_1724.isin(trio_ids), :]
		d= d.loc[:,~d.columns.str.contains('_y')]
		d= d.loc[:,~d.columns.str.contains('_x')]
		d.to_csv(output[0], sep= '\t', header= True, index= False)



rule get_hard_calls_missense:
	''
	input:
		'results/effect_origin/aux/ids/{sample}_toextract.txt',
                '/mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/imputed/all/vcf/2.vcf.gz'
	output:
		'results/UGT-missense/aux/GT/temp/{sample}_gt2'
	run:
		shell("bcftools query -S {input[0]} -r 2:234627536 -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' {input[1]} -o {output[0]}")

rule add_header_missense_hard_calls:
        ''
        input:
                'results/effect_origin/aux/ids/{sample}_toextract.txt',
                'results/UGT-missense/aux/GT/temp/{sample}_gt2'
        output: 
                'results/UGT-missense/GT/{sample}_GT2'
        run:
                cols= ['chr','pos','ref','eff'] + [line.strip() for line in open(input[0], 'r')]
                d= pd.read_csv(input[1], header= None, names= cols, sep= '\t')
                d.drop_duplicates(['chr', 'pos'], keep=False, inplace= True)
                d.to_csv(output[0], sep= '\t', header= True, index= False)

rule get_allele_missense_hard_calls:
        'Retrieve allele transmission from family trios (after phasing).'
        input:
                'results/UGT-missense/GT/fets_GT2',
                'results/UGT-missense/GT/moms_GT2',
                'results/effect_origin/aux/ids/parent_offspring_trios.txt',
                'results/UGT-missense/GT/dads_GT2'
        output:
                'results/UGT-missense/haplotypes/MT_PREG_ID',
                'results/UGT-missense/haplotypes/MnT_PREG_ID',
                'results/UGT-missense/haplotypes/PT_PREG_ID',
                'results/UGT-missense/haplotypes/PnT_PREG_ID'
        script:
                '../scripts/phase_by_transmission.py'

rule merge_missense_hard_calls_pheno:
        'Merge each haplotype and the pheno file.'
        input:
                'results/merge_data/delivery/jaundice-transmission.txt',
                'results/UGT-missense/haplotypes/MT_PREG_ID',
                'results/UGT-missense/haplotypes/MnT_PREG_ID',
                'results/UGT-missense/haplotypes/PT_PREG_ID',
                'results/UGT-missense/haplotypes/PnT_PREG_ID'
        output:
                'results/UGT-missense/delivery/jaundice-transmitted.txt'
        run:
                d= pd.read_csv(input[0], sep= '\t', header= 0)
                df_list= list()
                for i in range(1, len(input)):
                        x= pd.read_csv(input[i], sep= '\t', header= 0)
                        varnames= ('chr' + x.chr.apply(str) + '_' + x.pos.apply(str) + '_' + x.ref + '_' + x.eff).values.tolist()
                        x= pd.DataFrame(x.iloc[:, 4:].T)
                        haplo= input[i].split('/')[-1].replace('_PREG_ID', '')
                        x.columns= [i + '_' + haplo for i in varnames]
                        x['PREG_ID_1724']= x.index
                        df_list.append(x)
                x= reduce(lambda x, y: pd.merge(x, y, on = 'PREG_ID_1724', how = 'inner'), df_list)
                x['PREG_ID_1724']= x.PREG_ID_1724.apply(str)
                d['PREG_ID_1724']= d.PREG_ID_1724.apply(str)
                x= pd.merge(x, d, on= 'PREG_ID_1724')
                x.to_csv(output[0], sep= '\t', header= True, index= False)
