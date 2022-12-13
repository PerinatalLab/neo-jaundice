rule list_samples_plink:
	'Get list of non-related maternal samples for PLINK.'
	input:
		'results/effect_origin/delivery/jaundice.txt'
	output:
		'results/LD/aux/fets_samples.txt',
		'results/LD/aux/UGT_region.txt'
	run:
		d= pd.read_csv(input[0], sep= '\t', header= 0)
		d['IID']= d.Child
		d.to_csv(output[0], sep= '\t', header= False, index= False, columns = ['Child', 'IID'])
		x= pd.DataFrame({'CHR': [2], 'start': [234450000], 'end': [234750000], 'ID': ['id1']})
		x.to_csv(output[1], sep= '\t', header= False, index= False)

rule plink_LD_UGT:
	'Estimate R2 using plink.'
	input:
		'results/LD/aux/fets_samples.txt',
		'results/LD/aux/UGT_region.txt',
		expand('/mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/imputed/all/plink/2.{ext}', ext= ['bim', 'fam', 'bed'])
	output:
		'results/LD/delivery/UGT1A.ld'
	params:
		'/mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/imputed/all/plink/2',
		'results/LD/delivery/UGT1A'
	shell:
		"/home/pol.sole.navais/soft/plink --bfile {params[0]} --keep {input[0]} --extract range {input[1]} --r2 --ld-window-r2 0 --ld-window 7000 --out {params[1]}"
