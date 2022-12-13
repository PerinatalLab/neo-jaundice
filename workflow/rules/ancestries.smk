rule project_1kg:
	'Project 1000 Genomes samples onto the hapmap-moba pcs.'
	input:
		expand("/mnt/hdd/data/geno/references/1000g/phase3_chr{chr}.vcf.gz", chr=range(1,23)),
		loadings="hapmap_loadings", # from /mnt/archive/MOBAGENETICS/genotypes-base/aux/pca/hapmap
		meansd="hapmap_meansd"
	output:
		"1kg_all.bim",
		"1kg_projected"
	shell:
		'''
		# hapmap_* files taken from /mnt/archive/MOBAGENETICS/genotypes-base/aux/pca/hapmap/
		awk "NR>1{{print $1}}" {input.loadings} > hapmap_snps
		echo "" > filelist.txt
		for chr in {{1..22}}
		do
			plink --vcf /mnt/hdd/data/geno/references/1000g/phase3_chr$chr.vcf.gz \
			--extract hapmap_snps --make-bed \
			--out 1kg_chr$chr
			echo "1kg_chr$chr" >> filelist.txt
		done
		plink --merge-list filelist.txt --make-bed \
			--a1-allele {input.loadings} 2 1 \
			--out 1kg_all

		./flashpca_x86-64 --bfile 1kg_all --project \
			--inmeansd {input.meansd} --inload {input.loadings} \
			--outproj 1kg_projected

		'''
		
rule project_moba:
	'Project moba samples onto the hapmap-moba pcs.'
	input:
		"1kg_all.bim",
		expand("/mnt/archive/MOBAGENETICS/genotypes-base/imputed/all/plink/{chr}.bed", chr=range(1,23)),
		loadings="hapmap_loadings",
		meansd="hapmap_meansd",
		childlist="/mnt/archive2/p1724/v12/linkage_Child_PDB1724.csv"
	output:
		"mobagen_projected",
		"rs6755571.raw"
	shell:
		'''
		awk -F';' '{{print $3, $3}}' {input.childlist} > mobagen_clean_ids.txt
		awk '{{print $2}}' 1kg_all.bim > mobagen_pca_vars.txt

		echo "" > filelist.txt
		for chr in {{1..22}}
		do
			plink --bfile /mnt/archive/MOBAGENETICS/genotypes-base/imputed/all/plink/$chr \
			--extract mobagen_pca_vars.txt \
			--keep mobagen_clean_ids.txt \ 
			--make-bed --out mobagen_pca_chr$chr
			echo mobagen_pca_chr$chr >> filelist.txt
		done
		plink --merge-list filelist.txt --make-bed \
			--a1-allele {input.loadings} 2 1 \
			--out mobagen_pca_all

		./flashpca_x86-64 --bfile mobagen_pca_all --project \
			--inmeansd {input.meansd} --inload {input.loadings} \
			--outproj mobagen_projected

		# also extract the snp
		plink --bfile /mnt/archive/MOBAGENETICS/genotypes-base/imputed/all/plink/2 \
			--snp rs6755571 --recode A \
			--out rs6755571
		'''

rule analyze_ancestries:
	input:
		"/mnt/archive2/p1724/v12/PDB1724_MFR_541_v12.csv",
		"/mnt/archive2/p1724/v12/linkage_Child_PDB1724.csv",
		"/mnt/archive/MOBAGENETICS/genotypes-base/aux/flaglist-merged/mobagen-flaglist-n99259.txt",
		"1kg_projected",
		"mobagen_projected",
		"rs6755571.raw"
	output:
		"mafs_across_1kgpops.csv",
		"suppfig_pca.png"
	script:
		"../scripts/analyze-ancestries.R"

