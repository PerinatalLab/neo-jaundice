set -e
cd /mnt/hdd/data/liver-rnaseq/

# 0. download references
wget http://ftp.ensembl.org/pub/release-107/gtf/homo_sapiens/Homo_sapiens.GRCh38.107.gtf.gz
wget http://ftp.ensembl.org/pub/release-107/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.toplevel.fa.gz
# also install hisat2, download and untar their grch38 index

# 1. download data
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR613/004/ERR6130614/ERR6130614.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR986/000/ERR9867660/ERR9867660_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR986/000/ERR9867660/ERR9867660_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR986/003/ERR9867663/ERR9867663_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR986/003/ERR9867663/ERR9867663_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR986/004/ERR9867664/ERR9867664_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR986/004/ERR9867664/ERR9867664_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR986/001/ERR9867661/ERR9867661_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR986/001/ERR9867661/ERR9867661_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR613/005/ERR6130625/ERR6130625_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR613/005/ERR6130625/ERR6130625_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR613/007/ERR6130627/ERR6130627_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR613/007/ERR6130627/ERR6130627_2.fastq.gz

# 2. run fastqc
/mnt/hdd/soft/FastQC/fastqc ERR9867660_1.fastq.gz
/mnt/hdd/soft/FastQC/fastqc ERR9867660_2.fastq.gz
/mnt/hdd/soft/FastQC/fastqc ERR9867663_1.fastq.gz
/mnt/hdd/soft/FastQC/fastqc ERR9867663_2.fastq.gz
/mnt/hdd/soft/FastQC/fastqc ERR6130614.fastq.gz
cd /mnt/hdd/common/liver-rnaseq/qc/
mv /mnt/hdd/data/liver-rnaseq/*fastqc.html ./
mv /mnt/hdd/data/liver-rnaseq/*fastqc.zip ./
multiqc .

# adapter percentage <1%, so ignoring
# GC deviations etc at the start are, hopefully, due to start position bias of the library prep

# 3. map
# HISAT2: somehow broken? freezes on 1 thread?
#cd /mnt/hdd/common/liver-rnaseq/mapping/
#/mnt/hdd/soft/Hisat2/hisat2 -p 8 -x /mnt/quick/jjuod/grch38/genome \
#	-U /mnt/hdd/data/liver-rnaseq/ERR6130614.fastq.gz | \
#	samtools view -bh - | \
#	samtools sort -@2 - > ERR6130614.bam
#/mnt/hdd/soft/Hisat2/hisat2 -p 8 -x /mnt/quick/jjuod/grch38/genome \
#	-1 /mnt/hdd/data/liver-rnaseq/ERR9867660_1.fastq.gz \
#	-2 /mnt/hdd/data/liver-rnaseq/ERR9867660_2.fastq.gz | \
#	samtools view -bh - | \
#	samtools sort -@2 - > ERR9867660.bam
#/mnt/hdd/soft/Hisat2/hisat2 -p 8 -x /mnt/quick/jjuod/grch38/genome \
#	-1 /mnt/hdd/data/liver-rnaseq/ERR9867663_1.fastq.gz \
#	-2 /mnt/hdd/data/liver-rnaseq/ERR9867663_2.fastq.gz | \
#	samtools view -bh - | \
#	samtools sort -@2 - > ERR9867663.bam

# ALT:
# Create genome index w/ STAR:
cd /mnt/quick/jjuod/grch38
/mnt/hdd/soft/Star/STAR --runThreadN 10 \
	--runMode genomeGenerate --genomeDir starix/ \
	--genomeFastaFiles Homo_sapiens.GRCh38.dna.chromosome.2.fa \
	--sjdbGTFfile Homo_sapiens.GRCh38.107.gtf

# map reads w/ STAR:
/mnt/hdd/soft/Star/STAR --runThreadN 10 --genomeDir starix/ \
	--readFilesIn /mnt/hdd/data/liver-rnaseq/ERR6130614.fastq.gz \
	--readFilesCommand zcat \
	--outFileNamePrefix /mnt/hdd/common/liver-rnaseq/mapping/ERR6130614
# note: paired end reads in particular are really slow,
# (b/c mapping on one chr only) so reducing the seed numbers
/mnt/hdd/soft/Star/STAR --runThreadN 10 --genomeDir starix/ \
	--readFilesIn /mnt/hdd/data/liver-rnaseq/ERR9867660_1.fastq.gz \
	/mnt/hdd/data/liver-rnaseq/ERR9867660_2.fastq.gz \
	--readFilesCommand zcat \
	--outFileNamePrefix /mnt/hdd/common/liver-rnaseq/mapping/ERR9867660 \
	--seedPerWindowNmax 20 --seedPerReadNmax 200
/mnt/hdd/soft/Star/STAR --runThreadN 10 --genomeDir starix/ \
	--readFilesIn /mnt/hdd/data/liver-rnaseq/ERR9867663_1.fastq.gz \
	/mnt/hdd/data/liver-rnaseq/ERR9867663_2.fastq.gz \
	--readFilesCommand zcat \
	--outFileNamePrefix /mnt/hdd/common/liver-rnaseq/mapping/ERR9867663 \
	--seedPerWindowNmax 20 --seedPerReadNmax 200
/mnt/hdd/soft/Star/STAR --runThreadN 10 --genomeDir starix/ \
	--readFilesIn /mnt/hdd/data/liver-rnaseq/ERR9867664_1.fastq.gz \
	/mnt/hdd/data/liver-rnaseq/ERR9867664_2.fastq.gz \
	--readFilesCommand zcat \
	--outFileNamePrefix /mnt/hdd/common/liver-rnaseq/mapping/ERR9867664 \
	--seedPerWindowNmax 20 --seedPerReadNmax 200
/mnt/hdd/soft/Star/STAR --runThreadN 10 --genomeDir starix/ \
	--readFilesIn /mnt/hdd/data/liver-rnaseq/ERR9867661_1.fastq.gz \
	/mnt/hdd/data/liver-rnaseq/ERR9867661_2.fastq.gz \
	--readFilesCommand zcat \
	--outFileNamePrefix /mnt/hdd/common/liver-rnaseq/mapping/ERR9867661 \
	--seedPerWindowNmax 20 --seedPerReadNmax 200

/mnt/hdd/soft/Star/STAR --runThreadN 10 --genomeDir starix/ \
	--readFilesIn /mnt/hdd/data/liver-rnaseq/ERR6130625_1.fastq.gz \
	/mnt/hdd/data/liver-rnaseq/ERR6130625_2.fastq.gz \
	--readFilesCommand zcat \
	--outFileNamePrefix /mnt/hdd/common/liver-rnaseq/mapping/ERR6130625 \
	--seedPerWindowNmax 20 --seedPerReadNmax 200

/mnt/hdd/soft/Star/STAR --runThreadN 10 --genomeDir starix/ \
	--readFilesIn /mnt/hdd/data/liver-rnaseq/ERR6130627_1.fastq.gz \
	/mnt/hdd/data/liver-rnaseq/ERR6130627_2.fastq.gz \
	--readFilesCommand zcat \
	--outFileNamePrefix /mnt/hdd/common/liver-rnaseq/mapping/ERR6130627 \
	--seedPerWindowNmax 20 --seedPerReadNmax 200

# TODO should really find a paired variant for the first sample.

# 4. sort the aligned reads
# NOTE: paired ones must be sorted by name for featureCounts
cd /mnt/hdd/common/liver-rnaseq/mapping

samtools view -bh ERR6130614Aligned.out.sam | samtools sort -@6 - -o ERR6130614.bam
rm ERR6130614Aligned.out.sam 

samtools view -bh ERR9867660Aligned.out.sam | samtools sort -@6 -n - -o ERR9867660.bam
rm ERR9867660Aligned.out.sam 

samtools view -bh ERR9867663Aligned.out.sam | samtools sort -@6 -n - -o ERR9867663.bam
rm ERR9867663Aligned.out.sam 

samtools view -bh ERR9867664Aligned.out.sam | samtools sort -@6 -n - -o ERR9867664.bam
rm ERR9867664Aligned.out.sam 

samtools view -bh ERR9867661Aligned.out.sam | samtools sort -@6 -n - -o ERR9867661.bam
rm ERR9867661Aligned.out.sam 

samtools view -bh ERR6130625Aligned.out.sam | samtools sort -@6 -n - -o ERR6130625.bam
rm ERR6130625Aligned.out.sam 

samtools view -bh ERR6130627Aligned.out.sam | samtools sort -@6 -n - -o ERR6130627.bam
rm ERR6130627Aligned.out.sam 

# 5. extract total counts for normalization
awk '$0~/Uniquely mapped reads number/{print substr(FILENAME, 4, 7),$NF}' ERR*Log.final.out > total_counts.txt

