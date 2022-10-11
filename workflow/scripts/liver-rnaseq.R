# nj
library(dplyr)
library(tidyr)
library(ggplot2)

#BiocManager::install("Rsubread")
library(Rsubread)
setwd("~/Documents/results/nj/liver-rnaseq/")

annotsl = data.table::fread("Homo_sapiens.GRCh38.107_chr2.gtf")
colnames(annotsl) = c("chr", "source", "type", "start", "end", "V6", "V7", "V8", "description")
annotsl = annotsl[grep("UGT", annotsl$description),]
annotsl$gene_name = gsub('.*gene_name \"(.*?)\";.*', "\\1", annotsl$description)
annotsl$transcript_name = ifelse(annotsl$type=="exon", gsub('.*transcript_name \"(.*?)\";.*', "\\1", annotsl$description), NA)
annotsl$exon_number = ifelse(annotsl$type=="exon", gsub('.*exon_number \"(.*?)\";.*', "\\1", annotsl$description), NA)

total_counts = read.table("total_counts.txt")

# fc = featureCounts("ERR6130614.bam", annot.ext="Homo_sapiens.GRCh38.107_chr2.gtf",
#                    isGTFAnnotationFile=T, useMetaFeatures=F, GTF.featureType="exon",
#                    allowMultiOverlap=T) # very important - to allow the multiply-marked exons!

bams = list.files(pattern="*.bam")
bams = bams[bams!="ERR6130614.bam"]  # remove the single-end case
counts2 = tibble()
for(bam in bams){
  print(bam)
  fc = featureCounts(bam, annot.ext="Homo_sapiens.GRCh38.107_chr2.gtf",
                     isGTFAnnotationFile=T, useMetaFeatures=F, GTF.featureType="exon",
                     allowMultiOverlap=T, isPairedEnd=T)
  print(fc$stat)
  
  annots = fc$annotation
  counts = fc$counts
  counts = data.frame(cnt=counts[,1])
  counts = bind_cols(counts, annots)
  
  # features at same positions get same count entries
  counts = group_by(counts, Start, End) %>%
    summarize(nn=n(), cnt_min=min(cnt), cnt_max=max(cnt))
  counts$bam = substr(bam, 4, 10)
  stopifnot(all(counts$cnt_min==counts$cnt_max))
  
  totalcount = total_counts$V2[total_counts$V1==counts$bam[1]]
  
  # normalize to CPM (reads mapped to chr 2)
  counts$cnt_min = counts$cnt_min / totalcount * 1e6
  counts$cnt_max = counts$cnt_max / totalcount * 1e6 
  
  counts2 = bind_rows(counts2, counts)
}

counts2 = filter(annotsl, type=="exon") %>%
  left_join(counts2, by=c("start"="Start", "end"="End"))

filter(counts2, cnt_min>0) %>%
  ggplot(aes(x=start/1e3, xend=end/1e3, y=transcript_name, yend=transcript_name)) +
  geom_segment(aes(col=cnt_min), size=5) +
  geom_segment(data=filter(counts2, cnt_min==0), col="black", size=5) +
  scale_color_distiller(palette = "OrRd") + 
  facet_wrap(~bam) +
  theme_minimal() + xlab("kbp") + theme(panel.background = element_rect(fill="gray80"))

filter(counts2, cnt_min>0) %>%
  ggplot(aes(x=start/1e3, xend=end/1e3, y=transcript_name, yend=transcript_name)) +
  geom_segment(aes(col=cnt_min), size=5) +
  scale_color_distiller(palette = "OrRd") + 
  theme_minimal() + xlab("kbp") + theme(panel.background = element_rect(fill="gray80"))

filter(counts2, exon_number==1)

celltypes = data.frame(bam=as.character(c(9867660, 9867661,  9867663, 9867664,  6130625, 6130627)),
                       cells=c("forward progr. hepatoc. from stem cells",
                               "forward progr. hepatoc. from stem cells",
                               "fetal liver cell, direct", "fetal liver cell, direct",
                               "adult liver, cultured", "adult liver, cultured"))
counts3 = left_join(counts2, celltypes, by="bam") %>%
  group_by(cells, gene_name, transcript_name, exon_number, start, end) %>%
  summarise(n=n(), cnt_mean=mean(cnt_min), cnt_sd=sd(cnt_min))
filter(counts3, cnt_sd>0.7*cnt_mean) %>% View  # mostly similar across replicates
# except the first UGT1A1 exon in adults, hmm

filter(counts3, cnt_mean>0) %>%
  ggplot(aes(x=start/1e3, xend=end/1e3, y=transcript_name, yend=transcript_name)) +
  geom_segment(aes(col=cnt_mean), size=5) +
  geom_segment(data=filter(counts3, cnt_mean==0), col="black", size=5) +
  scale_color_distiller(palette = "OrRd") + 
  facet_wrap(~cells) +
  theme_minimal() + xlab("kbp") + theme(panel.background = element_rect(fill="gray80"))
