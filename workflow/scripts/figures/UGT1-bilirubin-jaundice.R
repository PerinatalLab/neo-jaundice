library(ggplot2)
library(dplyr)
library(tidyr)
library("ggh4x")
library(data.table)

setwd("/mnt/work/pol/neo-jaundice/")

# all in hg19

bed = read.table(gzfile(snakemake@input[[1]]), sep="\t", h=F, quote="")

colnames(bed) = c("CHR", "predictor", "type", "start", "end", "V6", "strand", "V8", "ANN")
bed = filter(bed, end>234.45e6, start<234.75e6)
nrow(bed)
bed = filter(bed, type %in% c("gene", "exon", "mRNA"))

# assign exons to genes
bed_exons = filter(bed, type=="exon")
bed_tran = filter(bed, type=="mRNA")
bed_genes = filter(bed, type=="gene")
bed_exons$transcript = gsub(".*transcript:(.*?);.*", "\\1", bed_exons$ANN)
bed_tran$transcript = gsub(".*transcript_id=(.*?);.*", "\\1", bed_tran$ANN)
bed_tran$gene = gsub(".*Parent=gene:(.*?);.*", "\\1", bed_tran$ANN)
bed_genes$gene = gsub(".*ID=gene:(.*?);.*", "\\1", bed_genes$ANN)
bed_genes$name = gsub(".*Name=(.*?);.*", "\\1", bed_genes$ANN)
bed_genes = filter(bed_genes, name!="AC114812.8")

bed_genes$y = seq_along(bed_genes$name)
bed_genes$y[bed_genes$name=="USP40"] = 2
bed_genes$y[bed_genes$name=="HJURP"] = 3
bed_genes$y[bed_genes$name=="MROH2A"] = 2
bed_genes$y = -5*bed_genes$y
bed_genes$y_label = ifelse(grepl('UGT', bed_genes$name), bed_genes$y -3, bed_genes$y)
bed_genes$x_label = ifelse(bed_genes$name=="UGT1A1", bed_genes$start, bed_genes$start)/1e6
bed_genes$x_label[bed_genes$name=="MROH2A"] = bed_genes$x_label[bed_genes$name=="MROH2A"]+ 0.02


# Note: exons from aberrant transcripts etc are dropped here
bed = inner_join(bed_exons[,c("type", "start", "end", "transcript")],
                 bed_tran[,c("transcript", "gene")], by="transcript")
bed = inner_join(bed, bed_genes[,c("gene", "y", "y_label")], by="gene")
bed$source = "Neonatal jaundice"
bed_genes$source = "Neonatal jaundice"

# load actual data

d= fread(snakemake@input[[2]])
d= filter(d, gene== 'ENSG00000241635', eqtl_data== 'GTEx_ge_colon_transverse')

link= fread(snakemake@input[[3]]) %>% select(ID, POS)

link= separate(link, ID, into= c('CHR', 'POS2', 'REF', 'EFF'), sep= ':')

d= inner_join(d, link, by= c('snp'= 'POS'))

d$POS2= as.numeric(d$POS2)
d= filter(d, POS2 > 234.45e6, POS2 < 234.75e6)

ld = read.table("resources/LD_variants_rs6755571.txt", h=T)
ld = separate(ld, Coord, c("CHR", "POS2"), sep=":")
ld$POS2 = as.numeric(ld$POS2)

pall = left_join(d, ld[,c("POS2", "Distance","R2")])

# plot
panel_labels = group_by(pall, source) %>% summarize(LOG10P=max(LOG10P)*0.9)
pos_break_fn = function(x) if(max(x)<10) { seq(0,10,2) } else { seq(0, max(x), 20)}

p1= pall %>%
  filter(!is.na(R2)) %>%
  ggplot(aes(x=POS2/1e6, y=SNP.PP.H4 * 1000)) +
    geom_point(size=0.8, aes(col=R2)) +
  geom_point(data=filter(pall, POS2==234627536), col="purple", pch=5, size=1.3) +
  geom_segment(data=bed, aes(x=start/1e6, xend=end/1e6, y=y, yend=y), size=3, col="darkblue") +
  geom_segment(data=bed_genes, aes(x=start/1e6, xend=end/1e6, y=y, yend=y), size=0.5, col="darkblue") +
  geom_text(data=filter(bed_genes, !name %in% c("UGT1A5", "UGT1A6", "UGT1A7", "UGT1A3", "UGT1A10")),
            aes(x=pmax(x_label, 234.44), label= name, y= y, vjust= 1.7-0.9 * grepl("UGT", name),
                hjust=-0.1+1.2*grepl("UGT", name)), size=3, col="grey30", fontface=3) +
  coord_cartesian(xlim=c(234.45, 234.72)) +
  scale_color_gradient(low="#5782AD", high="#ED1330", name=expression(R^2)) +
  scale_y_continuous(breaks=pos_break_fn, name= expression(-log[10]~p)) +
  theme_minimal() + 
  xlab("position, Mbp") +
  theme(panel.grid.major.x=element_blank(), 
        panel.grid.minor=element_blank(),
        panel.background = element_rect(fill=NA, colour="grey60"),
        axis.ticks = element_line(colour="grey30", linewidth = 0.1),
        strip.text = element_blank(),
        axis.line.x=element_line(colour="grey30", linewidth= 0.4))

ggsave("plot_eqtllocus.png", width= 8, height= 4)



ggsave(snakemake@output[[1]], plot= p1, width= 180, height= 90, units= 'mm', dpi= 300)

