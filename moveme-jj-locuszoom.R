library(ggplot2)
library(dplyr)
library(tidyr)
library(ggh4x)

setwd("~/Documents/results/nj/eqtls/")

# all in hg19

bed = read.table(gzfile("Homo_sapiens.GRCh37.87.chromosome.2.gff3.gz"), sep="\t", h=F, quote="")
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
bed_genes$y_label = ifelse(bed_genes$name=="UGT1A1", bed_genes$y, bed_genes$y-3)
bed_genes$x_label = ifelse(bed_genes$name=="UGT1A1", bed_genes$start+50e3, bed_genes$start)/1e6
bed_genes$x_label[bed_genes$name=="MROH2A"] = bed_genes$x_label[bed_genes$name=="MROH2A"]+0.02


# Note: exons from aberrant transcripts etc are dropped here
bed = inner_join(bed_exons[,c("type", "start", "end", "transcript")],
                 bed_tran[,c("transcript", "gene")], by="transcript")
bed = inner_join(bed, bed_genes[,c("gene", "y", "y_label")], by="gene")
bed$source = "GWAS"
bed_genes$source = "GWAS"

# load actual data
pgwas = data.table::fread("UGT1_fets.txt")
peqc = data.table::fread("UGT1A1_eQTL_colon.txt")
peql = data.table::fread("UGT1A1_eQTL_liver.txt")
colnames(peqc) = c("CHR", "POS", "RSID", "REF", "EFF", "LOG10P", "BETA", "SE", "EAF")
colnames(peql) = c("CHR", "POS", "RSID", "REF", "EFF", "LOG10P", "BETA", "SE", "EAF")
peql$EAF = as.numeric(peql$EAF)

nrow(peqc)
nrow(peql)
nrow(pgwas)

ld = read.table("LD_variants_rs6755571.txt", h=T)
ld = separate(ld, Coord, c("CHR", "POS"), sep=":")
ld$POS = as.numeric(ld$POS)

pall = bind_rows("GWAS"=pgwas, "eQTL_colon"=peqc, "eQTL_liver"=peql, .id="source") %>%
  filter(POS>234.45e6, POS<234.75e6)
pall = left_join(pall, ld[,c("POS", "Distance","R2")])

# plot
panel_labels = group_by(pall, source) %>% summarize(LOG10P=max(LOG10P)*0.9)
pos_break_fn = function(x) if(max(x)<10) { seq(0,10,2) } else { seq(0, max(x), 20) }
pall %>%
  filter(!is.na(R2)) %>%
  ggplot(aes(x=POS/1e6, y=LOG10P)) + 
  geom_point(size=0.8, aes(col=R2)) +
  geom_point(data=filter(pall, POS==234627536), col="purple", pch=5, size=1.3) +
  facet_grid2(source~., scales="free_y", axes="all", remove_labels="x") +
  geom_segment(data=bed, aes(x=start/1e6, xend=end/1e6, y=y, yend=y), size=3, col="darkblue") +
  geom_segment(data=bed_genes, aes(x=start/1e6, xend=end/1e6, y=y, yend=y), size=0.5, col="darkblue") +
  geom_text(data=panel_labels, aes(x=234.44, label=source), hjust=0, size=3.7) +
  geom_text(data=filter(bed_genes, !name %in% c("UGT1A5", "UGT1A6", "UGT1A7", "UGT1A3", "UGT1A10")),
            aes(x=pmax(x_label, 234.44), label=name, y=y_label, vjust=1.7-0.9*grepl("UGT", name),
            hjust=-0.1+1.5*grepl("UGT", name)), size=3, col="grey30", fontface=3) +
  geom_segment(data=filter(bed_genes, name %in% c("UGT1A8", "UGT1A9", "UGT1A4")),
               aes(x=x_label-0.007, xend=start/1e6-0.001, y=y-3, yend=y), size=0.3, col="grey30") +
  geom_segment(data=filter(bed_genes, name=="UGT1A1"),
               aes(x=x_label-0.029, xend=end/1e6+0.001, y=y-1, yend=y), size=0.3, col="grey30") +
  coord_cartesian(xlim=c(234.45, 234.72)) + 
  scale_color_gradient(low="#5782AD", high="#ED1330", name=expression(R^2)) + 
  scale_y_continuous(breaks=pos_break_fn, name=expression(-log[10]~p)) +
  force_panelsizes(rows=c(1,1,2)) +
  theme_minimal() + xlab("position, Mbp") +
  theme(panel.grid.major.x=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(fill=NA, colour="grey60"),
        axis.ticks = element_line(colour="grey30"),
        strip.text = element_blank(),
        axis.line.x=element_line(colour="grey30"))
ggsave("plot_eqtllocus.png", width=8, height=4)
