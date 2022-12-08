library(ggplot2)
library(dplyr)
library(tidyr)
library(ggh4x)
library(data.table)

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
bed_genes$y = -2*bed_genes$y
#bed_genes$y_label = ifelse(bed_genes$name=="UGT1A1", bed_genes$y, bed_genes$y-3)
bed_genes$y_label = ifelse(grepl('UGT', bed_genes$name), bed_genes$y -4, bed_genes$y)
#bed_genes$x_label = ifelse(bed_genes$name=="UGT1A1", bed_genes$start+50e3, bed_genes$start)/1e6
bed_genes$x_label = ifelse(bed_genes$name=="UGT1A1", bed_genes$start, bed_genes$start)/1e6
bed_genes$x_label[bed_genes$name=="MROH2A"] = bed_genes$x_label[bed_genes$name=="MROH2A"]

# Note: exons from aberrant transcripts etc are dropped here
bed = inner_join(bed_exons[,c("type", "start", "end", "transcript")],
                 bed_tran[,c("transcript", "gene")], by="transcript")
bed = inner_join(bed, bed_genes[,c("gene", "y", "y_label")], by="gene")
bed$source = "Neonatal\njaundice"
bed_genes$source = "Neonatal\njaundice"

# load actual data
pgwas = fread(snakemake@input[[2]])
pgwas= pgwas %>% filter(substr(MarkerName, 1, 2) == "2:")

pgwas = separate(pgwas, MarkerName, into= c('CHR', 'POS'), sep= ':')
names(pgwas)[11]= 'pvalue'
pgwas$LOG10P= -log10(pgwas$pvalue)

pgwas$CHR= as.numeric(pgwas$CHR)
pgwas$POS= as.numeric(pgwas$POS)

nrow(pgwas)

ld = read.table(snakemake@input[[3]], h=T)
ld= filter(ld, BP_A== 234627536 | BP_B == 234627536)
ld$BP_A= ifelse(ld$BP_A== 234627536, ld$BP_B, ld$BP_A)

ld= select(ld, CHR_A, BP_A, R2)

ld= rbind(ld, data.frame(CHR_A= 2, BP_A= 234627536, R2= 1))

names(ld)= c('CHR', 'POS', 'R2')
pall = bind_rows("Neonatal\njaundice"=pgwas, .id="source") %>%
  filter(POS > 234.45e6, POS < 234.75e6)

pall = left_join(pall, ld[, c("POS", "R2")])


# plot
panel_labels = group_by(pall, source) %>% summarize(LOG10P=max(LOG10P)*0.9)
pos_break_fn = function(x) if(max(x)<10) { seq(0,10,2) } else { seq(0, 15, 5) }

p1= pall %>%
filter(!is.na(R2)) %>%
  ggplot(aes(x= POS/1e6, y=LOG10P)) +
    geom_point(size= 0.6, aes(col= R2)) +
  geom_point(data= filter(pall, POS== 234627536), col="purple", pch=18, size= 3) +
  facet_grid2(source~., scales="free_y", axes="all", remove_labels="x") +
  geom_segment(data= bed, aes(x=start/1e6, xend=end/1e6, y=y, yend=y), size=3, col="darkblue") +
  geom_segment(data= bed_genes, aes(x=start/1e6, xend=end/1e6, y=y, yend=y), size=0.5, col="darkblue") +
  geom_text(data= panel_labels, aes(x=234.44, label=source), hjust=0, size=3.5) +
  geom_text(data= filter(bed_genes, !name %in% c("UGT1A5", "UGT1A6", "UGT1A7", "UGT1A3", "UGT1A10")),
            aes(x= pmax(x_label, 234.44), label= name, y= y, vjust= 1.7-0.9 * grepl("UGT", name),
                hjust= -0.1+1.2*grepl("UGT", name)), size= 2.5, col="grey30", fontface=3) +
  coord_cartesian(xlim= c(234.45, 234.72)) +
  scale_color_gradient(low= "#5782AD", high= "#ED1330", name= expression(R^2)) +
  scale_y_continuous(breaks= pos_break_fn, name= expression(-log[10]~pvalue)) +
  theme_minimal() + 
  xlab("Position, Mbp") +
  theme(text= element_text(family= "Roboto", size= 10),
	panel.grid.major.x=element_blank(), 
        panel.grid.minor=element_blank(),
        panel.background = element_rect(fill=NA, colour="grey60"),
        axis.ticks = element_line(colour="grey30", linewidth = 0.1),
        strip.text = element_blank(),
        axis.line.x=element_line(colour="grey30", linewidth= 0.4),
	legend.key.size= unit(4, 'mm'),
	legend.title = element_text(size= 6), #change legend title font size
        legend.text = element_text(size=6),
	plot.margin = unit(c(t= 0, r=0, b= 0, l=0), 'cm'),
	axis.text=element_text(colour="black"))

ggsave(snakemake@output[[1]], plot= p1, width= 180, height= 80, units= 'mm', dpi= 300)

