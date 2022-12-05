library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(ggrepel)
library(ggh4x)
library('showtext')

showtext_opts(dpi = 300)
showtext_auto(enable = TRUE)


colorBlindBlack8= c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


d= fread(snakemake@input[[1]])
pph= fread(snakemake@input[[2]])
eqtls= filter(pph, PP.H4.abf > 0.9) %>% pull(eqtl_data)
eqtls= c(eqtls, 'GTEx_ge_liver')

d= filter(d, gene== 'ENSG00000241635', eqtl_data %in% eqtls)

link= fread(snakemake@input[[3]]) %>% select(ID, POS, rsid)

link= separate(link, ID, into= c('CHR', 'POS2', 'REF', 'EFF'), sep= ':')

d= inner_join(d, link, by= c('snp'= 'POS'))

d$POS2= as.numeric(d$POS2)
d= filter(d, POS2 > 234.45e6, POS2 < 234.75e6)

ld = read.table(snakemake@input[[4]], h=T)
ld= filter(ld, BP_A== 234627536 | BP_B == 234627536)
ld$BP_A= ifelse(ld$BP_A== 234627536, ld$BP_B, ld$BP_A)

ld= select(ld, CHR_A, BP_A, R2)

ld= rbind(ld, data.frame(CHR_A= 2, BP_A= 234627536, R2= 1))

names(ld)= c('CHR', 'POS2', 'R2')

pall = left_join(d, ld[,c("POS2", "R2")]) %>% filter(!is.na(R2))

pall$z.df2= with(pall, ifelse(z.df1<0, z.df2* -1, z.df2))
pall$z.df1= with(pall, ifelse(z.df1<0, z.df1* -1, z.df1))

pall$label= ifelse(pall$SNP.PP.H4> 0.1 & pall$eqtl_data!= 'GTEx_ge_liver', paste0(pall$rsid, ' (', round(pall$SNP.PP.H4, 2), ')'), NA)

pall$eqtl_data= factor(pall$eqtl_data, levels= c("CEDAR_microarray_transverse_colon", "GTEx_ge_colon_transverse", "GTEx_ge_liver"),
                       labels= c("Tranverse colon\n(microarray)", "Tranverse colon\n(RNA-seq)", "Liver\n(RNA-seq)") )

p1= ggplot(pall, aes(z.df1, z.df2, colour= R2)) +
  geom_point(size= 0.5) +
  facet_grid(vars(eqtl_data)) +
  geom_point(data= filter(pall, POS2==234627536), col="purple", pch=18, size= 3) +
  geom_text_repel(aes(label= label), direction= 'y', hjust= 0, nudge_x= 1) +
  theme_minimal() +
  geom_hline(yintercept= 0) +
  geom_vline(xintercept= 0) +
  scale_color_gradient(low= "#5782AD", high= "#ED1330", name= expression(R^2) ) +
  xlab(expression(Z-score~on~neonatal~jaundice)) +
  xlim(0, 18) +
  ylab(expression(Z-score~on~italic(UGT1A1)~cis-eQTLS)) + 
  theme(text= element_text(family= "Roboto", size= 10, colour= 'black'),
        plot.margin = unit(c(t= 0, r=0, b= 0, l=0), 'cm'),
        legend.key.size= unit(4, 'mm'),
        legend.title = element_text(size= 6), #change legend title font size
        legend.text = element_text(size=6),
        axis.text=element_text(colour="black"))

ggsave(snakemake@output[[1]], plot= p1, width= 180, height= 180, units= 'mm', dpi= 300)
