library("dplyr")
library("tidyr")
library("cowplot")
library("ggrepel")
library("data.table")
library('showtext')
options(warn=-1)

colorBlindBlack8= c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

desat <- function(cols, sat=0.5) {
  X <- diag(c(1, sat, 1)) %*% rgb2hsv(col2rgb(cols))
  hsv(X[1,], X[2,], X[3,])
}

desat_colorBlindBlack8= desat(colorBlindBlack8, 0.5)

#d= fread('/mnt/work/pol/neo-jaundice/results/GWAS/delivery/MoBa-GWAS-jaundice-fets.txt.gz', h= T, select= c('ID', 'CHR', 'POS', 'LOG10P', 'nearestGene', 'rsid'))
d= fread(snakemake@input[[1]], h= T, select= c('ID', 'CHR', 'POS', 'LOG10P', 'nearestGene', 'rsid'))
d$pheno= 'Neonate'

d$GENE= ifelse(d$rsid== 'rs17868338', 'UGT1A*', ifelse(d$ID== '23:109792100:C:T', 'CHRDL1', ''))
d= filter(d, !duplicated(ID))

#x=  fread('/mnt/work/pol/neo-jaundice/results/GWAS/delivery/MoBa-GWAS-jaundice-moms.txt.gz', h= T, select= c('ID', 'CHR', 'POS', 'LOG10P', 'nearestGene', 'rsid'))
x= fread(snakemake@input[[2]], h= T, select= c('ID', 'CHR', 'POS', 'LOG10P', 'nearestGene', 'rsid'))
x$pheno= 'Mother'

x$GENE= ifelse(x$rsid == 'rs687621', 'ABO', ifelse(x$rsid== 'rs17868336', 'UGT1A*', ifelse(x$ID == '23:109840240:C:T', 'CHRDL1', '')))
x= filter(x, !duplicated(ID))

d= arrange(d, POS)

x= arrange(x, POS)


don= d %>%
    group_by(CHR)      %>%
    summarise(chr_len= max(POS)) %>%
    mutate(tot= cumsum(as.numeric(chr_len))-chr_len) %>% # Calculate cumulative position of each chromosome
    select(-chr_len) %>%
    left_join(d, ., by= 'CHR') %>%
    arrange(CHR, POS) %>% # Add a cumulative position of each SNP
    mutate(BPcum=POS+tot) %>%
         ungroup()

 don1= x %>%
    group_by(CHR)      %>%
    summarise(chr_len= max(POS)) %>%
    mutate(tot= cumsum(as.numeric(chr_len))-chr_len) %>% # Calculate cumulative position of each chromosome
    select(-chr_len) %>%
    left_join(x, ., by= 'CHR') %>%
    arrange(CHR, POS) %>% # Add a cumulative position of each SNP
    mutate(BPcum= POS+tot) %>%
         ungroup()

don= rbind(don, don1)

axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum)) / 2 )
  names(axisdf)= c('CHR', 'center')

HC= -log10(5*10**-8)

#font_add("arial", "arial.ttf", bold= 'arial_bold.ttf')

#font_add_google("Roboto", "Roboto")

showtext_opts(dpi = 300)
showtext_auto(enable = TRUE)

don$logpval= with(don, ifelse(pheno== 'Mother', -LOG10P, LOG10P))

don$CHR2= with(don, ifelse(CHR %% 2, 'odd', 'even'))
don$CHR2= paste0(don$CHR2, don$pheno, sep='_')

p1= ggplot(data= don, aes(x= BPcum, y= logpval, colour= CHR2)) +
   geom_hline(yintercept= 0, size= 0.25, colour= 'black') +
  geom_point(size= 0.07) + # Show all points
  theme_cowplot(font_size= 10) +
  scale_colour_manual(values= c(desat_colorBlindBlack8[6], desat_colorBlindBlack8[2], colorBlindBlack8[6], 
                                colorBlindBlack8[2]), guide= F) +
    scale_x_continuous(label = c(1:19, '', 21,'', 'X'), breaks= axisdf$center, expand= expansion(0)) +
  scale_y_continuous(breaks= seq(-10, 60, 10), labels= c(10, seq(0, 60, 10))) +
  ylab(expression(-log[10]~pvalue)) +
  xlab('Chromosome') +
  geom_hline(yintercept= c(HC, -HC), size= 0.2, linetype= 2, colour= '#878787') +
  coord_cartesian(clip = "off") +
  geom_text_repel(data= filter(don, GENE!= ''), aes(x= BPcum, y= logpval, label= GENE),
                  colour= c(colorBlindBlack8[2], colorBlindBlack8[2], colorBlindBlack8[6], colorBlindBlack8[6], colorBlindBlack8[6]),
                  alpha= 1,
                  size= 9/ .pt,
                  force_pull= 0, # do not pull toward data points
                  force= 0.1,
                  nudge_y      =  ifelse(filter(don, GENE!= '') %>% pull(logpval)>0, 1, -1),
                  direction    = "both",
                  hjust        = 1,
                  vjust=  0.5,
	            	  box.padding= 0,
  		            angle= 0,
                  segment.size = 0.1,
                  segment.square= TRUE,
                  segment.inflect= FALSE,
                  segment.colour= colorBlindBlack8[8],
                  segment.linetype = 4,
                  ylim = c(-Inf, 60),
                  xlim = c(-Inf, Inf),
	            	  fontface = "italic") +
  theme(legend.position= 'none',
	plot.margin = unit(c(t= 0, r=0, b= 0, l=0), 'cm'),
        text= element_text(family= "Roboto", size= 10),
	axis.line= element_line(size= 0.1),
	axis.ticks= element_line(size= 0.1))

save_plot(snakemake@output[[1]], plot= p1, base_height= 90, base_width= 180, units= 'mm', dpi= 300, type= 'cairo')



