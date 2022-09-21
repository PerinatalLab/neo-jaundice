library("dplyr")
library("tidyr")
library("cowplot")
library("ggrepel")
library("data.table")
library('showtext')
options(warn=-1)


d= fread(snakemake@input[[1]], h= T, select= c('ID', 'CHR', 'POS', 'pvalue', 'nearestGene'))
d$pheno= 'Neonate'

d$GENO= ifelse(d$rsid== 'rs17868338' 'UGT1A4', ifelse(d$ID= '23:109792100:C:T', 'RTL9', ''))
d= filter(d, !duplicated(ID))
x= fread(snakemake@input[[3]], h= T, select= c('ID', 'CHR', 'POS', 'pvalue', 'nearestGene'))
x$pheno= 'Mother'
x$GENO= ifelse(d$rsid == 'rs687621', 'ABO', ifelse(d$rsid== 'rs17868336', 'UGT1A4', ''))
x= filter(x, !duplicated(ID))


colorBlindBlack8= c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

d= arrange(d, POS)

x= arrange(d, POS)


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

font_add("arial", "arial.ttf", bold= 'arial_bold.ttf')

showtext_opts(dpi = 300)
showtext_auto(enable = TRUE)

don$logpval= with(don, ifelse(pheno== 'Mother', LOG10P, -LOG10P))

p1= ggplot(data= don, aes(x= BPcum, y= logpval, colour= pheno, alpha= factor(CHR %% 2))) +
  geom_point(size= 0.07) +   # Show all points
  theme_cowplot(font_size= 9) +
  scale_colour_manual(values= cols, guide= F) +
  scale_x_continuous(label = c(1:19, '', 21,'', 'X'), breaks= axisdf$center, expand= c(0.03, 0.03)) + # label = ifelse(axisdf$CHR== 23, 'X', axisdf$CHR)
  ylab('-log10(pvalue)') +
  xlab('Chromosome') +
  geom_hline(yintercept= 0, size= 0.25, colour= 'black') +
  geom_hline(yintercept= c(HC, -HC), size= 0.2, linetype= 2, colour= '#878787') +
  coord_cartesian(clip = "off") +
  geom_text_repel(data= filter(don, GENE!= ''), aes(x= BPcum, y= logpval, label= GENE),
                  size= 6/ .pt,
                  force_pull= 0, # do not pull toward data points
                  force= 0.1,
                  nudge_y      =  ifelse(filter(don, GENE!= '') %>% pull(logpval)>0, 1, -1), #43 - ((-log10(filter(don, GENE!= '')$pvalue))),
                  direction    = "both",
                  hjust        = 0,
                  vjust=  0.5,
		  box.padding= 0.1,
		  angle= 0,
                  segment.size = 0.1,
                  segment.square= TRUE,
                  segment.inflect= FALSE,
                  segment.colour= colorBlindBlack8[8],
                  segment.linetype = 4,
                  ylim = c(-Inf, 50),
                  xlim = c(-Inf, Inf)) +
  theme(legend.position= 'none',
	plot.margin = unit(c(t= 0, r=0, b= 0, l=0), 'cm'),
        text= element_text(family="arial", size= 9),
	axis.line= element_line(size= 0.1)) 

save_plot(snakemake@output[[1]], plot= p1, base_height= 90, base_width= 180, units= 'mm', dpi= 300)


