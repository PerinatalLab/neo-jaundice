library(data.table)
library(dplyr)
library(broom)
library(ggplot2)
library(showtext)
library(cowplot)
library(ggh4x)
library(ggrepel)

colorBlindBlack8= c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

showtext_opts(dpi = 300)
showtext_auto(enable = TRUE)

SNIN2 = "~/Documents/results/nj/eqtls/UGT1_fets.txt"
SNIN3 = "resources/bilirubin/delivery/total-bilirubin-GWAS.txt.gz"
SNIN4 = "results/LD/delivery/UGT1A.ld"


# -----------------------
# Plot locusZoom for fetal vs adult scores
# and frequency of jaundice and 95%CI by polygenic score


# raw gwas results:
resR = fread(snakemake@input[[1]])
resR = filter(resR, CHR== 2, POS > 234.45e6, POS < 234.75e6) %>% filter(!duplicated(ID))
minpos = min(resR$POS)
maxpos = max(resR$POS)

# adult summaries:
resA = fread(snakemake@input[[2]])
resA = filter(resA, CHROM==2, POS > minpos, POS < maxpos)

# convert string p-vals:
#pvals = strsplit(resA$pvalue, split="e")
#pvals = lapply(pvals, function(x) if(length(x)==1) log10(as.numeric(x)) 
#               else log10(as.numeric(x[[1]])) + as.numeric(x[[2]]))
resA$LOG10P= -1*(pnorm(-abs(resA$BETA / resA$SE),log.p = T)*1/log(10) + log10(2))
#resA$LOG10P = -as.numeric(pvals)

resA = filter(resA, POS %in% resR$POS)

topsnpA = resA$RSID[which.max(resA$LOG10P)]  # rs887829

res = bind_rows("Adult"=resA, "Neonatal"=resR, .id="FA")

# LD
ld = read.table(snakemake@input[[3]], h=T)
ld = filter(ld, SNP_A==topsnpA | SNP_B==topsnpA)
ld$POS = ifelse(ld$SNP_A==topsnpA, ld$BP_B, ld$BP_A)
ld = ld[,c("POS", "R2")]
ld = bind_rows(ld, data.frame(POS=resA$POS[resA$RSID==topsnpA], R2=1))

res = inner_join(res, ld, by="POS")

res = filter(res, !is.na(R2))
# plot
panel_labels = group_by(res, FA) %>% summarize(LOG10P=pmax(max(LOG10P),65)*0.95)
panel_scales = list(
  scale_y_continuous(),
  scale_y_continuous(limits=c(0, 65))
)

# snps showing LD w/ rs887829
p_FA = ggplot(res, aes(x=POS/1e6, y=LOG10P)) + 
  geom_point(aes(col=R2), size= 0.6) +
  geom_point(data=filter(res, POS==resA$POS[resA$RSID==topsnpA]),
             col="purple", pch=18, size=2.5) +
  facet_grid(FA~., scales="free_y") +
  geom_text(data=panel_labels, aes(label=FA, x=234.47), hjust=0) + 
  scale_color_gradient(low="#5782AD", high="#ED1330", name=expression(R^2)) +
  scale_x_continuous(expand = c(0,0), limits=c(234.46, 234.72)) +
  facetted_pos_scales(y=panel_scales) + 
  theme_bw() + xlab("Position, Mbp") + ylab(expression(-log[10]~p)) + 
  theme(text= element_text(family= "Arial", size= 10),
	panel.grid.major.x=element_blank(), 
        panel.grid.minor=element_blank(),
        panel.background = element_rect(fill=NA, colour="grey60"),
        axis.ticks = element_line(colour="grey30", linewidth = 0.1),
        strip.text = element_blank(),
        axis.line.x=element_line(colour="grey30", linewidth= 0.4),
	legend.key.size= unit(4, 'mm'),
        legend.title = element_text(size= 6), #change legend title font size
        legend.text = element_text(size=6),
	axis.text= element_text(color= 'black'))

# ---- PGS subfigure ----

d= fread(snakemake@input[[4]])

d$PGS_cat= ntile(d$fets_jaundice, 10)
d$PGS_cat_nochr2= ntile(d$fets_jaundice_nochr2, 10)

x= group_by(d, PGS_cat) %>%
  summarize(p= mean(jaundice, na.rm= T), lo95= binom.test(sum(jaundice), n())$conf.int[1], 
            up95= binom.test(sum(jaundice), n())$conf.int[2])

names(x)[1]= 'PGS'

x$PGS= factor(x$PGS)
x$chr= 'Full'

x2= group_by(d, PGS_cat_nochr2) %>%
  summarize(p= mean(jaundice, na.rm= T), lo95= binom.test(sum(jaundice), n())$conf.int[1], 
            up95= binom.test(sum(jaundice), n())$conf.int[2])

names(x2)[1]= 'PGS'
x2$chr= 'No UGT1A*'

x2$PGS= factor(x2$PGS)

x= rbind(x, x2)

p2= ggplot(data= x, aes(PGS, p*100, colour= chr, group= chr)) +
  geom_point(position= position_dodge2(width= 0.5)) +
  geom_linerange( aes(x= PGS, ymin= lo95*100, ymax= up95*100, colour= chr), alpha= 1, size= 0.7, 
                  position= position_dodge2(width= 0.5)) +
  scale_color_manual(values= colorBlindBlack8[c(6, 2)], name=NULL, labels= c('Full', expression(Excluding~italic('UGT1A*')~region))) +
  theme_cowplot(font_size= 10) + 
  xlab('Adult bilirubin polygenic\nscore decile') + 
  ylab('Jaundice prevalence, %') +
  scale_y_continuous(limits= c(0, 12), expand= expansion(add= c(0, 0.0002))) +
  theme(text= element_text(size= 10),
        axis.ticks.x= element_blank(),
        axis.line = element_line(color = "black", size = 0.2, lineend = "square"),
        axis.ticks.y = element_line(color = "black", size = 0.2),
        legend.position = c(0.1, 0.1),
	legend.text.align = 0,
	axis.text= element_text(color= 'black'))

p3= plot_grid(p_FA, p2, labels="AUTO", rel_widths = c(6,4), rel_heights= c(6, 4))

save_plot(snakemake@output[[1]], plot= p3, base_height= 80, base_width= 180, units= 'mm', dpi= 300)


# -----------------
# Plot density of polygenic score

p1= ggplot() +
  geom_density(data= d, aes(fets_jaundice), fill= colorBlindBlack8[6], size= 0.1, alpha= 0.6, color= 'black') +
  geom_density(data= d, aes(fets_jaundice_nochr2), fill= colorBlindBlack8[2], size= 0.1, alpha= 0.6, color= 'black') +
  theme_cowplot(font_size= 10) + 
  xlab('Polygenic score of adult bilirubin levels') + 
  ylab('Denstiy') +
  scale_y_continuous(limits= c(0, 5.5), expand= expansion(add= c(0, 0.0002))) +
  theme(axis.text.x= element_text(size= 8),
        axis.ticks.x= element_blank(),
        axis.line = element_line(color = "black", size = 0.2, lineend = "square"),
        axis.ticks.y = element_line(color = "black", size = 0.2))

save_plot(snakemake@output[[2]], plot= p1, base_height= 60, base_width= 110, units= 'mm', dpi= 300)
