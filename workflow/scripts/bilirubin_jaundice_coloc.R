library("dplyr")
library("knitr")
library("tidyr")
library(cowplot)
library(ggrepel)
library("data.table")
library('showtext')
library(ggtern)
options(warn=-1)

colorBlindBlack8= c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

font_add("arial", "arial.ttf", bold= 'arial_bold.ttf')

showtext_opts(dpi = 300)
showtext_auto(enable = TRUE)

d= fread(snakemake@input[[1]])

colT= colorBlindBlack8[4]
colR= colorBlindBlack8[1]
colL= colorBlindBlack8[2]

shbg$One_or_Other= shbg$PP.H0.abf + shbg$PP.H1.abf + shbg$PP.H2.abf
shbg$coloc= shbg$PP.H4.abf
shbg$shared_locus= shbg$PP.H3.abf

p1= ggtern(shbg, aes(One_or_Other, coloc, shared_locus)) +
geom_point(colour= colorBlindBlack8[8], fill= colorBlindBlack8[8], shape= 21) +
scale_alpha_continuous(range= c(0.6, 1), guide= F) +
scale_size_continuous(range= c(.001, 10), guide= F) +
theme_custom(tern.plot.background = NULL, tern.panel.background = 'white', col.T = colT, col.L = colL, col.R = colR, col.grid.minor = "white") +
Tarrowlab("Probability of shared causal variant") +
Larrowlab("Probability of locus not shared") +
Rarrowlab("Probability of shared locus (distinct causal variant)")  +
theme_showarrows()  +
theme_notitles() +
theme(text=element_text(family="arial", size= 10),
        tern.axis.arrow.T = element_blank(),
        tern.axis.arrow.L = element_blank(),
        tern.axis.arrow.R = element_blank(),
        tern.axis.text.T = element_text(color = colT),
        tern.axis.text.L = element_text(color = colL),
        tern.axis.text.R = element_text(color = colR),
        tern.axis.arrow.text.T = element_text(color = colT),
        plot.margin = margin(0, 0, 0, 0, "cm"),
        tern.axis.arrow.text.L = element_text(color = colL),
        tern.axis.arrow.text.R = element_text(color = colR),
        tern.panel.grid.major = element_line(linetype = 6, size = 0.3))
