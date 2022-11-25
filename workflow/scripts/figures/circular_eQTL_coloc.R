library(data.table)
library(dplyr)
library(broom)
library(ggplot2)
library(showtext)
library(cowplot)
library(ggrepel)

colorBlindBlack8= c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

d= fread(snakemake@input[[1]])

gene_names= data.frame(gene_id= c('ENSG00000167165','ENSG00000241119','ENSG00000241635',
                                  'ENSG00000242366','ENSG00000242515','ENSG00000244122','ENSG00000244474'), gene_symbol= c('UGT1A6', 'UGT1A9', 'UGT1A1', 'UGT1A8', 'UGT1A10', 'UGT1A7', 'UGT1A4'),

color_blind= colorBlindBlack8[2:8])

gene_names= filter(gene_names, gene_symbol == snakemake@wildcards[['UGT_genes']])

col_blind= unique(gene_names$color_blind)

d= inner_join(d, gene_names, by= c('gene'= 'gene_id'))

ugt1a1= d %>% arrange(PP.H4.abf)
ugt1a1= filter(ugt1a1, row_number() >= n() - 10)

ugt1a1$eqtl_data= gsub('.*ge_', '', ugt1a1$eqtl_data)
ugt1a1$eqtl_data= gsub('.*microarray_', '', ugt1a1$eqtl_data)
ugt1a1$eqtl_data= gsub('_', '\n', ugt1a1$eqtl_data)

ugt1a1$eqtl_data= factor(ugt1a1$eqtl_data, levels= unique(ugt1a1$eqtl_data))

p1= ggplot(ugt1a1) +
  # Make custom panel grid
  geom_hline(
    aes(yintercept = y), 
    data.frame(y = c(0:4) * 0.25),
    color = "lightgrey", size= 0.5
  ) + 
  # Add bars to represent the cumulative track lengths
  # str_wrap(region, 5) wraps the text so each line has at most 5 characters
  # (but it doesn't break long words!)
  geom_col(
    aes(
      x = eqtl_data,
      y = PP.H4.abf,
      fill = PP.H4.abf
    ),
    position = "dodge2",
    show.legend = F,
    alpha = .9
  ) +
  geom_segment(
    aes(
      x = eqtl_data,
      y = 0,
      xend = eqtl_data,
      yend = 1
    ),
    linetype = "dashed",
    color = "gray12"
  ) +
   # Make it circular!
  coord_polar(clip = 'off') +
    theme(
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    # Use gray text for the region names
    axis.text.x = element_text(color = "gray12", size = 7),
    # Move the legend to the bottom
    legend.position = "bottom",
    panel.background = element_rect(fill = "white", color = "white"),
    panel.grid = element_blank(),
    panel.grid.major.x = element_blank()) +
  scale_y_continuous(expand= c(0,0), limits = c(-1, 1), breaks = 0:4 * 0.25) +
  annotate(
    x = 0.75, 
    y = c(0:4 * 0.25), 
    label = c(0:4 * 0.25), 
    geom = "text", 
    color = "gray12",
    size= 2) +
  scale_fill_gradient2(low= colorBlindBlack8[1], high= col_blind) +
  annotate(
    x = 0, 
    y = -1, 
    label = snakemake@wildcards[['UGT_genes']], 
    geom = "text", 
    color = "gray12",
    fontface= 'italic', 
    size= 8/.pt)

save_plot(snakemake@output[[1]], p1, base_height= 60, base_width= 60, units= 'mm', dpi= 300)
