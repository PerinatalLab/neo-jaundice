library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)

# mfrfile = "/mnt/HARVEST/PDB1724_MFR_541_v12.csv"
# sentrixlinkfileF = "/mnt/HARVEST/linkage_Child_PDB1724.csv"
# genoflagfile = "/mnt/HARVEST/mobagen-flaglist-n99259.txt"
# pca1kgfile = "/mnt/HARVEST/pca/1kg_projected"
# pcamobafile = "/mnt/HARVEST/pca/mobagen_projected"
# snpfile = "/mnt/HARVEST/neoj/rs6755571.raw"

mfrfile = snakemake@input[[1]]
sentrixlinkfileF = snakemake@input[[2]]
genoflagfile = snakemake@input[[3]]
pca1kgfile = snakemake@input[[4]]
pcamobafile = snakemake@input[[5]]
snpfile = snakemake@input[[6]]
popfile = "../../resources/1kg_populations.txt"

# outresfile = "/mnt/HARVEST/neoj/mafs_across_1kgpops.csv"
outresfile = snakemake@output[[1]]
outplotfile = snakemake@output[[2]]

# ------------- SAMPLE PREP ---------------

## MFR
mfr = read.table(mfrfile, h=T, sep=";")
nrow(mfr) # 113k

# Singleton, live births
mfr = filter(mfr, is.na(DAAR) | DAAR!=FAAR,
             is.na(FLERFODSEL), is.na(MISD))
nrow(mfr)  # 104 k

mfr = mfr[, c("PREG_ID_1724", "ICTERUS")]
mfr$ICTERUS[is.na(mfr$ICTERUS)] = 0

## ------------------
# attaching genetics (no trio requirement)

# genotyping QC flags for the 30k data:
flags = read.table(genoflagfile, h=T)
nrow(flags)  # 99259

# bad/mixed up samples
flags = filter(flags, genotypesOK, phenoOK)
nrow(flags)  # 97679

# remove the case-control batch
flags = filter(flags, BATCH!="TED")
nrow(flags)

# sentrix-pregid converter, FETAL:
linkF = read.table(sentrixlinkfileF, h=T, sep=";")
linkF = linkF[,c("PREG_ID_1724", "SENTRIX_ID")]
linkF = filter(linkF, SENTRIX_ID %in% flags$IID)
linkF = filter(linkF, !duplicated(PREG_ID_1724))
nrow(linkF)

colnames(linkF)[ncol(linkF)] = "C_ID"

nrow(linkF) # 29321

mfr = inner_join(mfr, linkF[,c("PREG_ID_1724", "C_ID")])
mfr = mfr[,c("PREG_ID_1724", "C_ID", "ICTERUS")]
mfr$ICTERUS = mfr$ICTERUS + 1
anyDuplicated(mfr$PREG_ID_1724)
nrow(mfr)  # 28112


# ----------------- 
# 1kg pca

pca1kg = read.table(pca1kgfile, h=T)
pca1kgpops = read.table(popfile, sep="\t", h=T)
pca1kg = left_join(pca1kg, pca1kgpops[,c("Individual.ID", "Population")],
                   by=c("IID"="Individual.ID"))
pcamoba2 = read.table(pcamobafile, h=T)
pcamoba2$Population = "MoBa"

# cluster centers
popcentrs = pca1kg %>%
  gather(key="PC", value="val", PC1:PC5) %>%
  group_by(Population, PC) %>%
  summarize(centr = mean(val))

popnames = data.frame(matrix(c("African Caribbean in Barbados","ACB",
                    "African Ancestry in SW USA","ASW",
                    "Bengali in Bangladesh","BEB",
                    "British From England and Scotland","GBR",
                    "Chinese Dai in Xishuangbanna, China","CDX",
                    "Northern and Western European in Utah, USA","CEU",
                    "Colombian in Medellín, Colombia","CLM",
                    "Esan in Nigeria","ESN",
                    "Finnish in Finland","FIN",
                    "Gambian in Western Division – Mandinka","GWD",
                    "Gujarati Indians in Houston, Texas, USA","GIH",
                    "Han Chinese in Beijing, China","CHB",
                    "Han Chinese South","CHS",
                    "Iberian Populations in Spain","IBS",
                    "Indian Telugu in the UK","ITU",
                    "Japanese in Tokyo, Japan","JPT",
                    "Kinh in Ho Chi Minh City, Vietnam","KHV",
                    "Luhya in Webuye, Kenya","LWK",
                    "Mende in Sierra Leone","MSL",
                    "Mexican Ancestry in Los Angeles, USA","MXL",
                    "Peruvian in Lima Peru","PEL",
                    "Puerto Rican in Puerto Rico","PUR",
                    "Punjabi in Lahore, Pakistan","PJL",
                    "Sri Lankan Tamil in the UK","STU",
                    "Toscani in Italia","TSI",
                    "Yoruba in Ibadan, Nigeria","YRI"), ncol=2, byrow=T))
colnames(popnames)[1] = "Full.name"

# calc distance from pop centers
dists = pcamoba2[,1:5]
dists = dists %>%
  gather(key="PC", value="val", PC1:PC3) %>%
  inner_join(popcentrs, by="PC") %>%
  mutate(d=(val-centr)^2) %>%
  group_by(IID, Population) %>%
  summarize(d=sum(d))
dists = top_n(dists, 1, desc(d))
dists = dists[,c("IID", "Population")]
anyDuplicated(dists$IID)
table(dists$Population)

pca1kg %>%
  ggplot(aes(x=PC1, y=PC2, col=Population)) + 
  coord_fixed() +
  geom_point(size=1, alpha=0.5) +
  geom_point(data=spread(popcentrs, key="PC", value="centr"), pch=3, col="black") +
  geom_text(data=spread(popcentrs, key="PC", value="centr"), aes(label=Population), col="black", size=2)

pcamoba2 %>%
  inner_join(dists, by="IID") %>%
  ggplot(aes(x=PC1, y=PC2, col=Population.y)) + 
  coord_fixed() +
  scale_color_discrete(name="Inferred population") + 
  geom_point(size=1, alpha=0.5) +
  geom_point(data=spread(popcentrs, key="PC", value="centr"), pch=3, col="black") +
  geom_text_repel(data=spread(popcentrs, key="PC", value="centr"), aes(label=Population), col="black", size=2)
ggsave(outplotfile, width=6, height=6, units="cm")

# attach the MAFs and case status
snp = read.table(snpfile, h=T)

dists2 = dists %>%
  inner_join(snp[,c("IID", "rs6755571_A")], by="IID") %>%
  inner_join(mfr, by=c("IID"="C_ID"))
dists3 = filter(dists2, !is.na(ICTERUS)) %>%
  group_by(Population, ICTERUS) %>%
  summarize(n=n(), maf=mean(rs6755571_A)/2, .groups = "drop") %>%
  complete(Population, ICTERUS, fill=list(n=0, maf=NA))

dists3 = dists3 %>%
  mutate(ICTERUS=ifelse(ICTERUS==1, "contr", "case")) %>%
  left_join(popnames, by=c("Population"="X2")) %>%
  pivot_wider(names_from=ICTERUS, values_from=c(n, maf)) %>%
  arrange(desc(n_contr))

write.table(dists3, outresfile, sep="\t", col.names=T, row.names = F, quote = F)
