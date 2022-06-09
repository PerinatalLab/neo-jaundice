# Script to determine blood group based on two missense variants.

library(data.table)
library(dplyr)
library(tidyr)

format_haps= function(hap){
variants= paste(hap$chr, hap$pos, hap$ref, hap$eff, sep =':')
ids= names(hap)[5:ncol(hap)]
hap= as.data.frame(t(hap[, 5:ncol(hap)]))
names(hap)= variants
hap$PREG_ID= ids
return(hap)
}



d = fread(snakemake@input[[1]])

d= format_haps(d)

variants= names(d[, 1:2])

names(d)= c('V1', 'V2', 'IID')

d$V1= gsub('\\/', '\\|', d$V1)
d$V2= gsub('\\/', '\\|', d$V2)

d$V1= with(d, ifelse(V1== '0|0', 0, ifelse(V1== '1|0' | V1== '0|1', 1, ifelse(V1== '1|1', 2, NA))))
d$V2= with(d, ifelse(V2== '0|0', 0, ifelse(V2== '1|0' | V2== '0|1', 1, ifelse(V2== '1|1', 2, NA))))


d= mutate(d, ABO = ifelse(V2== 0, "O", NA)) %>% mutate(ABO= ifelse((V1== 2 & V2 == 2) | (V1== 1 & V2== 1), 'B', ABO)) %>% mutate(ABO= ifelse(V1== 0 & V2> 0, 'A', ABO)) %>% mutate(ABO= ifelse((V1== 1 & V2== 2), 'AB', ABO))

d= select(d, IID, ABO)

fwrite(d, file = snakemake@output[[1]], sep="\t")

