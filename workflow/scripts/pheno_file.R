library(data.table)
library(dplyr)
library(tidyr)

mfr= fread(snakemake@input[[1]])
ids= fread(snakemake@input[[2]])

names(ids)= c('PREG_ID', 'IID', 'BATCH', 'Role')

ids= filter(ids, BATCH != 'TED')

ids= ids[!duplicated(ids[,c('PREG_ID', 'IID')]), ]

flag= fread(snakemake@input[[3]])
flag= filter(flag, genotypesOK== T, phenoOK== T)

out= readLines(snakemake@input[[4]])


ids= filter(ids, !(IID %in% out), IID %in% flag$IID)

ids= spread(ids, key= Role, value= IID)

mfr= inner_join(mfr, ids, by= c('PREG_ID_315'= 'PREG_ID'))

mfr= filter(mfr, FLERFODSEL=='Enkeltfødsel', grepl('Levendefødt', DODKAT))

mfr$KJONN= with(mfr, ifelse(KJONN== 'Pike', 1, ifelse(is.na(KJONN), NA, 0)))

mfr$parity= with(mfr, ifelse(PARITET_5!= '0 (førstegangsfødende)', 1, ifelse(PARITET_5== '0 (førstegangsfødende)', 0, NA)))

mfr$jaundice= with(mfr, ifelse(ICTERUS== 'Ja', 1, ifelse(ICTERUS== 'Nei', 0, NA)))
mfr= arrange(mfr, desc(jaundice))

mfr$cohort= mfr$BATCH

moms= mfr[!duplicated(mfr$Mother, incomparables= NA), ]
fets= mfr[!duplicated(mfr$Child, incomparables= NA), ]
dads= mfr[!duplicated(mfr$Father, incomparables= NA), ]

moms= select(moms, Mother, jaundice, cohort, KJONN, parity) %>% filter(!is.na(Mother))

fets= select(fets, Child, jaundice, cohort, KJONN, parity) %>% filter(!is.na(Child))
dads= select(dads, Father, jaundice, cohort, KJONN, parity) %>% filter(!is.na(Father))

names(moms)[1]= 'IID'
names(fets)[1]= 'IID'
names(dads)[1]= 'IID'

fwrite(moms, snakemake@output[[1]], sep= '\t')
fwrite(fets, snakemake@output[[2]], sep= '\t')
fwrite(dads, snakemake@output[[3]], sep= '\t')
