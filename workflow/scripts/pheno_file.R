library(data.table)
library(dplyr)
library(tidyr)

mfr= fread(snakemake@input[[1]])
ids= fread(snakemake@input[[2]], header=F)
names(ids)= c('IID', 'BATCH', 'PREG_ID', 'ROLE')
ids= filter(ids, BATCH != 'TED')

ids= ids[!duplicated(ids[,c('PREG_ID', 'IID')]), ]

flag= fread(snakemake@input[[3]])
flag= filter(flag, genotypesOK== T, phenoOK== T)

out= readLines(snakemake@input[[4]])


ids= filter(ids, !(IID %in% out), IID %in% flag$IID)

ids= group_by(ids, PREG_ID, ROLE) %>% filter(row_number()== 1)

ids= spread(ids, key= ROLE, value= IID)

mfr= inner_join(mfr, ids, by= c('PREG_ID_1724'= 'PREG_ID'))

mfr= filter(mfr, is.na(FLERFODSEL)) # grepl('Levendefødt', DODKAT))
mfr= filter(mfr, is.na(DAAR) | DAAR != FAAR, is.na(MISD))
mfr$KJONN= with(mfr, ifelse(KJONN== 1, 1, ifelse(KJONN== 2, 0, NA)))

mfr$parity= with(mfr, ifelse(PARITET_5!= 0, 1, ifelse(PARITET_5== 0, 0, NA)))

mfr$jaundice= with(mfr, ifelse(!is.na(ICTERUS), 1, 0))

mfr= arrange(mfr, desc(jaundice))

mfr$cohort= mfr$BATCH

moms= mfr[!duplicated(mfr$Mother, incomparables= NA), ]
fets= mfr[!duplicated(mfr$Child, incomparables= NA), ]
dads= mfr[!duplicated(mfr$Father, incomparables= NA), ]

moms= filter(moms, !duplicated(PREG_ID_1724))
dads= filter(dads, !duplicated(PREG_ID_1724))
fets= filter(fets, !duplicated(PREG_ID_1724))

moms= select(moms, Mother, jaundice, cohort, KJONN, parity) %>% filter(!is.na(Mother))

fets= select(fets, Child, jaundice, cohort, KJONN, parity) %>% filter(!is.na(Child))
dads= select(dads, Father, jaundice, cohort, KJONN, parity) %>% filter(!is.na(Father))

names(moms)[1]= 'IID'
names(fets)[1]= 'IID'
names(dads)[1]= 'IID'

fwrite(moms, snakemake@output[[1]], sep= '\t')
fwrite(fets, snakemake@output[[2]], sep= '\t')
fwrite(dads, snakemake@output[[3]], sep= '\t')
