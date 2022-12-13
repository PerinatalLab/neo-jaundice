library(data.table)
library(dplyr)
library(tidyr)

mfr= fread(snakemake@input[[1]])
ids= fread(snakemake@input[[2]])


flag= fread(snakemake@input[[3]])
flag= filter(flag, genotypesOK== T, phenoOK== T)

out= readLines(snakemake@input[[4]])

ids= pivot_longer(ids, c('Child', 'Mother', 'Father'), names_to= 'ROLE', values_to= 'IID')
ids= ids[!duplicated(ids[,c('PREG_ID_1724', 'IID')]), ]

ids= filter(ids, (IID %in% out), IID %in% flag$IID)
ids= group_by(ids, PREG_ID_1724, ROLE) %>% filter(row_number()== 1)

ids= spread(ids, key= ROLE, value= IID)

mfr= inner_join(mfr, ids, by= 'PREG_ID_1724')

mfr= filter(mfr, is.na(FLERFODSEL)) # grepl('Levendefødt', DODKAT))
mfr= filter(mfr, is.na(DAAR) | DAAR != FAAR, is.na(MISD))
mfr$KJONN= with(mfr, ifelse(KJONN== 1, 1, ifelse(KJONN== 2, 0, NA)))

mfr$parity= with(mfr, ifelse(PARITET_5!= 0, 1, ifelse(PARITET_5== 0, 0, NA)))

mfr$jaundice= with(mfr, ifelse(!is.na(ICTERUS), 1, 0))

mfr= arrange(mfr, desc(jaundice))

moms= filter(mfr, !duplicated(PREG_ID_1724)) %>% select(Mother, jaundice, KJONN, parity, PREG_ID_1724, BATCH_moms) %>% filter(!is.na(Mother))
fets= filter(mfr, !duplicated(PREG_ID_1724)) %>% select(Child, jaundice, KJONN, parity, PREG_ID_1724, BATCH_fets) %>% filter(!is.na(Child))
dads= filter(mfr, !duplicated(PREG_ID_1724)) %>% select(Father, jaundice, KJONN, parity, PREG_ID_1724, BATCH_dads) %>% filter(!is.na(Father))

moms= moms[!duplicated(moms$Mother, incomparables= NA), ]
fets= fets[!duplicated(fets$Child, incomparables= NA), ]
dads= dads[!duplicated(dads$Father, incomparables= NA), ]

names(moms)[ncol(moms)]= 'cohort'
names(fets)[ncol(fets)]= 'cohort'
names(dads)[ncol(dads)]= 'cohort'

moms$cohort= ifelse(grepl('NORM', moms$cohort), 'NORMENT', moms$cohort)
dads$cohort= ifelse(grepl('NORM', dads$cohort), 'NORMENT', dads$cohort)
fets$cohort= ifelse(grepl('NORM', fets$cohort), 'NORMENT', fets$cohort)

names(moms)[1]= 'IID'
names(fets)[1]= 'IID'
names(dads)[1]= 'IID'

fwrite(moms, snakemake@output[[1]], sep= '\t')
fwrite(fets, snakemake@output[[2]], sep= '\t')
fwrite(dads, snakemake@output[[3]], sep= '\t')
