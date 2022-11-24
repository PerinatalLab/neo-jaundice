import pandas as pd
import numpy as np
import csv
import gzip
from functools import reduce
from scipy import stats as st
import time

##### load config and sample sheets #####

configfile: "config/config.yaml"

fam_ids = pd.read_csv(config["fam_member"], sep= "\t").set_index("fam_id", drop=False)
pheno_file= pd.read_csv(config["pheno"], sep= '\t')
CHROM= [i.strip() for i in open(config["CHROM"], 'r')]
autosomal_CHR= [i for i in CHROM if i!= 'X']
haplotypes= [i.strip() for i in open(config["haplotypes"], 'r')]
eQTL_catalogue_conditions= [i.strip() for i in open(config["eQTL_catalogue_conditions"], 'r')]

def selectUnrelated(input_kin, df, x):
        kin= pd.read_csv(input_kin, header= 0, sep= '\t')
        kin= kin.loc[kin.Kinship > 0.125, :]
        kin= kin.loc[kin.ID1.isin(x.values)]
        kin= kin.loc[kin.ID2.isin(x.values)]
        kin= kin.loc[:, ['ID1','ID2','Kinship']]
        kin_temp= kin.copy()
        kin_temp.columns= ['ID2', 'ID1', 'Kinship']
        kin_temp= kin_temp.append(kin)
        kin_temp['n']= kin_temp.groupby('ID1')['ID1'].transform('count')
        kin_temp['nn']= kin_temp.groupby('ID2')['ID2'].transform('count')
        kin_temp.sort_values(by=['n', 'nn'], inplace= True)
        to_keep= list()
        for i in range(0, len(kin_temp.index)):
                if kin_temp.iloc[i, 0] in kin_temp.iloc[0:i, 1].values:
                        kin_temp.iloc[i, 1]= "X"
                else:
                        to_keep.append(kin_temp.iloc[i, 0])
        to_remove= [i for i in kin_temp.ID1 if i not in to_keep]
        to_remove= list(set(to_remove))
        remove= pd.DataFrame({'FID': to_remove})
        remove['IID']= remove.FID
        return remove

