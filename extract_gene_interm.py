import os,sys
import json
import pandas as pd
import numpy as np


go = json.load(open('/Users/rongbinzheng/Documents/github/mycode/ChIPseq/ChIPdownstream/cistromego_py3/data/hg_kegg_symbol_dedup.json'))

rp = pd.read_csv('/Volumes/GoogleDrive/My Drive/Collaboration/Endocrine_BreastCancer/RawData/analysis/ELF3_TWF/allPeak_new/DESeq2/ELF3_Par_TWF20_TWF30_TWF20_vs_Par_UPpeak_go/TWF20_vs_Par_UPpeak_rp.txt',
    index_col = 0, sep = '\t')

enrich = pd.read_csv('/Volumes/GoogleDrive/My Drive/Collaboration/Endocrine_BreastCancer/RawData/analysis/ELF3_TWF/allPeak_new/DESeq2/ELF3_Par_TWF20_TWF30_TWF20_vs_Par_UPpeak_go/TWF20_vs_Par_UPpeak_kegg.txt',
    index_col = 0, sep = '\t')


terms = enrich[enrich['padj']<0.05]
res = {}
for Id in terms.index.tolist():
    name = terms.loc[Id, 'term']
    genes = go['ontology'][Id]
    res[Id+' : '+name] = rp[rp.index.isin(genes)]

res['hsa01100 : Metabolic pathways'].to_csv('hsa01100__Metabolic_pathways.tsv', sep = '\t')
res['hsa00062 : Fatty acid elongation'].to_csv('hsa00062__Fatty_acid_elongation.tsv', sep = '\t')
res['hsa05205 : Proteoglycans in cancer'].to_csv('hsa05205__Proteoglycans_in_cancer.tsv', sep = '\t')
res['hsa01212 : Fatty acid metabolism'].to_csv('hsa01212__Fatty_acid_metabolism.tsv', sep = '\t')


