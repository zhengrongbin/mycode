#/bin/python
### human and mouse homology gene
## download from http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt from MGI
import os,sys
import numpy as np
import pandas as pd

hm_path = 'HOM_MouseHumanSequence.rpt'

hm = pd.read_csv(hm_path, sep = '\t')

## make to a new table
IdAll = list(set(hm['DB Class Key'].tolist()))
hm_pairs = []
for Id in IdAll:
    tmp = hm[hm['DB Class Key'] == Id]
    gene_mouse = tmp[tmp['Common Organism Name'] == 'mouse, laboratory']['Symbol'].values.tolist() # usually mouse gene is one
    gene_human = tmp[tmp['Common Organism Name'] == 'human']['Symbol'].values.tolist()
    if len(gene_human) == 0 or len(gene_mouse) == 0:
        continue
    for i in gene_mouse:
        for j in gene_human:
            hm_pairs.append([i, j])

hm_pairs = pd.DataFrame(hm_pairs, columns = ['mouse_gene', 'human_gene'])
# hm_pairs.columns = ['mouse_gene', 'human_gene']
## remove human NA
hm_pairs = hm_pairs[(~pd.isna(hm_pairs['human_gene'])) & (hm_pairs['human_gene'] != '')]
hm_pairs.to_csv('human_mouse_homology_gene_pair.csv', index = None)
