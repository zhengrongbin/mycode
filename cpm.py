import os,sys
import pandas as pd
import numpy as np


exp_path = sys.argv[1]

dat = pd.read_csv(exp_path, sep = '\t', index_col = 0)
colsums = dat.sum() / 1e06
cpm = dat.apply(lambda row: row / colsums, axis = 1)
cpm.to_csv(exp_path+'.CPM.tsv', sep = '\t')
