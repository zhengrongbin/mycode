import os
import sys
import django
import cPickle as p
  ## mkvirtualenv dc2, and pip install -r requirements.txt
sys.path.append('/data/home/qqin/01_Projects/Programming/dc2/lib/python2.7/site-packages')
sys.path.append('/data/home/qqin/01_Projects/Programming/dc2')
sys.path.append('/data/home/qqin/01_Projects/Programming/dc2/dc2')
sys.path.append('/data/home/qqin/01_Projects/Programming/dc2/datacollection')
sys.path = sys.path[::-1]
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "dc2.settings")
django.setup()
from datacollection import models
import json as json
import re
import os,sys
import pandas as pd
import numpy as np
# import h5py
import math
# import pybedtools as pbed
import time

types = sys.argv[1] ## type of function, interval or similar
specie = sys.argv[2] ## species: hg38 or mm10
saveFile = sys.argv[3] # the file name for saving result

if types == 'interval':
    interval = sys.argv[4] # chrX:888888-9999999
    ntmp = saveFile+'.giggle'
    os.system("/data5/home/rongbin/software/giggle/bin/giggle search -i /data6/home/rongbin2/changxin/giggle_index/%stf_gindex_all -r %s > %s"%(specie, interval, ntmp))
    result_df = pd.read_csv(ntmp , sep="\t", index_col=False, header=None)
    result_df.index = [i.replace("_peaks.bed.gz", "").split(".")[-1] for i in result_df[0]]
    result_df[2] = [int(i.replace("overlaps:", "")) for i in result_df[2]]
    result_df = result_df[result_df[2]>0]
    result_df = result_df.sort_values(by=2, ascending=False)
    # npage = math.ceil(result_df.shape[0]/20.0)
    items = []
    for i in result_df.index:
        s = models.Samples.objects.filter(id=str(i))
        treat = s.values("species__name",
             "factor__name",
             "cell_line__name",
             "cell_type__name",
             "tissue_type__name",
             "disease_state__name",
             "unique_id")[0]
        tmp = {"sid":str(i),"specie":treat['species__name'], "uid": treat['unique_id'],
        "factor": treat['factor__name'], "tissue": treat['tissue_type__name'],
        "celltype": treat['cell_type__name'], "cellline": treat['cell_line__name'],
        "overlap": str(result_df.loc[i,2])}
        items.append(tmp)
    items = pd.DataFrame(items)
    items.to_csv(saveFile)

if types == 'similar': 
    tpeak = '1k'
    peak_path = sys.argv[4] # the path of peak file
    # sfile = open(peak_path, 'w')
    # sfile.write(peak_file.read())
    # sfile.close()
    os.system("/data5/home/rongbin/software/htslib-1.7/bgzip -c %s > %s.gz"%(peak_path, peak_path))
    os.system("/data5/home/rongbin/software/giggle/bin/giggle search -i /data6/home/rongbin2/changxin/giggle_index/%stf_gindex_%s/ -q %s.gz -s > %s.giggle"%(specie ,tpeak, peak_path, saveFile))
    result_df = pd.read_csv("%s.giggle"%saveFile, sep="\t", index_col=False)
    result_df.index = [i.replace("_peaks.bed.gz", "").split(".")[-1] for i in result_df["#file"]]
    result_df = result_df.sort_values(by="combo_score", ascending=False)#.head(num_neighbor)
    npage = math.ceil(result_df.shape[0]/20.0)
    items = []
    for i in result_df.index:
        s = models.Samples.objects.filter(id=str(i))
        treat = s.values("species__name",
                     "factor__name",
                     "cell_line__name",
                     "cell_type__name",
                     "tissue_type__name",
                     "disease_state__name",
                     "unique_id")[0]
        tmp = {"sid":str(i),"specie":treat['species__name'], "uid": treat['unique_id'],
                "factor": treat['factor__name'], "tissue": treat['tissue_type__name'],
                "celltype": treat['cell_type__name'], "cellline": treat['cell_line__name'],
                "distance": str(result_df.loc[i, 'combo_score'])}
        items.append(tmp)
    items = pd.DataFrame(items)
    items.to_csv(saveFile)








