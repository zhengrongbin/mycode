import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
import pandas as pd
import os, sys, subprocess

import cooler
import cooltools.lib.plotting
from packaging import version
if version.parse(cooltools.__version__) < version.parse('0.5.2'):
    raise AssertionError("tutorials rely on cooltools version 0.5.2 or higher,"+
                         "please check your cooltools version and update to the latest")
import cooltools
import bioframe
import pickle as pk

cool_path = sys.argv[1]
fa_path = sys.argv[2]
prefix = sys.argv[3]

clr = cooler.Cooler('%s::resolutions/100000'%cool_path)

if not os.path.exists('gc_cov_100kb.tsv'):
    bins = clr.bins()[:]
    genome = bioframe.load_fasta(fa_path);
    ## note the next command may require installing pysam
    gc_cov = bioframe.frac_gc(bins[['chrom', 'start', 'end']], genome)
    gc_cov.to_csv('gc_cov_100kb.tsv',index=False,sep='\t')

view_df = pd.DataFrame({'chrom': clr.chromnames,
                        'start': 0,
                        'end': clr.chromsizes.values,
                        'name': clr.chromnames}
                      )
cis_eigs = cooltools.eigs_cis(
                        clr,
                        gc_cov,
                        view_df=view_df,
                        n_eigs=3,
                        )
out = open(prefix+'.cis_eigs.pk', 'wb')
pk.dump(cis_eigs, out)
out.close()


# # cis_eigs[0] returns eigenvalues, here we focus on eigenvectors
# eigenvector_track = cis_eigs[1][['chrom','start','end','E1']]

