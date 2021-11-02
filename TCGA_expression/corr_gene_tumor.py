"""
#!/usr/bin/env python
#author          :Rongbin Zheng
#date            :2020-08-29
#python_version  :3.6 
"""

import pingouin as pg
import scipy.cluster.hierarchy as sch
import os
import sys
import math
import re
import traceback
import pandas as pd
import numpy as np
import collections
import pickle as pk
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
from scipy.stats.stats import pearsonr, spearmanr
from scipy.stats import wilcoxon, ranksums
from scipy.stats import stats
import common_function as cf
import pingouin as pg


plt.rcParams.update(plt.rcParamsDefault)
sns.set(style='ticks')
rc = {'axes.labelpad': 15, "axes.labelsize": 14,
      "figure.titleweight": "bold"}
plt.rcParams.update(**rc)


def info(mesg):
    # output message
    os.system('echo "++++ %s"' % mesg)


def partial_corr(C, method='pearson'):
    p = C.shape[1]
    index = C.columns.tolist()
    P_corr = np.zeros((p, p), dtype=np.float)
    P_pval = np.zeros((p, p), dtype=np.float)
    for i in range(p):
        P_corr[i, i] = 1
        P_pval[i, i] = 1
        for j in range(i+1, p):
            x, y = index[i], index[j]
            res = pg.partial_corr(data=C, x=x, y=y,
                                  covar=list(set(index) - set([x, y])),
                                  method=method)
            corr, pval = res['r'][0], res['p-val'][0]
            P_corr[i, j] = corr
            P_corr[j, i] = corr
            P_pval[i, j] = pval
            P_pval[j, i] = pval
    P_corr = pd.DataFrame(P_corr, index=index, columns=index)
    P_pval = pd.DataFrame(P_pval, index=index, columns=index)
    return (P_corr, P_pval)


def _scatter(exp_mat, gene1, gene2, regular_corr_mat, partial_corr_mat, pdf=False):
    """plot scatter plot to show gene correlation"""
    for cancer in regular_corr_mat.index.tolist():
        exp_mat_tmp = exp_mat[exp_mat.project_id == cancer]
        fig, ax = plt.subplots(figsize=(4, 4))
        sns.regplot(x=gene1, y=gene2, data=exp_mat_tmp, ax=ax,
                    scatter_kws={'s': 2, 'color': "navy", 'alpha': .8},
                    line_kws={"color": 'grey'})
        if partial_corr_mat.shape[0] == 0:
            rcorr, rpval = regular_corr_mat.loc[cancer,
                                                "corr"], regular_corr_mat.loc[cancer, "pval"]
            ax.set_title(cancer+'\nr = %s, p = %.2e' % (round(rcorr, 4), rpval),
                         fontdict={'fontweight': 'bold'})
        else:
            pcorr, ppval = partial_corr_mat.loc[cancer,
                                                "corr"], partial_corr_mat.loc[cancer, "pval"]
            ax.set_title(cancer+'\nPartial r = %s, p = %.2e' %
                         (round(pcorr, 4), ppval))
        ax.set(xlabel='%s expression level' %
               gene1, ylabel='%s expression level' % gene2)
        plt.tight_layout()
        pdf.savefig(fig) if pdf else None
        plt.close()


def _calculate_corr(gene1, gene2, adj_gene, exp_mat):
    """calculate correlation between genes or partial correlation"""
    # info('correlation analysis')
    regular_corr_mat = exp_mat.groupby('project_id')[[gene1, gene2]].apply(lambda df: pd.Series(pearsonr(df[gene1], df[gene2]),
                                                                                                index=['corr', 'pval']))
    if adj_gene:
        # calculate partial correlation by cancer types
        partial_corr_mat = exp_mat.groupby('project_id')[[gene1, gene2, adj_gene]].apply(
            lambda df: pd.Series([x.loc[gene1, gene2] for x in partial_corr(df)],
                                 index=['corr', 'pval']))
    else:
        partial_corr_mat = pd.DataFrame()
    # return two column dataframe which includes corr coefficient and p-value
    return(regular_corr_mat, partial_corr_mat)


def _corr_tumor(config):
    """
    correlation analysis between given two genes  
    """
    gene1, gene2, adj_gene = config['gene1'], config['gene2'], config['adj_gene']
    if gene1 == gene2:
        info('gene1 and gene2 are same, %s and %s' % (gene1, gene2))
        sys.exit(1)
    cancers, exp_mat = cf._test_data(config)
    # ignore zero value, since zero value will influence correlation analysis
    smallest = min([isinstance(x, np.float64)
                    for x in exp_mat.iloc[0, :].tolist()])
    info('smallest value %s' % smallest)
    filter_zero = (exp_mat[gene1] != smallest) & (exp_mat[gene2] != smallest)

    exp_mat = exp_mat.loc[filter_zero &
                          (exp_mat.sample_type == 'Tumor'), :]
    ct = exp_mat.groupby('project_id')['project_id'].count()
    # cancer type greater than 10
    info('too few non-zero sample: %s' %
         str(ct[ct < 10].index.tolist())) if ct[ct < 10].shape[0] > 0 else None
    ct = ct[ct >= 10]
    exp_mat = exp_mat[exp_mat.project_id.isin(ct.index)]
    # calculating
    regular_corr_mat, partial_corr_mat = _calculate_corr(
        gene1, gene2, adj_gene, exp_mat)
    info('output correlation result')
    regular_corr_mat.to_csv(
        config['prefix'] + '_%s_%s_pearson_correlation.csv' % (gene1, gene2))
    partial_corr_mat.to_csv(
        config['prefix'] + '_%s_%s_pearson_partial_correlation.csv' % (gene1, gene2)) if partial_corr_mat.shape[0] != 0 else None

    # plotting scatter
    pdf = PdfPages(config['prefix'] + '_%s_%s_scatter.pdf' % (gene1, gene2))
    _scatter(exp_mat, gene1, gene2, regular_corr_mat,
             partial_corr_mat, pdf)
    pdf.close()
