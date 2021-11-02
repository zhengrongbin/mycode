"""
#!/usr/bin/env python
#author          :Rongbin Zheng
#date            :2020-08-29
#python_version  :3.6 
"""

import scipy.cluster.hierarchy as sch
from scipy import linalg
from lifelines.statistics import logrank_test
from lifelines import KaplanMeierFitter, CoxPHFitter
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

plt.rcParams.update(plt.rcParamsDefault)
sns.set(style='ticks')
rc = {'axes.labelpad': 15, "axes.labelsize": 14,
      "figure.titleweight": "bold"}
plt.rcParams.update(**rc)


def info(mesg):
    # output message
    os.system('echo "++++ %s"' % mesg)


def _figsize_est(plot_df):
    """
    estimate boxplot fig size
    """
    width = 2.5 + (plot_df['project_id'].unique().shape[0] * 0.15)
    height = 4.5
    return(width, height)


def _set_star(x):
    if not x:
        return(None)
    if x < 0.001:
        star = '***'
    elif x < 0.01:
        star = '**'
    elif x < 0.05:
        star = '*'
    else:
        return (None)
    return(star)


def _get_pval(plot_df, gene):
    pval = plot_df.groupby('project_id').apply(lambda df: ranksums(df[df.sample_type == 'Tumor'][gene],
                                                                   df[df.sample_type == 'Normal'][gene])[1] if df[df.sample_type == 'Normal'].shape[0] > 0 else None)
    pval = pd.DataFrame(pval, columns=['pval'])
    pval['star'] = [_set_star(x) for x in pval['pval'].tolist()]
    return(pval)


def _get_delta_median(plot_df, gene):
    # take median value, and compare tumor vs normal
    med = plot_df.groupby(['sample_type', 'project_id'])[[gene]].median()
    med.reset_index(inplace=True)
    med = med.pivot_table(columns=['sample_type'], index=[
                          'project_id']).dropna()
    med.columns = ['Normal', 'Tumor']
    delta_med = med['Tumor'] - med['Normal']
    delta_med = delta_med.sort_values(
        ascending=False)  # order from high to low
    return (delta_med)


def _boxplot(plot_df, gene, pdf=False):
    """
    start to plot
    """
    figsize = _figsize_est(plot_df)  # estimate figsize
    delta_med = _get_delta_median(plot_df, gene)  # order from high to low
    pval = _get_pval(plot_df, gene)
    fig, ax = plt.subplots(figsize=figsize)
    bp = sns.boxplot(data=plot_df, x='project_id', y=gene, hue='sample_type',
                     ax=ax, fliersize=1, width=.7, palette={'Normal': 'grey',
                                                            'Tumor': 'red'},
                     boxprops={'edgecolor': 'black', 'linewidth': .7},
                     medianprops={'color': 'black', 'linewidth': 1},
                     whiskerprops={'linestyle': 'dashed', "linewidth": .7},
                     capprops={'linewidth': .7},
                     order=delta_med.index)
    # sns.despine(ax=ax, offset=10, trim=True)
    ax.set(xlabel='', ylabel='%s expression level' % gene)
    ax.tick_params(axis='x', which='major',  # labelsize=10,
                   pad=9, rotation=90)
    ax.legend(title='', fontsize=8, loc='lower left')
    # get lim
    max_y = plot_df[gene].max()
    min_y = plot_df[gene].min()
    extend = (max_y-min_y)/20
    ax.set_ylim(min_y-(extend/2), max_y+(extend*2))
    # color xlabel by cancer type
    n = 0
    for l in bp.get_xticklabels():
        c = l.get_text()
        if c in delta_med.index.tolist() and delta_med[c] > 0:
            l.set_color('darkred')
        elif c in delta_med.index.tolist() and delta_med[c] < 0:
            l.set_color('darkblue')
        else:
            pass
        star = pval.loc[c, 'star']
        ax.text(n - 0.15, max_y + extend,
                pval.loc[c, 'star'], fontsize=8, horizontalalignment='left') if star else None
        n += 1
    plt.tight_layout()
    pdf.savefig(fig) if pdf else None
    plt.close()
    return(pval)


def _metatstic_boxplot(config):
    # do boxplot for metatstic
    info('plot metastastic boxplot')
    gene_name = config['gene']
    cancers, exp_mat = cf._test_data(config)
    exp_mat = exp_mat[exp_mat.sampleType.isin(
        ['Primary Tumor', 'Metastatic', 'Solid Tissue Normal'])]
    exp_mat['sample_type'] = exp_mat['sampleType'].tolist()
    pdf_name = config['prefix'] + \
        '_%s_boxplot_Tumor_Normal_Metastastic.pdf' % gene_name
    pdf = PdfPages(pdf_name)

    figsize = _figsize_est(exp_mat)
    fig, ax = plt.subplots(figsize=figsize)
    bp = sns.boxplot(data=exp_mat, x='project_id', y=gene_name, hue='sample_type',
                     ax=ax, fliersize=1, width=.7, palette={'Primary Tumor': 'red',
                                                            'Metastatic': 'blue',
                                                            'Solid Tissue Normal': 'grey'},
                     boxprops={'edgecolor': 'black', 'linewidth': 1},
                     medianprops={'color': 'black', 'linewidth': 1.2})
    # sns.despine(ax=ax, offset=10, trim=True)
    ax.set(xlabel='', ylabel='%s expression level' % gene_name)
    ax.tick_params(axis='x', which='major',  # labelsize=10,
                   pad=9, rotation=90)
    ax.legend(title='', fontsize=8, loc='lower left')
    plt.tight_layout()
    pdf.savefig(fig) if pdf else None
    plt.close()
    pdf.close()


def _tumor_normal(config):
    """
    draw boxplot for given gene in tumor vs normal
    """
    gene_name = config['gene']
    cancers, exp_mat = cf._test_data(config)
    # pdf
    pdf_name = config['prefix']+'_%s_boxplot_Tumor_vs_Normal.pdf' % gene_name
    pdf = PdfPages(pdf_name)
    info('start to do boxplot between Tumor vs. Normal')
    # return p value between tumor and normal
    pval = _boxplot(exp_mat, gene_name, pdf)
    pdf.close()
    # write out pval
    pval.to_csv(config['prefix'] +
                '_%s_boxplot_Tumor_vs_Normal_pval.csv' % gene_name)


def _survial_analysis(gene_name, mat, title, pert=25, pdf=False):
    """
    do survival analysis based on the gene expression
    mat is the matrix, gene_name is one of the columns in matrix
    """
    top_cutoff, bottom_cutoff = np.percentile(
        mat[gene_name], 100 - pert), np.percentile(mat[gene_name], pert)
    top, bottom = mat[mat[gene_name] >
                      top_cutoff], mat[mat[gene_name] < bottom_cutoff]
    if bottom.shape[0] == 0 or top.shape[0] == 0:
        return
    T = top['OS.time'].astype('int')
    E = top['OS'].astype('int')
    T1 = bottom['OS.time'].astype('int')
    E1 = bottom['OS'].astype('int')
    try:
        p_logrank = logrank_test(
            T, T1, event_observed_A=E, event_observed_B=E1)
        p = p_logrank.p_value
        cph = CoxPHFitter()
        cph.fit(mat[[gene_name, 'OS.time', 'OS']], 'OS.time', event_col='OS')
        cox_p, cox_z, cox_hr = cph.summary['p'][0], cph.summary['z'][0], cph.hazard_ratios_[
            0]
    except:
        info('problem at fiting survival model')
        traceback.print_exc(sys.stderr)
        sys.exit(1)
    kmf = KaplanMeierFitter()
    fig, ax = plt.subplots(figsize=(5, 4.5))
    ax = kmf.fit(T, E).plot(ax=ax, color='red', label="%s High 25%% n=%s" % (gene_name, top.shape[0]),
                            ci_show=False, linewidth=1)
    ax = kmf.fit(T1, E1, label="%s Low 25%% n=%s" % (gene_name, bottom.shape[0])).plot(ax=ax,
                                                                                       color='blue', ci_show=False, linewidth=1)
    ax.set(xlabel='Days', ylabel='Overall Survival (%)',
           title=title)
    ax.set_ylim(-0.05, 1.05)
    ax.text(50, 0.3, "p value:", fontsize=11)
    ax.text(50, 0.24, "CoxPH = %s" % round(cox_p, 4), fontsize=11)
    ax.text(50, 0.17, "log rank = %s" % round(p, 4), fontsize=11)
    ax.text(50, 0.1, "log2 HR = %s" % round(np.log2(cox_hr), 4), fontsize=11)
    ax.legend()
    plt.tight_layout()
    pdf.savefig(fig) if pdf else None
    plt.close()


def _run_survival(config):
    """
    prepare survival analysis
    """
    gene_name = config['gene']
    cancers, df = cf._test_data(config)
    pdf = PdfPages(config['prefix']+'_%s_survial.pdf' % gene_name)
    for cancer in cancers:
        mat = df[df.project_id == cancer]
        mat = mat[mat.sample_type ==
                  'Tumor']  # only tumor samples will be used for survival analysis
        info('start to run survival analysis of %s' % cancer)
        # print(mat.columns)
        mat = cf._unique_patient(mat, gene_name)  # unique patient
        _survial_analysis(gene_name, mat, title=cancer, pdf=pdf)
    pdf.close()


def _dist_plot(gene_name, mat, title, pdf):
    """
    plot distribution of gene in tumor and normal sample
    """
    mat_tumor = mat[mat.sample_type == 'Tumor']
    mat_normal = mat[mat.sample_type == 'Normal']
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.hist(mat_tumor[gene_name], density=True, label='Tumor', alpha=.8)
    if mat_normal.shape[0] > 0:
        ax.hist(mat_normal[gene_name], density=True,
                label='Normal', alpha=.8)
    ax.legend()
    sns.despine(ax=ax, offset=5, trim=True)
    ax.set(xlabel='%s expression level' % gene_name, ylabel='Frequency')
    ax.set_title(label=title,
                 fontdict={'fontweight': 'bold', 'fontsize': 14})
    plt.tight_layout()
    pdf.savefig(fig) if pdf else None
    plt.close()


def _distribution(config):
    """
    distribution of the gene
    """
    gene_name = config['gene']
    cancers, df = cf._test_data(config)
    pdf = PdfPages(config['prefix']+'_%s_distribution.pdf' % gene_name)
    for cancer in cancers:
        mat = df[df.project_id == cancer]
        info('start to plot distribution of %s' % cancer)
        _dist_plot(gene_name, mat, cancer, pdf)
    pdf.close()
