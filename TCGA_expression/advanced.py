"""
#!/usr/bin/env python
#author          :Rongbin Zheng
#date            :2020-08-29
#python_version  :3.6 
"""

import os
import sys
import pandas as pd
import numpy as np
from scipy.stats import norm

import common_function as cf
from corr_gene_tumor import _calculate_corr


def info(mesg):
    # output message
    os.system('echo "++++ %s"' % mesg)


def _get_mat(gene, g, adj_gene, exp_mat, smallest):
    """get exp mat by getting only the needed columns, and remove zero samples"""
    if adj_gene:
        mat = exp_mat[[gene, g, adj_gene, 'project_id']]
    else:
        mat = exp_mat[[gene, g, 'project_id']]
    ## the smallest value can be 0 or log2(x), since the too many zero result in faked correlation coefficient
    filter_zero = (mat[gene] != smallest) & (mat[g] != smallest)
    mat = mat.loc[filter_zero, :]
    # count for each cancer type
    ct = mat.groupby('project_id')['project_id'].count()
    info('too few non-zero sample: %s, for gene %s' %
         (str(ct[ct < 10].index.tolist()), gene+'~'+g)) if ct[ct < 10].shape[0] > 0 else None
    # cancer type greater than 10
    ct = ct[ct >= 10]
    mat = mat[mat.project_id.isin(ct.index)]
    return (mat)


def _corr_all_gene(gene, protein_genes, adj_gene, exp_mat):
    """
    calculate regular correlation and partial correlation for given gene and protein genes
    """
    info('input matrix size: %s' % str(exp_mat.shape))
    smallest = min([x for x in exp_mat.iloc[0, :].tolist()
                    if isinstance(x, np.float64)])
    info('smallest value %s' % smallest)
    info('starting calculate correlation between given gene and protein coding genes')
    protein_genes = sorted(
        list(set(protein_genes) - set([gene])))  # avoid given gene
    regular_corr_all, regular_pval_all = pd.DataFrame(), pd.DataFrame()
    partial_corr_all, partial_pval_all = pd.DataFrame(), pd.DataFrame()
    # record sample count for calculating delta fisher z score
    sample_count = {}
    one_percent = len(protein_genes) / 100
    cnt = 0
    idx = []  # the gene name for indexing
    # np.random.choice(protein_genes, 500, replace=False).tolist():
    for g in protein_genes:
        cnt += 1
        if cnt % one_percent == 0:
            info('%s%%' % cnt / one_percent)
        # get matrix for needed genes
        mat = _get_mat(gene, g, adj_gene, exp_mat, smallest)
        if mat.shape[0] == 0:
            continue  # if no cancer type sample size greater than 10, skip this gene
        sample_count[g] = mat.groupby('project_id')['project_id'].count()
        # corr: return two column dataframe which includes corr coefficient and p-value, if not adj_gene, partial_corr will be empty dataframe
        regular_corr, partial_corr = _calculate_corr(gene, g, adj_gene, mat)
        # merge dataframe
        regular_corr_all = pd.concat(
            [regular_corr_all, regular_corr[['corr']].T])  # row = gene, col = cancer type
        regular_pval_all = pd.concat(
            [regular_pval_all, regular_corr[['pval']].T])  # row = gene, col = cancer type
        if adj_gene:
            partial_corr_all = pd.concat(
                [partial_corr_all, partial_corr[['corr']].T])  # row = gene, col = cancer type
            partial_pval_all = pd.concat(
                [partial_pval_all, partial_corr[['pval']].T])  # row = gene, col = cancer type
        idx.append(g)
    # add row names
    regular_corr_all.index = idx
    regular_pval_all.index = idx
    if adj_gene:
        partial_corr_all.index = idx
        partial_pval_all.index = idx
    return ({'regular_corr': regular_corr_all,
             "regular_pval": regular_pval_all,
             "partial_corr": partial_corr_all,
             "partial_pval": partial_pval_all,
             'sample_count': pd.DataFrame(sample_count).T})  # row = gene, col = cancer type


def _fisherz(r):
    # fisher z transfermation
    z = 0.5 * np.log((1+r)/(1-r))
    return(z)


def _delta_fisherz(z1, z2, n1, n2):
    dz = (z1 - z2) / np.sqrt((1/(n1-3)) + (1/(n2-3)))
    pval = 2 * (1 - norm.cdf(np.abs(dz)))
    pval = pd.Series(pval, index=dz.index)
    return(dz, pval)


def _run_delta_between_tumor_normal(tumor_corr, normal_corr, tumor_sample_count, normal_sample_count):
    """
    tumor_corr and normal_corr should be dataframe that rows = genes, and cols = cancers
    sample_count are the same
    """
    common_genes = np.intersect1d(tumor_corr.index, normal_corr.index)
    common_cancers = np.intersect1d(tumor_corr.columns, normal_corr.columns)
    if common_genes.shape[0] == 0 or common_cancers.shape[0] == 0:
        info(
            'no common genes between tumor corr and normal corr') if common_genes.shape[0] == 0 else None
        info(
            'no common cancer types between tumor corr and normal corr') if common_cancers.shape[0] == 0 else None
        return (None, None)
    # same order
    tumor_sample_count = tumor_sample_count.reindex(
        index=common_genes, columns=common_cancers)
    normal_sample_count = normal_sample_count.reindex(
        index=common_genes, columns=common_cancers)
    # fisher z transformation
    tumor_corr_fz = tumor_corr.reindex(
        index=common_genes, columns=common_cancers).apply(lambda col: _fisherz(col))
    normal_corr_fz = normal_corr.reindex(
        index=common_genes, columns=common_cancers).apply(lambda col: _fisherz(col))
    # calculate delta fisher z score
    delta_fz, delta_pval = pd.DataFrame(), pd.DataFrame()
    for cancer in common_cancers:
        dz, pval = _delta_fisherz(
            tumor_corr_fz[cancer], normal_corr_fz[cancer], tumor_sample_count[cancer], normal_sample_count[cancer])
        delta_fz = pd.concat(
            [delta_fz, pd.DataFrame(dz, columns=[cancer])], axis=1)
        delta_pval = pd.concat(
            [delta_pval, pd.DataFrame(pval, columns=[cancer])], axis=1)
    delta_fz.columns = common_cancers
    delta_pval.columns = common_cancers
    return(delta_fz, delta_pval)


def _new_ranking(delta_corr, tumor_corr, ascending, cancer):
    """get rank product for delta correlation between tumor vs normal and tumor correlation"""
    common_genes = np.intersect1d(delta_corr.index, tumor_corr.index)
    delta_corr = delta_corr[common_genes].rank(ascending=ascending)
    tumor_corr = tumor_corr[common_genes].rank(ascending=ascending)
    rank_product = np.sqrt(delta_corr * tumor_corr)
    rank_product = pd.DataFrame(
        rank_product, index=rank_product.index, columns=[cancer])
    return(rank_product)


def _prioritizing(tumor_corr, delta_corr):
    """prioritize gian pos corr in tumor vs normal to the top, and gain negative corr in tumor to the bottom
    tumor_corr: dataframe where rows = genes, col = cancers
    delta_corr: dataframe where rows = genes, col = cancers
    """
    # do ranking based on:
    # 1. positively correlated in tumor
    # 2. gained correlation in tumor vs normal
    # tumor positive correlated genes
    new_order_all = pd.DataFrame()
    for cancer in np.intersect1d(tumor_corr.columns, delta_corr.columns).tolist():
        tumor_corr_coef = tumor_corr.loc[:, cancer]
        delta_corr_value = delta_corr.loc[:, cancer]
        # common genes after NA removal
        common_genes = np.intersect1d(
            tumor_corr_coef.dropna().index, delta_corr_value.dropna().index)
        tumor_corr_coef = tumor_corr_coef[common_genes]
        delta_corr_value = delta_corr_value[common_genes]
        # postively correlated genes
        pos_genes = tumor_corr_coef[tumor_corr_coef > 0].index
        pos_rank_product = _new_ranking(
            delta_corr_value[pos_genes], tumor_corr_coef[pos_genes], ascending=True, cancer=cancer)
        # neg genes
        neg_genes = tumor_corr_coef[tumor_corr_coef <= 0].index
        neg_rank_product = _new_ranking(
            delta_corr_value[neg_genes], tumor_corr_coef[neg_genes], ascending=False, cancer=cancer)
        # new order: rank product of tumor correlation and delta correlation between tumor vs normal
        new_order = pd.concat([pos_rank_product, -neg_rank_product])
        # new_order[cancer] = np.log2(
        #     np.abs(new_order[cancer])) * (new_order[cancer]/np.abs(new_order[cancer]))
        new_order_all = pd.concat([new_order_all, new_order], axis=1)

    return(new_order_all)  # dataframe where rows = gene, col = ranking score


def _refer_normal_delta_corr(tumor_res, normal_res, config):
    """run delta corr between tumor and normal"""
    info('refer normal is requested, so calculating delta corr between tumor vs normal')
    gene, adj_gene = config['gene'], config['adj_gene']
    info('for regular correlation result')
    regular_delta_corr, regular_delta_pval = _run_delta_between_tumor_normal(
        tumor_res['regular_corr'], normal_res['regular_corr'], tumor_res['sample_count'], normal_res['sample_count'])
    if type(regular_delta_corr) == type(None) or type(regular_delta_pval) == type(None):
        info('stop delta correlation')
        return
    # ranking score
    regular_new_order = _prioritizing(
        tumor_res['regular_corr'], regular_delta_corr)
    # write out
    regular_delta_corr.columns = regular_delta_corr.columns + '_delta_corr'
    regular_delta_pval.columns = regular_delta_pval.columns + '_delta_corr_pval'
    regular_new_order.columns = regular_new_order.columns + '_ranking_score'
    res = pd.concat(
        [regular_new_order, regular_delta_corr, regular_delta_pval], axis=1, sort=True).sort_values(regular_new_order.columns[0])
    res.to_csv(
        config['prefix'] + '_%s_regular_corr_signature_tumor_vs_normal.csv' % gene)
    # for partial correlation
    if adj_gene:
        partial_delta_corr, partial_delta_pval = _run_delta_between_tumor_normal(
            tumor_res['partial_corr'], normal_res['partial_corr'], tumor_res['sample_count'], normal_res['sample_count'])
        if type(partial_delta_corr) == type(None) or type(partial_delta_pval) == type(None):
            info('stop delta correlation')
            return
        # ranking score
        partial_new_order = _prioritizing(
            tumor_res['partial_corr'], partial_delta_pval)
        # write out
        partial_delta_corr.columns = partial_delta_corr.columns + '_delta_corr'
        partial_delta_pval.columns = partial_delta_pval.columns + '_delta_corr_pval'
        partial_new_order.columns = partial_new_order.columns + '_ranking_score'
        res2 = pd.concat(
            [partial_new_order, partial_delta_corr, partial_delta_pval], axis=1, sort=True)
        res2.to_csv(
            config['prefix'] + '_%s_regular_corr_signature_tumor_vs_normal.csv' % gene)


def _output(config, res, label='tumor'):
    """write result"""
    sc = res["sample_count"].copy()
    sc.columns = sc.columns + '_sample_size'
    for key in ['regular', 'partial']:
        corr = res['%s_corr' % key].copy()
        corr.columns = corr.columns + '_corr'
        pval = res['%s_pval' % key].copy()
        pval.columns = pval.columns + '_pval'
        if corr.shape[0] != 0 and pval.shape[0] != 0:
            corr_res = pd.concat([corr, pval, sc], sort=True, axis=1)
            info('corr output of %s of %s' % (key, label))
            corr_res.to_csv(
                config["prefix"] + '_%s_corr_%s_%s.csv' % (key, config['gene'], label))
        else:
            info('empty dataframe %s in %s' % (key, label))


def _run_advanced(config):
    """
    1. get correlation of given gene to all other protein genes
    2. run GSEA analysis using the ranked gene list
    """
    gene, adj_gene, refer_normal, protein_gene = config['gene'], config[
        'adj_gene'], config['refer_normal'], config['protein_gene']
    cancers, exp_mat = cf._test_data(config)
    # do tumor:
    info('calculating all gene correlation in tumor')
    tumor_res = _corr_all_gene(
        gene, protein_gene, adj_gene, exp_mat[exp_mat.sample_type == 'Tumor'])
    # returned:
    # 'regular_corr': regular_corr_all,
    # "regular_pval": regular_pval_all,
    # "partial_corr": partial_corr_all,
    # "partial_pval": partial_pval_all,
    # 'sample_count': pd.DataFrame(sample_count).T
    # write out tumor correlation result
    _output(config, tumor_res, 'tumor')
    # do normal
    if refer_normal:
        normal_res = _corr_all_gene(
            gene, protein_gene, adj_gene, exp_mat[exp_mat.sample_type == 'Normal'])
        # write out  normal corr result
        _output(config, normal_res, 'normal')
        # calculate fisher z delta between tumor and normal
        if normal_res['regular_corr'].shape[0] != 0:
            _refer_normal_delta_corr(tumor_res, normal_res, config)
        else:
            info('normal correlation matrix is empty, maybe too few sample for normal in given cancer, skip delta corr')
