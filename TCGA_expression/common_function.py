"""
#!/usr/bin/env python
#author          :Rongbin Zheng
#date            :2020-08-29
#python_version  :3.6 
"""

import pandas as pd


def info(mesg):
    # output message
    os.system('echo "++++ %s"' % mesg)


def _filter_by_cancer(df, cancer):
    # focus on cancer types
    all_cancers = df.project_id.unique().tolist()
    if cancer != 'all':
        cancer = cancer.split(',') if ',' in cancer else [cancer]
        df = df[df.project_id.isin(cancer)]
        if df.shape[0] == 0:
            info('No cancer type selected with %s' % str(all_cancers))
            sys.exit(1)
    else:
        cancer = all_cancers
    return(df, cancer)


def _test_data(config):
    """
    test whether the matrix is correct
    """
    cancer, exp_mat = config['cancer'], config['exp_mat']
    # filter cancer types
    df, cancers = _filter_by_cancer(exp_mat, cancer)
    return (cancers, df)


def _unique_patient(mat, gene_name):
    """
    mean expression based on patient ID, take mean expression for same patient
    """
    mat_ann = mat.loc[:, ~mat.columns.isin([gene_name])]
    mat_ann = mat_ann[~mat_ann['_PATIENT'].duplicated()]
    gene_exp = mat[[gene_name, '_PATIENT']].groupby('_PATIENT').mean()
    res = pd.merge(gene_exp, mat_ann, left_index=True, right_on='_PATIENT')
    return(res)
