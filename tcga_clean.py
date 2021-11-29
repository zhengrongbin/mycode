import pandas as pd
import numpy as np
import os,sys
import collections 


tcga_dat_path = '/data/home/maj/applications/gdc-client/data/process_matrix/tcga_gdc.merged.fpkm_convert.to_log2tpm0.001.csv.gz'
clinic_path = '/data2/zhengrongbin/data/tcga_clear/clinical.tsv'
survival_path = '/data1/TCGA_xena/pancancer/GDC-PANCAN.survival.tsv'
purity_path = '/data2/zhengrongbin/data/tcga_clear/TCGA_mastercalls.abs_tables_JSedit.fixed.txt'
tmb_path = '/data2/zhengrongbin/data/tcga_clear/mutation/TCGA_sample_mutation_burden_mutect2.txt'
basic_pheo_path = '/data1/TCGA_xena/pancancer/GDC-PANCAN.basic_phenotype.tsv'
gene_ann_path = '/data2/zhengrongbin/data/annotation/gene_annotation_GENCODE.v22.csv'
subtype_path = '/data2/zhengrongbin/data/tcga_clear/TCGASubtype.20170308.tsv'

clinic = pd.read_csv(clinic_path, sep = '\t')
survival = pd.read_csv(survival_path, sep = '\t')
purity = pd.read_csv(purity_path, sep = '\t')
tmb = pd.read_csv(tmb_path, sep = '\t')
basic_pheno = pd.read_csv(basic_pheo_path, sep = '\t')
gene_ann = pd.read_csv(gene_ann_path)
subtype = pd.read_csv(subtype_path, sep = '\t')

##=============
tcga_dat = pd.read_csv(tcga_dat_path, compression = 'gzip', index_col = 0)
## only protein coding genes
protein_genes = gene_ann[gene_ann.gene_type == 'protein_coding'].gene_name
tcga_dat = tcga_dat.loc[tcga_dat.index.isin(protein_genes),:]
## convert log2(TPM+0.001) to log2(TPM + 1)
tcga_dat = tcga_dat.apply(lambda col: np.log2(np.exp2(col) - 0.001 + 1))
## ========= match sample annotation
#### correct sample ID for basic_pheno
basic_pheno['sample'] = [x[:-1] for x in basic_pheno['sample'].tolist()]
basic_pheno = basic_pheno[basic_pheno.program == 'TCGA']
basic_pheno = basic_pheno[~basic_pheno['sample'].duplicated()]
# basic_pheno.index = basic_pheno['sample'].tolist()

#### sample_type to Tumor and Normal
basic_pheno['sampleType'] = basic_pheno['sample_type'].tolist()
# ['Primary Tumor', 'Solid Tissue Normal', 'Metastatic',
#        'Recurrent Tumor',
#        'Primary Blood Derived Cancer - Peripheral Blood',
#        'Additional - New Primary', 'FFPE Scrolls',
#        'Additional Metastatic', 'Buccal Cell Normal',
#        'Bone Marrow Normal']
def _assign_new_sampleType(x):
    if x in ['Buccal Cell Normal', 'Bone Marrow Normal', 'Solid Tissue Normal']:
        return('Normal')
    else:
        return('Tumor')
basic_pheno['sample_type'] = [_assign_new_sampleType(x) for x in basic_pheno['sample_type'].tolist()]
## merge with other information
#### survival
survival['sample'] = [x[:-1] for x in survival['sample'].tolist()]
survival = survival[survival['sample'].str.startswith('TCGA')]
survival = survival[~survival['sample'].duplicated()]
##### merge survival
basic_pheno = pd.merge(basic_pheno[['sample', 'project_id', 'sample_type', 'sampleType', 'Age at Diagnosis in Years', 'Gender']],
                        survival, left_on = 'sample', right_on = 'sample')
basic_pheno.columns = ['sample', 'project_id', 'sample_type', 'sampleType', 'Age', 'Gender', 'OS', '_PATIENT', 'OS.time']
#### purity
basic_pheno = pd.merge(basic_pheno, purity[['array', 'purity']], left_on = 'sample', right_on = 'array', how = 'left').drop('array', axis = 1)
#### TMB
tmb['sample'] = [x[:-1] for x in tmb['sample'].tolist()]
tmb = tmb.groupby('sample')[['TMB']].mean()
basic_pheno = pd.merge(basic_pheno, tmb, left_on = 'sample', right_index = True, how = 'left')
#### stage
clinic = clinic[~clinic['case_submitter_id'].duplicated()]

def Normzlize_Stage(x):
    x = x.replace('Stage ', '')
    if x == '0':
        return('0')
    elif x in ['I', 'IA', 'IB', 'IC', 'IS']:
        return('I')
    elif x in ['II', 'IIA', 'IIB', 'IIC']:
        return('II')
    elif x in ['III', 'IIIA', 'IIIB', 'IIIC']:
        return('III')
    elif x in ['IV', 'IVA', 'IVB', 'IVC']:
        return('IV')
    elif x in ['X']:
        return('X')
    else:
        return(np.nan)
    
stage_index = pd.Series(clinic['ajcc_pathologic_stage'].tolist(), index = clinic['case_submitter_id'].tolist())
stage_index_new = stage_index.apply(lambda x: Normzlize_Stage(x))
basic_pheno['Stage'] = [stage_index_new[x] for x in basic_pheno['_PATIENT'].tolist()]

## subtype
subtype = subtype[subtype['Subtype_mRNA'] != 'Normal']
subtype['patient'] = [x[:12] if len(x)>12 else x for x in subtype['pan.samplesID'].tolist()]
###### checked that 11 patient with duplicated subtype annotation, but most are normal or different center, and subtype is the same, so it is not a big matter
subtype = subtype[~(subtype['Subtype_Selected'].str.endswith('.NA') | subtype['Subtype_Selected'].str.endswith('.-'))]
subtype = subtype[~subtype.patient.duplicated()]
subtype_patient = pd.Series(subtype['Subtype_Selected'].tolist(), index = subtype['patient'].tolist())
basic_pheno['Subtype'] = [subtype_patient[x] if x in subtype_patient.index.tolist() else np.nan for x in basic_pheno['_PATIENT'].tolist()]

## merge meta annotation to expression data
tcga_dat = tcga_dat.T
tcga_dat = pd.merge(tcga_dat, basic_pheno, left_index = True, right_on = 'sample')
tcga_dat.index = tcga_dat['sample'].tolist()
### update project_id to BRCA subtype and SKCM primary and metastatic
##### SKCM
skcm_subtypes = tcga_dat.loc[tcga_dat.project_id == 'TCGA-SKCM', 'sampleType']
skcm_subtypes.index = tcga_dat.loc[tcga_dat.project_id == 'TCGA-SKCM', 'sample']
skcm_subtypes = skcm_subtypes[skcm_subtypes != 'Solid Tissue Normal']
skcm_subtypes[skcm_subtypes=='Additional Metastatic'] = 'TCGA-SKCM-Metastatic'
skcm_subtypes[skcm_subtypes=='Metastatic'] = 'TCGA-SKCM-Metastatic'
skcm_subtypes[skcm_subtypes=='Primary Tumor'] = 'TCGA-SKCM-Primary'

###### add SKCM subtype
tcga_dat.loc[skcm_subtypes.index, 'project_id'] = skcm_subtypes.tolist()
######### SKCM normal
tcga_dat.loc[(tcga_dat.project_id == 'TCGA-SKCM'), 'project_id'] = 'TCGA-SKCM-Normal'

#### BRCA
brca_subtypes = tcga_dat.loc[(tcga_dat.project_id == 'TCGA-BRCA') & (tcga_dat.sample_type == 'Tumor'), 'Subtype']
brca_subtypes = brca_subtypes.str.replace('BRCA.', 'TCGA-BRCA-')
brca_subtypes = brca_subtypes.dropna() ## same BRCA sample were considered as Normal-like based on mRNA and the labelis NA
###### add BRCA subtype
tcga_dat.loc[brca_subtypes.index, 'project_id'] = brca_subtypes.tolist()
######### BRCA normal
tcga_dat.loc[(tcga_dat.project_id == 'TCGA-BRCA'), 'project_id'] = 'TCGA-BRCA-Normal'


tcga_dat.to_csv('TCGA_GDC_fpkm_to_log2tpm1_add_surivial_age_stage_subtype_purity_tmb.csv.gz', compression = 'gzip')



