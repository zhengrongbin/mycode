import pandas as pd
import numpy as np
from gtfparse import read_gtf

## ===== hg38
df = read_gtf("gencode.v38.annotation.gtf")

df_genes = df[df["feature"] == "gene"]
df_genes.to_csv('gencode.v38.annotation.gene_annotation.csv', index=None)

df_genes_protein = df_genes[df_genes.gene_type == 'protein_coding']
df_genes_protein.to_csv('gencode.v38.annotation.protein_coding.csv', index=None)

## CDS to gene length
cds = df[df['feature'] == 'CDS']
cds['location'] = cds['seqname']+'_'+cds['start'].astype('str')+'_'+cds['end'].astype('str')

cds_length = {}

def get_length(tmp):
    tmp = tmp[~tmp.location.duplicated()]
    tmp = tmp.sort_values('start')
    return(sum(abs(tmp['start'] - tmp['end'])))

gene_length = cds.groupby('gene_name')[['start', 'end', 'location']].apply(lambda df: get_length(df))
gene_length = pd.DataFrame(gene_length, columns = ['CDS_bp'])
gene_length.to_csv('gencode.v38.annotation.gene_cds_length.csv')

## ========= mm10

df = read_gtf("gencode.vM27.annotation.gtf")

df_genes = df[df["feature"] == "gene"]
df_genes.to_csv('gencode.vM27.annotation.gene_annotation.csv', index=None)

df_genes_protein = df_genes[df_genes.gene_type == 'protein_coding']
df_genes_protein.to_csv('gencode.vM27.annotation.protein_coding.csv', index=None)

## === CDS length
cds = df[df['feature'] == 'CDS']
cds['location'] = cds['seqname']+'_'+cds['start'].astype('str')+'_'+cds['end'].astype('str')

cds_length = {}

def get_length(tmp):
    tmp = tmp[~tmp.location.duplicated()]
    tmp = tmp.sort_values('start')
    return(sum(abs(tmp['start'] - tmp['end'])))

gene_length = cds.groupby('gene_name')[['start', 'end', 'location']].apply(lambda df: get_length(df))
gene_length = pd.DataFrame(gene_length, columns = ['CDS_bp'])
gene_length.to_csv('gencode.v38.annotation.gene_cds_length.csv')
