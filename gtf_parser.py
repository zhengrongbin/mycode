import pandas as pd
import numpy as np
from gtfparse import read_gtf

## ===== hg19
df = read_gtf("gencode.v19.annotation.gtf")

df_genes = df[df["feature"] == "gene"]
df_genes.to_csv('gencode.v19.annotation.gene_annotation.csv', index=None)

df_genes_protein = df_genes[df_genes.gene_type == 'protein_coding']
df_genes_protein.to_csv('gencode.v19.annotation.protein_coding.csv', index=None)

#### ==== transcript and promoters
df_transcript = df[(df.feature == 'transcript') & (df.gene_type == 'protein_coding')]
df_transcript_positive = df_transcript[df_transcript.strand == '+']
df_transcript_positive = df_transcript_positive[['seqname', 'start', 'end', 'gene_name', 'strand', 'transcript_id']]
df_transcript_positive['start'] = df_transcript_positive['start'].astype('int') - 5000
df_transcript_positive['end'] = df_transcript_positive['start'].astype('int') + 5000


df_transcript_negative = df_transcript[df_transcript.strand == '-']
df_transcript_negative = df_transcript_negative[['seqname', 'start', 'end', 'gene_name', 'strand', 'transcript_id']]
df_transcript_negative['start'] = df_transcript_negative['end'].astype('int') - 5000
df_transcript_negative['end'] = df_transcript_negative['end'].astype('int') + 5000

df_transcript_promoter_protein = pd.concat([df_transcript_positive, df_transcript_negative])
df_transcript_promoter_protein = df_transcript_promoter_protein.sort_values(['seqname', 'start'])
df_transcript_promoter_protein.to_csv('gencode.v19.annotation.protein_coding.transcript_promoter_UpDown5Kb.csv', 
    index = None, header = None, sep = '\t')

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
gene_length.to_csv('gencode.v19.annotation.gene_cds_length.csv')


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

df = read_gtf("gencode.vM23.annotation.gtf")

df_genes = df[df["feature"] == "gene"]
df_genes.to_csv('gencode.vM23.annotation.gene_annotation.csv', index=None)

df_genes_protein = df_genes[df_genes.gene_type == 'protein_coding']
df_genes_protein.to_csv('gencode.vM23.annotation.protein_coding.csv', index=None)

#### ==== transcript and promoters
df_transcript = df[(df.feature == 'transcript') & (df.gene_type == 'protein_coding')]
df_transcript_positive = df_transcript[df_transcript.strand == '+']
df_transcript_positive = df_transcript_positive[['seqname', 'start', 'end', 'gene_name', 'strand', 'transcript_id']]

df_transcript_negative = df_transcript[df_transcript.strand == '-']
df_transcript_negative = df_transcript_negative[['seqname', 'start', 'end', 'gene_name', 'strand', 'transcript_id']]

## TSS
df_transcript_positive['tss'] = df_transcript_positive['start'].astype('int').tolist()
df_transcript_negative['tss'] = df_transcript_negative['end'].astype('int').tolist()

df_transcript_positive['tss+1'] = df_transcript_positive['tss']+1
df_transcript_negative['tss+1'] = df_transcript_negative['tss']+1

df_tss_promoter_protein = pd.concat([df_transcript_positive, df_transcript_negative]).reindex(columns = ['seqname', 'tss', 'tss+1', 'gene_name', 'transcript_id', 'strand'])
df_tss_promoter_protein['gene_name'] = df_tss_promoter_protein['transcript_id']+':'+df_tss_promoter_protein['gene_name']
df_tss_promoter_protein = df_tss_promoter_protein.sort_values(['seqname', 'tss'])
df_tss_promoter_protein = df_tss_promoter_protein[(df_tss_promoter_protein['tss'] > 0)]
df_tss_promoter_protein['transcript_id'] = 0
df_tss_promoter_protein.to_csv('gencode.vM23.annotation.protein_coding.tss.csv', 
    index = None, header = None, sep = '\t')

## promoter 5Kb
# df_transcript_positive['tss'] = df_transcript_positive['start'].astype('int')
df_transcript_positive['start'] = df_transcript_positive['start'].astype('int') - 5000
df_transcript_positive['end'] = df_transcript_positive['start'].astype('int') + 5000

# df_transcript_negative['tss'] = df_transcript_negative['end'].astype('int')
df_transcript_negative['start'] = df_transcript_negative['end'].astype('int') - 5000
df_transcript_negative['end'] = df_transcript_negative['end'].astype('int') + 5000

df_transcript_promoter_protein = pd.concat([df_transcript_positive, df_transcript_negative])
df_transcript_promoter_protein = df_transcript_promoter_protein.sort_values(['seqname', 'start'])
df_transcript_promoter_protein = df_transcript_promoter_protein[(df_transcript_promoter_protein['start'] > 0) &
                                (df_transcript_promoter_protein['end'] > 0)]
df_transcript_promoter_protein.to_csv('gencode.vM23.annotation.protein_coding.transcript_promoter_UpDown5Kb.csv', 
    index = None, header = None, sep = '\t')


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
gene_length.to_csv('gencode.vM23.annotation.gene_cds_length.csv')




